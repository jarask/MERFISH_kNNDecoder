import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import scipy.io
from scipy.spatial.distance import hamming
from timeit import default_timer as timer
from tqdm import tqdm
import Codebook

pd.options.mode.chained_assignment = None  # default='warn'


def map_frames_to_bin(frames: set, barcode_length: int):
    """
    Maps a given set of frames into a barcode with a given length.
    The frames in which the spots are found will be set to 1. Otherwise there will be inserted a 0
    Note: the frames set should be modified to be zero-indexed!
    :param frames: Set of frames (Zero-indexed!)
    :param barcode_length: The length of the barcode as an integer
    :return: Binary string of the barcode
    """
    binary_string = ""
    for bit in range(barcode_length):
        if bit not in frames:
            binary_string += str(0)
        else:
            binary_string += str(1)
    return binary_string


class kNNDecoder:
    """
    Class for decoding a MERFISH experiment which have been processed by ThunderSTORM
    """
    def __init__(self, filename: str, codebook_name: str, barcode_length: int, scale_factor: int,
                 verbose: object = True) -> None:
        """

        :param filename:
        :param codebook_name:
        :param barcode_length:
        :param scale_factor:
        :param verbose:
        """
        self.filename = filename
        self.codebook_name = codebook_name
        self.barcode_length = barcode_length
        self.scale_factor = scale_factor
        self.verbose = verbose
        self.codebook = Codebook.Codebook(self.codebook_name)
        self.data = None
        self.run_time = None
        # start of execution time
        self.start_time = timer()
        try:
            # Load the file (currently only able to read from ThunderSTORM files)
            self.data = pd.read_csv(self.filename)
            self.data = self.data.loc[:, ['id', 'frame', 'x [nm]', 'y [nm]']]
            # self.data[['frame']] = self.data[['frame']] - 1  # Zero-indexing
            # Convert from nm to pixels
            self.data[['x [nm]', 'y [nm]']] = self.data[['x [nm]', 'y [nm]']] / self.scale_factor

            # For some reason, flip the x y coords...
            # TODO: Investigate the 'flipping' some more
            self.data[['x [nm]', 'y [nm]']] = self.data[['y [nm]', 'x [nm]']]

            self.data.frame = self.data.frame.astype(int)
            self.data.frame = self.data.frame - 1
        except IOError:
            IOError("Error: File does not appear to exist")

    def run_knn_decoding(self):
        """

        :return:
        """
        # TODO: Find out where this could be implemented better
        # self.codebook = self.Codebook(self.codebook_name)

        # Create a df of only XY coords for NearestNeighbors fit
        data_XY = self.data.loc[:, ['y [nm]', 'x [nm]']]
        # Nearest neighbors attempt
        knn = NearestNeighbors(n_neighbors=4, p=2, radius=0.5)  # Includes itself in n_neighbors
        knn.fit(data_XY)
        distance, indices = knn.kneighbors(data_XY)
        indices.sort()  # Should not be done when you look at the distances
        self.data['knn'] = indices.tolist()
        # self.data['knn_distances'] = distance.tolist()  # Can be used to optimize the distances

        # Validate knn groups
        kNNDecoder._knn_validation(self)

        # Once we have the confirmed knn_groups (and the others), we can add the frames of the neighbours to the df
        nbrs_list_master = []
        for nbrs in self.data.knn:
            nbrs_list = []
            for nbr in nbrs:
                nbrs_list.append(self.data.loc[nbr].frame)
            # nbrs_list.sort(key=float)  # Sorting the nbrs_list could make sense?
            nbrs_list_master.append(nbrs_list)
        self.data['knn_frames'] = nbrs_list_master

        # Using the frames of knn_frames, we can now assign the barcode to each knn_group
        barcode_set = []
        for frames_set in tqdm(self.data.knn_frames, desc='Assigning barcodes'):
            barcode = map_frames_to_bin(frames_set, self.barcode_length)
            barcode_set.append(barcode)
        self.data['barcodes'] = barcode_set

        # Correct the correctable barcodes
        kNNDecoder._knn_correction(self)

        # Add gene information to the points
        # self.data['nucleotide_ID'] = ''
        self.data['gene'] = ''

        # converts codebook list to array for easier handling
        row_id = 0
        for barcode in self.data.barcodes:
            for gene in range(len(self.codebook.genes)):
                if np.all(barcode == self.codebook.code_list[gene]):
                    # self.data.loc[row_id, 'nucleotide_ID'] = self.codebook.genes[gene]
                    # self.data.loc[row_id, 'gene'] = gene
                    self.data.loc[row_id, 'gene'] = self.codebook.genes[gene]
            row_id += 1

        # Set the 'unconfirmed' as NaN in gene and nucleotide_ID
        for point_index in tqdm(range(len(self.data)), desc='Set unconfirmed as NaN'):
            if self.data.loc[point_index, 'status'] == 'unconfirmed' or self.data.loc[point_index, 'gene'] == '':
                self.data.loc[point_index, 'status'] = 'unconfirmed'
                # self.data.loc[point_index, ('nucleotide_ID', 'gene')] = np.nan
                self.data.loc[point_index, 'gene'] = np.nan

        # Sort the data
        self.data = self.data.sort_values(by=['frame', 'x [nm]'], ascending=True)
        self.data = self.data.reset_index(drop=True)

        # stop timer
        self.run_time = timer() - self.start_time

        if self.verbose:
            print('A total of {} spots were confirmed as RNA spots.'.format(str(sum(self.data.status == 'confirmed'))))
            print('{} were corrected, and {} were unconfirmed.'.format(str(sum(self.data.status == 'corrected')),
                                                                       str(sum(self.data.status == 'unconfirmed'))))
            print(self.data.gene.value_counts().sort_index())

            print('Decoding done in {} sec.'.format(str(round(self.run_time, 2))))
            print('----------')

        return self.data

    def _knn_validation(self):
        """

        """
        # Set all points to 'unconfirmed' for knn validation
        self.data['status'] = 'unconfirmed'
        # kNN validation
        index_counter = 0
        for knn_group in tqdm(self.data.knn, desc='Validating data'):
            # The group should only be checked, if its status is unconfirmed
            if self.data.loc[index_counter].status == 'unconfirmed':
                knn_group_comparison = []
                # Get the knn_group of all members in the current knn_group
                for item in knn_group:
                    knn_group_comparison.append(self.data.loc[item].knn)
                # Check if the four knn_group are the same for the members
                if all(elem == knn_group_comparison[0] for elem in knn_group_comparison):
                    for item in knn_group:
                        self.data.loc[item, 'knn_set_group'] = int(4)
                        self.data.loc[item, 'status'] = 'confirmed'
            index_counter += 1

        # Get the unconfirmed points and see if any of them can match, and how many they are
        unconfirmed_ts = self.data[self.data.status == 'unconfirmed']
        knn_set_counter = []
        for knn_set in tqdm(unconfirmed_ts.knn, desc='Going through unconfirmed points'):
            knn_set_group = 0
            for other_knn_set in unconfirmed_ts.knn:
                if knn_set == other_knn_set:
                    knn_set_group += 1
            knn_set_counter.append(knn_set_group)
        unconfirmed_ts.loc[:, 'knn_set_group'] = knn_set_counter
        # Remove the confirmed indices in these knn_groups
        remove_indices = []
        for knn_set in unconfirmed_ts.knn:
            for knn_item in knn_set:
                if all(knn_item != unconfirmed_ts.index):
                    remove_indices.append(knn_item)
        for k in range(len(unconfirmed_ts.knn)):
            for r in remove_indices:
                if r in unconfirmed_ts.knn.iloc[k]:
                    unconfirmed_ts.knn.iloc[k].remove(r)

        # Update the original df
        self.data.update(unconfirmed_ts)

    def _knn_correction(self):
        """

        """
        # For the points with 3 knn-members, correct their barcode
        # TODO: There should be some kind of threshold for the distance?
        correctable_points = self.data[self.data.knn_set_group == 3]
        for cor_point_index in tqdm(range(len(correctable_points)), desc='Correcting data'):
            hamming_dist_list = []
            for i in range(len(self.codebook.code_list)):
                hamming_dist_list.append(hamming(
                    np.array(list(correctable_points.iloc[cor_point_index].barcodes)).astype(int),
                    np.array(list(self.codebook.code_list[i])).astype(int))
                                         * len(np.array(list(correctable_points.iloc[cor_point_index].barcodes))))
            index_min = np.argmin(hamming_dist_list)
            # print('The barcode had a hamming distance difference of ' + str(hamming_dist_list[index_min]))
            if hamming_dist_list[index_min] == 1:
                correctable_points.barcodes.iloc[cor_point_index] = self.codebook.code_list[index_min]
                correctable_points.status.iloc[cor_point_index] = 'corrected'

        # Update the original df
        self.data.update(correctable_points)
