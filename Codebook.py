"""

"""
import pandas as pd
import scipy.io


class Codebook:
    """
    Class for containing and handling the codebook.
    """

    def __init__(self, codebook_name: str):
        """
        :param codebook_name: Name of the codebook. Should be placed in the Codebooks/ folder
        """
        self.name = 'Codebooks/' + codebook_name
        # Attempt at a system for taking multiple types of codebooks
        if self.name.lower().endswith('.mat'):
            try:
                codebook = scipy.io.loadmat(self.name, squeeze_me=True)
                # flattens the mat-file to a list
                self.genes = [(k, v) for k, v in codebook.items()][3][1]
                self.code_list = self.genes['Code']
            except IOError:
                IOError("Error: File does not appear to exist")
        elif self.name.lower().endswith('.csv'):
            try:
                codebook = pd.read_csv(self.name, header=None, dtype=object)
                if codebook.shape[1] == 1:
                    if ';' in codebook[0][0]:
                        codebook = pd.read_csv(self.name, header=None, sep=';', dtype=object)
                    elif ',' in codebook[0][0]:
                        codebook = pd.read_csv(self.name, header=None, sep=',', dtype=object)
                elif codebook.shape[1] == 2:
                    pass
                else:
                    print('The CSV file does not contain comma or semi-colon as separator..')
                self.genes = codebook.loc[:, 0].to_list()
                self.code_list = codebook.loc[:, 1].to_list()
            except IOError:
                IOError("Error: File does not appear to exist")
