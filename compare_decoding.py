import numpy as np
import pandas as pd
import glob


def compare_decoding(amplitude, points_per_gene, bkg, sigma, search_tolerance):
    # set the merfish experiment parameters
    exp_name = 'a%s_p%s_bkg%s' % (amplitude, points_per_gene, bkg)

    # Get the number of detections from ThunderSTORM
    ts_output_df = pd.read_csv('data/p%s_sigma%s/%s/tsoutput.csv' % (points_per_gene, sigma, exp_name))
    len_ts_output = len(ts_output_df)
    del ts_output_df

    knn_df = pd.read_csv('data/p%s_sigma%s/%s/knn_decoding.csv' % (points_per_gene, sigma, exp_name), index_col=0)
    knn_df = knn_df.loc[knn_df.status.isin(['confirmed', 'corrected'])]
    knn_df = knn_df.reset_index()
    del knn_df['index']
    knn_df.knn = knn_df.knn.astype('str')
    knn_df = knn_df.iloc[knn_df[['knn', 'gene']].drop_duplicates().index]
    knn_df = knn_df[['x [nm]', 'y [nm]', 'gene']]
    knn_df.columns = ['y', 'x', 'target']
    data_list = {'knn': knn_df}

    for (index, key) in enumerate(data_list):
        if data_list[key].empty:
            data_list[key][['y', 'x', 'target']] = np.nan
        data_list[key] = data_list[key].sort_values(by='y')
        data_list[key]['used'] = False

    # Load in the ground truth
    gt_df = pd.read_csv('data/p%s_sigma%s/ground_truth.csv' % (points_per_gene, sigma), index_col=0)
    # Get the single features (their 1st occurrence) and sort by y
    gt_features = gt_df.iloc[gt_df[['y', 'x', 'target']].drop_duplicates().index]
    gt_features.sort_values(by='y', inplace=False)
    gt_features.reset_index(inplace=True)
    del gt_features['index']

    # Set up comparison dataframe
    compare_df = gt_features[['y', 'x', 'target', 'frame']]
    compare_df.reset_index(inplace=True)
    del compare_df['index']
    compare_df = compare_df.rename(columns={'y': 'y_gt', 'x': 'x_gt', 'target': 'target_gt'})
    additional_cols = []
    for analysis_name in data_list:
        additional_cols.append('y_%s' % analysis_name)
        additional_cols.append('x_%s' % analysis_name)
        additional_cols.append('target_%s' % analysis_name)
        additional_cols.append('target_match_%s' % analysis_name)
        additional_cols.append('classification_%s' % analysis_name)
    compare_df = compare_df.reindex(columns=compare_df.columns.to_list() + additional_cols)

    # Go through the dataframes and compare them to gt
    tolerance_radius = search_tolerance  # The radius of where to align coordinates
    for (index_df, key) in enumerate(data_list):
        # print('----- %s -----' % key)
        for (row_index, row_value) in data_list[key].iterrows():
            for (row_gt_index, row_gt_value) in gt_features.iterrows():
                if (((row_gt_value.y - tolerance_radius) <= row_value.y <= (row_gt_value.y + tolerance_radius)) and
                        ((row_gt_value.x - tolerance_radius) <= row_value.x <= (row_gt_value.x + tolerance_radius))):
                    # Only set the values in, if there isn't any values here yet - Duplicates will be handled later
                    if str(compare_df.loc[row_gt_index, 'target_%s' % key]) == 'nan':
                        compare_df.loc[row_gt_index, ('y_%s' % key, 'x_%s' % key, 'target_%s' % key)] = \
                            row_value.loc[['y', 'x', 'target']].tolist()
                    else:
                        # Put the values in the bottom of the df
                        compare_df.loc[len(compare_df), ('y_%s' % key, 'x_%s' % key, 'target_%s' % key)] = \
                            row_value.loc[['y', 'x', 'target']].tolist()
                    data_list[key].loc[row_index, 'used'] = True

            # If it doesn't match, then put it in the end of the df and mark as FP
            if not data_list[key].loc[row_index, 'used']:
                compare_df.loc[len(compare_df), ('y_%s' % key, 'x_%s' % key, 'target_%s' % key)] = \
                    row_value.loc[['y', 'x', 'target']].tolist()

    # Assign matching status to the ground truth aligned data
    for (index_df, key) in enumerate(data_list):
        for (row_index, row_value) in compare_df.iterrows():
            # Only do it for the ground truth values
            if row_index < len(gt_features):
                if compare_df.loc[row_index, 'target_gt'] == compare_df.loc[row_index, 'target_%s' % key]:
                    compare_df.loc[row_index, 'target_match_%s' % key] = 'match'
                elif str(compare_df.loc[row_index, 'target_%s' % key]) == 'nan':
                    compare_df.loc[row_index, 'target_match_%s' % key] = 'not_found'
                else:
                    compare_df.loc[row_index, 'target_match_%s' % key] = 'no_match'
            # Assign the 'out-ofbounds' as 'no_match'
            if (str(compare_df.loc[row_index, 'target_%s' % key]) != 'nan') and \
                    (str(compare_df.loc[row_index, 'target_match_%s' % key]) == 'nan'):
                compare_df.loc[row_index, 'target_match_%s' % key] = 'no_match'
            # Assign classification
            # TP = 'match'
            # FN = 'not_found'
            # FP = nan (the "out-of-bounds")
            if compare_df.loc[row_index, 'target_match_%s' % key] == 'match':
                compare_df.loc[row_index, 'classification_%s' % key] = 'TP'
            elif compare_df.loc[row_index, 'target_match_%s' % key] == 'not_found':
                compare_df.loc[row_index, 'classification_%s' % key] = 'FN'
            elif compare_df.loc[row_index, 'target_match_%s' % key] == 'no_match':
                compare_df.loc[row_index, 'classification_%s' % key] = 'FP'

    # Make a summarized df
    sum_df = pd.DataFrame([], [['Total features', 'TP', 'FN', 'FP', 'Recall', 'FDR', 'F1']],
                          [['gt'] + list(data_list.keys())])
    # Set the gt data first
    sum_df.loc[:, 'gt'] = (len(gt_features), 0, 0, 0, 0, 0, 0)
    sum_df.loc[:, 'ts'] = (len_ts_output, 0, 0, 0, 0, 0, 0)
    true_p = len(gt_features)
    for (index_df, key) in enumerate(data_list):
        len_df = len(data_list[key])
        num_tp = sum(compare_df.loc[:, 'classification_%s' % key] == 'TP')
        num_fn = sum(compare_df.loc[:, 'classification_%s' % key] == 'FN')
        num_fp = sum(compare_df.loc[:, 'classification_%s' % key] == 'FP')
        sum_df.loc[:, key] = (len_df, num_tp, num_fn, num_fp, 0, 0, 0)
        sum_df.loc['Recall', key] = num_tp / true_p
        if len_df == 0:
            sum_df.loc['FDR', key] = 0
        else:
            sum_df.loc['FDR', key] = num_fp / len_df
        sum_df.loc['F1', key] = (2 * num_tp) / (2 * num_tp + num_fp + num_fn)

    # Save the two dfs
    compare_df.to_csv('data/p%s_sigma%s/%s/comparison.csv' % (points_per_gene, sigma, exp_name))
    sum_df.to_csv('data/p%s_sigma%s/%s/summarized_data.csv' % (points_per_gene, sigma, exp_name))


if __name__ == '__main__':
    compare_decoding(amplitude=5, points_per_gene=1, bkg=0, sigma=2, search_tolerance=4)
