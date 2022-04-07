import knn_decoding
import compare_decoding
import pandas as pd
import numpy as np
import glob
import os
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt


if __name__ == '__main__':
    points_per_gene = 50
    sigma = 2
    dir_list = glob.glob('data/p%s_sigma%s/*' % (points_per_gene, sigma) + os.path.sep)

    data_exists = True

    if not data_exists:
        for (index, dir_name) in tqdm(enumerate(dir_list), total=len(dir_list), desc='Decoding and analyzing experiments'):
            short_dir_name = dir_name[15+(len(str(points_per_gene))-1):len(dir_name) - 1]
            short_dir_name = short_dir_name.split('_')
            amplitude = short_dir_name[0][1:]
            points_per_gene = short_dir_name[1][1:]
            if len(short_dir_name) == 4:
                short_dir_name[2] = '%s_%s' % (short_dir_name[2], short_dir_name[3])
            bkg = short_dir_name[2][3:]
            # Before decoding, check if the tsoutput file is empty or not
            ts_output = pd.read_csv('%s/tsoutput.csv' % dir_name)
            if len(ts_output) < 2:
                # Create an empty file as the df
                df = pd.DataFrame([], columns=('id', 'frame', 'x [nm]', 'y [nm]', 'knn', 'status',
                                               'knn_set_group', 'knn_frames', 'barcodes', 'gene'))
            else:
                # Run decoding
                knn_decoder = knn_decoding.kNNDecoder('%s/tsoutput.csv' % dir_name,
                                                      'Codebook.csv', 16, 1, verbose=False)
                df = knn_decoder.run_knn_decoding()
            df.to_csv('%s/knn_decoding.csv' % dir_name)

            # Compare decoding
            compare_decoding.compare_decoding(amplitude=amplitude, points_per_gene=points_per_gene,
                                              bkg=bkg, sigma=sigma, search_tolerance=3)

    # Read in all summarized csv files and save them to a dict
    sum_data_dict = {}
    for (index, dir_name) in enumerate(dir_list):
        sum_df = pd.read_csv('%s/summarized_data.csv' % dir_name, index_col=0)
        short_dir_name = dir_name[15 + (len(str(points_per_gene)) - 1):len(dir_name) - 1]
        short_dir_name = short_dir_name.split('_')
        amplitude = short_dir_name[0][1:]
        # if len(short_dir_name) == 4:
        #     short_dir_name[2] = 'bkgno readout \n noise'
        bkg = short_dir_name[2][3:]
        name = 'amp%s bkg%s' % (amplitude, bkg)
        sum_data_dict[name] = {'amp': amplitude, 'bkg': bkg, 'data': sum_df}

    # Go through each experiment and add them to another dict
    amp_dict = {}
    # Add gt and its total to the dict
    amp_dict['ground truth'] = {'Total targets': sum_data_dict[next(iter(sum_data_dict))]['data'].loc[:, 'gt'][0]}
    amp_list = ['1', '5', '10', '100']
    raw_analyses = ['knn']
    raw_df = pd.DataFrame([], columns=['Amplitude', 'Background level', 'Analysis',
                                       'Total features', 'Recall', 'FDR', 'F1'])
    ts_df = pd.DataFrame([], columns=['Amplitude', 'Background level', 'Particles detected'])
    for amp in amp_list:
        for (index, key) in enumerate(sum_data_dict):
            if sum_data_dict[key]['amp'] == amp:
                temp_df = sum_data_dict[key]['data']
                bkg = sum_data_dict[key]['bkg']
                for raw_col in raw_analyses:
                    raw_df.loc[len(raw_df)] = [amp, bkg, raw_col] + \
                                              list(temp_df.loc[['Total features', 'Recall', 'FDR', 'F1'], raw_col])
                    ts_df.loc[len(ts_df)] = [amp, bkg] + list(temp_df.loc[['Total features'], 'ts'])

    amp_dict['metadata'] = {'points_per_gene': points_per_gene,
                            'sigma': sigma,
                            'amp_list': amp_list,
                            'bkg_order': ['0', '1', '10', '50', '100', '1000']}

    # Plot data
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 7))
    sns.barplot(ax=ax1, x='Amplitude', y='Particles detected', hue='Background level', data=ts_df,
                hue_order=amp_dict['metadata']['bkg_order'], palette='viridis')
    sns.barplot(ax=ax2, x='Amplitude', y='F1', hue='Background level', data=raw_df,
                hue_order=amp_dict['metadata']['bkg_order'], palette='viridis')
    sns.barplot(ax=ax3, x='Amplitude', y='Recall', hue='Background level', data=raw_df,
                hue_order=amp_dict['metadata']['bkg_order'], palette='viridis')
    sns.barplot(ax=ax4, x='Amplitude', y='FDR', hue='Background level', data=raw_df,
                hue_order=amp_dict['metadata']['bkg_order'], palette='viridis')
    for ax in (ax1, ax2, ax3, ax4):
        ax.get_legend().remove()
    fig.suptitle('P=%s | Sigma=%s' % (points_per_gene, sigma), fontweight='bold')
    ax1.set_title('ThunderSTORM detections')
    for ax in (ax2, ax3, ax4):
        ax.set_title('kNN decoder performance')
        ax.set_ylim(0, 1)
    ax1.axhline(y=amp_dict['ground truth']['Total targets']*4, color='r', linestyle='--')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    ax2.legend(bbox_to_anchor=(0.72, 0.95), bbox_transform=plt.gcf().transFigure, borderaxespad=0.,
               ncol=6, prop={'size': 9}, title='Background level')
    fig.savefig('data/Figures/p%s_sigma%s.png' % (points_per_gene, sigma), dpi=1000)
    # plt.show()
    print('Figure saved!')
