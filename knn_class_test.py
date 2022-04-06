import knn_decoding
import numpy as np


ts_name = 'tsoutput_p1_sigma2_amp5_bkg0'
ts_path = 'Data/'
knn_decoder = knn_decoding.kNNDecoder('Data/p1_sigma2/a5_p1_bkg0/tsoutput.csv',
                                      'Codebook.csv', 16, 1, verbose=True)
df = knn_decoder.run_knn_decoding()
df.to_csv('Data/p1_sigma2/a5_p1_bkg0/knn_decoding.csv')
