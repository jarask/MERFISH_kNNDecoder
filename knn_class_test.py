import knn_decoding
import numpy as np

if __name__ == '__main__':
    knn_decoder = knn_decoding.kNNDecoder('data/p1_sigma2/a5_p1_bkg100/tsoutput.csv',
                                          'Codebook.csv', 16, 1, verbose=True)
    df = knn_decoder.run_knn_decoding()
    # df.to_csv('data/p1_sigma2/a5_p1_bkg100/knn_decoding.csv')
