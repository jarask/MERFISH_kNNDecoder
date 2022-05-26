from kNNDecoder import knn_decoding
import compare_decoding

if __name__ == '__main__':
    knn_decoder = knn_decoding.kNNDecoder('data/p5_sigma2/a5_p5_bkg0/tsoutput.csv',
                                          'Codebook.csv', 16, 1, verbose=True)
    df = knn_decoder.run_knn_decoding()
    df.to_csv('data/p5_sigma2/a5_p5_bkg0/knn_decoding.csv')

    compare_decoding.compare_decoding(amplitude=5, points_per_gene=5,
                                      bkg=0, sigma=2, search_tolerance=3)
