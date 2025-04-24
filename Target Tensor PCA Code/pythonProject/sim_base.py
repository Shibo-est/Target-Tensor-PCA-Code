import numpy as np
import pandas as pd
from functools import reduce
from scipy.stats import ortho_group
import copy
from itertools import product
import datetime
from TT_PCA_algo import *

np.set_printoptions(precision=8, suppress=True)

n_iter_total = 500
iter_start = 0
file_path = 'sim/results00_20240613.csv'
T_range = np.arange(40, 160, 20)
P_range = np.arange(5, 30, 5)
R_range = np.arange(2, 5)
gamma_range = np.arange(0.2, 1.1, 0.4)
starttime = datetime.datetime.now()
mse_results_all = pd.DataFrame(columns=['iter', 'P', 'R', 'T',
                                        'missing_mse', 'run_time'])
mse_results_all.to_csv(file_path, index=False)
for iter_no in range(iter_start, iter_start + n_iter_total):

    for (T, P, r, gamma) in product(T_range, P_range, R_range, gamma_range):
        N_Y = N_X = K = P
        for T in T_range:
            np.random.seed(iter_no + T)
            Y_full, X, G_true, Gamma_Y_true, W_true, Gamma_X_true = DGP(T, N_Y, N_X, K, r, gamma)
            Y = Y_full[:, 0, :]
            starttime = datetime.datetime.now()
            print('(iter, T, P, r, gamma):', (iter_no, T, P, r, gamma))
            Lambda_gamma_est, Lambda_X_est, W_est, G_est = TT_PCA(Y, X, gamma, r,
                                                                  Gamma_X_true, W_true)
            Y_est = np.zeros((T, K, N_Y))
            for i in range(r):
                Y_est += np.einsum('ij,k->ijk',
                                   np.outer(G_est[:, i], W_est[:, i]), Lambda_gamma_est[:N_Y, i])
            endtime = datetime.datetime.now()
            missing_mse = np.linalg.norm(Y_est[:, 1:, :] - Y_full[:, 1:, :])

            mse_results_all = pd.DataFrame({'iter': iter_no,
                                            'P': P,
                                            'R': r,
                                            'T': T,
                                            'missing_mse': missing_mse,
                                            'run_time': (endtime - starttime).seconds},
                                           index=[0])

            mse_results_all.to_csv(file_path, mode='a', index=False, header=False)
