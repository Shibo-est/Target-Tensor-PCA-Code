import numpy as np
import pandas as pd
from functools import reduce
from scipy.stats import ortho_group
import copy
from itertools import product
import datetime
from ast import literal_eval
from scipy import linalg
from scipy.stats import ortho_group


np.set_printoptions(precision=8, suppress=True)

def TT_PCA(Y, X, gamma, r): # 20240717
    N_Y = Y.shape[-1]
    (T, K, N_X) = X.shape
    
    Lambda_X_est = np.zeros((N_X, r))
    W_est = np.zeros((K, r))

    Z = np.column_stack((X.reshape((T, -1), order='F'), np.sqrt(gamma) * Y))
    Sigma_Z = np.cov(Z.T)
    eign_value, Lambda_gamma = np.linalg.eig(Sigma_Z)
    Lambda_gamma = np.real_if_close(Lambda_gamma, 1)
    eign_value = np.real_if_close(eign_value)
    Lambda_gamma_est = Lambda_gamma[:, np.flip(np.argsort(eign_value))[:r]]

    F = Z @ Lambda_gamma_est @ np.linalg.inv(Lambda_gamma_est.T @ Lambda_gamma_est)
    
    for i in range(r):
        H_iM = Lambda_gamma_est[N_Y:, i].reshape((K, N_X), order='F')
        Lambda_X_est
    

    update_magnitude = np.zeros(2)
    Lambda_X_est = Lambda_X_ini
    W_est = W_ini
    weight = np.zeros(r)
    for i in range(r):
        H_iM = Lambda_gamma_est[N_Y:, i].reshape((K, N_X), order='F')
        U, S, Vh = np.linalg.svd(H_iM)

        max_index = np.flip(np.argsort(eign_value))[0]
        Lambda_X_est[:, i] = Vh[max_index,:]
        W_est[:, i] = U[:, max_index]*S[max_index]

            
    Ux, Sx, Vhx  = np.linalg.svd(Lambda_X_est)
    Lambda_X_est = Ux[:, np.flip(np.argsort(Sx))[:r]]
    W_est = W_est@np.diag(1/np.flip(np.sort(Sx))[:r])@Vhx[np.flip(np.argsort(Sx))[:r],:]
    
    return Lambda_gamma_est, Lambda_X_est, W_est, F


def TT_PCA_ALS(Y, X, gamma, r, Lambda_X_ini, W_ini, max_iter=1000, stop_thresh=1e-1):
    N_Y = Y.shape[-1]
    (T, K, N_X) = X.shape

    Z = np.column_stack((X.reshape((T, -1), order='F'), np.sqrt(gamma) * Y))
    Sigma_Z = np.cov(Z.T)
    eign_value, Lambda_gamma = np.linalg.eig(Sigma_Z)
    Lambda_gamma = np.real_if_close(Lambda_gamma, 1)
    eign_value = np.real_if_close(eign_value)
    Lambda_gamma_est = Lambda_gamma[:, np.flip(np.argsort(eign_value))[:r]]

    F = Z @ Lambda_gamma_est @ np.linalg.inv(Lambda_gamma_est.T @ Lambda_gamma_est)

    update_magnitude = np.zeros(2)
    Lambda_X_est = Lambda_X_ini
    W_est = W_ini
    weight = np.zeros(r)
    for iter_num in range(max_iter):
        Lambda_X_pre = Lambda_X_est.copy()
        W_pre = W_est.copy()
        for i in range(r):
            H_iM = Lambda_gamma_est[N_Y:, i].reshape((K, N_X), order='F')

            Lambda_X_est[:, i] = W_est[:, i].T @ H_iM
            Lambda_X_est[:, i] = Lambda_X_est[:, i] / np.linalg.norm(Lambda_X_est[:, i])

            W_est[:, i] = H_iM @ Lambda_X_est[:, i]
            W_est[:, i] = W_est[:, i] / np.linalg.norm(W_est[:, i])

            weight[i] = W_est[:, i].T @ H_iM @ Lambda_X_est[:, i]

            update_magnitude[0] += (np.linalg.norm(Lambda_X_est[:, i] - Lambda_X_pre[:, i])) ** 2
            update_magnitude[1] += (np.linalg.norm(W_est[:, i] - W_pre[:, i])) ** 2

        if (update_magnitude[0] < stop_thresh) and (update_magnitude[1] < stop_thresh):
            W_est = W_est @ weight
            break

    return Lambda_gamma_est, Lambda_X_est, W_est, F


def DGP(T, N_Y, N_X, K, r, gamma, burn=200):
    Lamma_Y = ortho_group.rvs(dim=N_Y)[:, :r] / np.sqrt(gamma) * np.sqrt(N_Y + K * N_X)
    Lamma_X = ortho_group.rvs(dim=N_X)[:, :r] * np.sqrt(N_Y + K * N_X)
    W = np.random.rand(K, r) - 0.5
    W = W / np.linalg.norm(W, axis=0)

    phi_G = np.random.rand(r) - 0.5
    G = np.zeros((T + burn, r))
    for t in range(T + burn - 1):
        noise = np.random.normal(loc=0.0, scale=1, size=r)
        G[t + 1, :] = np.einsum('..., ...->...', phi_G, G[t, :]) + noise
    G = G[burn:, :]

    # Y = G @ Lamma_Y.T + np.random.normal(loc=0.0, scale=1, size=(T, N_Y))
    Y = np.zeros((T, K, N_Y))
    for i in range(r):
        Y += np.einsum('ij,k->ijk', np.outer(G[:, i], W[:, i]), Lamma_Y[:, i])

    Y += np.random.normal(loc=0.0, scale=1, size=(T, K, N_Y))

    X = np.zeros((T, K, N_X))
    for i in range(r):
        X += np.einsum('ij,k->ijk', np.outer(G[:, i], W[:, i]), Lamma_X[:, i])

    X += np.random.normal(loc=0.0, scale=1, size=(T, K, N_X))

    return Y, X, G, Lamma_Y, W, Lamma_X