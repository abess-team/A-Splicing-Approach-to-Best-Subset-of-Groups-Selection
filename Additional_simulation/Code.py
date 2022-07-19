import math
import random
import numpy as np
import scipy as sc
import math
from bnb import BNBTree
from sklearn.metrics import matthews_corrcoef
import abess
from time import time
from cd import cd_swaps

## Pre-processing and metrics
def gen_coef(k, i):
    random.seed(i)
    coef = np.random.randn(k + 1)
    coef = list(coef - np.mean(coef))
    del coef[0]
    return(coef)

def coeff(Tn, k, p, ind, i):
    beta = np.zeros(p)
    coef = list(map(gen_coef, [k]*Tn, [i]*Tn))
    index_tmp = list(
        map(lambda x: list(range(x * k , (x+1) * k )), ind))
    for i in range(Tn):
        beta[index_tmp[i]] = coef[i]
    return(beta)

def metrics(true, gsplicing, coef, intercept, beta):
    coef_full = np.append(coef, intercept)
    beta_full = np.append(beta, 0)
    err = np.linalg.norm(coef_full - beta_full)
    tp = np.sum(gsplicing & true)
    fp = np.sum(gsplicing & 1-true)
    fn = np.sum(1-gsplicing & true)
    tn = np.sum(1-gsplicing & 1-true)
    mcc = matthews_corrcoef(true, gsplicing)
    return [tp + fp, tp / (tp + fn), fp / (tn + fp), mcc, err]

def gen_synthetic(n, p, Tn, k, sigma_y, i):
    random.seed(i)
    ind = sorted(np.random.randint(1, p / k, Tn))
    beta = coeff(Tn, k, p, ind, i)
    R = np.array(np.random.randn(n * p).reshape(n, p))
    dim = int(p / k)
    sigma = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            if (i == j):
                sigma[i][j] = 1
            else:
                sigma[i][j] = 0.6
    Z = np.random.multivariate_normal(mean=[0] * dim, cov=sigma, size=n)
    z = list(map(lambda x: (
        Z[:, (math.floor((x - 1) / k))] + R[:, x]) / math.sqrt(2), list(range(p))))
    tmp = [i for j in z for i in j]
    X = np.array(tmp).reshape(n, p)
    pe = np.dot(X, beta)
    y = pe + [random.gauss(0, sigma_y) for _ in range(n)]
    return X, y, beta, ind

def cd(X, y, group_indices, lam_max, nlam):
    Params = {
        'NumIters': 1000, # maximum number of CD full cycles.
        'Tolerance': 1e-5 , # convergence tolerance.
        'Logging': False, # prints output while running
        'FeaturesOrder': 'Cyclic', # order in which CD updates the groups.
    }
    cd_swaps_object = cd_swaps(X, y, group_indices, Params)
    lam = [np.exp(x) for x in np.linspace(np.log(lam_max), np.log(0.005*lam_max), nlam)]
    n, p = np.shape(X)
    k = len(group_indices[1])
    coef = np.zeros((p, nlam))
    intercept = np.zeros(nlam)
    lost = np.zeros(nlam)
    non_group = np.zeros(nlam)
    gic = np.zeros(nlam)
    for i in range(nlam):
        if i == 0:
            coef[:, i] = cd_swaps_object.fit(lambda_0 = lam[i], lambda_2 = 0, max_swaps = 1, warm_start = np.zeros(p))[0]
        else:
            coef[:, i] = cd_swaps_object.fit(lambda_0 = lam[i], lambda_2 = 0, max_swaps = 1, warm_start = coef[:, i-1])[0]
        intercept[i] = np.mean(y) - np.dot(np.mean(X, axis = 0), coef[:, i])
        lost[i] = math.pow(np.linalg.norm(y - np.dot(X, coef[:, i])-intercept[i]), 2)/2/n
        non_group[i] = int(len(np.flatnonzero(coef[:, i]))/k)
        gic[i] = n*np.log(lost[i])+0.5*np.log(np.log(n))*np.log(p/k)*non_group[i]*k
    ind = np.argmin(gic)
    temp = np.flatnonzero(coef[:, ind])
    gr_ind = [int(temp[i*k]/k) for i in range(int(non_group[ind]))]
    return intercept[ind], coef[:, ind], gr_ind

def cd_bnb(X, y, group_indices, lam_max, nlam):
    Params = {
        'NumIters': 1000, # maximum number of CD full cycles.
        'Tolerance': 1e-5 , # convergence tolerance.
        'Logging': False, # prints output while running
        'FeaturesOrder': 'Cyclic', # order in which CD updates the groups.
    }
    cd_swaps_object = cd_swaps(X, y, group_indices, Params)
    lam = [np.exp(x) for x in np.linspace(np.log(lam_max), np.log(0.005*lam_max), nlam)]
    n, p = np.shape(X)
    k = len(group_indices[1])
    coef = np.zeros((p, nlam))
    intercept = np.zeros(nlam)
    lost = np.zeros(nlam)
    non_group = np.zeros(nlam)
    gic = np.zeros(nlam)
    for i in range(nlam):
        if i == 0:
            coef[:, i] = cd_swaps_object.fit(lambda_0 = lam[i], lambda_2 = 0, max_swaps = 1, warm_start = np.zeros(p))[0]
        else:
            coef[:, i] = cd_swaps_object.fit(lambda_0 = lam[i], lambda_2 = 0, max_swaps = 1, warm_start = coef[:, i-1])[0]
        intercept[i] = np.mean(y) - np.dot(np.mean(X, axis = 0), coef[:, i])
        lost[i] = math.pow(np.linalg.norm(y - np.dot(X, coef[:, i])-intercept[i]), 2)/2/n
        non_group[i] = int(len(np.flatnonzero(coef[:, i]))/k)
        gic[i] = n*np.log(lost[i])+0.5*np.log(np.log(n))*np.log(p/k)*non_group[i]*k
    ind = np.argmin(gic)
    temp = np.flatnonzero(coef[:, ind])
    gr_ind = [int(temp[i*k]/k) for i in range(int(non_group[ind]))]
    beta_seq = [np.linalg.norm(coef[range((i*k), (i*k+k)), ind]) for i in gr_ind]
    tree = BNBTree(X, y, group_indices)
    solver_output = tree.solve(lambda_0 = lam[ind], lambda_2 = 0, m = np.max(beta_seq), warm_start = coef[:, ind], time_limit = 500)
    coef_bnb = solver_output[1]
    intercept_bnb = np.mean(y) - np.dot(np.mean(X, axis = 0), coef_bnb)
    non_group_bnb = int(len(np.flatnonzero(coef_bnb))/k)
    temp = np.flatnonzero(coef_bnb)
    gr_ind_bnb = [int(temp[i*k]/k) for i in range(int(non_group_bnb))]
    return intercept_bnb, coef_bnb, gr_ind_bnb

n = 500 # sample size
p = 500  # dimensionality
k = 5  # group size
Tn = 10  # True model size
sigma_y = 4
lam_max = 5000 # the largest value of lambda for CD+LS 
nlam = 100 # the number of lambda
group_indices = [[k*i + j for j in range(k)] for i in range(int(p/k))]


## Simulation 
result_list = []
for i in range(1, 51):
    try:
        print(i) 
        X, y, beta, ind = gen_synthetic(n, p, Tn, k, sigma_y, i)
        group = np.repeat(list(range(int(p / k))), k)
        true = np.zeros_like(np.unique(group))
        true[ind] = 1
    
        # Sequence search (SGSplicing)
        t = time()
        model = abess.LinearRegression(ic_type="gic", ic_coef = 0.5)
        model.fit(X, y, group=group)
        t = time() - t
        nonzero = np.nonzero(model.coef_)[0]
        nonzero_group = np.unique(group[nonzero])
        gsplicing = np.zeros_like(true)
        gsplicing[nonzero_group] = 1
        res_seq = metrics(
            true,
            gsplicing,
            model.coef_,
            model.intercept_,
            beta) + [t]
    
        # Golden search (GGSplicing)
        t = time()
        model = abess.LinearRegression(ic_type="gic", path_type="gs", ic_coef = 0.5)
        model.fit(X, y, group=group)
        t = time() - t
        nonzero = np.nonzero(model.coef_)[0]
        nonzero_group = np.unique(group[nonzero])
        gsplicing = np.zeros_like(true)
        gsplicing[nonzero_group] = 1
        res_gs = metrics(
            true,
            gsplicing,
            model.coef_,
            model.intercept_,
            beta) + [t]
   
        # Approximate approach (CD+LS)
        t = time()
        result_cd = cd(X, y, group_indices, lam_max, nlam) 
        t = time() - t
        cd_seq = np.zeros_like(true)
        cd_seq[result_cd[2]] = 1
        res_cd = metrics(
            true,
            cd_seq,
            result_cd[1],
            result_cd[0],
            beta) + [t]
    
        # Exact MIO-based approach (BNB)
        t = time()
        result_bnb = cd_bnb(X, y, group_indices, lam_max, nlam) 
        t = time() - t
        bnb_seq = np.zeros_like(true)
        bnb_seq[result_bnb[2]] = 1
        res_bnb = metrics(
            true,
            bnb_seq,
            result_bnb[1],
            result_bnb[0],
            beta) + [t]
        
        result_list.append([res_seq, res_gs, res_cd, res_bnb])
    
    except ValueError:
        pass
    continue
