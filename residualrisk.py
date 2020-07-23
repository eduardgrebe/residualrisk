import math
import multiprocessing as mp
import numpy as np
import scipy.stats as stats
import scipy.integrate as si

def concentration(C0, doubling_time, t):
    concentration = C0 * 2 ** (t / doubling_time)
    return concentration

def prob_infectious_copies(n_copies, k):
    prob = 1.000000000001 - math.exp(-k * n_copies)
    return prob

def prob_infectious(t, 
                    C0, 
                    doubling_time, 
                    volume_transfused, 
                    k, 
                    copies_per_virion = 2
                   ):
    C = concentration(C0, doubling_time, t)
    n_copies = C * copies_per_virion * volume_transfused
    prob = prob_infectious_copies(n_copies, k)
    return prob
    
def prob_infectious_copies_wc(n_copies):
    if n_copies < 2:
        return 0.0
    elif n_copies >= 2:
        return 1.0
    
def prob_infectious_wc(t, 
                       C0, 
                       doubling_time, 
                       volume_transfused, 
                       copies_per_virion = 2
                      ):
    C = concentration(C0, doubling_time, t)
    n_copies = C * copies_per_virion * volume_transfused
    prob = prob_infectious_copies_wc(n_copies)
    return prob

def prob_pos_init(C,
                  doubling_time,
                  pool_size,
                  lod50,
                  lod95_lod50_ratio,
                  z
                 ):
    if (not isinstance(pool_size, int)) or pool_size < 1:
        raise Exception("pool_size must be an integer of at least 1")
    # C is in copies copies_per_virion * C when C in virions
    X = z * (math.log10(((C) / (pool_size * lod50))) / 
             math.log10(lod95_lod50_ratio))
    prob = stats.norm.cdf(X)
    return prob

def prob_neg_retest(C,
                    doubling_time,
                    pool_size,
                    lod50,
                    lod95_lod50_ratio,
                    retests,
                    z
                   ):
    if (not isinstance(pool_size, int)) or pool_size < 1:
        raise Exception("pool_size must be an integer of at least 1")
    if (not isinstance(retests, int)) or retests < 0:
        raise Exception("retests must be a positive integer")
    elif retests == 0:
        return 0
    elif retests >= 1:
        # C is in copies copies_per_virion * C when C in virions
        X = z * (math.log10(((C) / lod50)) / math.log10(lod95_lod50_ratio))
        #print(X)
        prob = (1 - stats.norm.cdf(X)) ** retests
        return prob

def prob_nondetection(t,
                      copies_per_virion,
                      C0,
                      doubling_time,
                      pool_size,
                      lod50,
                      lod95_lod50_ratio,
                      retests,
                      z = 1.6449
                     ):
    Cv = concentration(C0, doubling_time, t)
    Cc = copies_per_virion * Cv
    p_pos_init = prob_pos_init(Cc, doubling_time, pool_size, lod50, 
                               lod95_lod50_ratio, z)
    p_neg_retest = prob_neg_retest(Cc, doubling_time, pool_size, lod50, 
                                   lod95_lod50_ratio, retests, z)
    prob = 1 - p_pos_init * (1 - p_neg_retest)
    return prob

def prob_infectious_nondetection(t,
                                 copies_per_virion,
                                 C0,
                                 doubling_time,
                                 volume_transfused, 
                                 k,
                                 pool_size,
                                 lod50,
                                 lod95_lod50_ratio,
                                 retests,
                                 z = 1.6449
                                ):
    product = (prob_infectious(t, C0, doubling_time, volume_transfused, k) * 
            prob_nondetection(t, copies_per_virion, C0, doubling_time, 
                              pool_size,lod50,lod95_lod50_ratio,retests,z))
    return(product)

def risk_days(copies_per_virion,
              C0,
              doubling_time,
              volume_transfused, 
              k,
              pool_size,
              lod50,
              lod95_lod50_ratio,
              retests,
              z = 1.6449,
              limits = (-100,500)
             ):

    # Ideally we would integrate from -np.inf to np.inf, but that causes an 
    # overflow error, so we choose safe limits instead    
    rd = si.quad(prob_infectious_nondetection, limits[0], limits[1], 
              args=(copies_per_virion, C0, doubling_time, volume_transfused, k, 
                    pool_size, lod50, lod95_lod50_ratio, retests, z))[0]
    return rd

def prob_infectious_nondetection_wc(t,
                                 copies_per_virion,
                                 C0,
                                 doubling_time,
                                 volume_transfused, 
                                 pool_size,
                                 lod50,
                                 lod95_lod50_ratio,
                                 retests,
                                 z = 1.6449
                                ):
    product = (prob_infectious_wc(t, C0, doubling_time, volume_transfused) * 
            prob_nondetection(t, copies_per_virion, C0, doubling_time, 
                              pool_size,lod50,lod95_lod50_ratio,retests,z))
    return(product)

def risk_days_wc(copies_per_virion,
              C0,
              doubling_time,
              volume_transfused, 
              pool_size,
              lod50,
              lod95_lod50_ratio,
              retests,
              z = 1.6449,
              limits = (-20,100)
             ):
    # Ideally we would integrate from -np.inf to np.inf, but that causes an 
    # overflow error, so we choose safe limits instead
    rd = si.quad(prob_infectious_nondetection_wc, limits[0], limits[1], 
              args=(copies_per_virion, C0, doubling_time, volume_transfused, 
                    pool_size, lod50, lod95_lod50_ratio, retests, z))[0]
    return rd

def iwp_bs(k, 
           k_gamma_shape, 
           k_gamma_scale, 
           doubling_time, 
           doubling_time_norm_sd, 
           lod50, 
           lod50_sd, 
           lod95_lod50_ratio, 
           volume_transfused, 
           volume_transfused_range, 
           pool_size, 
           retests, 
           C0 = 0.00025, 
           copies_per_virion = 2, 
           alpha = 0.05, 
           n_bs = 10000, 
           seed = 126887
          ):
    iwp_pe = risk_days(copies_per_virion, C0, doubling_time, volume_transfused, 
                              k, pool_size, lod50, lod95_lod50_ratio, retests)
    np.random.seed(seed)
    doubling_time_draws = stats.truncnorm.rvs(0, 
                                              np.inf, 
                                              doubling_time, 
                                              doubling_time_norm_sd, 
                                              n_bs)
    k_draws = np.random.gamma(k_gamma_shape, k_gamma_scale, n_bs)
    lod50_draws = stats.truncnorm.rvs(0, np.inf, lod50, lod50_sd, n_bs)
    volume_transfused_draws = np.random.uniform(volume_transfused_range[0],
                                                volume_transfused_range[1],
                                                n_bs)
    iwp = []
    for i in range(n_bs):
        iwp.append(risk_days(copies_per_virion, C0, doubling_time_draws[i], 
                             volume_transfused_draws[i], k_draws[i], pool_size, 
                             lod50_draws[i], lod95_lod50_ratio, retests))
    iwp_ci = np.quantile(iwp, (alpha/2, 1 - alpha/2))
    return (iwp_pe, iwp_ci, iwp)

def residual_risk_iwp(iwp_pe, 
                      iwp_bs, 
                      incidence, 
                      incidence_norm_sd, 
                      per = 1e6,
                      alpha = 0.05,
                      seed = 126887):
    rr_pe = incidence * iwp_pe / 365.25 * per
    n_bs = len(iwp_bs)
    np.random.seed(seed)
    incidence_draws = stats.truncnorm.rvs(0, np.inf, incidence, 
                                          incidence_norm_sd, n_bs)
    rr = []
    for i in range(n_bs):
        rr.append(incidence_draws[i] * iwp_bs[i] / 365.25 * per)
    rr_ci = np.quantile(rr, (0.025,0.975))
    rr_se = np.std(rr)
    return (rr_pe, rr_ci, rr_se)
    
def residual_risk(k, 
                  k_gamma_shape, 
                  k_gamma_scale, 
                  doubling_time, 
                  doubling_time_norm_sd, 
                  lod50, 
                  lod50_sd, 
                  lod95_lod50_ratio, 
                  volume_transfused, 
                  volume_transfused_range, 
                  pool_size, 
                  retests, 
                  incidence, 
                  incidence_norm_sd, 
                  C0 = 0.00025, 
                  copies_per_virion = 2, 
                  per = 1e6, 
                  n_bs = 10000, 
                  seed = 126887
                 ):
    iwp_pe = risk_days(copies_per_virion, C0, doubling_time, volume_transfused, 
                              k, pool_size, lod50, lod95_lod50_ratio, retests)
    rr_pe = incidence * (iwp_pe / 365.25) * per
    if n_bs > 0:
        np.random.seed(seed)
        doubling_time_draws = stats.truncnorm.rvs(0, np.inf, doubling_time, 
                                                  doubling_time_norm_sd, n_bs)
        volume_transfused_draws = np.random.uniform(volume_transfused_range[0],
                                                    volume_transfused_range[1],
                                                    n_bs)
        k_draws = np.random.gamma(k_gamma_shape, k_gamma_scale, n_bs)
        lod50_draws = stats.truncnorm.rvs(0, np.inf, lod50, lod50_sd, n_bs)
        incidence_draws = stats.truncnorm.rvs(0, np.inf, incidence, 
                                              incidence_norm_sd, n_bs)
        iwp = []
        rr = []
        for i in range(n_bs):
            #prog = i/n_bs*100
            #print('Computing RR [%d%%]\r'%prog, end="")
            iwp.append(risk_days(copies_per_virion, C0, doubling_time_draws[i], 
                                 volume_transfused_draws[i], k_draws[i], 
                                 pool_size, lod50_draws[i], lod95_lod50_ratio, 
                                 retests))
            rr.append(incidence_draws[i] * iwp[i] / 365.25 * per)
        iwp_ci = np.quantile(iwp, (0.025,0.975))
        rr_ci = np.quantile(rr, (0.025,0.975))
        rr_se = np.std(rr)
        return (iwp_pe, iwp_ci, rr_pe, rr_ci, rr_se)
    else:
        return (iwp_pe, rr_pe)

def iwp_bs_par(k, 
               k_gamma_shape, 
               k_gamma_scale, 
               doubling_time, 
               doubling_time_norm_sd, 
               lod50, 
               lod50_sd, 
               lod95_lod50_ratio, 
               volume_transfused, 
               volume_transfused_range, 
               pool_size, 
               retests, 
               C0 = 0.00025, 
               copies_per_virion = 2, 
               alpha = 0.05,
               n_bs = 10000, 
               seed = 126887
               ):

    iwp_pe = risk_days(copies_per_virion, 
                       C0, 
                       doubling_time, 
                       volume_transfused, 
                       k, 
                       pool_size, 
                       lod50, 
                       lod95_lod50_ratio, 
                       retests)
    np.random.seed(seed)
    doubling_time_draws = stats.truncnorm.rvs(0, 
                                              np.inf, 
                                              doubling_time, 
                                              doubling_time_norm_sd, 
                                              n_bs)
    k_draws = np.random.gamma(k_gamma_shape, k_gamma_scale, n_bs)
    lod50_draws = stats.truncnorm.rvs(0, np.inf, lod50, lod50_sd, n_bs)
    volume_transfused_draws = np.random.uniform(volume_transfused_range[0],
                                                volume_transfused_range[1],
                                                n_bs)
    p = mp.Pool()
    result = [p.apply_async(risk_days,
                                args=(copies_per_virion,
                                  C0,
                                  doubling_time_draws[i],
                                  volume_transfused_draws[i],
                                  k_draws[i],
                                  pool_size,
                                  lod50_draws[i],
                                  lod95_lod50_ratio,
                                  retests)
                                  ) for i in range(n_bs)]
    iwp = [r.get() for r in result]
    p.close()
    p.join()
    iwp_ci = np.quantile(iwp, (alpha/2, 1 - alpha/2))
    return (iwp_pe, iwp_ci, iwp)

def count_cores():
    return mp.cpu_count()

if __name__ == '__main__':
    pass
