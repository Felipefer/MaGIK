import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt

def xi(m, a = 2.7):
    return m**(-a)

def Chi2_i(cor_iso, cor_obs, cor_err, mag_iso, mag_obs, mag_err):
    arg_cor = ((cor_obs-cor_iso)/cor_err)**2
    arg_mag = ((mag_obs-mag_iso)/mag_err)**2
    return(arg_cor+arg_mag)

def L_i(cor_iso, cor_obs, cor_err, mag_iso, mag_obs, mag_err):
    arg_err = 1/(cor_err*mag_err)
    Chi2_i_ = Chi2_i(cor_iso=cor_iso, cor_obs=cor_obs, cor_err=cor_err,
                     mag_iso=mag_iso, mag_obs=mag_obs, mag_err=mag_err)

    return(arg_err*np.exp(-Chi2_i_/2))

def Dm_i(m_iso):
    Dm = np.zeros(len(m_iso))
    Dm[1:-1] = m_iso[2:] - m_iso[:-2]
    Dm[0] = 2*(m_iso[1]-m_iso[0])
    Dm[-1] = 2*(m_iso[-1]-m_iso[-2])
    return Dm

def P_one(cor_iso, cor_obs, cor_err, mag_iso, mag_obs, mag_err, m_iso, a = 2.7):
    L_i_ = L_i(cor_iso=cor_iso, cor_obs=cor_obs, cor_err=cor_err,
               mag_iso=mag_iso, mag_obs=mag_obs, mag_err=mag_err)

    xi_i_ = xi(m_iso, a)
    Dm_i_ = Dm_i(m_iso)

    P_i_ = L_i_*xi_i_*Dm_i_
    P_i_ = P_i_.sum()
    return(P_i_)

def P_multiple(cor_iso, cor_obs, cor_err, mag_iso, mag_obs, mag_err, m_iso, a = 2.7):
    P_i = []
    for i in range(len(cor_obs)):
        P_i.append(P_one(cor_iso=cor_iso, cor_obs=cor_obs[i], cor_err=cor_err[i],
                         mag_iso=mag_iso, mag_obs=mag_obs[i], mag_err=mag_err[i],
                         m_iso=m_iso, a=a))

    P_i = np.array(P_i)
    return(P_i)

def running_mean_sdev(A, nline, ncol = None, mask = None):

    if mask is None:
        mask = np.empty(A.shape)
        mask[:] = 1

    if ncol is None:
        ncol = nline

    A_mean = np.zeros((A.shape[0], A.shape[1]))
    A_var = np.zeros((A.shape[0], A.shape[1]))

    for line in range(A.shape[0]):
        for col in range(A.shape[1]):
            line_min = line - nline if (line-nline) >= 0 else 0
            line_max = line + nline if (line+nline) <= A.shape[0] else A.shape[0]

            col_min = col - ncol if (col-ncol) >= 0 else 0
            col_max = col + ncol if (col+ncol) <= A.shape[1] else A.shape[1]

            A_ij = A[line_min:(line_max+1), col_min:(col_max+1)]
            # Aplicar mascara para desconsiderar vazios
            mask_ij = mask[line_min:(line_max+1), col_min:(col_max+1)]
            A_ij = A_ij[mask_ij == 1]

            if len(A_ij) == 0:
                A_mean[line, col] = np.nan
                A_var[line, col] = np.nan
            else:
                A_mean[line,col] = A_ij.mean()
                A_var[line,col] = A_ij.var()

    A_sdev = np.sqrt(A_var)

    return A_mean, A_sdev

def running_SN(density, background, nline, ncol = None):
    if ncol is None: ncol = nline

    SN = np.zeros((density.shape[0], density.shape[1]))

    for line in range(SN.shape[0]):
        for col in range(SN.shape[1]):
            line_min = line - nline if (line-nline) >= 0 else 0
            line_max = line + nline if (line+nline) <= SN.shape[0] else SN.shape[0]

            col_min = col - ncol if (col-ncol) >= 0 else 0
            col_max = col + ncol if (col+ncol) <= SN.shape[1] else SN.shape[1]

            background_k = background[line_min:(line_max+1), col_min:(col_max+1)]
            Nbins = background_k.shape[0] * background_k.shape[1]
            g = 1./Nbins
            
            sum_bg = background_k.sum()
            ruido = np.sqrt(background_k.var())

            SN[line,col] = (density[line,col] - g*sum_bg) / np.sqrt(density[line,col] + g**2 * sum_bg)

    return SN
#####################
# Teste             #
#####################

if 0:
    data = np.loadtxt('Isocrona_teste_select_star_v1.dat')
    B_iso = data[:,9]
    V_iso = data[:,10]

    m_iso = data[:,2]
    cor_iso = B_iso-V_iso
    mag_iso = V_iso

    #Simulate observed stars
    N_obs = 10000
    cor_obs = np.random.uniform(low = -0.5, high = 2, size = N_obs)
    mag_obs = np.random.uniform(low = -4, high = 14, size = N_obs)

    #mag_err = np.zeros(N_obs)
    #for i in range(N_obs): mag_err[i] = 0.1
    #cor_err = np.zeros(N_obs)
    #for i in range(N_obs): cor_err[i] = 0.015

    #def mag_err_func(mag):
    #    if mag <= 0:
    #        mag_err_i = 0.02
    #    elif mag > 0 and mag <= 4:
    #        mag_err_i = 0.04
    #    elif mag > 4 and mag <= 10:
    #        mag_err_i = 0.08
    #    else:
    #        mag_err_i = 0.15
    #
    #    return mag_err_i
    #
    #mag_err = np.zeros(N_obs)
    #for i in range(N_obs): mag_err[i] = mag_err_func(mag_obs[i])
    #cor_err = 0.15*mag_err

    mag_err = 0.01*(mag_obs+5)
    cor_err = 0.15*mag_err

    P = P_multiple(cor_iso=cor_iso, cor_obs=cor_obs, cor_err=cor_err,
                   mag_iso=mag_iso, mag_obs=mag_obs, mag_err=mag_err,
                   m_iso=m_iso)

    print P.min()
    #Remove P == 0
    P_min = (P[P != 0]).min()
    for i in range(len(P)):
        if P[i] == 0: P[i] = P_min

    plt.scatter(cor_obs, mag_obs, c = np.log(P), cmap = cm.jet, marker='o')
    plt.plot(cor_iso,mag_iso)
    plt.gca().invert_yaxis()
    plt.show()
