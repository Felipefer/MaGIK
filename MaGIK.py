import gera_probs2
import percentil_kernel_test2
from gera_probs2 import *
from percentil_kernel_test2 import *
import time


#fields = ['bootes_dois','hercules','kim_um','leo_t','pisces_dois','segue_dois', 'segue_tres', 'segue_um','willman_um'] 
#colors = ['r', 'y', 'g', 'c','b','m','w','r', 'k']#'w'
#l_field = [352.7,27.7,67.5,214.5,78.2,148.4,68.4,219.5,157.6]
#b_field = [67.8,35.9,-39,43.2,-47.1,-38.1,21.3,49.4,55.8]
#dists = [62000,140000,20000,,407000,183000,35000,17000,23000,38000]

#fields = ['bootes_um','canes_dois','canes_um','coma','leo_quatro','ursa_um','ursa_dois/',]
#colors = ['r', 'y', 'g', 'c','b','m']
#l_field = [209,193.3,201.6,185.7,172.2,157.6]
#b_field = [13.5,33.3,33,22.9,0,51]
#dists = [42000,150000,218000,44000,160000,100000,30000]



#window, bins, kernels, running_mean_size, log?
#def window(log):
#    if log: window = [1/6., 1/10.]
#    else: window = [1/2., 1/6.]
#    return window

window = [1/2., 1/6.]
#log = [True, False]
log = [False]
kernels = [ [[16, 14, 12, 10, 8, 6], [20,40,60,80,90]],
           [[25, 22, 18, 15, 10, 5], [20,40,60,80,90]],
           [[16, 12, 8, 4], [20,50,80]], 
           [[25, 18, 10, 5], [20,50,80]]] # [[Kernels], [percentils]]

running_mean_size = [15, 30, 60]

bins = [130, 156, 182]#, 520]

parameters = []

for i in range(len(log)):
    #for j in range(2):
    for j in range(len(window)):
        for k in range(len(kernels)):
            for w in range(len(running_mean_size)):
                for v in range(len(bins)):
                    #parameters.append([log[i], window(log[i])[j], kernels[k][0], kernels[k][1], running_mean_size[w], bins[v]])
                    parameters.append([log[i], window[j], kernels[k][0], kernels[k][1], running_mean_size[w], bins[v]])




#D_SEQ = [20000,20000,20000,20000,20000]
#COORDS1 = [90,147,90,146,90]
#COORDS2 = [-44,-44,-50,-44,-80]
#FIELD_PATH = ['field_cinco/','field_um/','field_dois/','field_tres/','field_quatro/'] 
#FIELD_BASE_NAME = ['field_cinco','field_um','field_dois','field_tres','field_quatro']


D_SEQ = [20000]
COORDS1 = [90]
COORDS2 = [-80]
FIELD_PATH = ['field_quatro/'] 
FIELD_BASE_NAME = ['field_quatro']

for i in range(len(FIELD_PATH)):
    ##################
    # Parameters  Iso#
    ############################################

    isoc_ages = [12e9]
    isoc_metal = [0.0001]#,0.005,0.015]#,0.02,0.03]
    d_seq = [D_SEQ[i]]#,60000,140000,183000,407000]
    isoc_path = 'isocs/'
    isoc_base_name = 'iso_magik'
    mag_r_lim = 22.6
    mag_i_lim = 22.6
    #coord_system = 'equatorial'
    coord_system = 'galactic'
    #####################
    # Parameters  fields#
    ############################################

    coords1 = [COORDS1[i]]
    coords2 = [COORDS2[i]]
    field_path = FIELD_PATH[i]
    field_base_name = FIELD_BASE_NAME[i]

    field_size = [1.4,1.4]
    N_bins_cut = 15

    #####################
    # Save to#
    ############################################
    save_path = FIELD_PATH[i]
    probs_path = FIELD_PATH[i]
    probs_base_name = 'probs'
    
    for j in range(len(parameters)):
        log_error = parameters[j][0]
        error_window = parameters[j][1]
        kernels = parameters[j][2]
        percentils = parameters[j][3]
        running_mean_size = parameters[j][4]
        bins = parameters[j][5]
        
        print('\n\nfield: {0}'.format(FIELD_BASE_NAME[i]))
        print('\n----- PARAMETERS ------')
        print('log_error = {0}'.format(parameters[j][0]))
        print('error_window_multiplier = {0}'.format(parameters[j][1]))
        print('kernels = {0}'.format(parameters[j][2]))
        print('percentils = {0}'.format(parameters[j][3]))
        print('running_mean_size = {0}'.format(parameters[j][4]))
        print('bins = {0}'.format(parameters[j][5]))

        startss = time.time()
        get_probs(  isoc_path=isoc_path, 
                        isoc_base_name=isoc_base_name, 
                        isoc_ages=isoc_ages, 
                        isoc_metal=isoc_metal, 
                        d_seq=d_seq, 
                        field_path=field_path, 
                        field_base_name=field_base_name, 
                        coords1=coords1, 
                        coords2=coords2, 
                        save_path=save_path,
                        coord_system=coord_system,
                        mag_r_lim = mag_r_lim,
                        mag_i_lim = mag_i_lim,
                        use_stages = True,
                        save_number = j,
                        log_error = log_error,
                        error_window = error_window)
        endss = time.time()
        print'probs',(endss - startss)
        startss = time.time()
        percentil_kernel(isoc_ages=isoc_ages, 
                        isoc_metal=isoc_metal, 
                        d_seq=d_seq, 
                        probs_path=probs_path, 
                        probs_base_name=probs_base_name, 
                        coords1=coords1, 
                        coords2=coords2, 
                        save_path=save_path,
                        field_size=field_size,
                        N_bins_cut=N_bins_cut,
                        bins=bins,
                        nline = running_mean_size,
                        kernels = kernels,
                        percentiles = percentils,
                        save_number = j)
        endss = time.time()
        print'probs',(endss - startss)

