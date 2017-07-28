# Importing libraries
from utils import *
from scipy.interpolate import interp1d
from time import time
import pyfits

# For the list of coordinates in :param coords1(RA or longitude) and :param coords1(DEC or latitude), this program loads the field data stored in the path
# :param field_path and which file names are in the form "{}_l{}_b{}.dat".format(:param field_base_name,
# :param gal_l[i], :param gal_b[i].
#
# The field files must be formatted as follows:
# col0:
# col1: object RA or longitude, expressed in degrees, ranging from 0 to 360.
# col2: object DEC or latitude, expressed in degrees, ranging from -90 to 90.
# col3, col4, col5: object magnitudes g, r, i, respectively.
# col6, col7, col8: object magnitude errors g_err, r_err, i_err, respectively.
#
# This program proceeds to calculate the probability of each star in each field to belong to the simple stellar
# populations defined by age, metallicity and distance given in :param isoc_ages, :param isoc_metal, :param d_seq.
# Where ages are given in years, metallicity in overall metal content fraction (Z), and distance in parsecs.
#
# Probabilities are calculated using a bayesian method similar to the one described by Jorgensen & Lindegren (2005).
# In our cases, observational quantities are magnitudes and colors. The probabilities are calculated both in the
# r--(g-r) space and the i--(g-i) space. For each age and metallicity provided, isochrones must be stored in the path
# :param isoc_path and the name of isochronal files must be in the format
# "{}_Z{}_t{}e9.dat".format(:param isoc_base_name, :param metal[i], :param age[i]/1e9).
#
# Isochrone files must be formatted as follows:
# col2: Stellar mass, in solar units.
# col9, col10, col11: magnitudes g, r, i.
# col-1: Evolutionary stage labeled 0=PMS, 1=MS, 2=SGB, 3=RGB, (4,5,6)=different stages of CHEB, 7=EAGB, 8=TPAGB.
#
# For each field represented by coordinated l and b, the probability of the stars to belong to each simple stellar
# populations are saved in the path :param save_path. The files containing the probability for each star to belong to a
# given stellar population of age t, # metallicity Z and distance d, are named
# "probs_l{}_b{}_age{}e9_Z{}_d{}.dat".format(l, b, t, Z, d) and structured as follows:
# col0: object RA or longitude, expressed in degrees, ranging from 0 to 360.
# col1: object DEC or latitude, expressed in degrees, ranging from -90 to 90.
# col2, col3: object magnitude r, and magnitude r error, respectively.
# col4, col5: object magnitude i, and magnitude i error, respectively.
# col6, col7: object colors g-r and g-i, respectively.
# col8: Probability that the star belong to the simple stellar population, calculated in the g--g-r space
# col9: Probability that the star belong to the simple stellar population, calculated in the g--g-i space

def get_name_iso(base_name,age,metal):
    """
    returns the name of an isochrone file for a given age and metallicity

    :param base_name: isochronal base name
    :type  base_name: string
    :param       age: isochronal age in yr
    :type        age: float
    :param     metal: isochronal metallicity
    :type      metal: float
    :return		: isochronal name
    :rtype		: string
    """
    return(base_name+'_Z{}_t{:.1f}e9.dat'.format(metal, age/1e9))


def get_name_field(base_name,l,b):
    """
    returns the name of a field file for given latitude and longitude

    :param base_name: field base name
    :type  base_name: string
    :param         l: galactic longitude in degrees
    :type	       l: float
    :param         b: galactic latitude in degrees
    :type	       b: float
    :return	        : field name
    :rtype          : string
    """
    return(base_name+'_l{0}_b{1}.fits'.format(l, b))

###################
## Parameters  Iso#
#############################################

#isoc_ages = [12e9]
#isoc_metal = [0.0001]#,0.005,0.015]#,0.02,0.03]
#d_seq = [150000]#,60000,140000,183000,407000]
#isoc_path = 'isocs/'
#isoc_base_name = 'iso_magik'

######################
## Parameters  fields#
#############################################

#gal_l = [193.8]
#gal_b = [33.8]
#field_path = 'canes_dois/'
#field_base_name = 'canes_dois'

######################
## Save to#
#############################################
#save_path = 'canes_dois/'
#save_SAB = ["SAB16_1.png","SAB16_2.png","SAB16_3.png",]
#coord_system = 'equatorial'

def get_probs(isoc_path, isoc_base_name, isoc_ages, isoc_metal, d_seq, field_path, field_base_name, coords1, coords2, save_path, mag_r_lim , mag_i_lim, coord_system,use_stages = False,save_number=0, log_error = False, error_window = 1/2.):
    """
    For each star in each field, calculate the probability of the star belonging to each simple stellar population defined 
    by each age, metallicity and distance given. Saves one file for each ssp containning the obtained probabilities for all
    stars.
    
    :param       isoc_path: path to the isochrone files
    :param  isoc_base_name: isochronal base name
    :param       isoc_ages: list of all ages to be considered. Ages in yr.
    :param      isoc_metal: list of all metallicities (Z) to be considered.
    :param           d_seq: list of all distances to be considered. Distances in parsecs.
    :param      field_path: path to the field files
    :param field_base_name: field base name
    :param           gal_l: list of galactic longitudes to be considered.
    :param           gal_b: list of galactic latitudes to be considered.
    :param       save_path: path to where result files are to be saved
    
    :type        isoc_path: string
    :type   isoc_base_name: string
    :type        isoc_ages: list of floats
    :type       isoc_metal: list of floats
    :type            d_seq: list of floats
    :type       field_path: string
    :type  field_base_name: string
    :type            gal_l: list of floats
    :type            gal_b: list of floats
    :type        save_path: string
    
    obs: stars with high magnitude error's are removed.
    """
    # for each given coord1(RA or longitude)
    SAB_counter = 0
    for coord1 in coords1:
        # for each given coord2(dec or latitude)
        for coord2 in coords2:
       
            # Load field data
            field_name = get_name_field(field_base_name, l = coord1, b = coord2)
            # TODO: dar um jeito de permitir que os argumentos de np.loadtxt sejam passados nos argumentos da funcao

            tgg = pyfits.open(field_path+field_name)
            field_data = tgg[1].data
            #print((tgg[1].columns))
            if coord_system == 'equatorial':
                obj_coord1 = field_data.field('ra')
                obj_coord2 = field_data.field('dec')
            if coord_system == 'galactic':
                obj_coord1 = field_data.field('l')
                obj_coord2 = field_data.field('b')

            obj_mag_g = field_data.field('g')
            obj_mag_r = field_data.field('r')
            obj_mag_i = field_data.field('i')
            obj_mag_g_err = field_data.field('psfMagErr_g')
            obj_mag_r_err = field_data.field('psfMagErr_r')
            obj_mag_i_err = field_data.field('psfMagErr_i')

            # Calculate colors from magnitudes
            obj_cor_gr = obj_mag_g - obj_mag_r
            obj_cor_gi = obj_mag_g - obj_mag_i
             
            # Adopting a color error
            if log_error:
                obj_mag_r_err = np.log(obj_mag_r_err) * error_window
                obj_mag_i_err = np.log(obj_mag_i_err) * error_window
            else:
                obj_mag_r_err = obj_mag_r_err * error_window
                obj_mag_i_err = obj_mag_i_err * error_window
            
            #obj_mag_r_err = obj_mag_r_err/6 # TODO: descobrir porque estamos usando isso
            #obj_mag_r_err = np.log(obj_mag_r_err)/10# TODO: descobrir porque estamos usando isso
            obj_cor_gr_err = 0.1*obj_mag_r_err
             
            #obj_mag_i_err = obj_mag_i_err/6 # TODO: descobrir porque estamos usando isso
            #obj_mag_i_err = np.log(obj_mag_i_err)/10 # TODO: descobrir porque estamos usando isso
            obj_cor_gi_err = 0.1*obj_mag_i_err
             
            # Removing stars with high magnitude error
            # TODO: fazer este passo diretamente em field_data para nao ter que escrever linha por linha
            obj_mag_filt_high_error = obj_mag_r_err < 0.2
            obj_coord1 = obj_coord1[obj_mag_filt_high_error]
            obj_coord2 = obj_coord2[obj_mag_filt_high_error]
            obj_mag_g = obj_mag_g[obj_mag_filt_high_error]
            obj_mag_r = obj_mag_r[obj_mag_filt_high_error]
            obj_mag_i = obj_mag_i[obj_mag_filt_high_error]
            obj_mag_g_err = obj_mag_g_err[obj_mag_filt_high_error]
            obj_mag_r_err = obj_mag_r_err[obj_mag_filt_high_error]
            obj_mag_i_err = obj_mag_i_err[obj_mag_filt_high_error]           
            obj_cor_gr = obj_cor_gr[obj_mag_filt_high_error]
            obj_cor_gi = obj_cor_gi[obj_mag_filt_high_error]
            obj_cor_gr_err = obj_cor_gr_err[obj_mag_filt_high_error]
            obj_cor_gi_err = obj_cor_gi_err[obj_mag_filt_high_error]
            
            
            for age in isoc_ages: # for each age
                for metal in isoc_metal: # for each metallicity
                    isoc_name = get_name_iso(isoc_base_name,age,metal)
                    isoc_data = np.loadtxt(isoc_path+isoc_name) # TODO: dar um jeito de poder passar argumentos para esta funcao
 
                    isoc_original_mass = isoc_data[:,2]
                    isoc_original_mag_g = isoc_data[:,9]
                    isoc_original_mag_r = isoc_data[:,10]
                    isoc_original_mag_i = isoc_data[:,11]
                    isoc_original_stage = isoc_data[:,-1]
                    isoc_original_cor_gr = isoc_original_mag_g - isoc_original_mag_r
                    isoc_original_cor_gi = isoc_original_mag_g - isoc_original_mag_i

                    for dist in d_seq: # for each distance
                        t0 = time()

                        if use_stages:
                            if dist <= 5000:
                                stages_to_use = [0,1]
                            elif dist > 5000 and dist <= 50000:
                                stages_to_use = [0, 1, 2]
                            elif dist > 50000:
                                stages_to_use = [0, 1, 2, 3, 4, 5, 6, 7]
                            
                            isoc_mag_r = isoc_original_mag_r[np.in1d(isoc_original_stage, stages_to_use)]
                            isoc_mag_i = isoc_original_mag_i[np.in1d(isoc_original_stage, stages_to_use)]
                            isoc_mag_g = isoc_original_mag_g[np.in1d(isoc_original_stage, stages_to_use)]
                            isoc_mass = isoc_original_mass[np.in1d(isoc_original_stage, stages_to_use)]
                            
                        else:
                            isoc_mag_r = isoc_original_mag_r
                            isoc_mag_i = isoc_original_mag_i
                            isoc_mag_g = isoc_original_mag_g
                            isoc_mass = isoc_original_mass

                        # Interpolating isochrone
                        isoc_mag_g_interp_fun = interp1d(isoc_mass, isoc_mag_g)
                        isoc_mag_r_interp_fun = interp1d(isoc_mass, isoc_mag_r)
                        isoc_mag_i_interp_fun = interp1d(isoc_mass, isoc_mag_i)

                        isoc_mass = np.linspace(isoc_mass[0], isoc_mass[-1], 1000)
                        isoc_mag_g = isoc_mag_g_interp_fun(isoc_mass)
                        isoc_mag_r = isoc_mag_r_interp_fun(isoc_mass)
                        isoc_mag_i = isoc_mag_i_interp_fun(isoc_mass)

                        # Calculating isochrone colors
                        isoc_cor_gr = isoc_mag_g - isoc_mag_r
                        isoc_cor_gi = isoc_mag_g - isoc_mag_i

                        #print "tempo de interpolacao {} sec.".format(time() - t0)

                        isoc_mag_r_dist = 5*(np.log10(dist)) - 5 + isoc_mag_r
                        isoc_mag_i_dist = 5*(np.log10(dist)) - 5 + isoc_mag_i

                        isoc_original_mag_r_dist = 5*(np.log10(dist)) - 5 + isoc_original_mag_r
                        isoc_original_mag_i_dist = 5*(np.log10(dist)) - 5 + isoc_original_mag_i
                        ###############################
                        #print(age, metal, dist)
                        # Obtain the probabilities of each star belonging to the ssp defined by age, metal and dist using the
                        # two different colors
                        
                        t0_probs = time()
                        Probs_gr_r = P_multiple(cor_iso = isoc_cor_gr,
                                                cor_obs = obj_cor_gr,
                                                cor_err = obj_cor_gr_err,
                                                mag_iso = isoc_mag_r_dist,
                                                mag_obs = obj_mag_r,
                                                mag_err = obj_mag_r_err,
                                                m_iso = isoc_mass)

                        Probs_gi_i = P_multiple(cor_iso = isoc_cor_gi,
                                                cor_obs = obj_cor_gi,
                                                cor_err = obj_cor_gi_err,
                                                mag_iso = isoc_mag_i_dist,
                                                mag_obs = obj_mag_i,
                                                mag_err = obj_mag_i_err,
                                                m_iso = isoc_mass)
                        
                        #print "tempo de calculo de probs: {} sec.".format(time()-t0)
                        # This solves a problem with the plot

                        Probs_gr_r_min = (Probs_gr_r[Probs_gr_r != 0]).min()
                        
                        for k in range(len(Probs_gr_r)):
                            if Probs_gr_r[k] == 0: Probs_gr_r[k] = Probs_gr_r_min
                            if obj_mag_r[k] >= mag_r_lim: Probs_gr_r[k] = Probs_gr_r_min
                            if Probs_gr_r[k] < 1e-300: Probs_gr_r[k] = 1e-300

                                

                        Probs_gi_i_min = (Probs_gi_i[Probs_gi_i != 0]).min()
                        for k in range(len(Probs_gi_i)):
                            if Probs_gi_i[k] == 0: Probs_gi_i[k] = Probs_gi_i_min
                            if obj_mag_i[k] >= mag_i_lim: Probs_gi_i[k] = Probs_gi_i_min
                            if Probs_gi_i[k] < 1e-300: Probs_gi_i[k] = 1e-300


                        # Normalizing the Probs
                        Probs_gr_r = Probs_gr_r/Probs_gr_r.max()
                        Probs_gi_i = Probs_gi_i/Probs_gi_i.max()

                        # Preparing data to be saved
                        save_data = np.zeros((len(Probs_gr_r), 10))
                        save_data[:,0] = np.round(obj_coord1,3)
                        save_data[:,1] = np.round(obj_coord2, 3)
                        save_data[:,2] = np.round(obj_mag_r, 2)
                        save_data[:,3] = np.round(obj_mag_r_err,2)
                        save_data[:,4] = np.round(obj_mag_i,2)
                        save_data[:,5] = np.round(obj_mag_i_err,2)
                        save_data[:,6] = np.round(obj_cor_gr,2)
                        save_data[:,7] = np.round(obj_cor_gi,2)
                        save_data[:,8] = Probs_gr_r
                        save_data[:,9] = Probs_gi_i

                        # Saving data
                        save_filename = save_path+'probs_l{0}_b{1}_age{2}e9_Z{3}_d{4}kpc_{5}.dat'.format(coord1,coord2,age/1e9,metal,dist/1000,save_number)
                        np.savetxt(save_filename, save_data, delimiter = ',')

                        print "tempo total: {} sec.".format(time()-t0)

                        # Plotting data
                        plt.rcParams['axes.linewidth'] = 5

                        fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(20,32))
                        filter1 = Probs_gr_r != Probs_gr_r.min()
                        filter2 = isoc_original_stage < 4

                        ax.scatter(obj_cor_gr[~filter1], obj_mag_r[~filter1], c = '#DDDDDD', cmap = cm.rainbow, edgecolor = 'none', s = 600, marker='o', alpha = 0.7)
                        ax.scatter(obj_cor_gr[filter1], obj_mag_r[filter1], c = np.log(Probs_gr_r[filter1]), cmap = cm.rainbow, edgecolor = 'none', s = 600, marker='o', alpha = 0.7)
                        ax.plot(isoc_original_cor_gr[filter2], isoc_original_mag_r_dist[filter2], 'w-', linewidth = 7)
                        ax.plot(isoc_original_cor_gr[~filter2][1:], isoc_original_mag_r_dist[~filter2][1:], 'w-', linewidth = 7)
                        #plt.plot(isoc_cor_gr, isoc_mag_r_dist, 'b-')
                        ax.set_xlim((-0.2, 1.5))
                        ax.set_ylim((23.5, 15))
                        color = '#FF2A00'
                        for label in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                            label.set_fontsize(70)
                            label.set_color("#FAFAFA")

                        ax.spines['bottom'].set_color(color)
                        ax.spines['top'].set_color(color)
                        ax.spines['right'].set_color(color)
                        ax.spines['left'].set_color(color)

                        ax.xaxis.label.set_fontsize(70)
                        ax.set_xlabel('g-r')
                        ax.xaxis.label.set_color("#FAFAFA")
                        ax.yaxis.label.set_color("#FAFAFA")
                        ax.yaxis.label.set_fontsize(70)
                        ax.set_ylabel('r')

                        ax.minorticks_on()
                        ax.tick_params(length = 25, width = 5, which='major', color = color)
                        ax.tick_params(length = 12.5, width = 3, which='minor', color = color)

                        save_SAB = ["SAB16_1.png","SAB16_2.png","SAB16_3.png",]
                        plt.savefig(save_SAB[SAB_counter], transparent = True)
                        SAB_counter += 1
                        plt.clf()
                        plt.close()

