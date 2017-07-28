import scipy.ndimage as ndimage
import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from utils import *
import scipy
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import matplotlib

# Use Green's theorem to compute the area
# enclosed by the given contour.
def area(vs):
    """
    Calculates the area enclosed by a given contour using Green's theorem.
    :param vs: list of vertices in the form [[x0,y0], [x1, y1], [x2, y2], ...].
    :type  vs: list.
    :return  : Area.
    :rtype   : float.
    """
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(x0*dy - y0*dx)
        x0 = x1
        y0 = y1
    return a

def get_maxima(x,y,h,neighborhood_size = 10, return_index = False):
    # algoritmo retirado de http://stackoverflow.com/questions/9111711/get-coordinates-of-local-maxima-in-2d-array-above-certain-value
    """
    Find the local maxima in a given array h(x, y).

    :param                 x: 1D array that defines the x values of the columns in the :param h array
    :type                  x: array_like
    :param                 y: 1D array that defines the y values of the rows in the :param h
    :type                  y: array_like
    :param                 h: 2D array containing the z value for each (x,y)
    :type                  h: array_like
    :param neighborhood_size: size parameter used in filters.maximum_filter
    :type  neighborhood_size: scalar or tuple
    :return                 : list of local maxima in the form [[x0,y0], [x1,y1], ...]
    :rtype                  : array

    if return index is set True, also returns two more 1D arrays, the first containing the column index, and the second
    the row index of the of the local maxima in h.
    """

    data_max = filters.maximum_filter(h, neighborhood_size)
    maxima = (h == data_max)
    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    
    X, Y = [], []
    for dy,dx in slices:
        X_center = (dx.start + dx.stop - 1)/2
        X.append(X_center)
        Y_center = (dy.start + dy.stop - 1)/2    
        Y.append(Y_center)
    
    maxima = []
    for i in range(len(X)):
        maxima.append([x[X[i]], y[Y[i]]])
    if return_index: return maxima, X, Y
    else: return maxima

def get_name_probs(base_name, coord1, coord2, age, metal, dist, save_number):
    """
    Returns th name of the probs file for coordinates gal_l, gal_b, and ssp defined by age, metal and dist
    :param base_name: probs base name file
    :param     gal_l: galactic longitude of the field in degrees
    :param     gal_b: galactic latitude of the field in degrees
    :param       age: age of the ssp in years
    :param     metal: metallicity (Z) of the ssp
    :param      dist: distance of the ssp in parsecs
    :return         : probs filename

    :type  base_name: string
    :type      gal_l: float
    :type      gal_b: float
    :type        age: float
    :type      metal: float
    :type       dist: float
    :rtype          : string
    """
    return (base_name+'_l{}_b{}_age{:.1f}e9_Z{}_d{:.0f}kpc_{}.dat'.format(coord1, coord2, age/1e9, metal, dist/1e3, save_number))

def caracterize_density_field(field,
                              ximage,
                              yimage,
                              bins,
                              N_fracs_H = 1000):

    Delta = field.max() - field.min() # Diferenca entre maior pico e vale mais fundo
    N_fracs_H = 1000

    Vertices_contorno = []

    # Obter maximos
    maxima, xindex, yindex = get_maxima(ximage, yimage, field, return_index = True)
    N_picos = len(maxima)
    V_picos = []

    for i in range(N_picos): 
        V_picos.append(field[yindex[i], xindex[i]])

    V_profundidade_pico = np.zeros(N_picos)
    Area_pico = np.zeros(N_picos)

    x_maxima = []
    y_maxima = []
    for i in range(N_picos):
        x_maxima.append(maxima[i][0])
        y_maxima.append(maxima[i][1])

    for i in range(N_picos):
        #print 'Determinando largura do pico {}'.format(i)

        Point = np.array([maxima[i]])
        Continue = True
        level = V_picos[i]

        V_profundidade_pico[i] = None
        Vertices_contorno.append(None)

        while Continue:
            level = level - Delta/N_fracs_H
            # Construir contornos
            contour = plt.contour(ximage[:bins], yimage[:bins], field, levels = np.array([level]))
            
            # Numero de contornos gerados para esse level
            N_contour = len(contour.collections[0].get_paths())

            got_level_value = False

            for j in range(N_contour):
                contour_vertices_x = contour.collections[0].get_paths()[j].vertices[:,0]
                contour_vertices_y = contour.collections[0].get_paths()[j].vertices[:,1]
                
                # Test if contour is closed
                closed = (contour_vertices_x[0] == contour_vertices_x[-1]) & (contour_vertices_y[0] == contour_vertices_y[-1])
                if closed:
                    vertices = contour.collections[0].get_paths()[j].vertices
                    Path = matplotlib.path.Path(vertices)

                    # Does path contains the point?
                    if Path.contains_points(Point):
                        # Is it the only point in the path?
                        other_points = []
                        for k in range(N_picos):
                            if k != i: other_points.append(maxima[k])
                        other_points = np.array(other_points)

                        if len(other_points) == 0:
                            V_profundidade_pico[i] = level
                            Vertices_contorno[i] = vertices
                            
                            got_level_value = True
                            contour_index = j
                            
                        else:
                            if not any(Path.contains_points(other_points)):
                                V_profundidade_pico[i] = level
                                Vertices_contorno[i] = vertices
                                
                                got_level_value = True
                                contour_index = j

            if not got_level_value:
                Continue = False
                if np.isnan(V_profundidade_pico[i]): 
                    Area_pico[i] = np.nan
                else:
                    contour = plt.contour(ximage[:bins], yimage[:bins], field, levels = [V_profundidade_pico[i]])
                    vertices = contour.collections[0].get_paths()[contour_index].vertices
                    Area_pico[i] = area(vertices)
    plt.clf()
    return_dict = {'N_picos': N_picos, 'Area_pico': Area_pico, 'field': field, 'ximage': ximage, 'yimage': yimage,
                   'V_picos': V_picos, 'V_profundidade_pico': V_profundidade_pico, 'Vertices_contorno': Vertices_contorno,
                   'x_maxima': x_maxima, 'y_maxima': y_maxima}
    
    return return_dict

def bg_kernel(probs_base_name, 
                 probs_path,
                 gal_l, 
                 gal_b, 
                 age, 
                 metal, 
                 dist,
                 color = 'gr',
                 kernel_bg = 20,
                 bins = 130,
                 field_size = [1.0,1.0],
                 return_xy = True,
                 save_number = 0):

    min_gal_l = gal_l
    max_gal_l = gal_l + field_size[0]
    min_gal_b = gal_b
    max_gal_b = gal_b + field_size[1]

    # load file containing isochronal likelihoods
    filename = probs_path + get_name_probs(probs_base_name, gal_l, gal_b, age, metal, dist, save_number)
    probs_data = np.loadtxt(filename, delimiter = ',')

    star_gal_l = probs_data[:,0]
    star_gal_b = probs_data[:,1]

    if color == 'gr':
        star_probs = probs_data[:,8]
    elif color == 'gi':
        star_probs = probs_data[:,9]
    elif color == 'both':
        star_probs = probs_data[:,9]*probs_data[:,8]

    ##### Delete-this #######
    #star_probs = probs_data[:,8]*probs_data[:,9]
    #########################

    # filter_bg is used to represent the background    
    filter_bg = star_probs == star_probs.min()
    #filter_bg = star_cor_gr > 0.8
    
    # Background
    hbg, xbg, ybg, pbg = plt.hist2d(x = star_gal_l[filter_bg],
                                    y = star_gal_b[filter_bg],
                                    range = np.array([[min_gal_l, max_gal_l], [min_gal_b, max_gal_b]]),
                                    bins = bins)

    density_field_bg = ndimage.gaussian_filter(hbg, kernel_bg)
    density_field_bg = np.transpose(density_field_bg)

    density_field_bg = density_field_bg / (bins*bins*density_field_bg.sum())

    if return_xy: return density_field_bg, star_gal_b, star_gal_l, xbg, ybg, filter_bg
    else: return density_field_bg, star_gal_b, star_gal_l


def field_kernel(probs_base_name, 
                 probs_path,
                 gal_l, 
                 gal_b, 
                 age, 
                 metal, 
                 dist,
                 color = 'gr',
                 percentiles = [20,40,60,80,90],
                 kernels =  [30,20,15,10,8,5],
                 kernel_bg = 40,
                 bins = 130,
                 field_size = [1.0,1.0],
                 return_xy = True,
                 save_number = 0):

    min_gal_l = gal_l
    max_gal_l = gal_l + field_size[0]
    min_gal_b = gal_b
    max_gal_b = gal_b + field_size[1]

    # load file containing isochronal likelihoods
    filename = probs_path + get_name_probs(probs_base_name, gal_l, gal_b, age, metal, dist, save_number)
    probs_data = np.loadtxt(filename, delimiter = ',')

    star_gal_l = probs_data[:,0]
    star_gal_b = probs_data[:,1]
    star_mag_r = probs_data[:,2]
    star_mag_r_err = probs_data[:,3]
    star_mag_i = probs_data[:,4]
    star_mag_i_err = probs_data[:,5]
    star_cor_gr = probs_data[:,6]
    star_cor_gi = probs_data[:,7]

    if color == 'gr':
        star_probs = probs_data[:,8]
    elif color == 'gi':
        star_probs = probs_data[:,9]
    elif color == 'both':
        star_probs = probs_data[:,9]*probs_data[:,8]

    ##### multi-filter #######
    ##### star_probs = probs_data[:,8]*probs_data[:,9]
    #########################

    # Create two filters to separate values for which star_probs is minimmum from other values
    # filter_min is used to represent the background    
    #filter_min = star_probs == star_probs.min()
    #filter_min = star_cor_gr > 0.8
    #filter_min == star_probs.min()
    filter_min = star_probs == star_probs.min()
    filter_else = star_probs != star_probs.min()
    

    # obtain percentils
    percentils = np.percentile(star_probs[filter_else], percentiles)

    # set kernels
    star_kernel = np.zeros(len(star_probs[filter_else]))
    for i in range(len(star_probs[filter_else])):
        if (star_probs[filter_else][i] <= percentils[0]):
            star_kernel[i] = kernels[0]

        elif (star_probs[filter_else][i] > percentils[-1]):
            star_kernel[i] = kernels[-1]

        else:
            for j in range(1, len(percentils)):
                if (star_probs[filter_else][i] > percentils[j-1]) and (star_probs[filter_else][i] <= percentils[j]):
                    star_kernel[i] = kernels[j]

    #print star_kernel
    # Background

    hbg, xbg, ybg, pbg = plt.hist2d(x = star_gal_l[filter_min],
                                    y = star_gal_b[filter_min],
                                    range = np.array([[min_gal_l, max_gal_l], [min_gal_b, max_gal_b]]),
                                    bins = bins)

    density_field_bg = ndimage.gaussian_filter(hbg, kernel_bg)
    
    # Density field
    density_field = density_field_bg
    
    for i in range(len(kernels)):
        kernel_filter = star_kernel == kernels[i]

        h, ximage, yimage, p = plt.hist2d(x = star_gal_l[filter_else][kernel_filter],
                                          y = star_gal_b[filter_else][kernel_filter],
                                          range = np.array([[min_gal_l, max_gal_l], [min_gal_b, max_gal_b]]),
                                          bins = bins)

        density_field_i = ndimage.gaussian_filter(h, kernels[i])
        density_field = density_field + density_field_i

    # Normalizar density_field
    density_field = density_field / (bins*bins*density_field.sum())
    density_field = np.transpose(density_field)

    if return_xy: return density_field, star_gal_b, star_gal_l, xbg, ybg, filter_min
    else: return density_field, star_gal_b, star_gal_l


def filter_picos_1degree(picos, coord1_center, coord2_center):
    picos['x_maxima'] = np.array(picos['x_maxima'])
    picos['y_maxima'] = np.array(picos['y_maxima'])
    
    filtro_x = (picos['x_maxima'] >= (coord1_center - 0.5)) & (picos['x_maxima'] < (coord1_center + 0.5))
    filtro_y = (picos['y_maxima'] >= (coord2_center - 0.5)) & (picos['y_maxima'] < (coord2_center + 0.5))
    filtro = filtro_x & filtro_y
    
    for key in ['Area_pico', 'V_picos', 'V_profundidade_pico', 'Vertices_contorno', 'x_maxima', 'y_maxima']:
        picos[key] = np.array(picos[key])
        picos[key] = (picos[key])[filtro]
    
    picos['N_picos'] = len(picos['V_picos'])
    return picos

def percentil_kernel(isoc_ages, isoc_metal, d_seq, probs_path, probs_base_name, coords1, coords2, save_path,field_size,N_bins_cut,bins,nline=30,save_number = 0, kernels = [16,14,12,10,8,6], percentiles = [20,40,60,80,90], kernel_bg = 20):
    
    for coord1 in coords1:
        for coord2 in coords2:
            coord1_center = coord1+field_size[0]/2.0
            coord2_center = coord2+field_size[1]/2.0

            line = 0
            Nlines = len(isoc_ages)*len(isoc_metal)*len(d_seq)
            results = np.zeros((Nlines, 17))

            for age in isoc_ages: # for each age
                for metal in isoc_metal: # for each metallicity  
                    for dist in d_seq:
                        print 'age:', age, ' metal:', metal, ' dist:', dist
                        results[line, 0] = coord1
                        results[line, 1] = coord2
                        results[line, 2] = age
                        results[line, 3] = metal
                        results[line, 4] = dist/1e3 # salvar distancia em kpc

                        density_field, star_gal_b, star_gal_l, ximage, yimage, filter_bg = field_kernel(probs_base_name = probs_base_name,
                                                     probs_path = probs_path,
                                                     gal_l = coord1,
                                                     gal_b = coord2,
                                                     age = age,
                                                     metal = metal,
                                                     dist = dist,
                                                     bins = bins,
                                                     kernels =  kernels,
                                                     percentiles = percentiles,
                                                     kernel_bg = kernel_bg,
                                                     field_size = field_size,
                                                     return_xy = True,
                                                     save_number = save_number)

                        Background, star_gal_b, star_gal_l, xbg, ybg, filter_bg = bg_kernel(probs_base_name = probs_base_name,
                                                                                            probs_path = probs_path,
                                                                                            gal_l = coord1,
                                                                                            gal_b = coord2,
                                                                                            age = age,
                                                                                            metal = metal,
                                                                                            dist = dist,
                                                                                            bins = bins,
                                                                                            kernel_bg = kernel_bg,
                                                                                            field_size = field_size,
                                                                                            return_xy = True,
                                                                                            save_number = save_number)

                        A_mean, A_sdev = running_mean_sdev(density_field, nline = nline)

                        Mean   = density_field-A_mean
                        Desvio = (density_field-A_mean)/A_sdev

                        SN = running_SN(density_field, Background, nline = nline)

                        ###########################################################################################################
                        field_picos = caracterize_density_field(field = density_field, 
                                                                ximage = ximage, 
                                                                yimage = yimage, 
                                                                bins = bins)
                        field_picos = filter_picos_1degree(field_picos, coord1_center, coord2_center)
                        
                        Mean_picos = caracterize_density_field(field = Mean, 
                                                                ximage = ximage, 
                                                                yimage = yimage, 
                                                                bins = bins)
                        Mean_picos = filter_picos_1degree(Mean_picos, coord1_center, coord2_center)


                        Desvio_picos = caracterize_density_field(field = Desvio, 
                                                                ximage = ximage, 
                                                                yimage = yimage, 
                                                                bins = bins)
                        Desvio_picos = filter_picos_1degree(Desvio_picos, coord1_center, coord2_center)

                        SN_picos = caracterize_density_field(field = SN, 
                                                             ximage = ximage, 
                                                             yimage = yimage, 
                                                             bins = bins)
                        SN_picos = filter_picos_1degree(SN_picos, coord1_center, coord2_center)

                        ###########################################################################################################

                        view_xmin = coord1
                        view_xmax = coord1 + field_size[0]
                        view_ymin = coord2
                        view_ymax = coord2 + field_size[1]

                        fig, ax = plt.subplots(nrows=3, ncols=3,figsize=(18,14))

                        ############################################################

                        ax[0,0].scatter(star_gal_l, star_gal_b, c = '#0000CC', s = 15, alpha = 0.4, marker='o', edgecolor='none')
                        ax[0,0].set_xlim([view_xmin, view_xmax])
                        ax[0,0].set_ylim([view_ymin, view_ymax])
                        ax[0,0].grid(True)
                        ax[0,0].minorticks_on()

                        ax[0,1].scatter(star_gal_l[filter_bg], star_gal_b[filter_bg], c = '#0000CC', s = 15, alpha = 0.4, marker='o', edgecolor='none')
                        ax[0,1].set_xlim([view_xmin, view_xmax])
                        ax[0,1].set_ylim([view_ymin, view_ymax])
                        ax[0,1].grid(True)
                        ax[0,1].minorticks_on()

                        filename = probs_path + get_name_probs(probs_base_name, coord1, coord2, age, metal, dist, save_number = save_number)
                        probs_data = np.loadtxt(filename, delimiter = ',')
                        star_probs_gr = probs_data[:,8]
                        
                        ax[0,2].scatter(star_gal_l[~filter_bg], star_gal_b[~filter_bg], c = star_probs_gr[~filter_bg], s = 15, alpha = 0.8, marker='o', edgecolor='none')
                        ax[0,2].set_xlim([view_xmin, view_xmax])
                        ax[0,2].set_ylim([view_ymin, view_ymax])
                        ax[0,2].grid(True)
                        ax[0,2].minorticks_on()

                        #population_filter_l = (star_gal_l >= 220.4) & (star_gal_l <= 220.6)
                        #population_filter_b = (star_gal_b >= 50.3) & (star_gal_b <= 50.5)
                        #population_filter = population_filter_l & population_filter_b
                        #print len(population_filter[population_filter])
                        #ax[0,2].scatter(star_gal_l[population_filter], star_gal_b[population_filter], c = "#00AA00", s = 15, alpha = 0.8, marker='o', edgecolor='none')
                        ############################################################

                        X, Y = np.meshgrid(ximage, yimage)

                        ax[1,0].pcolormesh(X, Y, density_field, cmap = plt.cm.cubehelix)
                        ax[1,0].set_xlim([view_xmin, view_xmax])
                        ax[1,0].autoscale(False)

                        ############################################################

                        dx = field_picos['ximage'][1]-field_picos['ximage'][0]
                        dy = field_picos['yimage'][1]-field_picos['yimage'][0]

                        ax[1,0].plot(field_picos['x_maxima']+dx, field_picos['y_maxima']+dy, 'ro')

                        for i in range(len(field_picos['Area_pico'])):
                            if not np.isnan(field_picos['Area_pico'][i]):
                                ax[1,0].text(field_picos['x_maxima'][i], field_picos['y_maxima'][i], "{:4.3f}".format(field_picos['Area_pico'][i]))
                            if field_picos['Vertices_contorno'][i] != None:
                                X_vs = []
                                Y_vs = []
                                for j in range(len(field_picos['Vertices_contorno'][i])):
                                    X_vs.append(field_picos['Vertices_contorno'][i][j][0])
                                    Y_vs.append(field_picos['Vertices_contorno'][i][j][1])
                                ax[1,0].plot(X_vs, Y_vs, 'y')

                        ax[1,0].pcolormesh(X, Y, density_field, cmap = plt.cm.cubehelix)
                        ax[1,0].set_xlim([view_xmin, view_xmax])
                        ax[1,0].set_ylim([view_ymin, view_ymax])
                        ax[1,0].autoscale(False)
                        #ax[1,0].ylabel("latitude", fontsize=12)
                        #ax[1,0].xlabel("longitude", fontsize=12)
                        #ax[1,0].title('Density')
                        ax[1,0].grid(True)
                        ax[1,0].minorticks_on()

                        ############################################################

                        dx = Mean_picos['ximage'][1]-Mean_picos['ximage'][0]
                        dy = Mean_picos['yimage'][1]-Mean_picos['yimage'][0]

                        ax[1,1].plot(Mean_picos['x_maxima']+dx, Mean_picos['y_maxima']+dy, 'ro')

                        for i in range(len(Mean_picos['Area_pico'])):
                            if not np.isnan(Mean_picos['Area_pico'][i]):
                                ax[1,1].text(Mean_picos['x_maxima'][i], Mean_picos['y_maxima'][i], "{:4.3f}".format(Mean_picos['Area_pico'][i]))
                            if Mean_picos['Vertices_contorno'][i] != None:
                                X_vs = []
                                Y_vs = []
                                for j in range(len(Mean_picos['Vertices_contorno'][i])):
                                    X_vs.append(Mean_picos['Vertices_contorno'][i][j][0])
                                    Y_vs.append(Mean_picos['Vertices_contorno'][i][j][1])
                                ax[1,1].plot(X_vs, Y_vs, 'y')

                        ax[1,1].pcolormesh(X, Y, Mean, cmap = plt.cm.cubehelix)
                        ax[1,1].set_xlim([view_xmin, view_xmax])
                        ax[1,1].set_ylim([view_ymin, view_ymax])
                        ax[1,1].autoscale(False)
                        ax[1,1].grid(True)
                        ax[1,1].minorticks_on()

                        ############################################################

                        dx = Desvio_picos['ximage'][1]-Desvio_picos['ximage'][0]
                        dy = Desvio_picos['yimage'][1]-Desvio_picos['yimage'][0]

                        ax[1,2].plot(Desvio_picos['x_maxima']+dx, Desvio_picos['y_maxima']+dy, 'ro')

                        for i in range(len(Desvio_picos['Area_pico'])):
                            if not np.isnan(Desvio_picos['Area_pico'][i]):
                                ax[1,2].text(Desvio_picos['x_maxima'][i], Desvio_picos['y_maxima'][i], "{:4.3f}".format(Desvio_picos['Area_pico'][i]))
                            if Desvio_picos['Vertices_contorno'][i] != None:
                                X_vs = []
                                Y_vs = []
                                for j in range(len(Desvio_picos['Vertices_contorno'][i])):
                                    X_vs.append(Desvio_picos['Vertices_contorno'][i][j][0])
                                    Y_vs.append(Desvio_picos['Vertices_contorno'][i][j][1])
                                ax[1,2].plot(X_vs, Y_vs, 'y')

                        ax[1,2].pcolormesh(X, Y, Desvio, cmap = plt.cm.cubehelix)
                        ax[1,2].set_xlim([view_xmin, view_xmax])
                        ax[1,2].set_ylim([view_ymin, view_ymax])
                        ax[1,2].autoscale(False)
                        ax[1,2].grid(True)
                        ax[1,2].minorticks_on()
                        ax[1,2].set_title("teste")

                        ###########################################################################################################

                        dx = SN_picos['ximage'][1]-SN_picos['ximage'][0]
                        dy = SN_picos['yimage'][1]-SN_picos['yimage'][0]

                        ax[2,0].plot(SN_picos['x_maxima']+dx, SN_picos['y_maxima']+dy, 'ro')

                        for i in range(len(SN_picos['Area_pico'])):
                            if not np.isnan(SN_picos['Area_pico'][i]):
                                ax[2,0].text(SN_picos['x_maxima'][i], SN_picos['y_maxima'][i], "{:4.3f}".format(SN_picos['Area_pico'][i]))
                            if SN_picos['Vertices_contorno'][i] != None:
                                X_vs = []
                                Y_vs = []
                                for j in range(len(SN_picos['Vertices_contorno'][i])):
                                    X_vs.append(SN_picos['Vertices_contorno'][i][j][0])
                                    Y_vs.append(SN_picos['Vertices_contorno'][i][j][1])
                                ax[2,0].plot(X_vs, Y_vs, 'y')

                        ax[2,0].pcolormesh(X, Y, SN, cmap = plt.cm.cubehelix)
                        ax[2,0].set_xlim([view_xmin, view_xmax])
                        ax[2,0].set_ylim([view_ymin, view_ymax])
                        ax[2,0].autoscale(False)
                        ax[2,0].grid(True)
                        ax[2,0].minorticks_on()

                        #######
                        ax[2,1].pcolormesh(X, Y, SN, cmap = plt.cm.cubehelix)
                        ax[2,1].set_xlim([view_xmin, view_xmax])
                        ax[2,1].autoscale(False)
                        ax[2,1].grid(True)
                        ax[2,1].minorticks_on()

                        ax[2,2].pcolormesh(X, Y, density_field-SN, cmap = plt.cm.cubehelix)
                        ax[2,2].set_xlim([view_xmin, view_xmax])
                        ax[2,2].autoscale(False)
                        ax[2,2].grid(True)
                        ax[2,2].minorticks_on()

                        ###########################################################################################################

                        filename = probs_path + get_name_probs(probs_base_name, coord1, coord2, age, metal, dist, save_number = save_number)
                        probs_data = np.loadtxt(filename, delimiter = ',')

                        star_mag_r = probs_data[:,2]
                        star_mag_i = probs_data[:,4]
                        star_cor_gr = probs_data[:,6]
                        star_cor_gi = probs_data[:,7]
                        star_probs_gr = probs_data[:,8]
                        star_probs_gi = probs_data[:,9]

#                        ax[2,1].scatter(star_cor_gr[filter_bg], star_mag_r[filter_bg], c = '#000088', s = 15, alpha = 0.2, marker='o',
#                                        edgecolor='none')

#                        vmin_gr = np.log(star_probs_gr).min()
#                        vmax_gr = np.log(star_probs_gr).max()

#                        ax[2,1].scatter(star_cor_gr[~filter_bg], star_mag_r[~filter_bg], c = star_probs_gr[~filter_bg], 
#                                        alpha = 0.7, cmap = cm.jet, marker='o', edgecolor='none', vmin = vmin_gr, vmax = vmax_gr,
#                                        s = 15)
                        ax[2,1].scatter(star_cor_gr, star_mag_r, c = star_probs_gr, alpha = 0.5, cmap = cm.jet, marker='o', 
                                                                edgecolor='none', s = 15)

                        #ax[2,1].scatter(star_cor_gr[population_filter], star_mag_r[population_filter], c = "#00AA00", alpha = 0.5, cmap = cm.jet, marker='o', 
                        #                                        edgecolor='none', s = 30)
                        

                        ax[2,1].set_xlim(-0.5, 1.4)
                        ax[2,1].set_ylim(13, 24)       
                        ax[2,1].invert_yaxis()  

                        ###########################################################################################################

                        ax[2,2].scatter(star_cor_gi, star_mag_i, c = star_probs_gi, alpha = 0.5, cmap = cm.jet, marker='o', 
                                        edgecolor='none', s = 15)

                        

                        ax[2,2].set_xlim(-0.8, 2.8)
                        ax[2,2].set_ylim(13, 24.5)       
                        ax[2,2].invert_yaxis()    

                        save_graph_name = save_path+'graph_l{0}_b{1}_age{2}e9_Z{3}_d{4}kpc_{5}.png'.format(coord1,
                                                                                                       coord2,
                                                                                                       age/1e9,
                                                                                                       metal,
                                                                                                       dist/1000,
                                                                                                       save_number)

                        plt.savefig(save_graph_name, format='png')
                        plt.cla()
                        plt.clf()
                        plt.close()

                        ###########################################################################################################
                        jf = bins-N_bins_cut
                        j0 = N_bins_cut

                        first_line = "{l} {b} {t} {Z} {d} ".format(  l        = coord1,
                                                                     b        = coord2,
                                                                     t        = age,
                                                                     Z        = metal,
                                                                     d        = dist) + \
                                     "{:10e} {:10e} ".format(density_field[j0:jf,j0:jf].min(),
                                                             density_field[j0:jf,j0:jf].max()) + \
                                     "{:10e} {:10e} ".format(Mean[j0:jf,j0:jf].min(),
                                                             Mean[j0:jf,j0:jf].max()) + \
                                     "{:10e} {:10e} ".format(Desvio[j0:jf,j0:jf].min(),
                                                             Desvio[j0:jf,j0:jf].max()) + \
                                     "{:10e} {:10e}\n".format(SN[j0:jf,j0:jf].min(),
                                                              SN[j0:jf,j0:jf].max())

                        save_file_name = save_path+'results_l{0}_b{1}_age{2}e9_Z{3}_d{4}kpc_{5}.dat'.format(coord1,
                                                                                                        coord2,
                                                                                                        age/1e9,
                                                                                                        metal,
                                                                                                        dist/1000,
                                                                                                        save_number)

                        with open(save_file_name, 'w') as f:
                            f.write("")

                        with open(save_file_name, 'a') as f:
                            f.write("#  l    b          age metal dist     dens_min     dens_max      mean_min     mean_max    desvio_min   desvio_max        SN_max       SN_min\n")
                            f.write(first_line)
                            f.write("#   l     b       V_pico       V_prof   Area   Field\n")
                            
                            for w in range(field_picos['N_picos']):
                                line_w = "{:.2f} {:.2f} {:10e} {:10e} {:6.4f} {:s}\n".format(field_picos['x_maxima'][w],
                                                                                             field_picos['y_maxima'][w],
                                                                                             field_picos['V_picos'][w],
                                                                                             field_picos['V_profundidade_pico'][w],
                                                                                             field_picos['Area_pico'][w],
                                                                                             'Density')
                                f.write(line_w)

                            for w in range(Mean_picos['N_picos']):
                                line_w = "{:.2f} {:.2f} {:10e} {:10e} {:6.4f} {:s}\n".format(Mean_picos['x_maxima'][w],
                                                                                             Mean_picos['y_maxima'][w],
                                                                                             Mean_picos['V_picos'][w],
                                                                                             Mean_picos['V_profundidade_pico'][w],
                                                                                             Mean_picos['Area_pico'][w],
                                                                                             'Mean')
                                f.write(line_w)

                            for w in range(Desvio_picos['N_picos']):
                                line_w = "{:.2f} {:.2f} {:10e} {:10e} {:6.4f} {:s}\n".format(Desvio_picos['x_maxima'][w],
                                                                                             Desvio_picos['y_maxima'][w],
                                                                                             Desvio_picos['V_picos'][w],
                                                                                             Desvio_picos['V_profundidade_pico'][w],
                                                                                             Desvio_picos['Area_pico'][w],
                                                                                             'Desvio')
                                f.write(line_w)

                            for w in range(SN_picos['N_picos']):
                                line_w = "{:.2f} {:.2f} {:10e} {:10e} {:6.4f} {:s}\n".format(SN_picos['x_maxima'][w],
                                                                                             SN_picos['y_maxima'][w],
                                                                                             SN_picos['V_picos'][w],
                                                                                             SN_picos['V_profundidade_pico'][w],
                                                                                             SN_picos['Area_pico'][w],
                                                                                             'SN')
                                f.write(line_w)













