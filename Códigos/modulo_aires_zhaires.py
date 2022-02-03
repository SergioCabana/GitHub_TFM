import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['k', 'royalblue', 'r', 'gold', 'limegreen', 'navy', 'crimson', 'turquoise', 'darkorange', 'darkgreen'])

def gcm2toh(t):
    ''' Transforma profundidad vertical t (g/cm2) a alturas en km
        Utiliza parametrizacion de Linsley
        h1km, h2km son inversas de dicha parametrizacion
    '''
    
    if t < 0.:
        raise TypeError('Negative atmospheric thickness')
    
    if t < 0.00128293:
        h2km = 1e4*(1.12829e-2-t)
        return h2km
    
    def h1km(t, params):
        a, b, c = params
        return c*np.log(b/(t-a))
    
    params = [ (-1.86556e2, 1.2227e3, 9.9419), (-9.4919e1, 1.1449e3, 8.7815), \
               (6.1289e-1, 1.3056e3, 6.3614),  (0.0, 5.4018e2, 7.7217) ]

    
    if t <= 3.0395:
        return h1km(t, params[3])
    
    elif t <= 271.6991:
        return h1km(t, params[2])
    
    elif t <= 631.1:
        return h1km(t, params[1])
    
    elif t <= 2004.7:
        return h1km(t, params[0])
    
    else:
        raise TypeError('Atmospheric thickness above 2004.647')
        

def htogcm2(h):
    ''' Transforma altura en km a profundidad vertical t (g/cm2)
        Utiliza parametrizacion de Linsley, t1 y t2
    '''
    
    if h < -5.801:
        raise TypeError('Altitude lower than -5.8 km')
    
    if h > 112.8:
        return 0.
    
    if h >= 100.0:
        t2 = 1.12829e-2-h/1e4
        return t2
    
    def t1(h, params):
        a, b, c = params
        return a+b*np.exp(-h/c)
    
    params = [ (-1.86556e2, 1.2227e3, 9.9419), (-9.4919e1, 1.1449e3, 8.7815), \
               (6.1289e-1, 1.3056e3, 6.3614),  (0.0, 5.4018e2, 7.7217) ]

    
    if h >= 40.0:
        return t1(h, params[3])
    
    elif h >= 10.0:
        return t1(h, params[2])
    
    elif h >= 4.0 :
        return t1(h, params[1])
    
    else:
        return t1(h, params[0])


def CherenkovRing(Xmax, ground, theta, UseSlant = True):
    ''' Devuelve parametros de la elipse Cherenkov a nivel del suelo
        Xmax -> Profundidad del maximo de la cascada en gcm2
                UseSlant True -> Este valor es medido a lo largo del eje de la cascada
                UseSlant False -> Este valor es medido en vertical
        ground -> Altitud del suelo a.s.l.
        theta -> Angulo cenital en deg
    '''
    R0 = 325.
    
    theta     = theta * np.pi / 180. # conversion a radianes
    costet    = np.cos(theta)
    
    Xmax_vert = Xmax 
    
    if UseSlant:
        Xmax_vert = Xmax * costet # profundidad vertical del maximo
    
    alt_max = gcm2toh(Xmax_vert)*1e3
    h_max   = alt_max-ground       # altura del maximo de la cascada (vertical)
    d       = h_max/costet         # distancia desde ground al maximo a lo largo del eje
    
    # indice de refraccion
    Rh = R0 * np.exp( -0.1218 * alt_max / 1e3 )
    nh = 1. + 1e-6 * Rh    # indice de refraccion a la altura del maximo
    
    costetc = 1./nh          # coseno del angulo cerenkov
    
    thetac  = np.arccos(costetc) 
    sintetc = np.sin(thetac)
    
    ########### CALCULO APROXIMADO, PRIMER OUTPUT #####################
    
    r1 = h_max * ( np.tan(theta+thetac) - np.tan(theta) )
    r2 = h_max * np.tan(thetac) / costet
    
    print('------ INPUTS -------\n')
    if UseSlant:
        print('X_max slanted (g/cm2): %.4f'%Xmax)
    print('X_max vert. (g/cm2): %.4f'%Xmax_vert)
    print('Angulo cenital (deg): %.4f'%( theta * 180/np.pi ))
    print('Altitud maximo (m): %.4f'%alt_max)
    print('Altitud sobre suelo del maximo (m): %.4f'%h_max)
    print('Indice de refraccion, nh: %.8f'%nh)
    print('Angulo Cerenkov (deg): %.4f\n'%( thetac * 180/np.pi ))
    print('--------- CALCULO APROXIMADO ----------\n')
    print('r1 (m): %.5f'%r1)
    print('r2 (m): %.5f\n'%r2)
    
    ############ CALCULO EXACTO ######################
    
    cosplus  = np.cos(theta + thetac)
    cosminus = np.cos(theta - thetac)
    a1       = d * sintetc / cosminus
    a2       = d * sintetc / cosplus
    a        = d * np.sin(2.*thetac) * costet / (2.*cosplus*cosminus)
    b        = d * sintetc * costet / np.sqrt(cosplus*cosminus)
    epsilon  = d * sintetc * sintetc * np.sin(theta) / (cosplus*cosminus)
    
    print('------ CALCULO EXACTO ----------\n')
    print('a1 (m): %.4f'%a1)
    print('a2 (m): %.4f'%a2)
    print('a (m): %.4f'%a)
    print('b (m): %.4f'%b)
    print('Ecc. : %.4f'%epsilon)
    print('D (m): %.4f'%d)
    print('4,5 * r2 (m) = %.4f'%(4.5*r2))
    
def readfile(filename):
    '''Lee el archivo de datos Aires que le demos, ademas busca el valor del ground
    
    '''
    file = open(filename, 'r')
    lineas = file.readlines()
    file.close()
    
    dataset = []
    for line in lineas:
        if line[0:5] == '# GRD':
            ground = line.split()
        elif line[0] != '#':
            dataset += [line.split()]
    
    datos = np.array(dataset, dtype = float)
    
    ground = float(ground[2])
    
    return datos[:,1:], ground

def Aires_Plot(graf, pcles, rootdir, const, extras, xscale='linear', yscale='linear', \
            xlims = [], ylims=[], legend = True, cols = 1, omitextra=False, \
            UG = False, slant = False):
    ''' 
        CADA TIPO DE GRAFICA TIENE UN NUMERO
        
        Longitudinal Development  = 0
        
        Energy Long. Development  = 1
        
        Lateral distribution      = 2
        
        Energy distrib. at ground = 3
        
        CADA PARTICULA TIENE UN NUMERO:
            
        gamma            = 0
        
        p                = 1
        
        e+e-             = 2
        
        mu+mu-           = 3
        
        pi+pi-           = 4
        
        K+K-             = 5
        
        All chgd. pcles. = 6
        
        All pcles.       = 7
        
        Si se quiere la grafica de Long. development para e+e-, mu+mu-, se indica:
            graf = 0
            pcles = [2, 3]
        
        Si ademas se quiere una seleccion mas estricta de archivos, se añaden en const
        los strings por los que se quiere buscar:
            const = ['9km']
    
        Si se quieren separar los archivos escogidos en dos criterios, se añaden en extra
            extra = ['0km', '9km']
    
        Si no se quiere separar, se indica con un string vacio:
            extra = ['']
    
        omitextra = True si se quiere que los extras no aparezcan en las leyendas
    
        UG = True si tenemos datos upgoing
        
        slant = True si tenemos profundidad slanted en el eje x
        
        cols -> columnas de la leyenda
    '''
    ext_LD  = [1001, 1022, 1205, 1207, 1211, 1213, 1291, 1293]
    ext_ELD = [500+ext for ext in ext_LD]
    ext_LatD = [1000+ext for ext in ext_LD]
    ext_EDG = [1500+ext for ext in ext_LD]

    ext = [ ['.t'+str(numero) for numero in lista] for lista in [ext_LD, ext_ELD, ext_LatD, ext_EDG] ]


    lbl = [r'$\gamma$', 'p', r'$e\pm$' , r'$\mu\pm$', r'$\pi\pm$', r'$K\pm$', 'All chgd. pcles.', 'All pcles.']
    xlabels = [r'$X_v$ (g/cm2)', r'$X_v$ (g/cm2)', 'Distance to core (m)', 'E (GeV)']
    ylabels = ['N', 'E (GeV)', 'N', 'N']
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    
    
    data = []

    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            for p in pcles:
                if file.endswith(ext[graf][p]) and all([c in file for c in const]): 
                    # si es grafica y particula adecuada, y cumple constraints
                    for extra in extras:
                        if extra == '' or extra in file:
                            data.append([subdir + os.sep + file, p, extra])

    if len(data) != len(extras)*len(pcles):
        raise TypeError('Faltan archivos necesarios para la grafica que se pide')
        
    orden = []
    for extra in extras:
        for serie in data:
            if extra in serie:
                orden.append(serie) # solo para que la leyenda salga en el orden de extras
                
    data = orden        

    for file, p, extra in data:
        values, grd = readfile(file)
        
        xvalues = values[:,0]
        yvalues = values[:,1]
        
        if UG and slant:
            if graf == 0 or graf == 1: # depth in x axis
                xvalues = grd - xvalues # depth traversed in upward direction
            
        if extra == '' or omitextra:
            ax.step(xvalues, yvalues, where = 'mid', label = lbl[p], linewidth = 2.0)
        else:
            ax.step(xvalues, yvalues, where = 'mid', label = lbl[p]+', '+extra, linewidth = 2.0)
        
    ax.set_xlabel(xlabels[graf], size = 12)
    ax.set_ylabel(ylabels[graf], size = 12) 
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
        
    if slant and UG:
        if graf == 0 or graf == 1: # depth in x axis
            ax.set_xlabel(r'$X_s$ (g/cm2) (upward)')
    elif slant:
        if graf == 0 or graf == 1: # depth in x axis
            ax.set_xlabel(r'$X_s$ (g/cm2)')
    elif UG:
        ax.invert_xaxis()
               
    if legend:
        ax.legend(loc = 'best', prop={'size':12}, ncol = cols, fancybox = True)
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
    
    return fig

def ZHAireS_Plot_t(graphs, antenas, file, xscale='linear', yscale='linear', xlims=[], \
                   ylims = [], legend = True):
    ''' Codigos de graficas en dominio temporal

    A vs. t, antenas especificadas  = 1
    
    Ax vs. t, antenas especificadas = 2
    
    Ay vs. t, antenas especificadas = 3
    
    Az vs. t, antenas especificadas = 4
    
    E vs. t, antenas especificadas  = 5
    
    Ex vs. t, antenas especificadas = 6
    
    Ey vs. t, antenas especificadas = 7
    
    Ez vs. t, antenas especificadas = 8
    
    Para representacion del maximo en tiempo de la magnitud representada en graf. i:
        
    Frente a coord x = i+10
    
    Frente a coord y = i+20
    
    Frente a coord z = i+30
    
    Por ejemplo, se quiere el maximo de la señal de E en tiempo respecto a coord y de antenas:
        5 (grafica de E) + 20 (coord y) = 25
    
    ESPECIFICACION DE GRAFICAS Y ANTENAS:
        El index de antenas empieza en 1 como en Aires
        
        Si queremos más de una grafica, damos: 
            graphs = [11, 13]
        Podemos especificar que antenas entran en cada grafica:
            antenas = [[1,2,3,4,5], [6,7,8,9,10]]
            
        Siempre una lista de antenas por grafica
        
        Si se quieren poner todas
            antenas = 'all'

'''
    data_time = np.loadtxt(file, comments = '#').T
    
    if len(antenas) != len(graphs):
        raise TypeError('longitud equivocada')
        
    if (not all([g>10 for g in graphs])) and (not all([g<10 for g in graphs])):
        raise TypeError('Todas las graficas deben ser del mismo tipo')
    
    plot_maxs_vs_coords = all([g>10 for g in graphs])
    
    n_ant = int(max(data_time[1])) #numero de antenas
    
    ant_coord = [] # aqui guardamos las coordenadas de las antenas
    
    i_ant = []
    for i in range(n_ant):
        index = list(data_time[1]).index(i+1)
        i_ant.append(index)
        ant_coord.append([data_time[2][index], data_time[3][index], data_time[4][index]])
        # con esto ya tenemos las coordenadas de las antenas
    
    i_ant.append(len(data_time.T))
    
    if antenas == 'all':
        antenas = [[i+1 for i in range(n_ant)] for _ in range(len(graphs))]
    
    graph_list = []
    '''
    Para poder hacer mas de una grafica en el mismo canvas, guardamos en una lista
    los puntos a representar:
        graph_list = [[xs_1, ys_1], [xs_2, ys_2], ...]
    '''

    if plot_maxs_vs_coords:
        
        for g, ant in list(zip(graphs, antenas)):
            coord = g // 10 -1
            g = g % 10
            x = np.array([ant_coord[a-1][coord] for a in ant])
            y = np.array([np.max(data_time[int(5+g)][i_ant[a-1]:i_ant[a]]) for a in ant])
            
            graph_list.append([x, y])
            
    else: 
        for i in range(len(graphs)):
            for a in antenas[i]:
                x = data_time[5][i_ant[a-1]:i_ant[a]]
                y = data_time[5+graphs[i]][i_ant[a-1]:i_ant[a]]
            
                graph_list.append([x, y])
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    magnitudes = [r'$A$ (V/m)', r'$A_x$ (V/m)', r'$A_y$ (V/m)', r'$A_z$ (V/m)', r'$E$ (V/m)', r'$E_x$ (V/m)', r'$E_y$ (V/m)', r'$E_z$ (V/m)']
    lbl_coords = ['N-S', 'E-W', 'Vertical']
    labels = []
    
    if plot_maxs_vs_coords:
        if all([((g%10-1)-(graphs[0]%10-1))==0 for g in graphs]):
            labels = [lbl_coords[g//10-1] for g in graphs]
            ax.set_ylabel(magnitudes[graphs[0]%10-1], size = 12)
        else:
            labels =  [magnitudes[g%10-1]+', '+lbl_coords[g//10-1] for g in graphs]
            
        for i in range(len(graphs)):#g, data in list(zip(graphs, graph_list)):
            ax.plot(graph_list[i][0], graph_list[i][1], 'o', label = labels[i])
            
        ax.set_xlabel('Dist. to shower core (m)', size = 12)
     
    else:
        if all([(g-graphs[0])==0 for g in graphs]):
            labels = ['Antena ' for _ in range(len(graphs))]
            ax.set_ylabel(magnitudes[graphs[0]-1], size = 12)
            
        else:
            labels =  [magnitudes[g-1]+', Antena ' for g in graphs]
        index = 0
        for i in range(len(graphs)):#g, data, a in list(zip(graphs, graph_list, antenas)):
            for j in range(len(antenas[i])):
                ax.plot(graph_list[index][0], graph_list[index][1], label = labels[i]+'%d'%antenas[i][j])
                index += 1
        ax.set_xlabel('t (ns)')
        
        
    if legend:
        ax.legend(loc='best', prop={'size':12})
    
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
        
    return fig

def ZHAireS_Plot_f(graphs, antenas, freq, file, xscale='linear', yscale='linear', xlims=[], \
                   ylims = [], legend = True, formato = '-'):
    ''' Codigos de graficas en dominio de frecuencias

    E vs. f, antenas especificadas  = 2
    
    Ex vs. f, antenas especificadas = 3
    
    Ey vs. f, antenas especificadas = 4
    
    Ez vs. f, antenas especificadas = 5
    
    Para representacion de la magnitud representada en graf. i, a frecuencia concreta f:
        
    Frente a coord x = i+10
    
    Frente a coord y = i+20
    
    Frente a coord z = i+30
    
    Por ejemplo, se quiere la señal de E a 300MHz respecto a coord y de antenas:
        2 (grafica de E) + 20 (coord y) = 22
    
    ESPECIFICACION DE GRAFICAS Y ANTENAS:
        El index de antenas empieza en 1 como en Aires
        
        Si queremos más de una grafica, damos: 
            graphs = [12, 13]
        Podemos especificar que antenas entran en cada grafica:
            antenas = [[1,2,3,4,5], [6,7,8,9,10]]
            
        Siempre una lista de antenas por grafica
        
        Si se quieren poner todas
            antenas = 'all'

    ESPECIFICACION DE FRECUENCIAS
    Si tenemos varias grapficas del tipo componente de fourier para magnityud x
     en funcion de coordenadas de antenas, damos una lista con las frecuencias concreta
     para cada grafico
     
        freq = [[i's donde i es el indice en la lista de frecuencias que sale], ...]
'''
    data_time = np.loadtxt(file, comments = '#').T
    
    if len(antenas) != len(graphs):
        raise TypeError('longitud equivocada de antenas')
        
    if (not all([g>10 for g in graphs])) and (not all([g<10 for g in graphs])):
        raise TypeError('Todas las graficas deben ser del mismo tipo')
    
    plot_maxs_vs_coords = all([g>10 for g in graphs])
    
    if plot_maxs_vs_coords and len(freq)!=len(graphs):
        raise TypeError('Indica la frecuencia en TODAS las graficas que quieras')
    
    n_ant = int(max(data_time[1])) #numero de antenas
    
    ant_coord = [] # aqui guardamos las coordenadas de las antenas
    
    i_ant = []
    for i in range(n_ant):
        index = list(data_time[1]).index(i+1)
        i_ant.append(index)
        ant_coord.append([data_time[2][index], data_time[3][index], data_time[4][index]])
        # con esto ya tenemos las coordenadas de las antenas
    
    i_ant.append(len(data_time.T))
    
    freq_list = [data_time[6][i] for i in range(i_ant[0], i_ant[1])]
    
    print('Frecuencias (MHz): ', freq_list)
    
    if antenas == 'all':
        antenas = [[i+1 for i in range(n_ant)] for _ in range(len(graphs))]
    
    
    graph_list = []
    '''
    Para poder hacer mas de una grafica en el mismo canvas, guardamos en una lista
    los puntos a representar:
        graph_list = [[xs_1, ys_1], [xs_2, ys_2], ...]
    '''
    if plot_maxs_vs_coords:
        
        for g, ant, f in list(zip(graphs, antenas, freq)):
            coord = g // 10 -1
            g = g % 10
            for frecuencia in f:
                x = np.array([ant_coord[a-1][coord] for a in ant])
                y = np.array([data_time[int(5+g)][i_ant[a-1]+frecuencia] for a in ant])
            
                graph_list.append([x, y])
            
    else: 
        for i in range(len(graphs)):
            for a in antenas[i]:
                x = data_time[6][i_ant[a-1]:i_ant[a]]
                y = data_time[5+graphs[i]][i_ant[a-1]:i_ant[a]]
            
                graph_list.append([x, y])
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    magnitudes = [r'$E$ (V/m/MHz)', r'$E_x$ (V/m/MHz)', r'$E_y$ (V/m/MHz)', r'$E_z$ (V/m/MHz)']
    lbl_coords = ['N-S', 'E-W', 'Vertical']
    labels = []
    
    if plot_maxs_vs_coords:
        if all([((g%10-1)-(graphs[0]%10-1))==0 for g in graphs]):
            
            labels = [lbl_coords[graphs[i]//10-1]+', '+str(freq_list[f])+' MHz'  for i in range(len(graphs)) for f in freq[i]]
            
            ax.set_ylabel(magnitudes[graphs[0]%10-2], size = 12)
        else:
            labels =  [magnitudes[graphs[i]%10-2]+', '+lbl_coords[graphs[i]//10-1]+', '+str(freq_list[f])+' MHz'  for i in range(len(graphs)) for f in freq[i]]
            
        for i in range(len(graph_list)):
            ax.plot(graph_list[i][0], graph_list[i][1], formato, label = labels[i])
            
        ax.set_xlabel('Dist. to shower core (m)', size = 12)
     
    else:
        if all([(g-graphs[0])==0 for g in graphs]):
            labels = ['Antena ' for _ in range(len(graphs))]
            ax.set_ylabel(magnitudes[graphs[0]-2], size = 12)
            
        else:
            labels =  [magnitudes[g-1]+', Antena ' for g in graphs]
            
        for i in range(len(graphs)):#g, data, a in list(zip(graphs, graph_list, antenas)):
            for j in range(len(antenas[i])):
                index = i*len(antenas[i])+j
                ax.plot(graph_list[index][0], graph_list[index][1], formato, label = labels[i]+'%d'%antenas[i][j])
        ax.set_xlabel('Freq (MHz)')
        
        
    if legend:
        ax.legend(loc='best', prop={'size':12})
    
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
        
    return fig

##########################################################################################

# rootdir = 'Carac_UG'
# graf  = 0
# pcles = [3]
# const = ['5km', '95deg', 'Xs']
# extras = ['1PeV', '10PeV', '100PeV', '1EeV', '10EeV']

# Aires_Plot(graf, pcles, rootdir, const, extras, yscale='log', cols=1, UG = True, slant=True)

###########################################################################################

# rootdir = 'Graph_paperANITA'

# graphs = [1, 8]

# antenas = [[i+1 for i in range(10, 15)], [32]]

# file = rootdir+'/timefresnel-only300MHz.dat'

# ZHAireS_Plot_t(graphs, antenas, file)

###########################################################################################

# rootdir = 'Graph_paperANITA'
# file = rootdir+'/freqfresnel-variousfreqs.dat'

# graphs = [24]

# antenas = [[i+1 for i in range(50, 100)]]

# freq = [[0, 1, 2, 4]]

# ZHAireS_Plot_f(graphs, antenas, freq, file)