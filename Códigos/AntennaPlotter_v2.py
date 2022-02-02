import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['k', 'royalblue', 'r', 'gold', 'limegreen', 'navy', 'crimson', 'turquoise', 'darkorange', 'darkgreen'])

''' PROGRAMA PARA HACER GRAFICAS DE CAMPO ELECTRICO QUE SALEN DE ZHAIRES
    DOS DOMINIOS: TIEMPO Y FRECUENCIA
    
    FORMATO DE ARCHIVOS EN DOMINIO TEMPORAL:
        # -----------Time Fresnel file----------
        # Columns:
        #        1 - Shower Number
        #        2 - Antenna Number
        #        3 - Antenna X (m)
        #        4 - Antenna Y (m)
        #        5 - Antenna Z (m)
        #        6 - Time (ns)
        #        7 - |A| (V/M)
        #        8 -  Ax  (V/M)
        #        9 -  Ay  (V/M)
        #        10 - Az  (V/M)
        #        11 -|E| (V/M)
        #        12 - Ex  (V/M)
        #        13 - Ey  (V/M)
        #        14 - Ez  (V/M)
        #
        
    FORMATO DE ARCHIVOS EN FRECUENCIAS:
        # -----------Freq Fresnel file----------
        # Columns:
        #        1 - Shower Number
        #        2 - Antenna Number
        #        3 - Antenna X (m)
        #        4 - Antenna Y (m)
        #        5 - Antenna Z (m)
        #        6 - Frequency #
        #        7 - Frequency (MHz)
        #        8 - |E| (V/M MHz)
        #        9 -  Ex  (V/M MHz)
        #        10 - Ey  (V/M MHz)
        #        11 - Ez  (V/M MHz)
        
    EL PROGRAMA ESTA DISEÑADO PARA ARCHIVOS CON UNA UNICA SHOWER.
    LAS GRAFICAS IMPLEMENTADAS SON:
        Dominio t: E, E_xyz, A, A_xyz vs. t en antena concreta
                   Maximo de la señal (E, E_xyz, A, A_xyz) en t vs. coord xyz de antenas
        
        Dominio f: E, E_xyz vs. freq en antena concreta
                   Para freq concreta: E, E_xyz vs coord xyz de antenas
                   
'''
#%%
################ GRAFICAS EN DOMINIO TEMPORAL #############################

''' Codigos de graficas

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
rootdir = 'Graph_paperANITA'

graphs = [1, 8]

antenas = [[i+1 for i in range(10, 15)], [32]]

# ESCALAS EN LA GRAFICA: 'linear' ó 'log'
xscale = 'linear'
yscale = 'linear'

# LIMITES EN LA GRAFICA
xlims = []
ylims = []

#LEYENDA EN LA GRAFICA
legend = True
#%%
file_time = rootdir+'/timefresnel-only300MHz.dat'
data_time = np.loadtxt(file_time, comments = '#').T

#%%
################################## CAJA NEGRA ############################

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
#%%
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


fig = plt.figure(1)
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
    
#%%

################ GRAFICAS EN DOMINIO DE FRECUENCIAS #############################

''' Codigos de graficas

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
rootdir = 'Graph_paperANITA'
file_time = rootdir+'/freqfresnel-variousfreqs.dat'

graphs = [24]

antenas = [[i+1 for i in range(50, 100)]]

freq = [[0, 1, 2, 4]]

# ESCALAS EN LA GRAFICA: 'linear' ó 'log'
xscale = 'linear'
yscale = 'linear'

# LIMITES EN LA GRAFICA
xlims = []
ylims = []

#LEYENDA EN LA GRAFICA
legend = True

# FORMATO DE LA GRAFICA
formato = '-' # '-', ':', '--'


data_time = np.loadtxt(file_time, comments = '#').T

#%%
################################## CAJA NEGRA ############################

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


#%%

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


fig = plt.figure(1)
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
    