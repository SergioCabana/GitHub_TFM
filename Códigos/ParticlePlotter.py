import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['k', 'royalblue', 'r', 'gold', 'limegreen', 'navy', 'crimson', 'turquoise', 'darkorange', 'darkgreen'])



''' ESTE PROGRAMA BARRE LOS ARCHIVOS DENTRO DE UN DETERMINADO DIRECTORIO
    ESCOGIENDO QUE FICHEROS SE CONSIDERAN, LOS GRAFICA
'''

''' EXTENSIONES DE LOS ARCHIVOS
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
    K+K-'            = 5
    All chgd. pcles. = 6
    All pcles.       = 7
    
# Si se quiere la grafica de Long. development para e+e-, mu+mu-, se indica:
    graf = 0
    pcles = [2, 3]
    
# Si ademas se quiere una seleccion mas estricta de archivos, se añaden en constraint
los strings por los que se quiere buscar:
constraint = '9km'

Si se quieren separar los archivos escogidos en dos criterios, se añaden en extra
extra = ['0km', '9km']

# Si no se quiere separar, se indica con un string vacio:
extra = ['']

# omitextra = True si se quiere que los extras no aparezcan en las leyendas

# UseSlant = True si se squiere añadir un eje extra con profundidad slanted
# En ese caso hay que indicar el angulo theta de la cascada
'''
################################### INPUTS AQUI ############################################
# DIRECTORIO
rootdir = 'Carac_UG'

# GRAFICA QUE SE QUIERE
graf  = 0
pcles = [2]
constraint = ['5km', '95deg']
extras = ['1PeV', '10PeV', '100PeV', '1EeV', '10EeV']

# ESCALAS EN LA GRAFICA: 'linear' ó 'log'
xscale = 'linear'
yscale = 'log'

# LIMITES EN LA GRAFICA
xlims = [0, 1050]
ylims = []

#LEYENDA EN LA GRAFICA
legend = True
cols = 1
omitextra = False

# EJE X INVERTIDO, IMPORTANTE PARA GRAFICAS UPGOING
invert_x = True

UseSlant = False
theta    = 123

#################################### CAJA NEGRA #######################################

ext_LD  = [1001, 1022, 1205, 1207, 1211, 1213, 1291, 1293]
ext_ELD = [500+ext for ext in ext_LD]
ext_LatD = [1000+ext for ext in ext_LD]
ext_EDG = [1500+ext for ext in ext_LD]

ext = [ ['.t'+str(numero) for numero in lista] for lista in [ext_LD, ext_ELD, ext_LatD, ext_EDG] ]

data = []
# (perdon por el spaghetti)
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        for p in pcles:
            if file.endswith(ext[graf][p]) and all([c in file for c in constraint]): 
                # si es grafica y particula adecuada
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
    
# data es una lista con los nombres de los archivos que nos interesan, junto con que 
# particula es cada uno y que extra llevan
# solo hay que hacer las graficas ahora

lbl = [r'$\gamma$', 'p', r'$e\pm$' , r'$\mu\pm$', r'$\pi\pm$', r'$K\pm$', 'All chgd. pcles.', 'All pcles.']
xlabels = [r'$X_v$ (g/cm2)', r'$X_v$ (g/cm2)', 'Distance to core (m)', 'E (GeV)']
ylabels = ['N', 'E (GeV)', 'N', 'N']


fig = plt.figure(1)
ax = fig.add_subplot(111)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
for file, p, extra in data:
    datos = np.loadtxt(file, comments = '#').T
    xvalues = datos[1]
    yvalues = datos[2]
    if extra == '' or omitextra:
        ax.step(xvalues, yvalues, where = 'mid', label = lbl[p], linewidth = 2.0)
    else:
        ax.step(xvalues, yvalues, where = 'mid', label = lbl[p]+', '+extra, linewidth = 2.0)

ax.set_xlabel(xlabels[graf], size = 12)
ax.set_ylabel(ylabels[graf], size = 12) 
ax.set_xscale(xscale)
ax.set_yscale(yscale)

if UseSlant: # Eje extra con profundidad slanted
    if graf == 0 or graf == 1:
        ax2 = ax.twiny()
        ax2.set_xlabel(r'$X_s$ (g/cm2)')
        ax2.set_xlim([xlim/np.cos(theta) for xlim in xlims])
    else:
        raise TypeError('Eje x no es profundidad vertical en esta grafica')
    
if len(xlims)>0:
    ax.set_xlim(xlims)
if len(ylims)>0:
    ax.set_ylim(ylims)
    
if invert_x:
    if UseSlant:
        ax2.invert_xaxis()
    ax.invert_xaxis()
    
if legend:
    ax.legend(loc = 'best', prop={'size':12}, ncol = cols, fancybox = True)
    
    

