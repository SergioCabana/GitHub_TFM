import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import pandas            as pd
import os

from scipy.fftpack       import fft, fftfreq, fftshift


colores = ['k', 'royalblue', 'r', 'gold', 'limegreen', 'navy', 'crimson', 'turquoise', 'darkorange', 'darkgreen']
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colores)

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
    
def Long(h, theta):
    ''' Convierte la altura h (km) de un punto del eje de la cascada a la 
        distancia en km medida a lo largo del eje de dicho punto
        hasta el nivel del suelo
        
        theta es el angulo cenital de la cascada a nivel del suelo
    '''
    
    RT = 6.3781e3 # radio de la tierra en km
    
    d  = np.sqrt((RT+h)**2-RT**2*np.sin(theta)**2)-RT*np.cos(theta)
    
    return d


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
        if line[0:5] == '# GRD' or line[0:5] == '#  GT':
            ground = line.split()
        elif line[0] != '#':
            dataset += [line.split()]
    
    datos = np.array(dataset, dtype = float)
    
    ground = float(ground[2])
    
    return datos[:,1:], ground

def cos_localtheta(height, theta):
    ''' cosine of zenital angle at height height
    '''
    RT = 6370
    return np.sqrt(1. - (RT**2 * np.sin(theta)**2/(RT+height)**2))


def Xs_of_injh(inj_h, theta, step):
    ''' Returns the (approximate) value of slanted depth, measured along shower axis,
        of the injection point (height of the first interacion), given the 
        zenital angle of the upgoing shower and the injection height in km
        
        Default consideration: Xs(h -> infty) = 0, as in Aires
        
        step: 
            Size of height increment. As I don't know the functional form
            of the density with h, I approximate the integral
            
            Danger when inj_h is very big, depends on htogcm2 function !
            if h > 112.8 km, then Xv = 0
    '''
    if inj_h > 112.8:
        raise TypeError('Height cannot be more than 112.8 km')
    
    heights = np.arange(inj_h, 112.8+step, step)
    
    Xs = 0
    
    for i in range(len(heights)-1):
        Xv_0, Xv_1 = htogcm2(heights[i]), htogcm2(heights[i+1])
        Xs += (Xv_0-Xv_1)/cos_localtheta(heights[i+1], theta)
        
    return Xs


def ZHAireS_file_mixer(file, rootdir, domain = 'f'):
    ''' Combina archivos de output de ZHAireS en uno unico, en el caso de que 
        hayan hecho falta varias runs con varios conjuntos de antenas
        
        files: lista de archivos que se quieren unir
        rootdir : directorio general, donde se guardara el nuevo archivo
        domain  : tipo de outut (tiempo 't', frecuencia 'f')
    '''
    
    datasets = []
    n_ant = 0
    
    for i in range(len(file)):
        try:
            data = np.loadtxt(file[i], comments = '#').T

        except ValueError:
            if rootdir == '':
                print('Error leyendo, prueba a poner rootdir')
            else:
                print('Error leyendo, estamos en ello ...')
                
                name = 'timefresnel-rootnew_'+str(i+1)+'.dat' if domain == 't' else 'freqfresnel-rootnew_'+str(i+1)+'.dat'
                fileout = rootdir + name
                
                f = open(file[i],'r')
                filedata = f.read()
                f.close()
            
                newdata = filedata.replace("-"," -")
            
                f = open(fileout,'w')
                f.write(newdata)
                f.close()
                print('Archivo corregido: ', fileout)
                
            data = pd.read_csv(fileout, comment = '#', delimiter = '\s+').to_numpy().T
        
        data[1] = data[1]+n_ant
        
        n_ant = np.max(data[1])
        
        datasets.append(data)
        
        
    new_file = 'timefresnel-root_mix.dat' if domain == 't' else 'freqfresnel-root_mix.dat'
    
    output   = rootdir + new_file 

    
    if domain == 'f':
        
        header = ' -----------Freq Fresnel file----------\n'
        header += ' Archivos combinados: \n'
        for archivo in file:
            header += ' '+archivo+'\n'
        header += ' \n'
        header += ' Columns: \n'
        header += '        1 - Shower Number\n'
        header += '        2 - Antenna Number\n'
        header += '        3 - Antenna X (m)\n'
        header += '        4 - Antenna Y (m)\n'
        header += '        5 - Antenna Z (m)\n'
        header += '        6 - Frequency #\n'
        header += '        7 - Frequency (MHz)\n'
        header += '        8 - |E| (V/M MHz)\n'
        header += '        9 -  Ex  (V/M MHz)\n'
        header += '        10 - Ey  (V/M MHz)\n'
        header += '        11 - Ez  (V/M MHz)\n'
        header += '\n'
        header += ' new shower    1\n'
        header += '\n' 
        
    elif domain == 't':
        
        header = ' -----------Time Fresnel file----------\n'
        header += ' Archivos combinados: \n'
        for archivo in file:
            header += ' '+archivo+'\n'
        header += ' \n'
        header += 'Columns: \n'
        header += '        1 - Shower Number\n'
        header += '        2 - Antenna Number\n'
        header += '        3 - Antenna X (m)\n'
        header += '        4 - Antenna Y (m)\n'
        header += '        5 - Antenna Z (m)\n'
        header += '        6 - Time (ns)\n'
        header += '        7 - |A| (V/M)\n'
        header += '        8 -  Ax  (V/M)\n'
        header += '        9 -  Ay  (V/M)\n'
        header += '        10 - Az  (V/M)\n'
        header += '        11 -|E| (V/M)\n'
        header += '        12 - Ex  (V/M)\n'
        header += '        13 - Ey  (V/M)\n'
        header += '        14 - Ez  (V/M)\n'
        header += '\n'
        header += ' new shower    1\n'
        header += '\n' 
        
    full_data = np.hstack(tuple(d for d in datasets))
    
    np.savetxt(output, full_data.T, fmt = '%.3e', header = header)
    
    return output


def Aires_Plot(graf, pcles, rootdir, const, extras, xscale='linear', yscale='linear', \
            xlims = [], ylims=[], legend = True, loc_leg = 'best', cols = 1, omitextra=False, \
            UG = False, slant = False, DistAlongAxis = False, firstin = True, \
            ang=[], inj_h = [], step = .005, export_graph = False):
    ''' 
        CADA TIPO DE GRAFICA TIENE UN NUMERO
        
        Longitudinal Development     = 0
        
        Energy Long. Development     = 1
        
        Egy depos., Long Development = 2
        
        Lateral distribution         = 3
        
        Energy distrib. at ground    = 4
        
        Energy Long. Dev (per part) = 5
        
        CADA PARTICULA TIENE UN NUMERO:

        gamma            = 0
        
        p                = 1
        
        e+e-             = 2
        
        mu+mu-           = 3
        
        pi+pi-           = 4
        
        K+K-             = 5
        
        Ot. neut. pcles  = 6
        
        All chgd. pcles. = 7
        
        All pcles.       = 8
        
        Si se quiere la grafica de Long. development para e+e-, mu+mu-, se indica:
            graf = 0
            pcles = [2, 3]
        
        Si ademas se quiere una seleccion mas estricta de archivos, se a??aden en const
        los strings por los que se quiere buscar:
            const = ['9km']
    
        Si se quieren separar los archivos escogidos en dos criterios, se a??aden en extra
            extra = ['0km', '9km']
    
        Si no se quiere separar, se indica con un string vacio:
            extra = ['']
            
        Con minotacion de archivos, para separar por particula primaria hay que
        indicarlo como _p_ por ejemplo
    
        omitextra = True 
            Si se quiere que los extras no aparezcan en las leyendas
        
        cols : columnas de la leyenda
    
        UG = True 
            Si tenemos datos de cascadas upgoing, con Xv
        
        slant = True
            Si tenemos datos slanted en el eje x
        
        DistAlongAxis = True 
            Transforma los valores de X_v en eje x a distancia
            en km a lo largo del eje. Incompatible con datos slanted.
            
        firstin: Las distancias y profundidades slanted se miden desde la primera
            interaccion (altura de inyeccion)
            
        Si tenemos DistAlongAxis ?? slant, necesitamos dar los angulos cenitales
        y alturas de inyeccion asociados a cada grafico

        ang   : list angulos cenitales en deg (relevante para slant y DistAlongAxis)
        inj_h : list Altura de inyeccion en km (relevante para slant y DistAlongAxis)
        
        
        step: paso de integracion usado en Xs_of_injh
        
        export_graph: devuelve lista de valores a graficar
            [[label, valx, valy], ...]
    '''
    ext_LD   = [1001, 1022, 1205, 1207, 1211, 1213, 1092, 1291, 1293]
    ext_ELD  = [500+ext for ext in ext_LD]
    ext_EDLD = [7801, 0000, 7905, 7907, 0000, 0000, 7892, 7991, 7993]
    ext_LatD = [1000+ext for ext in ext_LD]
    ext_EDG  = [1500+ext for ext in ext_LD]

    ext = [ ['.t'+str(numero) for numero in lista] for lista in [ext_LD, ext_ELD, ext_EDLD, ext_LatD, ext_EDG] ]


    lbl     = [r'$\gamma$', 'p', r'$e\pm$' , r'$\mu\pm$', r'$\pi\pm$', r'$K\pm$', 'Other neutral pcles.', 'All chgd. pcles.', 'All pcles.']
    xlabels = [r'$X_v$ [$g/cm^2$]', r'$X_v$ [$g/cm^2$]', r'$X_v$ [$g/cm^2$]', 'Distance to core [m]', 'E [GeV]', r'$X_v$ [$g/cm^2$]']
    ylabels = ['N', 'E [GeV]', 'E [GeV]', 'N', 'N', 'E per part. [GeV]']
    
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    
    export = []
    data   = []
    
    if graf != 5:
        
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
        for p in pcles:
            for extra in extras:
                for serie in data:
                    if p in serie and extra in serie:
                        orden.append(serie) 
                        # solo para que la leyenda salga en el orden de extras
                    
        data = orden
        
    else:
        data_n = []
        data_E = []
        for subdir, dirs, files in os.walk(rootdir):
            for file in files:
                for p in pcles:
                    if file.endswith(ext[0][p]) and all([c in file for c in const]): 
                        # si es grafica y particula adecuada, y cumple constraints (numero)
                        for extra in extras:
                            if extra == '' or extra in file:
                                data_n.append([subdir + os.sep + file, p, extra])
                    elif file.endswith(ext[1][p]) and all([c in file for c in const]): 
                        # si es grafica y particula adecuada, y cumple constraints (energia)
                        for extra in extras:
                            if extra == '' or extra in file:
                                data_E.append([subdir + os.sep + file, p, extra])
        
        if len(data_n) != len(extras)*len(pcles) or len(data_n) != len(data_E):
            raise TypeError('Faltan archivos necesarios para la grafica que se pide')
        
        orden_n, orden_E = [], []
        for p in pcles:
            for extra in extras:
                for serie in data_n:
                    if p in serie and extra in serie:
                        orden_n.append(serie) 
                for serie in data_E:
                    if p in serie and extra in serie:
                        orden_E.append(serie) 
                        # solo para que la leyenda salga en el orden de extras
                    
        data_n, data_E = orden_n, orden_E

    ang = [t*np.pi/180. for t in ang]
    
    angle_index  = 0
    h_index = 0
    
    if graf != 5:
        for file, p, extra in data:
            values, grd = readfile(file)
            
            xvalues = values[:,0]
            yvalues = values[:,1]
        
            if any([graf == i for i in [0,1,2]]): # depth in x axis
                if DistAlongAxis:
                    if slant or 'Xs' in const:
                        raise TypeError('DistAlongAxis trabaja con datos Xv originales')
                    else:
                        if firstin:
                            xvalues = [Long(inj_h[h_index], ang[angle_index]) - Long(gcm2toh(xv), ang[angle_index]) for xv in xvalues]
                            if UG:
                                xvalues = [-x for x in xvalues]
                        else: 
                            xvalues = [Long(gcm2toh(xv), ang[angle_index]) for xv in xvalues]
                        
                        angle_index  += 1
                        h_index += 1
                        
                        
                elif slant:
                    if firstin:
                        xvalues =  xvalues - Xs_of_injh(inj_h[h_index], ang[angle_index], step) if inj_h[h_index] != 0  else xvalues - grd
                        if UG:
                            xvalues = [-x for x in xvalues]
                    else:
                        xvalues = grd-xvalues if UG else xvalues
    
                    # depth traversed in upward direction, starting from first interaction
                    # consistency check
                    print('Check: GRD (Aires): %.2f  GRD (Xs_of_injh): %.2f'%(grd, Xs_of_injh(0, ang[angle_index], step)))
                    angle_index  += 1
                    h_index += 1
                    
            if extra == '' or omitextra:
                ax.step(xvalues, yvalues, where = 'mid', label = lbl[p], linewidth = 2.0)
                export.append([lbl[p], xvalues, yvalues])
                
            else:
                extra = 'Primary $'+extra[1:-1]+'$' if extra.startswith('_') and extra.endswith('_') else extra
                # esta linea es solo por si tenemos un extra particula, que se da como _p_ con mi notacion
                extra = extra[1:] if extra.startswith('_') else extra
                ax.step(xvalues, yvalues, where = 'mid', label = lbl[p]+', '+extra, linewidth = 2.0)
                export.append([lbl[p]+', '+extra, xvalues, yvalues])
            
    else: #graf = 5
        for dataset_n, dataset_E in list(zip(data_n, data_E)):
             
            file_n, p, extra = dataset_n
            file_E, _, _     = dataset_E
            
            values_n, grd = readfile(file_n)
            values_E, grd = readfile(file_E)
            
            xvalues = values_n[:,0]
            yvalues = np.array([a/b if b!=0 else 0 for a, b in list(zip(values_E[:,1], values_n[:,1]))])
        
            if DistAlongAxis:
                if slant or 'Xs' in const:
                    raise TypeError('DistAlongAxis trabaja con datos Xv originales')
                else:
                    if firstin:
                        xvalues = [Long(inj_h[h_index], ang[angle_index]) - Long(gcm2toh(xv), ang[angle_index]) for xv in xvalues]
                        if UG:
                            xvalues = [-x for x in xvalues]
                    else: 
                        xvalues = [Long(gcm2toh(xv), ang[angle_index]) for xv in xvalues]
                    
                    angle_index  += 1
                    h_index += 1
                        
                        
            elif slant:
                if firstin:
                    xvalues =  xvalues - Xs_of_injh(inj_h[h_index], ang[angle_index], step) if inj_h[h_index] != 0  else xvalues - grd
                    if UG:
                        xvalues = [-x for x in xvalues]
                else:
                    xvalues = grd-xvalues if UG else xvalues

                # depth traversed in upward direction, starting from first interaction
                # consistency check
                print('Check: GRD (Aires): %.2f  GRD (Xs_of_injh): %.2f'%(grd, Xs_of_injh(0, ang[angle_index], step)))
                angle_index  += 1
                h_index += 1
                    
            if extra == '' or omitextra:
                ax.step(xvalues, yvalues, where = 'mid', label = lbl[p], linewidth = 2.0)
                export.append([lbl[p], xvalues, yvalues])
                
            else:
                extra = 'Primary $'+extra[1:-1]+'$' if extra.startswith('_') and extra.endswith('_') else extra
                # esta linea es solo por si tenemos un extra particula, que se da como _p_ con mi notacion
                extra = extra[1:] if extra.startswith('_') else extra
                ax.step(xvalues, yvalues, where = 'mid', label = lbl[p]+', '+extra, linewidth = 2.0)
                export.append([lbl[p]+', '+extra, xvalues, yvalues])
            
    ax.set_xlabel(xlabels[graf], size = 12)
    ax.set_ylabel(ylabels[graf], size = 12) 
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
        
    if any([graf == i for i in [0,1,2,5]]): # depth in x axis
        if UG:
            if DistAlongAxis:
                if firstin:
                    ax.set_xlabel(r'Dist. along axis (upward, from first interaction) [km]')
                else:
                    ax.set_xlabel(r'Dist. along axis (upward) [km]')
            elif slant:
                if firstin:
                    ax.set_xlabel(r'$X_s$ (upward, from first interaction) [$g/cm^2$]')
                else:
                    ax.set_xlabel(r'$X_s$ (upward) [$g/cm^2$]')
            else:
                ax.invert_xaxis()
                
        else:
            if DistAlongAxis:
                if firstin:
                    ax.set_xlabel(r'Dist. along axis (from first interaction) [km]')
                else:
                    ax.set_xlabel(r'Dist. along axis [km]')
            elif slant:
                if firstin:
                    ax.set_xlabel(r'$X_s$ (from first interaction) [$g/cm^2$]')
                else:
                    ax.set_xlabel(r'$X_s$ [$g/cm^2$]')
                
            
              
    if legend:
        ax.legend(loc = loc_leg, prop={'size':12}, ncol = cols, fancybox = True)
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
    
    if export_graph:
        return export
    
def plot_shower_dev(data, thetas, max_height = 100, \
                    num_scale = 7, sep = 100):
    
    ''' Grafica del desarrollo de la cascada, (histogramas en angulo simplemente)
        
        Parametros:
            data: output de AiresPlot (export_graph = True)
                Si no es output, algo como [[leyenda, data_x, data_y], ...]
            theta: angulos de las cascadas (deg)
            max_height: limite de altura en la grafica
            num_scale: orden de magnitud en que medimos numero de particulas (ajustar a mano)
            sep: separacion entre graficas. Ajustar a mano
            
        Las entradas salen de AiresPlot (DistAlongAxis = True, firstin = False)
            
    '''
   
    plt.figure()
    plt.yticks(size = 12)
    n_show = len(data)
    
    if len(thetas) != n_show:
        raise TypeError('Faltan (o sobran) angulos para las cascadas')
    
    thetas = np.array([t*np.pi/180 for t in thetas])
    
    index = np.argsort(thetas)

    offset  = 1
    c_index = 1
    
    for i in index:
        d, N = np.array(data[i][1]), np.array(data[i][2])/10**num_scale
        
        xpos =  np.cos(-thetas[i]+np.pi/2)*d-np.sin(-thetas[i]+np.pi/2)*N
        ypos =  np.sin(-thetas[i]+np.pi/2)*d+np.cos(-thetas[i]+np.pi/2)*N
        
        xneg = np.cos(-thetas[i]+np.pi/2)*d+np.sin(-thetas[i]+np.pi/2)*N
        yneg = np.sin(-thetas[i]+np.pi/2)*d-np.cos(-thetas[i]+np.pi/2)*N
        
        plt.plot(xpos+offset, ypos, color = 'k')
        plt.plot(np.array(xneg)+offset, np.array(yneg), color = 'k')
        
        plt.fill(np.append(xpos+offset, xneg[::-1]+offset), np.append(ypos, yneg[::-1]), color = colores[c_index], label = data[i][0])
        offset += sep
        c_index += 1
        
    alturas = np.linspace(0, max_height, 1000)
    long_axis = np.linspace(plt.gca().get_xlim()[0], plt.gca().get_xlim()[1])
    for i in range(len(alturas)-1):
        plt.fill_between(long_axis, alturas[i], alturas[i+1], alpha = np.exp(-alturas[i]/8.31), color = 'mediumblue')
        
    plt.ylim(0, max_height)
    plt.ylabel('h [km]', size = 12)
    plt.xticks([])
    plt.legend(loc = 'best', fontsize = 12)


def ZHAireS_Plot_t(graphs, antenas, file, xscale='linear', yscale='linear', xlims=[], \
                   ylims = [], legend = True, formato = '-', labels = [], size = 12):
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
    
    Por ejemplo, se quiere el maximo de la se??al de E en tiempo respecto a coord y de antenas:
        5 (grafica de E) + 20 (coord y) = 25
    
    ESPECIFICACION DE GRAFICAS Y ANTENAS:
        El index de antenas empieza en 1 como en Aires
        
        Si queremos m??s de una grafica, damos: 
            graphs = [11, 13]
        Podemos especificar que antenas entran en cada grafica:
            antenas = [[1,2,3,4,5], [6,7,8,9,10]]
            
        Siempre una lista de antenas por grafica
        
        Si se quieren poner todas
            antenas = 'all'
            
    labels : lista de strings, labels que queremos explicitamente
    size   : tama??o de letra

'''
    data_time = np.loadtxt(file, comments = '#').T
    
    if len(antenas) != len(graphs):
        raise TypeError('longitud equivocada')
        
    if (not all([g>10 for g in graphs])) and (not all([g<10 for g in graphs])):
        raise TypeError('Todas las graficas deben ser del mismo tipo')
    
    plot_maxs_vs_coords = all([g>10 for g in graphs])
    
    n_ant               = int(max(data_time[1])) #numero de antenas
    
    ant_coord           = [] # aqui guardamos las coordenadas de las antenas
    
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
            coord = g // 10 - 1 
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
    ax  = fig.add_subplot(111)
    plt.xticks(fontsize = size)
    plt.yticks(fontsize = size)
    
    magnitudes = [r'$A$ (V/m)', r'$A_x$ (V/m)', r'$A_y$ (V/m)', r'$A_z$ (V/m)', r'$E$ (V/m)', r'$E_x$ (V/m)', r'$E_y$ (V/m)', r'$E_z$ (V/m)']
    lbl_coords = ['N-S', 'E-W', 'Vertical']
    autolabel = False
    
    if plot_maxs_vs_coords:
        if all([((g%10-1)-(graphs[0]%10-1))==0 for g in graphs]):
            if len(labels)<1:
                labels = [lbl_coords[g//10-1] for g in graphs]
            ax.set_ylabel(magnitudes[graphs[0]%10-1], size = size)
        else:
            if len(labels)<1:
                labels =  [magnitudes[g%10-1]+', '+lbl_coords[g//10-1] for g in graphs]
            
        for i in range(len(graphs)):#g, data in list(zip(graphs, graph_list)):
            ax.plot(graph_list[i][0], graph_list[i][1], formato, label = labels[i])
            
        ax.set_xlabel('Dist. a core (m)', size = size)
     
    else:
        if all([(g-graphs[0])==0 for g in graphs]):
            if len(labels)<1:
                autolabel = True
                labels = ['Antena ' for _ in range(len(graphs))]
            ax.set_ylabel(magnitudes[graphs[0]-1], size = size)
            
        else:
            if len(labels)<1:
                autolabel = True
                labels =  [magnitudes[g-1]+', Antena ' for g in graphs]
        index = 0
        for i in range(len(graphs)):#g, data, a in list(zip(graphs, graph_list, antenas)):
            for j in range(len(antenas[i])):
                label = labels[i]+'%d'%antenas[i][j] if autolabel else labels[index]
                ax.plot(graph_list[index][0], graph_list[index][1], formato, label = label)
                index += 1
        ax.set_xlabel('t (ns)')
        
        
    if legend:
        ax.legend(loc='best', prop={'size':size})
    
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
        
    return fig

def ZHAireS_Plot_f(graphs, antenas, freq, file, xscale='linear', yscale='linear', xlims=[], \
                   ylims = [], legend = True, formato = '-', labels = [], size = 12):
    ''' Codigos de graficas en dominio de frecuencias

    E vs. f, antenas especificadas  = 2
    
    Ex vs. f, antenas especificadas = 3
    
    Ey vs. f, antenas especificadas = 4
    
    Ez vs. f, antenas especificadas = 5
    
    Para representacion de la magnitud representada en graf. i, a frecuencia concreta f:
        
    Frente a coord x = i+10
    
    Frente a coord y = i+20
    
    Frente a coord z = i+30
    
    Por ejemplo, se quiere la se??al de E a 300MHz respecto a coord y de antenas:
        2 (grafica de E) + 20 (coord y) = 22
    
    ESPECIFICACION DE GRAFICAS Y ANTENAS:
        El index de antenas empieza en 1 como en Aires
        
        Si queremos m??s de una grafica, damos: 
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
                 
    labels : lista de strings, labels que queremos explicitamente
    size   : tama??o de letra
'''
    data_time = np.loadtxt(file, comments = '#').T
    
    if len(antenas) != len(graphs):
        raise TypeError('longitud equivocada de antenas')
        
    if (not all([g>10 for g in graphs])) and (not all([g<10 for g in graphs])):
        raise TypeError('Todas las graficas deben ser del mismo tipo')
    
    plot_maxs_vs_coords = all([g>10 for g in graphs])
    
    if plot_maxs_vs_coords and len(freq)!=len(graphs):
        raise TypeError('Indica la frecuencia en TODAS las graficas que quieras')
    
    n_ant     = int(max(data_time[1])) #numero de antenas
    
    ant_coord = [] # aqui guardamos las coordenadas de las antenas
    
    i_ant     = []
    
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
            coord = (g // 10) -1
            g     = g % 10
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
    ax  = fig.add_subplot(111)
    plt.xticks(fontsize = size)
    plt.yticks(fontsize = size)
    magnitudes = [r'$E$ (V/m/MHz)', r'$E_x$ (V/m/MHz)', r'$E_y$ (V/m/MHz)', r'$E_z$ (V/m/MHz)']
    lbl_coords = ['N-S', 'E-W', 'Vertical']
    autolabel = False
    
    if plot_maxs_vs_coords:
        if all([((g%10-1)-(graphs[0]%10-1))==0 for g in graphs]):
            if len(labels) < 1:
                labels = [lbl_coords[graphs[i]//10-1]+', '+str(freq_list[f])+' MHz'  for i in range(len(graphs)) for f in freq[i]]
            
            ax.set_ylabel(magnitudes[graphs[0]%10-2], size = size)
        else:
            if len(labels) < 1:
                labels =  [magnitudes[graphs[i]%10-2]+', '+lbl_coords[graphs[i]//10-1]+', '+str(freq_list[f])+' MHz'  for i in range(len(graphs)) for f in freq[i]]
            
        for i in range(len(graph_list)):
            ax.plot(graph_list[i][0], graph_list[i][1], formato, label = labels[i])
            
        ax.set_xlabel('Dist. a core (m)', size = size)
     
    else:
        if all([(g-graphs[0])==0 for g in graphs]):
            if len(labels) < 1:
                autolabel = True
                labels = ['Antena ' for _ in range(len(graphs))]
            ax.set_ylabel(magnitudes[graphs[0]-2], size = size)
            
        else:
            if len(labels)<1:
                autolabel = True
                labels =  [magnitudes[g-1]+', Antena ' for g in graphs]
            
        for i in range(len(graphs)):#g, data, a in list(zip(graphs, graph_list, antenas)):
            for j in range(len(antenas[i])):
                index = i*len(antenas[i])+j
                label = labels[i]+'%d'%antenas[i][j] if autolabel else labels[index]
                ax.plot(graph_list[index][0], graph_list[index][1], formato, label = label)
        ax.set_xlabel('Freq (MHz)')
        
        
    if legend:
        ax.legend(loc='best', prop={'size':size})
    
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
        
    return fig


def FFT(file, graphs, antenas, xscale='log', yscale='log', xlims=[], \
        ylims = [], legend = True, formato = '-', labels = [], size = 12):

    ''' Tranformada de Fourier rapida de output de ZHAireS en tiempo
    
    Codigos de graficas en dominio temporal

    A vs. t, antenas especificadas  = 1
    
    Ax vs. t, antenas especificadas = 2
    
    Ay vs. t, antenas especificadas = 3
    
    Az vs. t, antenas especificadas = 4
    
    E vs. t, antenas especificadas  = 5
    
    Ex vs. t, antenas especificadas = 6
    
    Ey vs. t, antenas especificadas = 7
    
    Ez vs. t, antenas especificadas = 8
    
    
    ESPECIFICACION DE GRAFICAS Y ANTENAS:
        El index de antenas empieza en 1 como en Aires
        
        Si queremos m??s de una grafica, damos: 
            graphs = [11, 13]
        Podemos especificar que antenas entran en cada grafica:
            antenas = [[1,2,3,4,5], [6,7,8,9,10]]
            
        Siempre una lista de antenas por grafica
        
        Si se quieren poner todas
            antenas = 'all'
            
    labels: lista de str, labels explicitos para cada grafico  
    size  : tama??o de letra     
    '''
    
    data_time = np.loadtxt(file, comments = '#').T
    
    if len(antenas) != len(graphs):
        raise TypeError('longitud equivocada')
        
    n_ant               = int(max(data_time[1])) #numero de antenas
    
    ant_coord           = [] # aqui guardamos las coordenadas de las antenas
    
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
    
    for i in range(len(graphs)):
        for a in antenas[i]:
            x = data_time[5][i_ant[a-1]:i_ant[a]]
            y = data_time[5+graphs[i]][i_ant[a-1]:i_ant[a]]
            
            N = len(y) # number of signal points
            
            dT = (x[1]-x[0])*1e-9 #time bin in ns
            
            xf = np.abs(fftfreq(N, dT)*1e-6) #frequencies in MHz
            
            yf = np.abs(fft(y))/N # FFT of data
            
            graph_list.append([xf, yf])
            
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    plt.xticks(fontsize = size)
    plt.yticks(fontsize = size)
    
    magnitudes = [r'$A$ (V/m/MHz)', r'$A_x$ (V/m/MHz)', r'$A_y$ (V/m/MHz)', r'$A_z$ (V/m/MHz)', r'$E$ (V/m/MHz)', r'$E_x$ (V/m/MHz)', r'$E_y$ (V/m/MHz)', r'$E_z$ (V/m/MHz)']

    autolabel = False
    if all([(g-graphs[0])==0 for g in graphs]):
        if len(labels) < 1:
            autolabel = True
            labels = ['Antena ' for _ in range(len(graphs))]
        ax.set_ylabel(magnitudes[graphs[0]-1], size = size)
        
    else:
        if len(labels)<1:
            autolabel = True
            labels =  [magnitudes[g-1]+', Antena ' for g in graphs]
        
    index = 0
    for i in range(len(graphs)):#g, data, a in list(zip(graphs, graph_list, antenas)):
        for j in range(len(antenas[i])):
            label = labels[i]+'%d'%antenas[i][j] if autolabel else labels[index]
            ax.plot(graph_list[index][0], graph_list[index][1], formato, linewidth = 2., label = label)
            index += 1
    ax.set_xlabel('f (MHz)', size = size)
        
        
    if legend:
        ax.legend(loc='best', prop={'size':size})
    
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    
    if len(xlims)>0:
        ax.set_xlim(xlims)
    if len(ylims)>0:
        ax.set_ylim(ylims)
        
    plt.tight_layout()
    return fig


def ZHAireS_arrayplot(graph, file, domain = 'f', lims = [],  d3 = False, rootdir='',\
                      export = True, xoffset=0, yoffset=0, height = '', draw_ant = False, aspect = 'auto'):
    ''' Grafica de campo electrico en un array 2d de antenas
        En dominio t, se representa el maximo alcanzado para la variable 
        pedida en cada antena
        
        C??digos de graficas en dominio de frecuencia
        
        E  = 1
        
        Ex = 2
        
        Ey = 3
        
        Ez = 4
        
        Codigos de graficas en dominio temporal
        
        A  = 1
        
        Ax = 2
        
        Ay = 3
        
        Az = 4
        
        E  = 5
        
        Ex = 6
        
        Ey = 7
        
        Ez = 8
        
        ----------- INPUTS -------------
        
        graph: numero que indica la grafica que queremos, segun lo anterior
        
        file: archivo de datos. Output estandar de ZHAireS o ZHAireS_file_mixer
        
        domain: datos en dominio de frecuencia ('f') o tiempo ('t')
            Si se escoge 'f', debe introducirse manualmente la frecuencia deseada
            
        lims: lista, limites para la escala de colores (eje z)
        
        d3: bool, plot 2d con mapa de calor o superficie 3d
        
        rootdir: directorio con archivos, necesario en caso de error leyendo
        
        box_size: tama??o de cuadraditos, ajustar si sale mal
        
        export: True si devuelve datos np.array([xs, ys, zs])
        
        xoffset, yoffset: para centrar la grafica, coordenadas del "centro" de cascada [km]
        
        height: altura a mano
        
        draw_ant: indicar explicitamente posicion de antenas en mapa 2d
        
        aspect: proporcion de los ejes, defalut 'auto' (normal en matplotlib)
            poniendo 'equal' las escalas en ejes x e y se hacen iguales
    '''
    try:
        data = np.loadtxt(file, comments = '#').T

    except ValueError:
        if rootdir == '':
            print('Error leyendo, prueba a poner rootdir')
        else:
            print('Error leyendo, estamos en ello ...')
            fileout = rootdir+'timefresnel-rootnew.dat' if domain == 't' else rootdir+'freqfresnel-rootnew.dat'
            f = open(file,'r')
            filedata = f.read()
            f.close()
        
            newdata = filedata.replace("-"," -")
        
            f = open(fileout,'w')
            f.write(newdata)
            f.close()
            print('Archivo corregido: ', fileout)
            
        data = pd.read_csv(fileout, comment = '#', delimiter = '\s+').to_numpy().T
    
    n_ant     = int(max(data[1])) #numero de antenas
    
    ant_coord = [] # aqui guardamos las coordenadas de las antenas
    
    i_ant     = []
    
    for i in range(n_ant):
        index = list(data[1]).index(i+1)
        i_ant.append(index)
        ant_coord.append([data[2][index], data[3][index], data[4][index]])
        # con esto ya tenemos las coordenadas de las antenas
    
    i_ant.append(len(data.T))


    height = round(100 - ant_coord[0][2]/1000) if height == '' else height  #altura
    
    if domain == 'f':
        labels = [r'$E$ (V/m/MHz)', r'$E_x$ (V/m/MHz)', r'$E_y$ (V/m/MHz)', r'$E_z$ (V/m/MHz)']
        freq_list = [data[6][i] for i in range(i_ant[0], i_ant[1])]
        
        print('Frecuencias (MHz): ', freq_list)
        
        f = int(input('Escoge frecuencia (indice empezando en 0): '))
        
        g = graph+6
        
        values = np.array([data[int(g)][i_ant[a]+f] for a in range(n_ant)]).reshape(n_ant,1)
    else:
        
        labels = [r'$A$ (V/m)', r'$A_x$ (V/m)', r'$A_y$ (V/m)', r'$A_z$ (V/m)', r'$E$ (V/m)', r'$E_x$ (V/m)', r'$E_y$ (V/m)', r'$E_z$ (V/m)']
    
        g = graph+5
        
        values = np.array([np.max(data[int(g)][i_ant[a]:i_ant[a+1]]) for a in range(n_ant)]).reshape(n_ant,1)
        
    puntos = np.hstack((np.array(ant_coord)[:,:2], values))
    
    x, y, z = puntos.T[0]/1000-xoffset, puntos.T[1]/1000-yoffset, puntos.T[2]
    
    fig = plt.figure(figsize=(8,6))
    
    if d3:
        ax = fig.add_subplot(111, projection = '3d')
        ax.plot_trisurf(x,y,z, cmap = 'gnuplot', alpha = .9)
        ax.set_xlabel('x [km]', size = 14)
        ax.set_ylabel('y [km]', size = 14)
        ax.set_zlabel(labels[graph-1], size = 14)
        
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.zaxis.set_tick_params(labelsize=12)
        
        title = 'h = '+ str(height)+' km, '+ str(freq_list[f])+' MHz' if domain =='f' else 'h = '+ str(height)+' km'
        ax.set_title(title, size = 14)
        
        if len(lims) > 0:
            ax.set_zlim(lims)
            
    else:
        
        df = pd.DataFrame.from_dict(np.array([x, y, z]).T)
        
        df.columns = ['x [km]', 'y [km]', labels[graph -1]]
        
        pivotted = df.pivot('y [km]','x [km]',labels[graph -1]).to_numpy()
        
        ny, nx = pivotted.shape
        
        x_min, x_max, y_min, y_max =  min(x), max(x), min(y), max(y)
        
        deltax, deltay = (x_max-x_min)/(nx), (y_max-y_min)/(ny)
        
        extent = [x_min-deltax/2, x_max+deltax/2, y_min-deltay/2, y_max+deltay/2]

        ax = fig.add_subplot(111)
        
        if len(lims) > 0:
            sc = ax.imshow(pivotted , origin = 'lower', cmap = 'gnuplot', aspect = 'auto', extent = extent, vmin = lims[0], vmax = lims[1])
        else:
            sc = ax.imshow(pivotted , origin = 'lower', cmap = 'gnuplot', aspect = 'auto', extent = extent)
            
        cbar = fig.colorbar(sc)
        cbar.set_label(labels[graph-1], size = 14, rotation = 270, labelpad = 30)
        cbar.ax.tick_params(labelsize=12)
        
        ax.set_xlabel('x [km]', size = 14)
        ax.set_ylabel('y [km]', size = 14)
        
        ax.xaxis.set_tick_params(labelsize=12)
        ax.yaxis.set_tick_params(labelsize=12)
        
        title = 'h = '+ str(height)+' km, '+ str(freq_list[f])+' MHz' if domain =='f' else 'h = '+ str(height)+' km'
        ax.set_title(title, size = 14)
        
    plt.tight_layout()
    plt.axis(aspect)
    
    if draw_ant:
        ax.scatter(x, y, marker = '1', color = 'k')
    
    if export:
        return np.array([x, y, z])
    
    
################################# EJEMPLOS Aires_Plot ##################################
# ---------------------------------------------------------------------------------------

# rootdir = 'Carac_UG'
# graf  = 0
# pcles = [2]
# const = ['1EeV', '5km', 'Xv']
# extras = ['95deg', '100deg', '105deg', '110deg', '115deg', '120deg', '125deg', '130deg']
# xlim=[0, 500]
# #como voy a usar DistAlongAxis, necesito dar angulos y inj_h en el orden que se van haciendo las graficas

# thetas = [(85-5*i)*np.pi/180 for i in range(10)] # tengo 10 graficas a 95 grados
# inj_h = [5 for _ in range(10)]

# maz.Aires_Plot(graf, pcles, rootdir, const, extras, yscale='log', xlims = xlim, cols=2, \
#                UG = True, DistAlongAxis=True, ang=thetas, inj_h = inj_h)

# -------------------------------------------------------------------------------

# rootdir = 'Carac_UG'
# graf  = 0
# pcles = [2]
# const = ['1EeV', '5km', 'Xs']
# extras = ['95deg', '100deg', '105deg', '110deg', '115deg', '120deg', '125deg', '130deg']
# xlim=[0, 500]
# #como voy a usar slant, necesito dar angulos en el orden que se van haciendo las graficas

# thetas = [(85-5*i)*np.pi/180 for i in range(10)] # tengo 10 graficas a 95 grados
# inj_h = [5 for _ in range(10)]

# maz.Aires_Plot(graf, pcles, rootdir, const, extras, yscale='log', xlims = xlim, cols=2, 
#               UG = True, slant=True, ang=thetas, inj_h = inj_h)

# -------------------------------------------------------------------------------

################################## Ejemplos ZHAireS_plot_t #######################################

# rootdir = 'Graph_paperANITA'

# graphs = [1, 8]

# antenas = [[i+1 for i in range(10, 15)], [32]]

# file = rootdir+'/timefresnel-only300MHz.dat'

# ZHAireS_Plot_t(graphs, antenas, file)

################################## Ejemplos ZHAireS_plot_f ###################################

# rootdir = 'Graph_paperANITA'
# file = rootdir+'/freqfresnel-variousfreqs.dat'

# graphs = [24]

# antenas = [[i+1 for i in range(50, 100)]]

# freq = [[0, 1, 2, 4]]

# ZHAireS_Plot_f(graphs, antenas, freq, file)