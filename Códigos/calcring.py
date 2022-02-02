import numpy as np

''' CALCULO DE ANILLO CERENKOV A NIVEL DEL SUELO

Inputs y parametros:
----------------------
UseSlant=True => el programa usa xmaxslant para los cálculos.  (false => usa xmaxvert)

xmaxslant=700;    // only used if UseSlant is true
xmaxvert= 304.0;  // Xmax in vertical g/cm2 -> Solo si UseSlant=False
groundalt=1400.;  // ground altitude in meters
thetadeg=70.;     // zenith angle in degrees
 
Modelo de índice de refraccion -> ZHAireS default

Output:
----------------------
Inputs $$$$$$$$$$$$$$$$$
XmaxSlant=700 g/cm2
XmaxVert=239.414 g/cm2 Theta=70 deg
altxmax=10806.7 m, hxmax=9406.66 m, nh=1.00009, thetac=0.756386 deg

Approximate calculation $$$$$$$$$$$$$$$$$$$$$$$$
r1=1101.6 m r2=363.103 m  ---> 
eje mayor y menor de la elipse 
(asume centro de la elipse es el core -> buena aproximacion para angulos pequeños)

Exact calculation $$$$$$$$$$$$$$$$$$$$$$$$
a1=1024.48 m  a2=1101.6 m  a=1063.04 m --> Eje mayor (a1 -> lado "early", a2 -> lado "late", a -> media) 
(Esta asimetria aparece porque el core no está en el centro de la elipse, esta mas proximo al lado "early")
b=363.342 m                            --> Minor axis
Ecc.=38.5593                        --> Ellipse eccentricity 
D=27503.2 m                            --> Distancia del core a Xmax a lo largo del eje
4.5*r2=1633.96 m                       --> 
Estimacion de distancia hasta 4.5x el anillo cherenkov a lo largo de la direccion del eje menor
'''

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
    
    

    
    
    
