Carac_UG:

Caracterizacion de cascadas atmosféricas hacia arriba. Lista de resultados:
	
	Primary p, 1EeV, inyeccion a 5km, variando angulos (95-130 deg) : desarrollo de electrones y muones vs X_v; dist y X_s (desde primera interaccion) 
	Primary p, 1EeV, 95deg, variando altura de inyeccion(0, 9km)    : desarrollo de electrones y muones vs X_v; dist y X_s (desde grd)
	Primary p, 9deg, inyeccion a 5km, variando energia (1PeV-10EeV) : desarrollo de electrones y muones vs X_v, dist y X_s (desde grd)

	Primary p, 1EeV, inyeccion a 0km, variando angulos (95-130deg)  : desarrollo de electrones y muones vs dist y X_s (desde primera interaccion)
	Primary e, 1EeV, inyeccion a 0km, variando angulos (95-130deg)  : desarrollo de electrones y muones vs dist y X_s (desde primera interaccion)

	Comparativa de desarrollos de electrones y muones vs dist y X_s (desde primera interaccion) segun primario p, e. E=1EeV, h=0km para 100 y 120 deg

Codigos:

Programas en Python para tratamiento de datos, graficas y geometria de cascadas:

	calcring, AntennaPlotter, ParticlePlotter: programas iniciales para calculos del cono Cherenkov y graficas. Integrados en modulo_aires_zhaires

	modulo_aires_zhaires: Incorpora las siguientes funciones:
		gcm2toh, htogcm2 : Conversion de profundidad vertical a altura y viceversa
		Long             : Conversion de alturas a distancia hasta el suelo, medida a lo largo del eje de la cascada
		CherenkovRing    : Parametros de la elipse que fomra el cono Cherenkov a nivel del suelo
		readfile         : Lectura de outputs de Aires
		cos_localtheta   : Coseno del angulo cenital local a una altura h de una cascada a angulo theta
		Xs_of_injh       : Profundidad slanted de la altura de inyeccion, a lo largo del eje de la cascada (calculo aproximado)
		Aires_Plot       : Graficas de outputs de Aires
		ZHAireS_Plot_t   : Graficas de outputs de ZHAireS en dominio temporal
		ZHAireS_Plot_f   : Graficas de outputs de ZHAireS en dominio de frecuencias

Documento: 

Archivo .TeX, .pdf y .bib del TFM
		


