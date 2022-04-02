Carac_UG:

Caracterizacion de cascadas atmosféricas hacia arriba. Lista de resultados:

	Cascadas con primario p - primera caract:
	
	Primary p, E = 1EeV,         h = 5km,      theta = (95-130 deg)  : desarrollo de electrones y muones vs X_v; dist y X_s (desde primera interaccion) 
	Primary p, E = 1EeV,         h = (0, 9km), theta = 95deg         : desarrollo de electrones y muones vs X_v; dist y X_s (desde grd)
	Primary p, E = (1PeV-10EeV), h = 5km,      theta = 95deg         : desarrollo de electrones y muones vs X_v, dist y X_s (desde grd)
	
	Cascadas con primario p-e, 0km:

	Primary p, E = 1EeV, h = 0km, theta = (95-130deg)   : desarrollo de electrones, muones, piones, all chged vs Xv, dist y X_s (desde primera interaccion)
	Primary e, E = 1EeV, h = 0km, theta = (95-130deg)   : desarrollo de electrones, muones, piones, all chged vs Xv, dist y X_s (desde primera interaccion)

	Primary (p,e), E = 1EeV, h = 0km, theta = (100,120deg)  : desarrollo de electrones, muones, piones vs dist y X_s (desde primera interaccion)
	Primary (p,e), E = 1EeV, h = 0km, theta = (100,120deg)  : energia depositada por electrones, muones, all chgd vs dist y X_s (desde primera interaccion)

Comp_UGDG:

Comparacion del desarrollo de cascadas atmosféricas hacia arriba y hacia abajo, para mismo ángulo. Lista de resultados:
	
	Primary p, E = 1EeV, h = 0 (100) km, theta = 130 (50) deg  : desarrollo de electrones, muones y piones vs X_s atravesada (escala lin y log)
	Primary p, E = 1EeV, h = 0 (100) km, theta = 180 (0)  deg  : desarrollo de electrones, muones y piones vs X_s atravesada (escala lin y log)

	Primary p, E = 1EeV, h = 0 (100) km, theta = 130 (50) deg  : energía depositada por muones vs X_s atravesada (escala lin y log)
	Primary p, E = 1EeV, h = 0 (100) km, theta = 180 (0)  deg  : energía depositada por muones vs X_s atravesada (escala lin y log)

	Primary p, E = 1EeV, h = 0 (100) km, theta = 130 (50) deg  : desarrollo de electrones, muones y piones en la atmósfera (rep. gráfica)	

Codigos:

Programas en Python para tratamiento de datos, graficas y geometria de cascadas:

	calcring, AntennaPlotter, ParticlePlotter: programas iniciales para calculos del cono Cherenkov y graficas. Integrados en modulo_aires_zhaires

	modulo_aires_zhaires: Incorpora las siguientes funciones:
		gcm2toh, htogcm2 : Conversion de profundidad vertical a altura y viceversa
		Long             : Conversion de alturas a distancia hasta el suelo, medida a lo largo del eje de la cascada
		CherenkovRing    : Parametros de la elipse que forma el cono Cherenkov a nivel del suelo
		readfile         : Lectura de outputs de Aires
		cos_localtheta   : Coseno del angulo cenital local a una altura h de una cascada a angulo theta
		Xs_of_injh       : Profundidad slanted de la altura de inyeccion, a lo largo del eje de la cascada (calculo aproximado)
		Aires_Plot       : Graficas de outputs de Aires
		plot_shower_dev  : Representacion del desarrollo longitudinal de cascadas. Utiliza outputs de AiresPlot
		ZHAireS_Plot_t   : Graficas de outputs de ZHAireS en dominio temporal
		ZHAireS_Plot_f   : Graficas de outputs de ZHAireS en dominio de frecuencias

Documento: 

Archivo .TeX, .pdf y .bib del TFM, y carpeta con graficos y diagramas
		


