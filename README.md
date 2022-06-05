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

Comparacion del desarrollo de cascadas atmosféricas hacia arriba y hacia abajo, para mismo ángulo:

	*(Cascada extraña: Resultados con posible fluctuacion por semilla del random)*
		
	10cascadas: Gráficos de desarrollos promediando las simulaciones de 10 cascadas. Resultados:

	Primary p, E = 1EeV, h = 0 (100) km, theta = 0 (180) deg : desarrollo de e, mu, pi cargados, k cargados y otras particulas neutras vs. Xs
	
	Primary p, E = 1EeV, h = 0 (100) km, theta = 0 (180) deg : comparativas de desarrollo hacia arriba y abajo para e, mu, pi, k, y demas particulas neutras vs Xs

	Ídem para 50-130 deg y 85-95 deg; y para desarrollo de energia y desarrollo de energia promedio

	ComparativaPrimarios: Desarrollos longitudinales hacia arriba y abajo, para primario proton o electron. Comparativas para e, mu, pi, demas particulas neutras vs Xs
	

	* los archivos que comienzan por Desarrollo... son una representacion grafica del desarrollo en la atmosfera de cada componente de la cascada. 

Radio:

Simulaciones de la emision en radio, para una cascada vertical hacia abajo. Señal en tiempo y transformada de Fourier


Radio_UG:

Caracterizacion de la emision en radio en cascadas hacia arriba:
	
	Run1_prueba: Primary p, E = 1EeV, theta = 0deg, B = 50uT, I = 0deg, grd = 0km, invertphi = false

	Vertical   : Resultados para cascadas verticales (20, 36 km)
		     Primary p, E = 1EeV, theta = 0deg, B = 50uT, I = 0deg, grd = 0km, invertphi = false (usando uprimary)
	45deg      : Resultado para cascada inclinada a 45deg, 36km 
		     Primary p, E = 1EeV, theta = 45deg, B = 50uT, I = 0deg, grd = 0km, invertphi = false (usando uprimary)

Codigos:

Programas en Python para tratamiento de datos, graficas y geometria de cascadas:

	calcring, AntennaPlotter, ParticlePlotter: programas iniciales para calculos del cono Cherenkov y graficas. Integrados en modulo_aires_zhaires

	modulo_aires_zhaires: Incorpora las siguientes funciones:
		gcm2toh, htogcm2   : Conversion de profundidad vertical a altura y viceversa
		Long               : Conversion de alturas a distancia hasta el suelo, medida a lo largo del eje de la cascada
		CherenkovRing      : Parametros de la elipse que forma el cono Cherenkov a nivel del suelo
		readfile           : Lectura de outputs de Aires
		cos_localtheta     : Coseno del angulo cenital local a una altura h de una cascada a angulo theta
		Xs_of_injh         : Profundidad slanted de la altura de inyeccion, a lo largo del eje de la cascada (calculo aproximado)
		ZHAireS_file_mixer : Combina outputs de ZHAireS en unico archivo, cuando son necesarias mas de una run (>200 antenas)
		Aires_Plot         : Graficas de outputs de Aires
		plot_shower_dev    : Representacion del desarrollo longitudinal de cascadas. Utiliza outputs de AiresPlot
		ZHAireS_Plot_t     : Graficas de outputs de ZHAireS en dominio temporal
		ZHAireS_Plot_f     : Graficas de outputs de ZHAireS en dominio de frecuencias
		FFT                : Transformada de Fourier rapida de outputs en tiempo de ZHAireS.
		ZHAireS_arrayplot  : Graficas de campo electrico (output ZHAireS) en un array 2d de antenas a cierta altura, en dominio t o f

Documento: 

Archivo .TeX, .pdf y .bib del TFM, y carpeta con graficos y diagramas
		


