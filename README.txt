Primera version: 2/2/2022

Códigos Python:
	calcring2.py -> Parámetros de la elipse que fomra el cono Cerenkov a nivel del suelo. Incluye las funciones gcm2toh y htogcm2
	AntennaPlotter, ParticlePlotter -> Códigos para hacer gráficas de manera rápida

Caract_UG:
	Caracterizacion de cascadas hacia arriba. Gráficas de desarrollos longitudinales variando altura, ángulo y energía del primario
	Por ahora graficas frente a Xv, para electrones y muones

3/2/2022:
	Añadidas gráficas de desarrollos longitudinales de e, mu, en función de X_slanted, a lo largo del eje, para caracterizar cascadas UG
	Nuevo código -> modulo_aires_zhaires , incluye funciones para convertir gcm2 a h y vversa, funciones para plots, calculos del cono cherenkov

5/2/2022:
	Añadidas graficas de desarrollos longitudinales de e, mu en funcion de la distancia a lo largo del eje (km) en Caract_UG
	Nueva funcion en modulo_aires_zhaires: Long(h, theta), para la altura h de un pto. de la cascada a angulo theta calcula la 
		longitud medida a lo largo del eje