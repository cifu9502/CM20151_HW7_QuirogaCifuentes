all:  times2.csv tiempos.csv lineal.dat exponencial.dat tangente.dat  bestmodels.txt

#Primero se ve si existe el archivo times2.csv, si no, se aplica convert.sh para pulir el archivo times.csv  y guardarlo en times2.csv
times2.csv: convert.sh 
	chmod u+x convert.sh&
	./convert.sh

#Se leen las fechas dates.R aplicando times2.csv, lo cual produce el documento tiempos.csv
tiempos.csv: times2.csv  dates.R
	Rscript dates.R 

#Usa tiempos.csv y los datos dados para hallar la mejor aproximacion lineal
lineal.dat: 2alin.py  tiempos.csv
	python 2alin.py
#Usa tiempos.csv y los datos dados para hallar la mejor aproximacion lineal con gaussiana
exponencial.dat: 2aexp.py  tiempos.csv
	python 2aexp.py

#Usa tiempos.csv y los datos dados para hallar la mejor aproximacion lineal con tangente
tangente.dat: 2atan.py  tiempos.csv
	python 2atan.py

#Lee los archivos creados anteriormente y produce las graficas de los mejores fits al igual que el archivo bestmodels.txt con los mejores modelos
bestmodels.txt : tiempos.csv 2atan.py 2alin.py 2aexp.py 2b.py 
	rm *.png
	python 2b.py
