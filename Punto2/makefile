all:  times2.csv tiempos.csv lineal.dat exponencial.dat tangente.dat

times2.csv: convert.sh 
	chmod u+x convert.sh&
	./convert.sh

tiempos.csv: times2.csv  dates.R
	Rscript dates.R 


lineal.dat: 2alin.py  tiempos.csv
	python 2alin.py

exponencial.dat: 2aexp.py  tiempos.csv
	python 2aexp.py

tangente.dat: 2atan.py  tiempos.csv
	python 2atan.py


