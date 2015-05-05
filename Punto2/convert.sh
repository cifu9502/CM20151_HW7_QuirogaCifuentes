#pule el archivo times.csv y lo guarda en times2.csv para que pueda ser leido 
#por dates.R
cat times.csv | sed 's\_TAI\\' | sed 's/\_/ /g' | sed s/\\./-/g > times2.csv


