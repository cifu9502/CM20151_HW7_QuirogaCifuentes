cat times.csv | sed 's\_TAI\\' | sed 's/\_/ /g' | sed s/\\./-/g > times2.csv


