library('lubridate')
destfile <- "times.csv"
dates <- read.csv(destfile, stringsAsFactors = FALSE)
date <- as.POSIXct(dates[,2], tz = "GMT")
options(digits = 12)
date <- decimal_date(date) 
write.csv(date, file = "tiempos.csv")
