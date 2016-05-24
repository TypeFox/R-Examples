Rain <- list()
Rain$rain <- as.matrix(read.csv("Rain_rain.csv", header=F)) + 1
Rain$months <- as.matrix(read.csv("Rain_months.csv", header=F))
save(Rain, file="Rain.RData")
