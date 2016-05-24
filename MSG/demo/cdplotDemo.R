# original data
fail = factor(c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1,2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1), levels = 1:2,labels = c("no", "yes"))
temperature = c(53, 57, 58, 63, 66, 67, 67, 67, 68,69, 70, 70, 70, 70, 72, 73, 75, 75, 76, 76, 78, 79, 81)
par(mar=c(4,4,.5,2))
cdplot(fail~temperature, col = c("lightblue", "red"))
points(temperature, c(0.25, 0.75)[as.integer(fail)],col = "blue", bg = "yellow", pch = 21)
