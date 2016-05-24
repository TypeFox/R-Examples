#fivesec <- read.csv(file="data/fivesec.txt", head=TRUE, sep=",")
#save(fivesec, file="data/fivesec.rda")

data(fivesec)
fivesec$x = filter(fivesec$x, rep(1/3, 3))
fivesec$y = filter(fivesec$y, rep(1/3, 3))

# data is in seconds, 200Hz
# dispersion, e.g. 2 cm at a viewing distance of 80 cm is ~1.4 deg
# minimal fixation duration 0.1 s is 20 samples
fixations = emov.idt(fivesec$time, fivesec$x, fivesec$y, 2, 20)

plot(fivesec$x, fivesec$y, xlim=c(-20,30), ylim=c(-12,-6), type="l", 
     xlab="Horizontal position (deg) ", 
     ylab="Vertical position (deg)")
par(new=T)
lines(fixations$x, fixations$y, type="b", col="red")

