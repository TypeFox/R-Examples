# A. Screen devices
dev.new(); plot(1:3)
windows(width = 6, height = 6); plot(10:15)
windows(width = 5, height = 5, bg = 'pink'); plot(rnorm(10))
dev.list()  # three screen devices
dev.cur()
dev.size()

# B. Four file devices; different file sizes
setwd("C:/aErer")
pdf(file = "test.pdf", width = 6, height = 4)
set.seed(5); plot(rnorm(1000))
dev.list()  # three screen devices + one file device
dev.off()

win.metafile(file = "test.wmf", width = 6, height = 4)
set.seed(5); plot(rnorm(1000))
dev.off()

png(file = "test.png", width = 6, height = 4, units = "in", res = 600)
set.seed(5); plot(rnorm(1000))
dev.off()

jpeg(file = "test.jpeg", width = 6, height = 4, units = "in", res=600)
set.seed(5); plot(rnorm(1000))
dev.off()
graphics.off()  # close all graphics devices

# Comparison: sizes of the graphs saved on the local drive
# test.pdf = 12 kb,  test.wmf  = 264 kb 
# test.png = 144 kb, test.jpeg = 797 kb 

# C1. Switching between screen and file devices by commenting in/out
# win.graph(width = 6, height = 6); bringToTop(stay = TRUE)  # screen
png("test2.png", width = 6, height = 4, units = "in", res = 600)  # file
plot(rnorm(1000))
# ... more ...
# ... plotting ...
# ... statements ...
dev.off()  # file 

# C2. Switching between screen and file devices by recordPlot()
win.graph(width = 6, height = 4); bringToTop(stay = TRUE)  # screen
plot(rnorm(1000))     # more plotting statements
out <- recordPlot()   # record at the end
str(out); class(out)  # "recordedplot" class

png("test2.png", width = 6, height = 4, units = "in", res = 600)  # file
replayPlot(out)
dev.off() 


