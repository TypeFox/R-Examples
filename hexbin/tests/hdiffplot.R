library(hexbin)

if(R.version$major != "1" || as.numeric(R.version$minor) >= 7)
    RNGversion("1.6")
set.seed(213)
x1 <- rnorm(10000)
y1 <- rnorm(10000)

x2 <- rnorm(10000,mean = .3)
y2 <- rnorm(10000,mean = .3)

rx <- range(x1,x2)
ry <- range(y1,y2)

str(bin1 <- hexbin(x1,y1, xbnds = rx, ybnds = ry))
str(bin2 <- hexbin(x2,y2, xbnds = rx, ybnds = ry))

str(erode(bin1))

str(smbin1 <- smooth.hexbin(bin1))
(smbin2 <- smooth.hexbin(bin2))

str(erodebin1 <- erode.hexbin(smbin1))
(erodebin2 <- erode.hexbin(smbin2))

if(FALSE)## does not work -- what funny stuff is hdiffplot() doing???
    par(mfrow = c(2,1))

if(exists("hdiffplot", mode="function")) { ## not yet in new hexbin
hdiffplot(bin1,bin2, main = "Original N(0,*) Random bins")

hdiffplot(smbin1,smbin2, main = "smooth.hexbin() smoothed bins")

plot.new()
hdiffplot(erodebin1,erodebin2, main = "erode.hexbin()d smoothed bins")
}# not yet
