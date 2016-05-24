# This demo illustrates parallel coordinate-type displays for data with mixed categorical continuous variables. The displays are similar to Parallel Sets. http://www.vrvis.at/via/research/parsets/index.html
# For the example here to work, you must first load the following functions.

require(PairViz)
library(alr3)

source(system.file("demo","demo_fns.R", package = "PairViz"))


#---------------------------------
# Now the example

data(donner)


d <- na.omit(donner)
d <- d[,-4]
d <- d[,c(2,3,4,1)]

colvar <- 2 # any of the categorical variables 1:3 will work here
d <- d[order(d[,colvar]),]

cols <- colour_var(d[,colvar])
catvars <- 1:3

dspread <- factor_spreadout(d[,catvars])
ds <- d
ds[,catvars] <- dspread[["data"]]



o <- hpaths(1:ncol(d),matrix=FALSE)




dev.new(width=5,height=2.5)
par(mar=c(2,1,2,1))
par( cex.axis=.8,cex=.8)
catpcp(ds,order=o,col=cols,lwd=2,pcpbars=dspread$bars,barvars=catvars,main="Donner Data",pcpbars.labels=T)


# benefits: handles any number of categorical, and continuous
# see al pairwise relationships- eg, single and hired people are in older age groups and these are almost all male (1 exception).

#--------------------------


y <- as.data.frame(as.table(HairEyeColor))

colvar <- 3 # any of 1:3 will do
y <- y[order(y[,colvar]),] # ensures that cases are ordered by colour within each factor level
ylong <- apply(y[,-4],2, function(x) rep(x,times=y[,4]))


cols <- colour_var(ylong[,colvar])



ds <- factor_spreadout(ylong)


o <- c(1,2,3,1) 

dev.new(width=5,height=2.5)
par(mar=c(2,1,2,1))
par( cex.axis=.8,cex=.8)



catpcp(ds$data,order=o,col=cols,pcpbars=ds$bars,pcpbars.labels=T,main="Hair Eye data")

