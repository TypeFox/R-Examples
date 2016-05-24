opar <- par(ask = dev.interactive(orNone = TRUE))
#get(getOption("device"))(width=12,height=12)
#oldpar<-par(mfrow=c(2,4))
lshorth(runif(25))
lshorth(runif(50))
lshorth(runif(100))
lshorth(runif(200))
lshorth(rnorm(25))
lshorth(rnorm(50))
lshorth(rnorm(100))
lshorth(rnorm(200))
#par(oldpar)

#dev.off()
#get(getOption("device"))(width=6,height=6)

# note: these are time series data
lshorth(AirPassengers)

# note: these are time series data
lshorth(LakeHuron)

# note: these are time series data
lshorth(Nile)


# note: these are time series data
work <- diff(WWWusage)
oldpar <- par(mfrow = c(2,1))
lshorth(WWWusage)
lshorth(work)
par(oldpar)

oldpar <- par(mfrow = c(2,1))
lshorth(faithful[,1])
lshorth(faithful[,2])
par(oldpar)

lshorth(rivers)

# note: these are time series data
oldpar <- par(mfrow = c(2,1))
lshorth(sunspot.month)
lshorth(sunspot.year)
par(oldpar)

par(opar)