

#########  get some data 
data(KH)


##### informaiton on the data structure:
names(KH)

#####  seismic wiggles are in JSTR


dt = KH$dt[1]

y =  KH$JSTR[[1]]

x =  seq(from=0, by=dt, length=length(y))

##########  plot the data:
plot(x,y, type='l', xlab='s', ylab='amp', main="Vertical Comp Seismic Volcano Explosion")

readline("To Continue Hit Enter Key\n")


################    filter the data

#######  use a band pass, butterworth filter from 10 s to 1 Hz

fy = butfilt(y, .1, 1, dt, "BP", "BU" )


plot(x,fy, type='l', xlab='s', ylab='amp', main="Bandpass filter  10 s to 1 Hz" )

readline("To Continue Hit Enter Key\n")

################ make a plot of both traces on the same figure:
PLOT.MATN(cbind(y, fy), tim=x, dt=dt , notes=c("raw", "filtered" ) , COL=c('brown', 'blue'))

readline("To Continue Hit Enter Key\n")


################ zoom in on part of the plot
## could use locator here as in win=locator(2); win=win$x

win =c(100, 200)


PLOT.MATN(cbind(y, fy), tim=x, dt=dt , WIN=win, notes=c("raw", "filtered" ) , COL=c('brown', 'blue'))

readline("To Continue Hit Enter Key\n")


##############################################
########   repeat  filtering example with pre and post-tapering
###  remove mean of data
By = y  - mean(y)
####
#### taper the data
By = applytaper(By, p=0.05)
flow =1/100
fhi = 1
fy = butfilt(By, flow, fhi, dt, "BP", "BU" )
fy = applytaper(fy, p=0.05)

flab =paste(sep=" ", "BP filtered:", flow, "to", fhi, "Hz" )

  
PLOT.MATN(cbind(y, fy), tim=x, dt=dt , notes=c("raw", flab) , COL=c('brown', 'blue'))


readline("To Continue Hit Enter Key\n")

##############################################
##############################################   suite of filters:
##############################################

###  set up filters:
fl=rep(1/100, 5)
fh=1/c(1,2,5,10,20)



###  show a filter spread using the above definition:

FILT.spread(x, y, dt, fl = fl, fh = fh, sfact = 1, WIN = NULL, PLOT = TRUE, TIT = NULL, TAPER = 0.05, POSTTAPER=0.1)


readline("To Continue Hit Enter Key\n")


###  zoom in on filter spread
FILT.spread(x, y, dt, fl = fl, fh = fh, sfact = 1, WIN = win, PLOT = TRUE, TIT = NULL, TAPER = 0.05)


readline("To Continue Hit Enter Key\n")


##############################################
