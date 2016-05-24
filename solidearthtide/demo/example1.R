# Example 1
# Time series at a single location 
# Calculate two weeks of vertical displacement in Austin, TX (30.25 N,97.75 W) at hourly resolution

# original MATLAB code:
# Chad Greene, 2015
# Solid Earth Tide Predictions
# http://www.mathworks.com/matlabcentral/fileexchange/50327-sold-earth-tide-predictions/content/earthtide

# define latitude and longitude
latIn <- 30.25
lonIn <- -97.75
# define start and end dates
a<-Datenum(as.Date("may 27, 1984", "%b %d, %Y"))
b<-Datenum(as.Date("june 10, 1984", "%b %d, %Y"))
# generate a sequence of time values at hourly resolution
tseq<-seq(a,b,by=1/24)
tseq2<-UTC2GPS(tseq,TRUE)# generate a vector of GPS time values
# call Earthtide
up = Earthtide(30.25,-97.75,tseq2)
# plot the result
plot(1:337,up*100,type='l',xlab='',ylab='vertical displacement (cm)' , xaxt = "n")
axis(1, at=c(0,150,350), labels=c('5/27','6/03','6/10'))
