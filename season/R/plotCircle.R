# plotCircle.R
# Function to plot monthly rates using a grey shaded circular "country"
# Assumes estimates are in numerical order, 1=Jan, 2=Feb, etc
# April 2009

plotCircle<-function(months,dp=1, ...){
if ((length(months)==12)==FALSE) {stop('Length of monthly estimates must be 12')}
op <- par(no.readonly = TRUE) # the whole list of settable par's.
# create grey scale, standardise results to [0.2,1]
# add small number to prevent very dark colours
stan0to1<-(months-min(months))/(max(months)-min(months))
stan1to0<-abs(stan0to1-1) # reverse so that darker shades are higher
stan<-(stan1to0*0.7)+0.3 # [0.3,1]
dens=gray(stan)
plot(c(0),c(0),bty='n',xlab='',main='',ylab='',type='p',xlim=c(-1,1),ylim=c(-1,1),
 col='white',yaxt='n',xaxt='n', ...) # blank plot

# important constants
bins<-12
clockstart=pi/2 # default clock start at 12 o'clock
half<- 2*pi/(bins*2) # for moving text half-way round
start<-pi/4
scale<-1
detail<-50
# Polygons for each month
for (m in 1:12){
   poly<-matrix(data=NA,nrow=((detail+1)*2)+3,ncol=2)
   index<-0
   index<-index+1
   mstart<-(2*pi*m/bins)+start
   poly[index,1]<- -1*scale*cos(mstart);
   poly[index,2]<-1*scale*sin(mstart);
   index<-index+1
   poly[index,1]<- -0.7*scale*cos(mstart);
   poly[index,2]<-0.7*scale*sin(mstart);
   for (i in 0:detail){
      index<-index+1
      x<-(mstart+(2*pi*i/(detail*bins)))
      poly[index,1]<- -0.7*scale*cos(x);
      poly[index,2]<-0.7*scale*sin(x);
   }
   index<-index+1
   poly[index,1]<- -1*scale*cos(x);
   poly[index,2]<-1*scale*sin(x);
   for (i in detail:0){
      index<-index+1
      x<-(mstart+(2*pi*i/(detail*bins)))
      poly[index,1]<- -1*scale*cos(x);
      poly[index,2]<-1*scale*sin(x);
   }
   polygon(poly,col=dens[m],border='black')
}
# Month labels
for (j in 1:bins){
   x<- -0.85*cos((2*pi*j/bins)+start+half)
   y<-0.85*sin((2*pi*j/bins)+start+half)
   label<-month.abb[j]
   text(x,y,label)
}
## Add a colour bar to the centre of the plot
# total bar height = 1
bartoplot<-seq(1,0.3,-0.1)
dens<-gray(bartoplot)
width<-0.1
height<-0.05
nbars<-length(bartoplot)
# unstandardise results for text label
unstan<-(((bartoplot-0.3)/0.7)*(max(months)-min(months)))+min(months)
for (i in 1:nbars){
   x<-c(-width,width,width,-width,-width)-width
   yref<-((i-1)-(nbars/2))/(nbars+2)
   y<-c(yref-height,yref-height,yref+height,yref+height,yref-height)
   polygon(x,y,col=dens[i],border=NA,lwd=1)
# text label
   clabel2<-formatC(unstan[nbars-i+1], format="f", digits=dp) # convert to character
   text(2*width,yref,clabel2)
}
par(op) # restore graphic settings
} # end of function
