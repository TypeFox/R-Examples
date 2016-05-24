crosssection <- function (humi,wind,temperature,plot.temp=TRUE,plot.wind=TRUE,
colors=c("brown", "yellow","green"),ylab_tics=c(0,0.2,0.4,0.6,0.8,0.9),
ylab=c(1000,800,600,400,200,100),h_limit=54,h_step=3,p_nr=10) {

humi2=as.matrix(humi)
if (plot.wind==TRUE)
{
wind2=as.matrix(wind)
}
if (plot.temp==TRUE)
{
temperature2=as.matrix(temperature)
}

rgb.palette <- colorRampPalette(colors, space = "rgb")
filled.contour(humi2,plot.axes={axis(2,ylab_tics,labels=ylab)},col=rgb.palette(20))
par(new=TRUE)
par(mar=c(4, 3, 2, 7))

if (plot.temp==TRUE)
{
contour(temperature2,axes=FALSE)
}

if (plot.wind==TRUE)
{
for (i in seq(1,p_nr))
{
        for(j in seq(1,h_limit/h_step))
                {
                station.symbol((j-1)*h_step/h_limit,(i-1)/p_nr,
speed=sqrt((wind2[j,i]^2+wind2[j,p_nr+i]^2)),
direction=180+atan2(wind2[j,i],wind2[j,p_nr+i])*57.295,
fill=0,cex=1,circle=FALSE)
                }
}
}

}

temperatures <- function (temperature.data,plot.dewt=TRUE,plot.surf=TRUE,
plot.min_max=TRUE) {

Temperature <- NULL
minT <- NULL
maxT <- NULL
Tdew <- NULL
surf.temp <- NULL

tt=ggplot(data=temperature.data, aes(x=time, y=Temperature, ymin=minT, ymax=maxT),col="red") + geom_line(col="red") + geom_text(aes(label=Temperature),size=3,vjust=-2,col="red")

if (plot.dewt==TRUE)
{
tt=tt+geom_line(data=temperature.data, aes(x=time, y=Tdew),col="blue")
}
if (plot.surf==TRUE)
{
tt=tt+geom_line(data=temperature.data, aes(x=time, y=surf.temp),col="black")
}
if (plot.min_max==TRUE)
{
tt=tt+geom_line(data=temperature.data, aes(x=time, y=maxT),lwd=0.1) + geom_line(data=temperature.data, aes(x=time, y=minT),lwd=0.1) + geom_ribbon(alpha=0.1)
}

tt
}

