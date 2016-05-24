rgl.sphsun=function (Ydate = c(3, 21), radius = 1, col='yellow', type='s', sunrad=0.02, addeclip=TRUE, addsun=TRUE) 
{
    if(addsun){
        if(Ydate[1]=='get'){Ydate=c(format(Sys.Date(), "%m"),format(Sys.Date(), "%d"));Ydate=as.numeric(Ydate)}
        year = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        yearfrac = sum(year[1:Ydate[1]]) - Ydate[1] + Ydate[2]
        yearfrac = yearfrac - 108
        if (yearfrac < 0) {
            yearfrac = yearfrac + 365
        }
        yearfrac = (yearfrac * 2 * pi/365)*180/pi
        sunloc = cbind(yearfrac,sin(yearfrac*pi/180)*23.4,radius)
        plot3d(sph2car(sunloc,deg=TRUE),type=type,col=col,radius=sunrad,add=TRUE)
    }
    if(addeclip){rgl.sphcirc(CrossEq = 0, PeakDec = 23.4,col=col,radius=radius)}
}

