rgl.sphgrid <-
function (radius = 1, col.long='red', col.lat='blue', deggap = 15, longtype='H', add=FALSE, radaxis=TRUE, radlab='Radius') 
{
    if(add==F){open3d()}
    for(lat in seq(-90,90,by=deggap)){
    if(lat==0){col.grid='grey50'}else{col.grid='grey'}
    plot3d(sph2car(long=seq(0,360,len=100),lat=lat,radius=radius,deg=T),col=col.grid,add=T,type='l')
    }
    for(long in seq(0,360-deggap,by=deggap)){
    if(long==0){col.grid='grey50'}else{col.grid='grey'}
    plot3d(sph2car(long=long,lat=seq(-90,90,len=100),radius=radius,deg=T),col=col.grid,add=T,type='l')
    }
    if(longtype=='H'){scale=15}
    if(longtype=='D'){scale=1}
    rgl.sphtext(long=0,lat=seq(-90,90,by=deggap), radius=radius, text = seq(-90,90,by=deggap), deg=TRUE, col = col.lat)
    rgl.sphtext(long=seq(0,360-deggap,by=deggap),lat=0, radius=radius, text = seq(0,360-deggap,by=deggap)/scale, deg=TRUE, col = col.long)

    if(radaxis){
        radpretty=pretty(c(0,radius))
        radpretty=radpretty[radpretty<=radius]
        lines3d(c(0,0),c(0,max(radpretty)),c(0,0),col='grey50')
        for(i in 1:length(radpretty)){
            lines3d(c(0,0),c(radpretty[i],radpretty[i]),c(0,0,radius/50),col='grey50')
            text3d(0,radpretty[i],radius/15,radpretty[i],col='darkgreen')
        }
        text3d(0,radius/2,-radius/25,radlab)
    }
}

