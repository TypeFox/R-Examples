choropleth.spplot <-
function(sp,dem="P0010001",cuts=list("quantile",seq(0, 1, 0.25)),color=list(fun="hsv",attr=list(h = c(.4,.5,.6,.7), s = .6, v = .6, alpha=1)),main=NULL,sub="Quantiles (equal frequency)",legend=list(pos="bottomleft",title="Population Count"),border="transparent",object=NULL,...){

### Sp.layout information)
bb<-bbox(sp)

narrow <- list("SpatialPolygonsRescale", layout.north.arrow(),offset=c(bb[1,2]-.5,bb[2,2]-.45),scale = .4)


### Quantiles of lnden
cuts<-round(do.call(cuts[[1]],list(x=sp[[dem]],probs=cuts[[2]])),digits = 6)
### Color scheme for the quantiles
colForRegions<-do.call(color$fun,color$attr)

### help(spplot)
if(is.null(object)){
spplot(sp[,dem], scales=list(draw=TRUE), xlab=expression("Longitude"^o), ylab=expression("Latitude"^o), cex=1, col= border, font=2,sp.layout=list(narrow),colorkey=list(labels=list(at=cuts)),at=cuts,main=main, sub=sub,col.regions=colForRegions,...)
}else{
spplot(sp[,dem], scales=list(draw=TRUE), xlab=expression("Longitude"^o), ylab=expression("Latitude"^o), cex=1, col= border, font=2,sp.layout=c(list(narrow),object),colorkey=list(labels=list(at=cuts)),at=cuts,main=main, sub=sub,col.regions=colForRegions,...)
}

}

