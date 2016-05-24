#' Make an animation of the model output in relation to bare soil.
#' 
#' It requires ImageMagick installed in the system to work.
#'
#' @return a GIF animated file
#'
#' @keywords animation
#'
#' @export
#' 
#' @examples
#' ## Not run animate_bs()

animate_bs<-function(){
#require(raster)
out<-get('out')
nr<-get('nr')
nc<-get('nc')
NN<-get('NN')

brks <- seq(0, 30, by=1) 
nb <- length(brks)-1 

#require(animation)

animation::saveGIF({
for (i in seq(1, dim(out)[1], by = 1)){
 i2<-seq(1984.25,2008,0.25)
 par(mfrow=c(2,1))
 par(mar=c(2,1,2,0.8))
 image(raster::raster(matrix(nrow = nr, ncol = nc, out[i, (3*NN+2):(4*NN+1)])),breaks=brks, col=rev(topo.colors(nb)), lab.breaks=brks, zlim=c(0,30),main="Bare Soil")
 par(mar=c(13,1,3,1.5))
 print(barplot(i2[i],col="black",horiz=T,xlim=c(1984,2008),axes=F,cex.sub=1.4, main="Time",cex.main=1.5))
 par(las=2)
 axis(1,at=c(1984,1992,1995,1997,2001,2008))
 }},movie.name='movie_bs.gif')

}

### PLOT GOOGLEVIS CHART

## Data preparation
#out2<-as.data.frame(out)
#b<-1:96
#DatosApilados <- stack(out2[, c('Salt marsh','Salt steppe','Reed beds','Bare soil')])
#tiempo<-rep(b,4)
#DatosApilados2<-DatosApilados
#DatosApilados2[,1]<-DatosApilados$ind
#DatosApilados2[,2]<-tiempo
#DatosApilados2[,3]<-DatosApilados$values
#names(DatosApilados2) <- c("factor", "tiempo","valor")

## Plot
#require(googleVis)
#Motion=gvisMotionChart(DatosApilados2, idvar="factor", timevar="tiempo", options=list(height=350, width=400))
## Display chart
#plot(Motion) 
## Create Google Gadget
##cat(createGoogleGadget(Motion), file="motionchart.xml")

###
