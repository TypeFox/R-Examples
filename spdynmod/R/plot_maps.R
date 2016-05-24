#' Plot abundance maps of plant communities.
#' 
#' Plot abundance maps of plant communities in a given year. 
#'
#' @param year year to plot (from 1984 to 2008)
#'
#' @return by default plots the final plant communities map (year = 2008).
#'
#' @keywords plot
#'
#' @export
#' 
#' @examples
#' ## Not run plot_maps()

plot_maps<-function(year = 2008) { 

out<-get('out')
nr<-get('nr')
nc<-get('nc')
NN<-get('NN')

i <- (year-1984)*4

if(i==0){i = 1}

out[out<0]<-0
out[out>25]<-25

print(paste('year = ',year))

a0<-raster(matrix(nrow = nr, ncol = nc, out[i, 2:(NN+1)]))

b0<-raster(matrix(nrow = nr, ncol = nc, out[i, (NN+2):(2*NN+1)]))

c0<-raster(matrix(nrow = nr, ncol = nc, out[i, (2*NN+2):(3*NN+1)]))

d0<-raster(matrix(nrow = nr, ncol = nc, out[i, (3*NN+2):(4*NN+1)]))

volc <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#own<-c("white","orange","red")
#own2<-c('#ffeda0','#feb24c','#f03b20')
own2<-c('#fed976','#feb24c','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')# '#d9f0a3'
map0<-stack(a0,b0,c0,d0)
names(map0)<-c('Salt marsh','Salt steppe','Reed beds','Bare soil')
#map<-spplot(map0,col.regions = terrain.colors(25))
map<-spplot(map0,zlim=c(0,25),col.regions = colorRampPalette(own2)(25))
map
}

