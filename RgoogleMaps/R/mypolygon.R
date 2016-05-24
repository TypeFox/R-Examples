`mypolygon` <-structure(function#simple wrapper function to plot colored polygons
###same as \link{polygon}, execept the value for color is taken from the 1st element of the exra column 'col'
(
  x, ##<< matrix containing columns X,Y,col
  ...##<< extra arguments passed to  \link{polygon}
){
  polygon(x[,c("X","Y")],col=x[1,"col"],...)  
})


