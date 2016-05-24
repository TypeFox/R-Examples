kFcheck <-
function( mat.output, name.iso,
           plot.lines=TRUE, title1=' ')
{
### mat.output is part of the outcome of kf.sort
### isopter is the position of the isopter in mat.output
### isopter == 3 is III4e, 4 is I4e, and 5 is I2e
### 
warning("This function only runs after kf.sort has been run")

arcsM<- rad(seq( 0,  345, by=15)) ### angles for centra lines of arcs

ind.isopter<- ifelse(name.iso=='III4e', 3,
              ifelse(name.iso=='I4e', 4, 5))

xy<- cbind(
     mat.output[, ind.isopter]*cos(arcsM),
     mat.output[, ind.isopter]*sin(arcsM))

set.template()

title(main=title1)
if(!plot.lines) 
{points(xy[,1], xy[,2])
id<- mat.output[, 1] 
print(paste('There are ',length(unique(id)),' datasets', sep=''),quote=FALSE)
}
else
{
k0<-0
id<- mat.output[, 1] 
print(paste('There are ',length(unique(id)),' datasets', sep=''),quote=FALSE)

for (i in unique( id))  
    {
        k0<-k0+1
        b0<- i==id
        x0<- c(xy[b0,1], xy[b0,1][1])
        y0<- c(xy[b0,2], xy[b0,2][1])
        lines(x0,y0, col=k0, lwd=1)
     }
}
invisible(NULL)
}
