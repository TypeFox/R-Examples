delete.points <-
function(outer.iso)
{
### this programme deletes practice and error points selected using locator()
### outer.iso is a matrix containing practice and error points
### output is a matrix excluding the identified points called outer.iso too

set.template()
title(main='Delete unwanted isopter points')
points(outer.iso[,1], outer.iso[,2], pch=19, col='black')
n.pts<- dim( outer.iso)[1]
s1<- 1:(n.pts-1)


deleted.index<- identify( outer.iso[,1], outer.iso[,2], pos=TRUE)$ind
### extracts indices of manually selected points
if(length(deleted.index >0)) outer.iso<- outer.iso[-deleted.index, ]
dimnames( outer.iso)[[1]]<- 1:dim( outer.iso)[1]

invisible(outer.iso)
}
