stop.identify <-
function( xy)
{
#### xy is an n x 2 matrix with positions of points
#### ordered.list is a vector containing an ordered list of  identified points

set.template()
title(main='Join isopter points in order starting at 180 degrees working clockwise')

points( xy, pch=19)
ordered.list<-NA

stop1<- FALSE
for (i in 1:1000)
{
while(!stop1)
{
aux<- identify( xy, n=1)
stop1<- length( aux)==0
ordered.list<- c( ordered.list, aux)
}
}
ordered.list[-1]
}
