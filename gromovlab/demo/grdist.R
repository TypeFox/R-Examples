#various examples
library("phangorn")
set.seed(666)
tr1<-rtree(n=10)
tr2<-rtree(n=10)
tr3<-rtree(n=10)
gromovdist(tr1,tr2,"l1")
gromovdist(tr1,tr2,"l2")
gromovdist(tr1,tr2,"linf")
gromovdist(tr1,tr2,"lp")
gromovdist(tr1,tr2,"lp",p=4)
gromovdist(tr1,tr2,"lp",p=1.5)

#list
gromovdist(list(tr1,tr2,tr3))


#igraph
library("igraph")
cp1<-graph.formula("5"-"root"-"4x" -"4" -"4x"-"3x"-"3"-"3x"-"2x"-"2"-"2x"-"1" )
plot(cp1)
cp2<-graph.formula("2"-"root"-"4x" -"4" -"4x"-"3x"-"3"-"3x"-"2x"-"5"-"2x"-"1" )
plot(cp2)
gromovdist(cp1,cp2)
gromovdist(cp1,cp2,leavesonly=FALSE)

#matrix

dm1<-matrix(c(0,1,2,1,0,1,2,1,0),nr=3)
dm2<-matrix(c(0,2,1,2,0,1,1,1,0),nr=3)

gromovdist(dm1,dm2)


#distances
library("cluster")
abc<-data.frame(a=rnorm(n=10),b=rnorm(n=10),c=rnorm(n=10))
gromovdist(lapply(c("euclidean", "manhattan", "gower"),FUN=function(t)daisy(t(abc),metric=t)))





