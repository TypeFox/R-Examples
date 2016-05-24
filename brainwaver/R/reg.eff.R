reg.eff<-function(n.nodes.rand,n.edges.rand,dist.mat)
{

eff<-0
loc<-0
cost<-0

z<-.C("Rregsim",
	as.integer(n.nodes.rand),
	as.integer(n.edges.rand),
	mat=integer(n.nodes.rand*n.nodes.rand), PACKAGE="brainwaver")

mat<-z$mat
mat<-matrix(mat,n.nodes.rand,n.nodes.rand)
eff<-global.efficiency(mat,dist.mat)

loc<-local.efficiency(mat,dist.mat)
cost<-global.cost(mat,dist.mat)


list(eff=eff$eff,loc=loc$eff,cost=cost)
}






