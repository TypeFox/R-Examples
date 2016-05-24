rand.eff<-function(nsim,n.nodes.rand,n.edges.rand,dist.mat,dat="reduced")
{

if(dat=="all"){
reg.eff<-rep(0,n.nodes.rand)
loc.eff<-rep(0,n.nodes.rand)
}

eff<-0
loc<-0
cost<-0
for(sim in 1:nsim){


z<-.C("Rrandsim",
	as.integer(n.nodes.rand),
	as.integer(n.edges.rand),
	mat=integer(n.nodes.rand*n.nodes.rand), PACKAGE="brainwaver")

mat<-z$mat
mat<-matrix(mat,n.nodes.rand,n.nodes.rand)

tmp.eff<-global.efficiency(mat,dist.mat)

eff<-tmp.eff$eff+eff

tmp.loc<-local.efficiency(mat,dist.mat)

tmp.cost<-global.cost(mat,dist.mat)

cost<-tmp.cost+cost
loc<-tmp.loc$eff+loc

if(dat=="all"){

loc.eff<-loc.eff+tmp.loc$loc.eff

}

}

eff<-eff/nsim
loc<-loc/nsim
cost<-cost/nsim


list(eff=eff,loc=loc,cost=cost)
}


