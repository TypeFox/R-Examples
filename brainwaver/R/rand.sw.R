rand.sw<-function(nsim,n.nodes.rand,n.edges.rand,dist.mat,dat="reduced")
{

if(dat=="all"){
in.degree.all<-rep(0,n.nodes.rand)
Lp.all<-rep(0,n.nodes.rand)
Cp.all<-rep(0,n.nodes.rand)
}

in.degree<-0
Lp.rand<-0
Cp.rand<-0
countLp<-0
countCp<-0

for(sim in 1:nsim){


z<-.C("Rrandsim",
	as.integer(n.nodes.rand),
	as.integer(n.edges.rand),
	mat=integer(n.nodes.rand*n.nodes.rand), PACKAGE="brainwaver")

mat<-z$mat
mat<-matrix(mat,n.nodes.rand,n.nodes.rand)

z<-.C("Rkfun",
	as.integer(n.nodes.rand),
	as.integer(mat),
	deg=double(n.nodes.rand), PACKAGE="brainwaver")

if(dat=="all") in.degree.all<-in.degree.all+z$deg

in.degree<-in.degree+mean(z$deg)			

z<-.C("Rlpfun",
	as.integer(n.nodes.rand),
	as.integer(mat),
	as.double(dist.mat),
	lp=double(n.nodes.rand), PACKAGE="brainwaver")

lp<-z$lp

if(dat=="all") Lp.all<-lp+Lp.all

lp.true<-lp[lp>-1]
countLp <- countLp+length(lp.true)
Lp.rand<-Lp.rand+sum(lp.true)


z<-.C("Rcpfun",
	as.integer(n.nodes.rand),
	as.integer(mat),
	cp=double(n.nodes.rand), PACKAGE="brainwaver")


cp<-z$cp

if(dat=="all") Cp.all<-cp+Cp.all


cp.true<-cp[cp>-1]
countCp <- countCp+length(cp.true)
Cp.rand<-Cp.rand+sum(cp.true)

}



in.degree<-in.degree/(sim)
if(dat=="all") in.degree.all<-in.degree.all/(sim)

if(countLp!=0){ 
Lp.rand<- Lp.rand/countLp
if(dat=="all") Lp.all<- Lp.all/sim
}else{
Lp.rand<- NaN
}

if(countCp!=0){ 
Cp.rand<- Cp.rand/countCp
if(dat=="all") Cp.all<- Cp.all/sim
}else{
Cp.rand<- NaN
}

if(dat=="all"){
list(in.degree=in.degree,Lp.rand=Lp.rand,Cp.rand=Cp.rand,in.degree.all=in.degree.all,Lp.all=Lp.all,Cp.all=Cp.all)
}else{ 
list(in.degree=in.degree,Lp.rand=Lp.rand,Cp.rand=Cp.rand)
}

}








