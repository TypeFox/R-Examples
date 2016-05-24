equadist.rand.sw<-function(nsim,dat="reduced",dist.mat,degree.dist)
{
n.nodes.rand<-length(degree.dist)
n.edges.rand<-sum(degree.dist)/2

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
countSim<-0

for(sim in 1:nsim){

z<-sim.equadist(degree.dist)

if(sum(z)==sum(degree.dist)){
countSim<-countSim+1
mat<-z

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

}

if(countSim==0) stop("No simulation matched the given degree distribution, sorry!")
in.degree<-in.degree/(countSim)
if(dat=="all") in.degree.all<-in.degree.all/(countSim)


if(countLp!=0){ 
Lp.rand<- Lp.rand/countLp
if(dat=="all") Lp.all<- Lp.all/countSim
}else{
Lp.rand<- NaN
}

if(countCp!=0){ 
Cp.rand<- Cp.rand/countCp
if(dat=="all") Cp.all<- Cp.all/countSim
}else{
Cp.rand<- NaN
}

if(dat=="all"){
list(in.degree=in.degree,Lp.rand=Lp.rand,Cp.rand=Cp.rand,in.degree.all=in.degree.all,Lp.all=Lp.all,Cp.all=Cp.all)
}else{ 
list(in.degree=in.degree,Lp.rand=Lp.rand,Cp.rand=Cp.rand)
}



}









