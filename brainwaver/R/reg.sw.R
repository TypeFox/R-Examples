reg.sw<-function(n.nodes.rand,n.edges.rand,dist.mat)
{

in.degree<-0
Lp.reg<-0
Cp.reg<-0
countLp<-0
countCp<-0

z<-.C("Rregsim",
	as.integer(n.nodes.rand),
	as.integer(n.edges.rand),
	mat=integer(n.nodes.rand*n.nodes.rand), PACKAGE="brainwaver")

mat<-z$mat
mat<-matrix(mat,n.nodes.rand,n.nodes.rand)

z<-.C("Rkfun",
	as.integer(n.nodes.rand),
	as.integer(mat),
	deg=double(n.nodes.rand), PACKAGE="brainwaver")

in.degree<-in.degree+mean(z$deg)			

z<-.C("Rlpfun",
	as.integer(n.nodes.rand),
	as.integer(mat),
	as.double(dist.mat),
	lp=double(n.nodes.rand), PACKAGE="brainwaver")

lp<-z$lp

lp.true<-lp[lp>-1]
countLp <- countLp+length(lp.true)
Lp.reg<-Lp.reg+sum(lp.true)


z<-.C("Rcpfun",
	as.integer(n.nodes.rand),
	as.integer(mat),
	cp=double(n.nodes.rand), PACKAGE="brainwaver")


cp<-z$cp

cp.true<-cp[cp>-1]
countCp <- countCp+length(cp.true)
Cp.reg<-Cp.reg+sum(cp.true)


if(countLp!=0){ 
Lp.reg<- Lp.reg/countLp
}else{
Lp.reg<- NaN
}

if(countCp!=0){ 
Cp.reg<- Cp.reg/countCp
}else{
Cp.reg<- NaN
}

list(in.degree=in.degree,Lp.reg=Lp.reg,Cp.reg=Cp.reg)


}








