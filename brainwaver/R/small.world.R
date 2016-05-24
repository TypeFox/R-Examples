small.world <- function(wave.adj.mat,dat="reduced",distance="norm",coord=0, export.data=FALSE)
{
n.regions<-dim(wave.adj.mat)[1]

if(distance=="euclid"){

euclid <- 2 * dist(t(coord), method = "euclidean")

x.euclid <- as.matrix(euclid)
}

# For the computation of Cp and Lp

if(distance=="euclid"){
dist.mat<-x.euclid
}else{
dist.mat<-matrix(1,n.regions,n.regions)
}


mat<-wave.adj.mat

# Test : diagonal of mat. 

if(sum(diag(mat)>0)>0) stop("terms on the diagonal of the adjacency matrix are not all equal to 0")

z<-.C("Rkfun",
	as.integer(n.regions),
	as.integer(mat),
	deg=double(n.regions), PACKAGE="brainwaver")
deg<-z$deg

if(dat=="all") in.degree<-deg

in.degree.mean<-mean(deg)


z<-.C("Rlpfun",
	as.integer(n.regions),
	as.integer(mat),
	as.double(dist.mat),
	lp=double(n.regions), PACKAGE="brainwaver")

lp<-z$lp

if(dat=="all") Lp<-lp

lp<-lp[lp>-1]

mlp=mean(lp)

if(mlp=='NaN'){
Lp.mean<- NaN
}else{
Lp.mean<-mlp
}

z<-.C("Rcpfun",
	as.integer(n.regions),
	as.integer(mat),
	cp=double(n.regions), PACKAGE="brainwaver")

cp<-z$cp



if(dat=="all") Cp<-cp

cp<-cp[cp>-1]

mcp<-mean(cp)

if(mcp=='NaN'){
Cp.mean<- NaN
}else{
Cp.mean<-mcp
}

z<-.C("Rconnex_components",
	as.integer(n.regions),
	as.integer(mat),
	acc=double(n.regions+3), PACKAGE="brainwaver")

size.large.connex<-z$acc[n.regions+2]





#small-world parameters


if(export.data==TRUE){

if(dat == "all"){
write.table(in.degree,"in_degree.txt")
write.table(Lp,"Lp_txt")
write.table(Cp,"Cp_txt")
}

write.table(in.degree.mean,"in_degree_mean.txt")
write.table(Lp.mean,"Lp_mean.txt")
write.table(Cp.mean,"Cp_mean.txt")
}

if(dat=="all"){ 
list(in.degree.mean=in.degree.mean,in.degree=in.degree, Cp=Cp, Cp.mean=Cp.mean,Lp.mean=Lp.mean, Lp=Lp, size.large.connex=size.large.connex)
}else{
list(in.degree.mean=in.degree.mean,Cp.mean=Cp.mean,Lp.mean=Lp.mean,size.large.connex=size.large.connex)
}

}

