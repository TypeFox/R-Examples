fit_t_comp<-function(phylo,data,model=c("MC","DDexp","DDlin"),par=NULL,geography.object=NULL,method="Nelder-Mead",bounds=NULL){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(is.null(names(data))){stop("data missing taxa names")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}
if(is.null(par)){par<-c(log(var(data)/max(nodeHeights(phylo))),0)}
if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
        bounds$lower = -Inf
        bounds$upper = Inf
    }

if(is.null(geography.object)){
	if(model=="MC"){
		opt<-optim(par,likelihood_t_MC,phylo=phylo,data=data,method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])
		S = opt$par[2]
		V = .VCV.rescale(phylo,sig2,0,S)
		data<-as.matrix(data[rownames(V)])
		IV<-solve(V)
		I<-matrix(rep(1,length(phylo$tip.label)))
		z0<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = sig2, S = S, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par,likelihood_t_DD,phylo=phylo,data=data,model="DDexp",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])
		r = opt$par[2]
		V = .vcv.rescale.DDexp(phylo,sig2,r)
		data<-as.matrix(data[rownames(V)])
		IV<-solve(V)
		I<-matrix(rep(1,length(phylo$tip.label)))
		z0<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = sig2, r = r, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		opt<-optim(par,likelihood_t_DD,phylo=phylo,data=data,model="DDlin",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])
		b = opt$par[2]
		V = .vcv.rescale.DDlin(phylo,sig2,b)
		data<-as.matrix(data[rownames(V)])
		IV<-solve(V)
		I<-matrix(rep(1,length(phylo$tip.label)))
		z0<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = sig2, b = b, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object)){
#check to make sure length matches length of nodeDiff
	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(model=="MC"){
		opt<-optim(par,likelihood_t_MC_geog,phylo=phylo,geo.object=geography.object,data=data,method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])
		S = -abs(opt$par[2])
		V = .VCV.rescale.geog(phylo,sig2,0,S,geography.object)
		data<-as.matrix(data[rownames(V)])
		IV<-solve(V)
		I<-matrix(rep(1,length(phylo$tip.label)))
		z0<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = sig2, S = S, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		opt<-optim(par,likelihood_t_DD_geog,phylo=phylo,geo.object=geography.object,data=data,model="DDexp",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])
		r = opt$par[2]
		V = .VCV.rescale.DDexp_geog(phylo,sig2,r,geography.object,check=FALSE)
		data<-as.matrix(data[rownames(V)])
		IV<-solve(V)
		I<-matrix(rep(1,length(phylo$tip.label)))
		z0<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = sig2, r = r, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		geography.matrix<-geography.object$geography.object
		maxN<-max(vapply(geography.matrix,function(x)max(rowSums(x)),1))
		opt<-optim(par,likelihood_t_DD_geog,phylo=phylo,geo.object=geography.object,data=data,model="DDlin",maxN=maxN,method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])
		b = opt$par[2]
		V = .VCV.rescale.DDlin_geog(phylo,sig2,b,geography.object,check=FALSE)
		data<-as.matrix(data[rownames(V)])
		IV<-solve(V)
		I<-matrix(rep(1,length(phylo$tip.label)))
		z0<-solve(t(I)%*%IV%*%I)%*%t(I)%*%IV%*%data[,1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.paramaters = 3, sig2 = sig2, b = b, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}
}