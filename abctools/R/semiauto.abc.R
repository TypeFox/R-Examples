semiauto.abc<-function(obs, param,sumstats,obspar=NULL,abcmethod=abc,saprop=.5,abcprop=.5,overlap=FALSE,satr=list(),plot=FALSE,verbose=TRUE,do.err=FALSE,final.dens=FALSE,errfn=rsse,...){

if(!is.matrix(obs)) {
    obs <- matrix(obs, nrow=1)
}
if(is.data.frame(obs)){
    obs <- as.matrix(obs)
}
if(!is.matrix(param)|is.data.frame(param)){
        param<-as.matrix(param)
}
if(!is.matrix(sumstats)|is.data.frame(sumstats)){
        sumstats<-as.matrix(sumstats)
}
if(!is.null(obspar)|is.data.frame(obspar)){
        if(!is.matrix(obspar)){
                obspar<-matrix(obspar,byrow=T,ncol=ncol(param))
        }
        if(nrow(obs)!=nrow(obspar)){
                stop("Please supply observed statistics and observed parameter matrices with the same number of rows!\n")
        }
}

q<- (!is.null(obspar))&(do.err)

if(!q){
        do.err<-FALSE
}

##SPLIT INTO GROUPS
nsims<-nrow(sumstats)
size1<-floor(nsims*saprop)
size2<-floor(nsims*abcprop)

if(overlap){
	# allow for "overlapping" indices of sims for construction / abc.

	tobuild <- sample(1:nsims,size1)
	forabc <- sample(1:nsims,size2)
}
else{
	tobuild <- sample(1:nsims,size1)

	nobuild <- setdiff(1:nsims,tobuild)
	forabc <- sample(nobuild,size2)
}

##SEMI-AUTOMATIC ABC

if(verbose){
	cat("Doing statistics regression with sample size:",size1,"\n")
}

sa.param <- param[tobuild,]
abc.param <- param[forabc,]
if(!is.matrix(sa.param)){
        sa.param<-as.matrix(sa.param)
}
if(!is.matrix(abc.param)){
        abc.param<-as.matrix(abc.param)
}


sumstats.tr<-sa.ss.tr<-obs.tr<-NULL

# check ss transformations:

#if(!is.list(satr)){
#	satr<-as.list(satr)
#}

ntr<-length(satr)

fncheck<-function(s,f,i){

	nrs<-nrow(s)

	te<-try(f(s),silent=TRUE)

	if(class(te)=="try-error"){
		stop(paste("satr function",i,"not valid!!\n"))
	}
	else{
		out<-as.matrix(te)
	}

	if((nrow(out)%%nrs)!=0){
		stop(paste("satr function",i,"does not give valid output!!\n"))
	}

return(matrix(out,nrow=nrs))

}


#fns<-sapply(satr,is.function)

#if(any(!fns)){
#	notfns<-which(!fns)
#	if(any(!fns)){
#		stop(paste("statistics transformations are not all valid.  Please check transformations: ",notfns,"!\n",sep=""))
#	}
#}

if(ntr==0){
	# if no transformations are given, use identity (do nothing).
	sumstats.tr <-sumstats
	obs.tr <-obs
}
else{	# cbind all transformed statistics together.  There is probably a more efficient way of doing this.
	for(i in 1:ntr){
		# first do simulated statistics:
		#trss<-eval(satr[[i]](sumstats))
		trss<-fncheck(sumstats,satr[[i]],i)
		#trss<-matrix(trss,nrow=nsims)
		sumstats.tr<-cbind(sumstats.tr,trss)
		# now do observed statistics:
		# if a non-valid function is given, we shouldn't get to here...
		trss<-eval(satr[[i]](obs))
		trss<-matrix(trss,nrow=nrow(obs))
		obs.tr<-cbind(obs.tr,trss)
	}
}

sa <- saABC(sa.param,sumstats.tr[tobuild,],plot=plot)

B <- sa$B
B[is.na(B)] <- 0 ##NAs may exist due to collinearity of sumstats.tr[tobuild,] but can safely be set to zero
ss.sa <- sumstats.tr %*% t(B)
obs.sa <-  obs.tr %*% t(B)

if(verbose){
	cat("Doing ABC with sample size:",size2,"\n")
}

argl <- list(...)
    targind <- match(names(argl), "tol")
    targind <- which(!is.na(targind))
    margind <- match(names(argl), "method")
    margind <- which(!is.na(margind))
    if ((length(targind) == 0) & identical(abcmethod, abc)) {
        argl$tol <- 0.01
    }
    if ((length(margind) == 0) & identical(abcmethod, abc)) {
        argl$method <- "rejection"
    }
    argl$target=obs.sa
    argl$param=param[forabc,]
    argl$sumstat=ss.sa[forabc,]
    abcout.sa <-do.call(abcmethod,argl)
	
    #abcout.sa <- abcmethod(obs.sa, param[forabc,], ss.sa[forabc,], ...)
if(is.null(abcout.sa$adj.values)){
	vals<-abcout.sa$unadj.values
}
else{
        vals<-abcout.sa$adj.values
}

# return useful things:

l<-list()

if(final.dens){
        l$post.sample<-vals
}

if(do.err){
	err<-errfn(vals,obspar,apply(param,2,var))
        l$err<-err
}

sainfo<-list(saprop=saprop,abcprop=abcprop,overlap=overlap,satr=satr)
l$sainfo<-sainfo

return(l)

}
