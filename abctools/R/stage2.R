stage2 <-function (obs,param, sumstats, obspar=NULL, init.best, dsets = 100, sumsubs = 1:ncol(sumstats), limit = length(sumsubs), do.only = NULL, do.err=FALSE,final.dens=FALSE,...) {

if(!is.matrix(obs)|is.data.frame(obs)){
        obs<-matrix(obs,nrow=1)
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

if (!length(colnames(param))) {
        colnames(param) <- paste("P", 1:ncol(param), sep = "")
}
if (!length(colnames(sumstats))) {
        colnames(sumstats) <- paste("C", 1:ncol(sumstats), sep = "")
}

    argl <- list(...)
    aargind <- match(names(argl),"abcmethod") 
    aargind<-which(!is.na(aargind))
    if (length(aargind)==0){
        abcmethod<-abc
    }
    else{
	abcmethod<-eval(argl[[aargind]])
    }

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
    argl$param=param

sumstats<-sumstats[,sumsubs]
obs<-obs[,sumsubs]


ndatasets <- nrow(obs)

err <- best<-vals<-NULL

nr<-nrow(param)
npar<-ncol(param)

eps2 <- dsets/nr
 cat("initial best subset is:", init.best, "\n")

cm <- combmat(length(sumsubs), limit)
if (is.null(do.only)) {
        do.only <- 1:nrow(cm)
}
if (nrow(cm) < init.best) {
        stop("value of init.best is too big for value of limit!")
}
if (length(init.best) == 1) {
        init.best <- which(cm[init.best, ] == 1)
}
l.true <- abc(obs[init.best], param, sumstats[, init.best], tol = eps2,method="rejection")
closest <- which(l.true$region)
obss <- sumstats[closest,]
obst <- param[closest,]

tmp<-selectsumm(obss,param, sumstats, obst, sumsubs = sumsubs,limit = limit,do.only = do.only, do.crit = FALSE, do.err=TRUE,final.dens=FALSE, ...) 

ave<- rowMeans(tmp$err)
best<-which.min(ave)  
besti <- do.only[best]

Ifin <- matrix(which(cm[besti, ] == 1),nrow=1)
rownames(Ifin) <- besti

if (final.dens) {
        cat("getting final posterior sample...\n")
	argl$target=obs[Ifin]
	argl$sumstat=sumstats[,Ifin]
	valsI <-do.call(abcmethod,argl)
        if(is.null(valsI$adj.values)){
                post.sample<-valsI$unadj.values
        }
        else{
                post.sample<-valsI$adj.values
        }
}

l<-list()

l$best<-Ifin
l$closest<-closest
l$posssubs<-do.only
if(do.err){
	l$err<-ave
}

if (final.dens) {
        l$post.sample <- post.sample
}
l$sumsubs<-sumsubs


return(l)

}

