mlm.em <- function(y,subj,pred,xcol,start,maxits=200,epsil=.0001){
        # Calculates ML estimates for the multivariate linear model with fixed-effects only
        #
        if(any(is.na(pred)))
                stop("missing values in pred not allowed")
        # change y and pred to matrices, if necessary
        if(is.vector(y)) y <- matrix(y,ncol=1)
        if(is.vector(pred)) pred <- matrix(pred,ncol=1)
        m <- as.integer(length(table(subj)))
        ntot <- as.integer(nrow(y))
        nmax <- as.integer(max(table(subj)))
        r <- as.integer(ncol(y))
        p <- length(xcol)
        xcol <- as.integer(xcol)
        pcol <- as.integer(ncol(pred))
		storage.mode(subj)<-"integer"
        #
        # for now starting values are assumed to be supplied
        {if(missing(start)){
                beta <- matrix(0,p,r)
                sigma <- matrix(0,r,r)
                eps <- matrix(0,ntot,r)
                sflag <- as.integer(0)}
         else{
                beta <- start$beta
                sigma <- start$sigma
                eps <- matrix(0,ntot,r)
                sflag <- as.integer(1)
                storage.mode(eps) <- "double"
                storage.mode(beta) <- "double"
                storage.mode(sigma) <- "double"}}
        cat(" ###  performing mle for mult. linear model with NA values ### ")
        now <- proc.time()
        #
        # create rmat, npatt and patt to keep track of missingness patterns
        rmat <- 1-1*is.na(y)
        storage.mode(rmat) <- "integer"
        revcpatt <- rep("",ntot)
        for(i in 1:r) revcpatt <- paste(as.character(rmat[,i]),revcpatt,sep="")
        nulpat0 <- ""
        nulpat2 <- ""
        for(i in 1:r){
           nulpat0 <- paste(nulpat0,"0",sep="")
           nulpat2 <- paste(nulpat2,"2",sep="")}
        revcpatt[revcpatt==nulpat0] <- nulpat2
        tmp <- rev(table(revcpatt))
        npatt <- length(tmp)
        if(any(revcpatt==nulpat2)) npatt <- npatt-1
        ww <- !duplicated(revcpatt)
        upatt <- revcpatt[ww]
        rmat <- rmat[ww,]
        if(r==1) rmat <- matrix(rmat,ncol=1)
        ww <- rev(order(upatt))
        upatt <- upatt[ww]
        rmat <- matrix(rmat,ncol=r,nrow=length(rev(order(upatt))))
        rmat <- rmat[ww,]
        if(r==1) rmat <- matrix(rmat,ncol=1)
        if(any(upatt==nulpat2)){
                rmat <- rmat[-1,]
                upatt <- upatt[-1]}
        patt <- integer(ntot)
        patt[revcpatt==nulpat2] <- 0
        for(i in 1:npatt) patt[revcpatt==upatt[i]] <- i
        storage.mode(npatt) <- "integer"
        storage.mode(rmat) <- "integer"
        storage.mode(patt) <- "integer"
        iposn <- as.integer(1:ntot)
        ww <- order(patt)
        iposn <- iposn[ww]
        pstfin <- matrix(0,npatt,2)
        {if(any(patt==0)){
                sst <- tmp[1]+1
                for(i in 2:(npatt+1)){
                        pstfin[i-1,1] <- sst
                        pstfin[i-1,2] <- sst+tmp[i]-1
                        sst <- sst+tmp[i]}}
        else{
                sst <- 1
                for(i in 1:npatt){
                        pstfin[i,1] <- sst
                        pstfin[i,2] <- sst+tmp[i]-1
                        sst <- sst+tmp[i]}}}
        storage.mode(pstfin) <- "integer"
        #
        storage.mode(y) <- "double"
        y[is.na(y)] <- -999.99
        storage.mode(pred) <- "double"
        #####
        tmp <- .Fortran("mlem",ntot,m,r,p,as.integer(subj),ist=integer(m),
                ifin=integer(m),nmax,iposn=iposn,npatt=npatt,pstfin=pstfin,
                patt=patt,nstar=integer(1),nstari=integer(m),rmat,pcol,xcol,
                pred,wxbeta=matrix(0,ntot,r),wkrrpt=array(0,c(r,r,npatt)),y,
                #
                ey=matrix(0,ntot,r),eyyt=matrix(0,r*nmax,r*nmax),
                eyxyxt=matrix(0,r*nmax,r*nmax),iter=integer(1),msg=integer(1),
                sigma=sigma,beta=beta,xtx=array(0,c(p,p,m)),
                xtw=matrix(0,p*r,nmax*r),xtwx=matrix(0,p*r,p*r),
                xtwy=numeric(p*r),xtwxinv=matrix(0,p*r,p*r),
                wkrr1=matrix(0,r,r),wkrr2=matrix(0,r,r),
                wkeyxyxt=matrix(0,r*nmax,r*nmax),
                wkqnm1=matrix(0,r*nmax,r*nmax),
                cvgd=integer(1),obeta=matrix(0,p,r),osigma=matrix(0,r,r),
                bigm=matrix(0,r,r),maxits=as.integer(maxits),
                llovec=numeric(as.integer(maxits)),
                epsil=epsil,sflag=sflag,eps=eps,wkpr=matrix(0,p,r),
                wkpp=matrix(0,p,p),xtxinv=matrix(0,p,p),PACKAGE="mlmmm")
        clock <- proc.time()-now
        cat("\n")
        iter <- tmp$iter
        msg <- tmp$msg
        {if(msg==1)
         warning("xtx is not full rank, failed for calculating beta <- (0)")
        else if(msg==4)
          warning("GLS failed for start vals, xtwx not full rank")
        else if(msg==6)
          warning("Value of sigma became non-pos.def. during iterations")}
        llovec <- tmp$llovec[1:iter]
        converged <- tmp$cvgd==as.integer(1)
        if(!converged) warning(paste("did not converge by",
           format(iter),"iterations"))
        #
        list(beta=tmp$beta,sigma=tmp$sigma,xtwxinv=tmp$xtwxinv,converged=converged,iter=iter,
		     npatt=npatt,pstfin=pstfin,iposn=iposn,patt=patt,rmat=rmat,
             logoll=llovec,clock=clock)}
