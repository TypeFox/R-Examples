####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Likelihood.r
### Description:
### This file contains a set of procedures
### for maximum likelihood fitting of
### random fields.
### Last change: 28/03/2013.
####################################################

### Procedures are in alphabetical order.

### Optim call for log-likelihood maximization
Likelihood <- function(corrmodel,data,fixed,flagcor,flagnuis,grid,lower,model,namescorr,
                       namesnuis,namesparam,numcoord,numpairs,numparamcor,numrep,numtime,
                       optimizer,param,setup,spacetime,varest,taper,type,upper)
{
    ### START Defining the objective functions
    # Restricted log-likelihood for multivariate normal density:
    LogNormDenRestr <- function(const,cova,ident,dimat,nuisance,setup,stdata)
    {
        llik <- -1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['nugget']+nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        # Cholesky decomposition of the covariance matrix:
        cholvarcov <- try(chol(varcov),silent=TRUE)
        if(!is.matrix(cholvarcov)) return(llik)
        # Computes the determinat of the covariance matrix:
        detvarcov <- sum(log(diag(cholvarcov)))
        ivarcov <- chol2inv(cholvarcov)
        sumvarcov <- sum(ivarcov)
        p <- ivarcov-array(rowSums(ivarcov),c(dimat,1))%*%colSums(ivarcov)/sumvarcov
        llik <- -0.5*(const+2*detvarcov+log(sumvarcov)+crossprod(t(crossprod(stdata,p)),stdata))
        return(llik)
    }
        # Tapering log-likelihood for multivariate normal density:
    LogNormDenTap <- function(const,cova,ident,dimat,nuisance,setup,stdata)
    {
        lliktap <- -1.0e8
        # Computes the vector of the correlations:
        cova[cova==(nuisance['sill'])] <- nuisance['sill']+nuisance['nugget']
        varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
        cholvctap <- try(spam::chol.spam(varcovtap),silent=TRUE)
        if(class(cholvctap)=="try-error") return(lliktap)
        logdet <- c(determinant(cholvctap)$modulus)
        inv <- spam::solve.spam(cholvctap)
        slot(varcovtap,"entries") <- inv[setup$idx]*setup$taps
        lliktap= -0.5*(const+2*logdet+drop(t(stdata)%*%varcovtap%*%stdata))
        return(lliktap)
    }
    # Standard log-likelihood function for full normal density:
    LogNormDenStand <- function(const,cova,ident,dimat,nuisance,setup,stdata)
    {
        llik <- -1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['nugget']+nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        # Cholesky decomposition of the covariance matrix:
        cholvarcov <- try(chol(varcov),silent=TRUE)
        if(!is.matrix(cholvarcov)) return(llik)
        # Computes the determinat of the covariance matrix:
        detvarcov <- sum(log(diag(cholvarcov)))
        llik <- -0.5*(const+2*detvarcov+sum(stdata*backsolve(cholvarcov,
                 forwardsolve(cholvarcov,stdata,transpose=TRUE,upper.tri=TRUE))))
        return(llik)
    }
    ### END Defining the objective functions
    # Call to the objective functions:
    loglik <- function(const,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,model,namescorr,namesnuis,param,setup)
      {
        loglik <- -1.0e8
        # Set the parameter vector:
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        # Standardizes the data:
        stdata <- data-nuisance['mean']
        # Computes the vector of the correlations:
        cc=.C(corrmat,corr=corr,as.integer(corrmodel),as.double(nuisance),
           as.double(paramcorr),PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
        corr<-cc$corr
        if(corr[1]==-2) return(loglik)
        # Computes the correlation matrix:
        cova <- corr*nuisance['sill']
        # Computes the log-likelihood
        loglik <- sum(apply(stdata,1,fname,const=const,cova=cova,dimat=dimat,
                            ident=ident,nuisance=nuisance,setup=setup))
        return(loglik)
      }

          ### END Defining the objective functions
    # Call to the objective functions:
    loglik_tap_comp <- function(const,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,model,namescorr,namesnuis,param,setup)
      {
        loglik <- -1.0e8
        # Set the parameter vector:
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        # Standardizes the data:
        stdata <- data-nuisance['mean']
        # Computes the vector of the correlations:
        cc=.C(corrmat,corr=corr,as.integer(corrmodel),as.double(nuisance),
           as.double(paramcorr),PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
        corr<-cc$corr
        if(corr[1]==-2) return(loglik)
        # Computes the correlation matrix:
        cova <- corr*nuisance['sill']
        # Computes the log-likelihood

        loglik <-  LogNormDenTap(const=const,cova=cova,dimat=dimat,
                            ident=ident,nuisance=nuisance,setup=setup,stdata=stdata)
        sumloglik=sumloglik+loglik


        return(sumloglik)
      }

    ### START the main code of the function:
    dimat <- numcoord*numtime# set the data format
    numpairstot <- dimat*(dimat-1)/2
    const<-dimat*log(2*pi)# set the likelihood constant
    corr<-double(numpairstot)# initialize the correlation
    corrmat<-"CorrelationMat"# set the type of correlation matrix
    ident <- diag(dimat)# set the identity matrix
    # settings for space-time:
    if(spacetime){  data=matrix(c(data),numrep,dimat,byrow=T)
                    corrmat<-"CorrelationMat_st"}
    # detects the type of likelihood:
    if(type==3){
        fname<-"LogNormDenRestr"# set the name of object function
        const <- const-log(2*pi)}
    if(type==4) fname <- 'LogNormDenStand'
    if(type==5){
        corrmat<-"CorrelationMat_tap"
        if(spacetime) corrmat<-"CorrelationMat_st_tap"
        fname <- 'LogNormDenTap'
        corr <- double(numpairs)
        tapcorr <- double(numpairs)
        tapmod <- setup$tapmodel
        tp=.C(corrmat, tapcorr=tapcorr,as.integer(tapmod),as.double(c(0,0,1)),
           as.double(1),PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
        setup$taps<-tp$tapcorr
        }
    # Optimize the log-likelihood:
    if(optimizer=='L-BFGS-B')
      Likelihood <- optim(param,loglik,const=const,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(fnscale=-1,factr=1,
                          pgtol=1e-14,maxit=1e8),data=data,dimat=dimat,fixed=fixed,
                          fname=fname,grid=grid,ident=ident,lower=lower,
                          method=optimizer,model=model,namescorr=namescorr,
                          namesnuis=namesnuis,upper=upper,setup=setup)
    else
      Likelihood <- optim(param,loglik,const=const,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(fnscale=-1,factr=1,
                          pgtol=1e-14,maxit=1e8),data=data,dimat=dimat,
                          fixed=fixed,fname=fname,grid=grid,ident=ident,
                          method=optimizer,model=model,namescorr=namescorr,
                          namesnuis=namesnuis,setup=setup)
    # Set useful quantities:
    param<-Likelihood$par# set the parameter vector
    numparam<-length(param)# parameter size
    Likelihood$clic <- NULL
    if(type==3||type==4) Likelihood$clic <- -2*(Likelihood$value-numparam)# penalty in the AIC case
    ### Some checks of the output from the optimization procedure:
    if(Likelihood$convergence == 0)
      Likelihood$convergence <- 'Successful'
    else
      if(Likelihood$convergence == 1)
        Likelihood$convergence <- 'Iteration limit reached'
      else
        Likelihood$convergence <- "Optimization may have failed"
    ### START Computing the asymptotic variance-covariance matrices:
    if(varest){
        gnames <- namesparam
        if(flagnuis[1]) {# cheks if the mean is included
            numparam<-numparam-1
            gnames <- gnames[gnames!="mean"]}
        param<-c(param,fixed)# update with the fixed
        paramcorr<-param[namescorr]# correlation components
        numparamcorr<-length(paramcorr)# correlation size
        nuisance<-param[namesnuis]# nuisance components
        eps<-(.Machine$double.eps)^(1/3)
        numfish<-numparam*(numparam-1)/2+numparam# set variance-covariance matrix size
        # computing the off diagonal elements of the correlation matrix
        cc=.C(corrmat,corr=corr,as.integer(corrmodel),as.double(nuisance),as.double(paramcorr),
           PACKAGE='CompRandFld',DUP=TRUE,NAOK=TRUE)
        corr<-cc$corr
        varian<-corr*nuisance['sill']# computes covariance components
        if(!spacetime) dname<-"DCorrelationMat"
        else dname<-"DCorrelationMat_st"
        numparamcorr<-length(paramcorr[flagcor==1])# set the effective number of the corr param
        namescorr<-namescorr[flagcor==1]# set the effective names of the corr param
        # ML and REML cases:
        if(type==3||type==4){
            # Computing variance-covariance matrix of the random field:
            varcov<-(nuisance['nugget']+nuisance['sill'])*ident# computes variance components
            varcov[lower.tri(varcov)]<-varian
            varcov<-t(varcov)
            varcov[lower.tri(varcov)]<-varian
            cholcov<-chol(varcov)
            invar<-chol2inv(cholcov) # inverse of varcov matrix
            fish<-double(numfish)
            # Restricted likelihood case
            if(type==3) P<-invar-array(rowSums(invar),c(dimat,1))%*%colSums(invar)/sum(invar)
            # set array of gradient matrices
            gradient<-array(0,dim=c(dimat,numparam,dimat))# vector derivatives
            colnames(gradient) <- gnames
            dcorr <- double(numpairstot*numparamcorr)
            # correlation gradient vector
            dc=.C(dname,as.integer(corrmodel),dcorr=dcorr,as.double(eps),
               as.integer(flagcor),as.integer(numparamcorr),
               as.double(paramcorr),corr,PACKAGE='CompRandFld',
               DUP=TRUE,NAOK=TRUE)
            dcorr<-dc$dcorr
            dim(dcorr) <- c(numpairstot,numparamcorr)
            # Computing the gradient matrices:
            if(flagnuis[2]) gradient[,namesnuis[2],] <- ident
            # gradient matrix, derivatives with respect to the sill (correlation matrix)
            if(flagnuis[3]){
                R<-ident
                R[lower.tri(R)]<-corr
                R<-t(R)
                R[lower.tri(R)]<-corr
                gradient[,namesnuis[3],] <- R}
            # gradient matrix, derivatives with respect to the correlation parameters
            for(i in 1:numparamcorr){
                grad<-matrix(0,dimat,dimat)# set gradient matrix
                grad[lower.tri(grad)]<-dcorr[,i]
                grad<-t(grad)
                grad[lower.tri(grad)]<-dcorr[,i]
                if(flagcor[namescorr][i])  gradient[,namescorr[i],] <- grad  }
            i<-1
            # Computing the gradient matrices:
            k<-1
            # Computing off diagonal elements of the asymptotic Fisher/Godambe information matrix:
            for(i in 1:numparam)
                for(j in i:numparam){
                    if(type==3) fish[k]<-sum(diag(P%*%gradient[,i,]%*%P%*%gradient[,j,]))/2# REML case
                    if(type==4) fish[k]<-sum(diag(invar%*%gradient[,i,]%*%invar%*%gradient[,j,]))/2# ML case
                    k<-k+1}
            # Building Fisher/Godambe information matrix:
            fisher<-diag(0,numparam)# full and restricted likelihood cases
            fisher[lower.tri(fisher,diag=TRUE)] <- fish
            fisher<- t(fisher)
            fisher[lower.tri(fisher,diag=TRUE)] <- fish
            # Adding mean to the asymptotic Fisher/Godambe information matrix:
            if(flagnuis[1]){
                if(type==4) fishmean<-sum(invar)
                zeros<-rep(0,numparam)
                fisher<-rbind(c(fishmean,zeros),cbind(zeros,fisher))}
            invfisher <- try(solve(fisher),silent=TRUE)
            if(!is.matrix(invfisher)) invfisher<-NULL
            Likelihood$sensmat <- NULL
            Likelihood$varimat <- NULL
        }
        # Computing of the asymptotic variance covariance matrix: case TAP
        if(type==5){
            varcov <- varian
            # define the variance-covariance vector
            varcov[varcov==(nuisance['sill'])] <- nuisance['nugget']+nuisance['sill']
            # define the sparse variance-covariance matrix
            spamvar <- new("spam",entries=varcov,colindices=setup$ja,rowpointers=setup$ia,
                           dimension=as.integer(rep(dimat,2)))
            varcov <- as.matrix(spamvar)
            # define tha variance-covariance tapered matrix
            covtap <- new("spam",entries=varcov[setup$idx]*setup$taps,colindices=setup$ja,
                          rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
            cholcovtap <- chol(covtap) # Chol. decomp. of the covariance tapered matrix
            invtap<-as.matrix(spam::solve.spam(cholcovtap)) # inverse of the Chol. decomp. cov tap
            covtap <- as.matrix(covtap)
            HH <- double(numfish)# upper triang matrix (with diag) of sensitivity matrix
            JJ <- double(numfish)# upper triang matrix (with diag) of variability matrix
            # computing the matrix of derivatives
            gradient <- array(0,c(numpairs,numparam))# vector derivatives
            colnames(gradient) <- gnames
            if(!spacetime) dname<-"DCorrelationMat_tap"
            else dname<-"DCorrelationMat_st_tap"
            dcorr <- double(numpairs*numparamcorr)
            # correlation gradient vector
            dc=.C(dname,as.integer(corrmodel),dcorr=dcorr,as.double(eps),
               as.integer(flagcor),as.integer(numparamcorr),
               as.double(paramcorr),corr,PACKAGE='CompRandFld',
               DUP=TRUE,NAOK=TRUE)
            dcorr<-dc$dcorr
            dim(dcorr) <- c(numpairs,numparamcorr)
            if(flagnuis[2]) gradient[,namesnuis[2]] <- ident[setup$idx]
            if(flagnuis[3]) gradient[,namesnuis[3]] <- corr
            gradient[,namescorr] <- dcorr
            # computing matrix derivatives
            k <- 1
            # computing uppper triangular of the Fisher matrix
            H<-diag(0,numparam)# sensitivity matrix
            J<-diag(0,numparam)# variability matrix
            for(i in 1:numparam)
                for(j in i:numparam){
                    gradtapI <- gradtapJ <- bigI <- bigJ <- spamvar
                    slot(gradtapI,"entries") <- gradient[,i]*setup$taps
                    slot(gradtapJ,"entries") <- gradient[,j]*setup$taps
                    HH[k] <- -0.5*sum(diag(gradtapI%*%invtap%*%gradtapJ%*%invtap))
                    slot(bigI, "entries") <- (invtap%*%gradtapI%*%invtap)[setup$idx]*setup$taps
                    slot(bigJ, "entries") <- (invtap%*%gradtapJ%*%invtap)[setup$idx]*setup$taps
                    JJ[k] <- 0.5*sum(diag(bigI%*%varcov%*%bigJ%*%varcov))
                    k <- k+1}
            # Building Godambe information  matrix
            H[lower.tri(H,diag=TRUE)] <- HH
            H<- t(H)
            H[lower.tri(H,diag=TRUE)] <- HH
            J[lower.tri(J,diag=TRUE)] <- JJ
            J<- t(J)
            J[lower.tri(J,diag=TRUE)] <- JJ
            # Adding the mean parameter
            if(flagnuis[1]){
                slot(spamvar,"entries") <- invtap[setup$idx]*setup$taps
                zeros <- rep(0,numparam)
                fishmH <- sum(spamvar)
                H <- rbind(c(fishmH,zeros),cbind(zeros,H))
                fishmJ <- sum(spamvar%*%varcov%*%spamvar)
                J <- rbind(c(fishmJ,zeros),cbind(zeros,J))}
            invH <- try(solve(H),silent = TRUE)
              Likelihood$sensmat <- H
              Likelihood$varimat <- J
            if(!is.matrix(invH) || !is.matrix(H)){
                invfisher<-NULL
                Likelihood$clic <- NULL}
            else{
                invfisher<-invH%*%J%*%invH
                Likelihood$clic <- -2*(Likelihood$value-sum(diag(J%*%invH)))# penalty in the TIC case
                }
        }
        ### END Computing the asymptotic variance-covariance matrices
        Likelihood$varcov <- invfisher

        #Checks if the resulting variance and covariance matrix:
        if(is.null(Likelihood$varcov)){
            warning("Asymptotic information matrix is singular")
            Likelihood$varcov <- 'none'
            Likelihood$stderr <- 'none'}
        else{
            dimnames(Likelihood$varcov)<-list(namesparam,namesparam)
            Likelihood$stderr<-diag(Likelihood$varcov)
        if(any(Likelihood$stderr < 0)) Likelihood$stderr <- 'none'
        else{
            Likelihood$stderr<-sqrt(Likelihood$stderr)
            names(Likelihood$stderr)<-namesparam}}}
    ### END the main code of the function:
    return(Likelihood)
  }
