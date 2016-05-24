####################################################
### Authors: Simone Padoan and Moreno Bevilacqua.
### Emails: simone.padoan@unibocconi.it,
### moreno.bevilacqua@uv.cl
### Institutions: Department of Decision Sciences,
### University Bocconi of Milan and
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: CompositeLikelihood.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 28/03/2013.
####################################################


### Procedures are in alphabetical order.

### Optim call for Composite log-likelihood maximization

CompLikelihood <- function(coordx, coordy, corrmodel, data, distance, flagcorr, flagnuis, fixed, grid,
                           likelihood, lower, model, namescorr, namesnuis, namesparam,
                           numparam, numparamcorr, optimizer, param, spacetime, threshold, type,
                           upper, varest, vartype, winconst, winstp)
  {
    ### Define the object function:
    comploglik <- function(corrmodel, data, fixed, fun, namescorr, namesnuis, param, threshold)
      {
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        result <- .C(fun, as.integer(corrmodel), as.double(data), as.double(nuisance),
                     as.double(paramcorr), as.double(threshold), res=double(1),
                     PACKAGE='CompRandFld', DUP=TRUE , NAOK=TRUE)$res
        return(result)
      }
    fname <- NULL
    hessian <- FALSE
    if(all(model==1,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss'
    if(all(model==1,likelihood==3,type==1)) fname <- 'Comp_Diff_Gauss'
    if(all(model==1,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss'
    if(all(model==2,likelihood==1,type==2)) fname <- 'Comp_Cond_BinGauss'
    if(all(model==2,likelihood==3,type==1)) fname <- 'Comp_Diff_BinGauss'
    if(all(model==2,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinGauss'
                                             if(varest & vartype==2) hessian <- TRUE}
    if(all(model==3,likelihood==3,type==2)){ fname <- 'Comp_Brow_Resn'
                                             if(varest & vartype==2) hessian <- TRUE}
    if(all(model==4,likelihood==3,type==2)){ fname <- 'Comp_Ext_Gauss'
                                             if(varest & vartype==2) hessian <- TRUE}
    if(all(model==5,likelihood==3,type==2)){ fname <- 'Comp_Ext_T'
                                             if(varest & vartype==2) hessian <- TRUE}

    if(spacetime) fname <- paste(fname,"_st",sep="")
    if(optimizer=='L-BFGS-B')
      CompLikelihood <- optim(param, comploglik, corrmodel=corrmodel, control=list(fnscale=-1,
                              factr=1, pgtol=1e-14, maxit=1e8), data=data, fixed=fixed,
                              fun=fname, hessian=hessian, lower=lower, method=optimizer,
                              namescorr=namescorr, namesnuis=namesnuis, threshold=threshold,
                              upper=upper)
    else
      CompLikelihood <- optim(param, comploglik, corrmodel=corrmodel, control=list(fnscale=-1,
                              reltol=1e-14, maxit=1e8), data=data, fixed=fixed, fun=fname,
                              hessian=hessian, method=optimizer, threshold=threshold,
                              namescorr=namescorr, namesnuis=namesnuis)
    # check the optimisation outcome
    if(CompLikelihood$convergence == 0)
      CompLikelihood$convergence <- 'Successful'
    else
      if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
      else
        CompLikelihood$convergence <- "Optimization may have failed"

    # Compute the nugget in the binary case:
    if(model==2) fixed["nugget"] <- 1-CompLikelihood$par["sill"]
        ### Computation of the variance-covariance matrix:
        if(varest)
          {
            # The sensitivity (H) and the variability (J) matrices
            # are estimated by the sample estimators contro-parts:
            dimmat <- numparam^2
            dmat <- numparam*(numparam+1)/2
            eps <- (.Machine$double.eps)^(1/3)
            param <- c(CompLikelihood$par, fixed)
            score <- double(numparam)
            paramcorr <- param[namescorr]
            nuisance <- param[namesnuis]
            if(hessian) sensmat <- -as.double(CompLikelihood$hessian[lower.tri(CompLikelihood$hessian, diag=TRUE)])
            else sensmat <- double(dmat)
            varimat <- double(dmat)
            # Set the window parameter:
            GOD=.C('GodambeMat',as.double(coordx),as.double(coordy),as.integer(corrmodel),
               as.double(data),as.integer(distance),as.double(eps),as.integer(flagcorr),
               as.integer(flagnuis),as.integer(grid),as.integer(likelihood),as.integer(model),
               as.integer(numparam),as.integer(numparamcorr),as.double(paramcorr),as.double(nuisance),
               score=score,sensmat=sensmat,as.integer(spacetime),as.double(threshold),as.integer(type),
               varimat=varimat,as.integer(vartype),as.double(winconst),as.double(winstp),
               PACKAGE='CompRandFld', DUP=TRUE, NAOK=TRUE)
            score<-GOD$score
            sensmat<-GOD$sensmat
            varimat<-GOD$varimat
            # Set score vectore:
            CompLikelihood$winconst<-winconst
            CompLikelihood$winstp<-winstp
            CompLikelihood$score <- score
            # Set sensitivity matrix:
            CompLikelihood$sensmat <- matrix(rep(0,dimmat),ncol=numparam)
            # Set variability matrix:
            CompLikelihood$varimat <- matrix(rep(0,dimmat),ncol=numparam)
            namesgod <- c(namesnuis[as.logical(flagnuis)], namescorr[as.logical(flagcorr)])
            names(CompLikelihood$score) <- namesgod
            CompLikelihood$score <- CompLikelihood$score[namesparam]
            dimnames(CompLikelihood$sensmat) <- list(namesgod, namesgod)
            dimnames(CompLikelihood$varimat) <- list(namesgod, namesgod)
            if(numparam>1){
              CompLikelihood$sensmat[lower.tri(CompLikelihood$sensmat, diag=TRUE)] <- sensmat
              CompLikelihood$sensmat <- t(CompLikelihood$sensmat)
              CompLikelihood$sensmat[lower.tri(CompLikelihood$sensmat, diag=TRUE)] <- sensmat
              CompLikelihood$varimat[lower.tri(CompLikelihood$varimat, diag=TRUE)] <- varimat
              CompLikelihood$varimat <- t(CompLikelihood$varimat)
              CompLikelihood$varimat[lower.tri(CompLikelihood$varimat, diag=TRUE)] <- varimat
              CompLikelihood$sensmat <- CompLikelihood$sensmat[namesparam, namesparam]
              CompLikelihood$varimat <- CompLikelihood$varimat[namesparam, namesparam]}
            else {CompLikelihood$sensmat[1,1] <- sensmat
                  CompLikelihood$varimat[1,1] <- varimat}

            isensmat <- try(solve(CompLikelihood$sensmat), silent = TRUE)
            if(!is.matrix(isensmat) || !is.matrix(CompLikelihood$varimat))
              {
                warning("observed information matrix is singular")
                CompLikelihood$varcov <- 'none'
                CompLikelihood$stderr <- 'none'
              }
            else
              {
                penalty <- CompLikelihood$varimat %*% isensmat
                CompLikelihood$clic <- -2 * (CompLikelihood$value - sum(diag(penalty)))
                CompLikelihood$varcov <- isensmat %*% penalty
                CompLikelihood$stderr <- diag(CompLikelihood$varcov)
                if(any(CompLikelihood$stderr < 0))
                  CompLikelihood$stderr <- 'none'
                else
                  CompLikelihood$stderr <- sqrt(CompLikelihood$stderr)
              }
          }
    return(CompLikelihood)
  }

