##Function written by Xiang Li (last update 5/18/11),  Saonli Basu (last update 7/18/12), Rob Kirkpatrick (last update May 2013).
###############################################################################
fgls <- function(
  fixed, data=parent.frame(), tlist=tlist, sizelist=sizelist, med=c("UN","VC"),
  vmat=NULL, start=NULL, theta=NULL, drop=NULL, get.hessian=FALSE, optim.method="BFGS",
  control=list(), weights=NULL, sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL){
  
  if(length(med)>1){med <- med[1]}
  if( !(med %in% c("UN","VC")) ){
    stop("Argument 'med' must be a character string, and either 'UN' or 'VC'.")
  }
  
  #Idiot checks pertaining to family sizes: ##################################
  if( !is.null(sizeLab) ){
    if( !is.character(sizeLab) ){
      warning("Value provided for 'sizeLab' is not a character string, and will be ignored.")
    }
    else{
      warning("Use of argument 'sizeLab' is optional, and may be eliminated in future package versions.")
      if(nchar(sizeLab)>max(sizelist)){
        warning("Value provided for 'sizeLab' specifies a larger maximum family size than appears in the data.")
      }
      if(nchar(sizeLab)<max(sizelist)){
        warning("Value provided for 'sizeLab' specifies a smaller maximum family size than appears in the data.")
      }
    }}
  if( !is.null(Mz) | !is.null(Bo) | !is.null(Ad) | !is.null(Mix) | !is.null(indobs) ){
    warning("Use of arguments 'Mz', 'Bo', 'Ad', 'Mix', and 'indobs' is deprecated; their values are ignored.")
  }
  #End of these idiot checks. #############################################
  
  #Create objects needed in the rest of the function: ########################################################
#   m <- match.call(expand.dots=FALSE)
#   temp <- c("", "data", "weights")
#   m <- m[ match(temp, names(m), nomatch=0)]
#   temp.fixed <- fixed
#   m$formula <- fixed
#   m[[1]] <- as.name("model.frame")
#   m <- eval(m, sys.parent())\
  call <- match.call()
  M <- model.frame(formula=fixed,data=data,na.action=NULL) #<--Use na.action=NULL here...
  na.rows <- NULL
  if(any(rowSums(is.na(M))>0)){ #<--...because we need to know if other objects need to be adjusted.
    na.rows <- which(rowSums(is.na(M))>0)
    weights <- weights[-na.rows]
  }
  M <- model.frame(formula=fixed,data=data,na.action=na.omit) #<--Redefine M.
  #Terms <- terms(fixed)
  X <- model.matrix(terms(M), data=data) #<--Rows with NA now omitted.
  coef.names <- dimnames(X)[[2]] #<--RMK May'13: Used later so the names of the regression coefficients correspond to variables' names.
  dimnames(X) <- NULL                            
  Y <- model.extract(M, "response") #<--Rows with NA now omitted.
  Y <- as.vector(Y)
  if("ID" %in% names(data)){id <-data$ID}
  else{id <- seq_len(length(Y))}
  n <- length(Y)
  #weights <- model.extract(m, 'weights')
  
  #RMK May'13: I decided to allow NAs in the variables provided in the formula
  #(there shouldn't be any when fgls() is called internally by gls.batch()).
  #Participants with NAs get trimmed out above,
  #and the corresponding rows and columns of the cov matrix get trimmed out later.
  
  #RMK May'13--generic start values:
  inivar <- summary(lm(Y~X,weights=weights))$sigma^2
  if(med=="VC"){defaultstart <- c(0.3*inivar,0.3*inivar,0.4*inivar)}
  else{defaultstart <- c(rep(0.46,8),rep(inivar,4))}
  #Done creating these objects. ##########################################################
  
  #Idiot checks pertaining to previously estimated residual covariance matrix: ###########################
  if( !is.null(vmat) & !is.null(theta) ){stop("At least one of arguments 'vmat' and 'theta' must be NULL.")}
  if( !is.null(theta) & (med=="UN" & length(theta)!=12) ){
    theta <- NULL
    warning("Non-NULL value to argument 'theta' must be a numerical vector of length 12 (NA's are accepted) when med='UN'; forcing 'theta' to NULL instead.")
  }
  if( !is.null(theta) & (med=="VC" & length(theta)!=3) ){
    theta <- NULL
    warning("Non-NULL value to argument 'theta' must be a numerical vector of length 3 (NA's are accepted) when med='VC'; forcing 'theta' to NULL instead.")
  }
  #Done with these idiot checks. ##########################################################
    
  ### produce the residual covariance matrix: ###############################
  if(is.null(vmat)){
    if( !is.null(theta) ){#<--RMK May'13: If parameters were supplied as theta.
      if(med=="VC"){blocks <- getblocks.ACE(theta=theta,tlist=tlist,sizelist=sizelist)}
      else{blocks <- getblocks(theta=theta,tlist=tlist,sizelist=sizelist,
                               pad=FALSE)} #<--RMK May'13: pad=F, matrix built from theta AS-IS.
      tkmat <- bdsmatrix(sizelist,blocks=blocks$blocks)
      estimates <- NA
      hessian.out <- iter <- NULL
    }
    else{
      #RMK May'13--Identify residual-covariance parameters that need to be dropped: ###################################
      if(med=="VC"){
        if(all(sizelist==1)){drop <- c(drop,1,2)} #drop A and C
        #(remember that C, as defined here, is identified by the spousal covariance alone)
        if( #When would A be unidentified, but C identified?
          !any(tlist %in% c("ccmf","bbmf","bamf","cmf","bmf","ccm","ccf","bam","baf")) & #When none of these are present,
            ( any(c("aa","ba","mf","af","am","aaf","aam","aamf") %in% tlist) + #and when fewer than two of these 3 sets are represented
                any(c("bb","bm","bf","cm","cf") %in% tlist) + ("cc" %in% tlist) < 2) ){drop <- c(drop,1)}
        if(length(drop)>0){
          drop <- as.integer(drop)
          drop <- sort(unique(drop[drop>0 & drop<4]))
        }
        if(any(drop==3)){drop <- drop[drop!=3]} #<--Parameter 3, "E", cannot be dropped.
        if(sum(drop)>=3){
          warning("Either no compatible family structure is identifiable from 'tlist', or too many parameters are supplied to 'drop'.")
          drop <- c(1:2)
        }
      }
      else{ #"UN"
        if( !any(tlist %in% c("ccmf","bbmf","aamf","bamf","cmf","bmf","amf","mf")) ){drop <- c(drop,1)} #r.mf
        if( !any(tlist %in% c("ccmf","bbmf","bamf","cmf","bmf","ccm","bbm","bam","cm","bm")) ){drop <- c(drop,2)}#r.bm
        if( !any(tlist %in% c("ccmf","bbmf","bamf","cmf","bmf","ccf","bbf","baf","cf","bf")) ){drop <- c(drop,3)}#r.bf
        if( !any(tlist %in% c("ccmf","ccm","ccf","cc")) ){drop <- c(drop,4)}#r.cc
        if( !any(tlist %in% c("bbmf","bbm","bbf","bb")) ){drop <- c(drop,5)}#r.bb
        if( !any(tlist %in% c("aamf","bamf","amf","aam","bam","am")) ){drop <- c(drop,6)}#r.am
        if( !any(tlist %in% c("aamf","bamf","amf","aaf","baf","af")) ){drop <- c(drop,7)}#r.am
        if( !any(tlist %in% c("aamf","bamf","aam","bam","aaf","baf","aa","ba")) ){drop <- c(drop,8)}#r.aa
        if( all(tlist %in% c("mf","m","f","INDPT")) ){drop <- c(drop,9)}#v.o
        if( all(tlist %in% c("ccf","bbf","aaf","baf","cf","bf","af","cc","bb","aa","ba","c","b","a","f","INDPT")) ){
          drop <- c(drop,10)#v.m
        }
        if( all(tlist %in% c("ccm","bbm","aam","bam","cm","bm","am","cc","bb","aa","ba","c","b","a","m","INDPT")) ){
          drop <- c(drop,11)#v.f
        }
        if( !any(tlist=="INDPT") ){drop <- c(drop,12)}
        if(length(drop)>0){
          drop <- as.integer(drop)
          drop <- sort(unique(drop[drop>0 & drop<13]))
        }
        if(sum(drop)>=78){
          warning("Either no compatible family structure is identifiable from 'tlist', or too many parameters are supplied to 'drop'.")
          drop <- c(1:11)
        }
      }
      #Done identifying drops. #########################################
      
      ########################### Optimization #################################
      #RMK May'13--first, check start values:
      if(is.null(start)){start <- defaultstart}
      if(med=="VC"){
        if(length(start)!=3){
          warning("start values must be a vector of length 3 (NA's are accepted) when med='VC'; using generic start values instead.")
          start <- defaultstart
        }
        start[ which(is.na(start) | (!is.na(start) & start<0)) ] <- 
          defaultstart[ which(is.na(start) | (!is.na(start) & start<0)) ]
      }
      else{ #"UN"
        if(length(start)!=12){
          warning("start values must be a vector of length 12 (NA's are accepted) when med='UN'; using generic start values instead.")
          start <- defaultstart
        }
        if( length(which( !is.na(start[1:8]) & abs(start[1:8])>1 ))>0 ){
          start[which( !is.na(start[1:8]) & abs(start[1:8])>1 )] <- 0.46
        }
        if( length(which( !is.na(start[9:12]) & start[9:12]<0 ))>0 ){
          start[8+which( !is.na(start[9:12]) & start[9:12]<0 )] <- inivar
        }
        if(any(is.na(start))){
          start[which(is.na(start))] <- defaultstart[which(is.na(start))]
        }
      }

      #RMK May'13--actual optimization:
      first_try <- TRUE
      nfit <- try(optim.logfun(start=start, drop=drop, med=med, force.PD=FALSE,
                               X=X, Y=Y, inivar=inivar, tlist=tlist, sizelist=sizelist, 
                               id=id, weights=weights, na.rows=na.rows, optim.method=optim.method, 
                               get.hessian=get.hessian, control=control),silent=TRUE)
      if(class(nfit)=="try-error"){ #<--RMK May'13: If initial optimization attempt fails.
        first_try <- FALSE
        print("Initial optimization attempt failed; attempting once more.")
        nfit <- optim.logfun(start=start, drop=drop, med=med, 
                             force.PD=TRUE, #<--RMK May'13: Force positive-definiteness on the second try only.
                             X=X, Y=Y, inivar=inivar, tlist=tlist, 
                             sizelist=sizelist, id=id, weights=weights, na.rows=na.rows,
                             optim.method=optim.method, get.hessian=get.hessian, control=control)
      }
      #RMK May'13--make some output objects:
      mssage <- ifelse(test=is.null(nfit$message),yes=" ",no=nfit$message)
      iter <- data.frame(iterations=as.vector(nfit$counts[1]),convergence=nfit$convergence,
                         message=mssage,first_try=first_try, stringsAsFactors=F)
      if(med=="VC"){
        estimates <- rep(NA,3)
        names(estimates) <- c("A","C","E")
      }
      else{
        estimates <- rep(NA,12)
        names(estimates) <- c("cor(m,f)","cor(c/b,m)","cor(c/b,f)","cor(c,c)","cor(b,b)","cor(a,m)","cor(a,f)","cor(a,a)","var(O)","var(m)","var(f)","var(ind)")
      }
      if(is.null(drop)){
        estimates[1:length(estimates)] <- as.vector(nfit$par) #Preserves previously-given names of estimates.
        print(estimates)
      }
      else{
        estimates[-drop] <- as.vector(nfit$par)
        print(estimates[-drop])
      }
      if(get.hessian==T){
        hessian.out <- nfit$hessian
        if(is.null(drop)){dimnames(hessian.out) <- list(names(estimates), names(estimates))}
        else{dimnames(hessian.out) <- list(names(estimates)[-drop], names(estimates)[-drop])}
        #print(hessian.out)
      } 
      else{hessian.out <- NULL}
      
      ########################### construct tkmat matrix from blocks ############
      if(med=="VC"){blocks <- getblocks.ACE(theta=estimates,tlist=tlist,sizelist=sizelist)}
      else{blocks <- getblocks(theta=estimates,force.PD=TRUE,tlist=tlist,sizelist=sizelist,pad=TRUE,inivar=inivar)}
      tkmat <- bdsmatrix(sizelist,blocks=blocks$blocks)
  }}
  #RMK May'13--if vmat was non-NULL:
  else{                     #Saonli's suggested change, 11/19/12
    if(is.character(vmat)){
      vmat.temp <- read.csv(vmat,header=T,colClasses="numeric") #<--Read in blocks from file.
      tkmat <- bdsmatrix(sizelist,vmat.temp[,1])
      rm(vmat.temp)
      estimates <- NA
      hessian.out <- iter <- NULL
    }
    else{
      tkmat <- vmat
      estimates <- NA
      hessian.out <- iter <- NULL
    }}
  ### done making residual covariance matrix ########################################################
  
  if(length(na.rows)>0){tkmat <- tkmat[-na.rows,-na.rows]}
  if(dim(tkmat)[1]!=length(Y)){stop("Dimensions of 'vmat' not equal to number of observations in data.")}
  list.vmat<-listbdsmatrix(tkmat,diag=T,id=F)
  vmat1<-sparseMatrix(list.vmat[,2],list.vmat[,1],x=list.vmat[,3],symmetric=T)
  vmat.Inv<-as(solve(vmat1,full=T),"sparseMatrix")
  vmat.Inv<-forceSymmetric(vmat.Inv)
  if(!is.null(weights)){ #<--RMK May'13: what to do if user supplies case weights.
    vmat.Inv <- as(vmat.Inv %*% diag(weights,nrow=dim(vmat.Inv)[1]),"sparseMatrix")
    vmat.Inv <- forceSymmetric(vmat.Inv)
  }
  gkmat<-as(chol(vmat.Inv),"sparseMatrix")
  
  Lambda<-1/diag(gkmat)
  newz <-gkmat%*%as.matrix(X)
  newy <- gkmat%*%Y
  lvd <- sum(log(Lambda))
  lfit <- lm(newy[,1]~0+as.matrix(newz))
  n <- length(newy[,1])
  loglik <- sum(lfit$residuals^2)/2 + lvd
  
  ls <- summary(lfit)
  rownames(ls$coefficients) <- coef.names
  fitted <- c(X %*% lfit$coef)  #fitted, on the original scale
  residuals <- Y - fitted
  coef.variance <- ls$cov.unscaled * ls$sigma^2
  dimnames(coef.variance) <- list(coef.names,coef.names)
  #RMK May'13--R^2:
  Ybar <- ( t(gkmat%*%matrix(1,nrow=n,ncol=1)) %*% (gkmat%*%matrix(Y,ncol=1)) /
    t(gkmat%*%matrix(1,nrow=n,ncol=1)) %*% (gkmat%*%matrix(1,nrow=n,ncol=1))
  )[1,1]
  Rsqd <- 1 - (
    ( t(gkmat%*%matrix(residuals,ncol=1)) %*% (gkmat%*%matrix(residuals,ncol=1)) ) /
      ( t(gkmat %*% matrix(Y-Ybar,ncol=1)) %*% (gkmat%*%matrix(Y-Ybar,ncol=1)) )
  )
  
  fcoef <- lfit$coef
  names(fcoef) <- coef.names
  call$fixed <- fixed
  fit <- list(
    ctable = ls$coefficients,    
    Rsqd=Rsqd[1,1],
    estimates = estimates,
    drop = drop,
    iter = iter,
    loglik = loglik,
    sigma = tkmat,
    hessian = hessian.out,
    n=n,
    df.residual = lfit$df.residual,
    residuals= residuals,
    fitted.values= fitted,
    #coefficients=fcoef,#list(fixed=fcoef),
    variance= coef.variance,
    call = call)
  #na.action <- attr(m, "na.action")
  #if (length(na.action)) fit$na.action <- na.action
  
  oldClass(fit) <- c('fgls')
  return(fit)
}
