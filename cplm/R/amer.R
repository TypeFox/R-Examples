#######################################################
# These functions are copied from the retired R
# package 'amer', written by Fabian Scheipl (fabian.scheipl@googlemail.com)
#######################################################
# 
# I changed the original expand.call in amer to remove the use of .Internal
expand.call <-
  function(call = sys.call(sys.parent(1)))
  {
    #given args:
    ans <- as.list(call)
    # ans1 <- ans[[1]]
    # ans <- lapply(ans[-1], eval, envir = sys.frame(sys.parent(2)))
    # ans <- c(ans1, ans)
    
    #possible args:
    frmls <- formals(safeDeparse(ans[[1]]))
    #remove formal args with no presets:
    frmls <- frmls[!sapply(frmls, is.symbol)]
    
    add <- which(!(names(frmls) %in% names(ans)))
    return(as.call(c(ans, frmls[add])))
  }


indsF <-
  function(m, fct, fctterm){
    # add assign-like info to fctterm:
    # which penalization/ranef groups and coefficients (fixed/random) belong to which function 
    # also include info on global intercept and by-level intercepts  
    ranefinds <- reinds(m@Gp)
    
    indIntercept <- ifelse("(Intercept)" %in% names(fixef(m)), 1, 0)
    
    for(i in 1:length(fctterm)){
      if(length(fct[[i]]$Z) == 1){
        
        attr(fctterm[[i]], "indGrp") <- match(names(fct)[i], colnames(m@flist)) 
        if(eval(attr(fct[[i]], "call")$allPen)) {
          #add pen. group(s) with grouping factor u.x.by
          indUGrp <- match(sub("f.", "u.", names(fct)[i]), colnames(m@flist))
          attr(fctterm[[i]], "indGrp") <-  c(attr(fctterm[[i]], "indGrp"), which(attr(m@flist, "assign")==indUGrp))
        }
        attr(fctterm[[i]], "indPen") <- unlist(ranefinds[attr(fctterm[[i]], "indGrp")])
        
        if(!(eval(attr(fct[[i]], "call")$allPen)||ncol(fct[[i]]$X[[1]])==0)){
          attr(fctterm[[i]], "indUnpen") <-  sapply(paste("^",colnames(fct[[i]]$X[[1]]),"$",sep=""),
                                                    grep, x=names(m@fixef))
          names(attr(fctterm[[i]], "indUnpen")) <- colnames(fct[[i]]$X[[1]])
        } else attr(fctterm[[i]], "indUnpen") <- 0
        
        attr(fctterm[[i]], "indConst") <- indIntercept
        
        attr(fctterm[[i]], "indGrp") <- list(attr(fctterm[[i]], "indGrp"))
        attr(fctterm[[i]], "indPen") <- list(attr(fctterm[[i]], "indPen"))
        attr(fctterm[[i]], "indUnpen") <- list(attr(fctterm[[i]], "indUnpen"))
        attr(fctterm[[i]], "indConst") <- list(attr(fctterm[[i]], "indConst"))
      } else {
        by <- eval(attr(fct[[i]],"call")$by, m@frame)
        attr(fctterm[[i]], "indGrp") <- vector(mode="list", length=nlevels(by))
        attr(fctterm[[i]], "indPen") <-	vector(mode="list", length=nlevels(by))
        attr(fctterm[[i]], "indUnpen") <- vector(mode="list", length=nlevels(by))
        attr(fctterm[[i]], "indConst") <- vector(mode="list", length=nlevels(by))
        for(j in 1:nlevels(by)){
          attr(fctterm[[i]], "indGrp")[[j]] <- grep(paste("^",paste(names(fct)[i],".",names(fct[[i]]$Z)[j],sep=""), "$", sep=""), colnames(m@flist))
          attr(fctterm[[i]], "indPen")[[j]] <- ranefinds[[attr(fctterm[[i]], "indGrp")[[j]]]]
          if( ncol(fct[[i]]$X[[j]]) == 0){
            attr(fctterm[[i]], "indUnpen")[[j]] <-	 0
          } else {
            attr(fctterm[[i]], "indUnpen")[[j]] <- sapply(
              paste("^",colnames(fct[[i]]$X[[j]]),"$",sep=""), grep, x=names(m@fixef))
            names(attr(fctterm[[i]], "indUnpen")[[j]]) <- colnames(fct[[i]]$X[[j]])
          }	
          #add by-level intercept:
          indBy <- grep(paste("^",safeDeparse(attr(fct[[i]],"call")$by), levels(by)[j],"$", sep=""), names(m@fixef))
          indBy <- indBy[!(indBy %in% attr(fctterm[[i]], "indUnpen")[[j]])]
          attr(fctterm[[i]], "indConst")[[j]] <- c(indIntercept, indBy)
        }
      }
    }
    return(fctterm)
  }


safeDeparse <- function(expr){
  ret <- paste(deparse(expr), collapse="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}


expandBasis <-
  function(basis, by, varying, bySetToZero = T){
    #multiply X, Z with varying
    #set up lists of design matrices if by-variable and/or allPen are given:   
    allPen <- eval(attr(basis, "call")$allPen, parent.frame())
    X <- basis$X
    Z <- basis$Z
    
    
    xName <- safeDeparse(attr(basis, "call")$x)
    if(!is.null(varying)){
      xName <- paste(xName, "X", safeDeparse(attr(basis, "call")$varying), sep="")
      X <- cbind(varying, X * varying)
      Z <- Z * varying
    }	
    
    
    if(!is.null(by)){
      X.o <- X
      Z.o <- Z
      byName <- safeDeparse(attr(basis, "call")$by)
      if(!allPen){
        basis$X <- basis$Z <- vector(mode="list", nlevels(by))
        for(i in 1:nlevels(by)){
          if(bySetToZero){
            keep <- 1*(by == levels(by)[i])
          } else keep <- rep(1, length(by))	
          #set X,Z partially to zero for each level of by-variable
          if(NCOL(X.o)){ 
            basis$X[[i]] <- X.o * keep
            #naming scheme: x.fx.bylevel1.fx1, x.fx.bylevel1.fx2, ...,x.fx.bylevel2.fx1, or xXvarying.fx.bylevel1.fx1
            
            colnames(basis$X[[i]]) <- paste(xName,".",byName,levels(by)[i],
                                            paste(".fx",1:NCOL(basis$X[[i]]),sep=""), sep="")
          }
          basis$Z[[i]] <- Z.o * keep
        }
        #naming scheme: bylevel1, bylevel2, ...
        names(basis$X) <- names(basis$Z) <- paste(byName, levels(by), sep="")
      } else {
        basis$X <- basis$Z <- vector(mode="list", 1)
        by <- C(by[, drop=TRUE], contr.treatment) #make sure treatment contrasts are used, unused levels dropped
        #basis$Z[[1]] <- Matrix(model.matrix(~ 0 + Z.o:by)) #TODO: can this be done without intermediate dense matrix?
        basis$Z[[1]] <- model.matrix(~ 0 + Z.o:by) #FIXME: ?constructing directly as sparse breaks transposing Z in subAZ?
        #cbind Z set partially to zero for each level of by-variable:
        ## for(i in 1:nlevels(by[, drop=TRUE])) {
        ##     if(bySetToZero){
        ##         keep <- (by == levels(by[, drop=TRUE])[i])
        ##     } else keep <- 1
        ##     basis$Z[[1]] <- cBind(basis$Z[[1]], Z.o * keep)
        ##     #basis$Z[[1]] <- Matrix(model.matrix(~ 0 + Z.o:by[, drop=TRUE]))
        ## }	
        basis$X[[1]] <- X.o
        if(NCOL(X.o)) colnames(basis$X[[1]]) <- paste(xName,".",byName,paste(".fx",1:NCOL(X),sep=""), sep="")
        #naming scheme: u.x.by (also: name of duplicated by-variable in expandMf, subFcts)
        names(basis$X) <- paste("u", xName, byName, sep=".")
      }
    } else {
      if(NCOL(X)) colnames(X) <- paste(xName,paste(".fx",1:NCOL(X),sep=""), sep="")
      basis$X <- list(X)
      basis$Z <- list(Z)
    }  
    
    return(basis)
  }



subFcts <-
  function(rhs, fctterm, fct, fr)
    # replace formula parts for smooth functions with  xi + (xi^2+ )... + xi^dimUnpen + (1|fcti) or
    # by*(xi + xi^2+ ... + xi^dimUnpen) + (1|fcti.1) + ... + (1|fcti.N) for by-variable with N levels  
  {
    for(i in 1:length(fct)){
      by <- eval(attr(fct[[i]],"call")$by, fr)
      allPen <- eval(attr(fct[[i]],"call")$allPen)
      diag <- eval(attr(fct[[i]],"call")$diag)
      
      replacement <-  
        if(is.null(by)){
          # 1 + x.fx1 + x.fx2+ ... + (1|f.x)
          paste(ifelse(ncol(fct[[i]]$X[[1]])!=0,
                       paste(as.vector(sapply(fct[[i]]$X,colnames)),collapse=" + "),
                       "1"),
                " + (1|",names(fct)[i],")",sep="")		
        } else {
          if(allPen){ 
            if(!diag){
              # add correlated random effects for normally unpenalized part of basis grouped according to by and fake random intercept
              # (1 + x.fx1 + x.fx2+ ...|u.x.by)  + (1|f.x.by)
              paste(
                paste(paste("(1",
                            paste(as.vector(sapply(fct[[i]]$X,colnames)), collapse="+"),
                            sep="+"),
                      "|", 
                      names(fct[[i]]$X),")",
                      sep=""),
                paste("(1|",names(fct)[i],")",sep="", collapse=" + "),
                sep =" + ")
            } else {
              # add independent random effects for normally unpenalized part of basis grouped according to by and fake random intercept
              # (1|u.x.by) + x.fx1|u.x.by) + x.fx2|u.x.by) + ...  + (1|f.x.by)
              paste(
                paste(c("(1", paste("(0+", as.vector(sapply(fct[[i]]$X,colnames)),sep="")),"|", names(fct[[i]]$X),")",sep="",collapse=" + "),
                paste("(1|",names(fct)[i],")",sep="", collapse=" + "),
                sep =" + ")
            }
            
          } else { 
            # add fixed effect for unpenalized part of basis + fake random intercept for each by-level
            # by + x.fx1.BYlevel1 + x.fx2.BYlevel1 +...+ (1|f.x.BYlevel1) + ... + x.fx1.BYlevelD + x.fx2.BYlevelD +... + (1|f.x.BYlevelD)
            paste(#deparse(attr(fct[[i]],"call")$by), 
              ifelse(ncol(fct[[i]]$X[[1]])!=0,
                     paste(as.vector(sapply(fct[[i]]$X,colnames)),collapse=" + "),
                     "1"),
              paste("(1|",names(fct)[i],".",names(fct[[i]]$Z),")",sep="", collapse=" + "),
              sep =" + ")
          }
        }
      rhs <- sub(safeDeparse(fctterm[[i]]), replacement, rhs, fixed=T)
    }
    return(rhs)
  }


expandMf <-
  function(fr, fct)
    # cbind model frame with design matrices for unpenalized&penalized parts of the smooth fcts.  
  {
    for(i in 1:length(fct)){
      #matrix with all unpenalized terms for fct
      newX <- do.call(cBind, fct[[i]]$X)
      
      #factor variables with no. of levels = no. of penalized basis fcts  
      #newFact <-   replicate(length(fct[[i]]$Z), rep(1:ncol(fct[[i]]$Z[[1]]), length=nrow(fct[[i]]$Z[[1]])))
      newFact <- data.frame(factor(rep(1:ncol(fct[[i]]$Z[[1]]), length=nrow(fct[[i]]$Z[[1]]))))
      if(length(fct[[i]]$Z) > 1){
        for(j in 2:length(fct[[i]]$Z)){
          newFact <- cbind(newFact, factor(rep(1:ncol(fct[[i]]$Z[[1]]), length=nrow(fct[[i]]$Z[[1]]))))
        }
      }
      
      colnames(newFact) <- if(length(fct[[i]]$Z) == 1){ 
        names(fct)[i]
      } else {
        paste(names(fct)[i],".",names(fct[[i]]$Z),sep="")
      }
      
      if(eval(attr(fct[[i]],"call")$allPen)){
        # duplicate grouping factor for allPen-function groups so that assignment (which entries in ranef belong to which penalization 
        # group) can be reconstructed from the fitted model object m if there is another random effect associated with the by-variable.
        # will need this for predict etc.. since attr(m@flist,"assign") only works the other way around....
        newFact <- cBind(newFact, eval(attr(fct[[i]],"call")$by, fr))
        colnames(newFact)[ncol(newFact)] <- names(fct[[i]]$X)
      } 
      
      fr <- cBind(cBind(fr, newX),newFact)
    }
    return(fr)
  }


subAZ <-
  function(m, fct)
    # replace design matrices for fake factors with designs for penalized spline basis
  {
    for(i in 1:length(fct)){
      if(length(fct[[i]]$Z) == 1){
        ind <- which(names(m$FL$fl)[attr(m$FL$fl, "assign")]==names(fct[i])) 
        Zt <- as(t(fct[[i]]$Z[[1]]), "sparseMatrix")
        m$FL$trms[[ind]]$A <- m$FL$trms[[ind]]$Zt <- Zt
        dimnames(m$FL$trms[[ind]]$ST) <- list(safeDeparse(attr(fct[[i]], "call")[[1]]), safeDeparse(attr(fct[[i]], "call")[[1]])) 
      } else {
        for(j in 1:length(fct[[i]]$Z)){
          ind <- grep(paste("^",paste(names(fct[i]),names(fct[[i]]$Z)[j],sep='.'), "$", sep=""), names(m$FL$fl)[attr(m$FL$fl, "assign")])
          Zt <- as(t(fct[[i]]$Z[[j]]), "sparseMatrix")
          m$FL$trms[[ind]]$A <- m$FL$trms[[ind]]$Zt <- Zt
          dimnames(m$FL$trms[[ind]]$ST) <- list(safeDeparse(attr(fct[[i]], "call")[[1]]), safeDeparse(attr(fct[[i]], "call")[[1]]))
        }
      }
    }
    return(m)
  }


set.mfrow <-
  function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE)
    # taken from coda 
  {
    mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
      if (Nchains == 2) {
        switch(min(Nparms, 5), c(1, 2), c(2, 2), c(3, 2), 
               c(4, 2), c(3, 2))
      }
      else if (Nchains == 3) {
        switch(min(Nparms, 5), c(2, 2), c(2, 3), c(3, 3), 
               c(2, 3), c(3, 3))
      }
      else if (Nchains == 4) {
        if (Nparms == 1) 
          c(2, 2)
        else c(4, 2)
      }
      else if (any(Nchains == c(5, 6, 10, 11, 12))) 
        c(3, 2)
      else if (any(Nchains == c(7, 8, 9)) || Nchains >= 13) 
        c(3, 3)
    }
    else {
      if (nplots == 1) {
        mfrow <- switch(min(Nparms, 13), c(1, 1), c(1, 2), 
                        c(2, 2), c(2, 2), c(3, 2), c(3, 2), c(3, 3), 
                        c(3, 3), c(3, 3), c(3, 2), c(3, 2), c(3, 2), 
                        c(3, 3))
      }
      else {
        mfrow <- switch(min(Nparms, 13), c(1, 2), c(2, 2), 
                        c(3, 2), c(4, 2), c(3, 2), c(3, 2), c(4, 2), 
                        c(4, 2), c(4, 2), c(3, 2), c(3, 2), c(3, 2), 
                        c(4, 2))
      }
    }
    return(mfrow)
  }


tp <-
  function(x, degree=1, k = 15, by=NULL, allPen = FALSE, varying = NULL, diag=FALSE,
           knots= quantile(x, probs = (1:(k - degree))/(k - degree  + 1)), centerscale=NULL, scaledknots=FALSE)
  {
    call <- as.list(expand.call(match.call()))
    
    stopifnot(is.numeric(x), is.factor(by)||is.null(by), is.numeric(varying)||is.null(varying), degree >= 0)
    
    degree <- as.integer(degree); call$degree <- degree
    
    knots <- eval(knots)
    if(is.null(centerscale)){
      x <- scale(x)
      #make sure prediction data uses the same center/scale as the data used to fit the model:
      call$centerscale <- c(attr(x, "scaled:center"),attr(x, "scaled:scale"))
      x <- as.vector(x)
    } else x <- (x - centerscale[1])/centerscale[2]
    if(!scaledknots){
      knots <- (knots - call$centerscale[1])/call$centerscale[2]
      call$scaledknots <- TRUE
    }	
    
    if(length(unique(knots))!=length(knots)) warning("duplicate knots detected and removed.")
    knots <- sort(unique(knots))
    
    call$knots <- knots
    if(k != length(knots)+ degree){
      k <- length(knots) + degree; call$k <- k
      warning("set k to ", k," to conform with given knots and degree.")
    }
    if((knots[1]<min(x)||(knots[k-degree]>max(x)))) warning("knots outside range of variable.")
    
    if(is.null(by) && allPen) stop("allPen = TRUE only makes sense for smooths with a by-variable.")
    
    #design for unpenalised part: global polynomial trends (no intercept)
    if(degree>0){	
      X <- outer(x, 1:degree, "^")#poly(x, degree)
      #colnames(X) <- paste(as.character(call$x),".fx",1:NCOL(X),sep="")
    } else{
      X <- matrix(nrow=length(x), ncol=0)
    }	
    #design for penalised part: 
    Z <- outer(x,knots,"-")^degree*outer(x,knots,">")
    
    
    res <- list(X=X, Z=Z, knots=knots)
    attr(res, "call") <- as.call(call)
    
    return(res)
  }

bsp <- function(x, k=15, spline.degree = 3, diff.ord = 2, 
                knots=NULL, by=NULL, allPen = FALSE, varying = NULL, diag=FALSE)
{
  
  call <- as.list(expand.call(match.call()))
  stopifnot(diff.ord>=0, spline.degree>=0, is.numeric(x), k > spline.degree)
  
  if(is.null(call$knots)){
    knots.no <- k - spline.degree + 1 
    #generate a B-Spline-Matrix with equidistant knots (adapted from code by Thomas Kneib):
    n<-length(x)
    xl<-min(x)
    xr<-max(x)
    xmin<-xl-(xr-xl)/100
    xmax<-xr+(xr-xl)/100
    dx<-(xmax-xmin)/(knots.no-1)
    knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
    call$knots <- knots
  }
  
  
  B<-spline.des(knots,x,spline.degree+1)$design
  
  
  #generate Penalization-Matrix
  D<-diag(k)
  if((d<-min(diff.ord,spline.degree))<diff.ord) warning(paste("order of differences > degree of splines:\n new order of differences=",d,"\n"))
  if(d > 0) for(i in 1:d) D<-diff(D)
  call$diff.ord <- d
  
  #reparametrization: X = unpenalized part, Z =penalized part
  X <-rep(1,k)
  if(diff.ord>1) {for(i in 2:diff.ord) X <-cbind(X, knots[1:k]^(i-1)) }
  
  X <- (B%*%X)[,-1, drop=FALSE]
  Z  <- B%*%t(D)%*%solve(tcrossprod(D))
  
  res <- list(X=unname(X), Z=unname(Z))
  attr(res, "call") <- as.call(call)	
  return(res)
}

