#utils::globalVariables(c("dummy.data.frame",".kolmogorov.dependence.nptest"))
i<-permSpace<-testType<-statTest<-return.permIDs<-P<-idClust<-test <-j <- otherParams<- perms <- NULL

#####trace of a matrix
.tr <- function(sigma) sum(diag(sigma))

.getW <- function(data){
  if(!is.null(data$W)) return(data)
  data$W=array(,dim(data$Y))
  dimnames(data$W)=dimnames(data$Y)
  if(!is.null(data$covs))
  for(j in 1:nrow(data$covs)) {
    data$W[j,]=1/sqrt(diag(data$Su) + diag(matrix(data$covs[j,,])))
  }	else
    for(j in 1:nrow(data$se)) {
      data$W[j,]=1/sqrt(diag(data$Su) + diag(matrix(data$se[j,]^2)))
    }  
  data
}

#################### widelly taken from gt{globaltest}
.getXY <- function(Y,X,Z,data,rotationTest,dummyfy=NULL,statTest,Strata=NULL){
  call <- match.call()
  if(is.null(dummyfy)) {
    dummyfy=list(X=TRUE,Y=TRUE) 
  } else{
    DUM=list(X=TRUE,Y=TRUE)
    DUM[names(dummyfy)]=dummyfy
    dummyfy=DUM
    rm(DUM)
  }
  if(!is.function(statTest)){
    if(statTest%in%c("Wilcoxon", "Kruskal-Wallis", "rank")) dummyfy$X=TRUE
    if(statTest%in%c("AD", "Kolmogorov-Smirnov")) {dummyfy$Y=TRUE ; Yordinal=TRUE} else  Yordinal=FALSE
    if(statTest%in%c("chisq","chisq.separated")) {
      dummyfy=list(X=TRUE,Y=TRUE) ; 
      Yordinal <- Xordinal <-FALSE
      forceFactor=TRUE} else  forceFactor=FALSE
    oldRefCat=options()$ref.cat
    options(exclude.ref.cat=TRUE)
    if(statTest%in%c("AD", "Kolmogorov-Smirnov","chisq","chisq.separated")) {
      options(ref.cat=NULL)
      options(exclude.ref.cat=FALSE) }
  } else { #statTest is a function
    forceFactor <- Yordinal <- Xordinal <-FALSE
    oldRefCat=options()$ref.cat
  }
  
  # data default
  # if (missing(data) || is.null(data))
  # if(is.data.frame(Y) | (is.matrix(Y)))
  # data <- Y else data <- NULL
  # if (is.matrix(data))  
  # data <- data.frame(data)
  if(missing(data)) data=NULL
  if(!is.null(data) && is.matrix(data))  
    data <- data.frame(data)
  
  
  if (missing(X) || is.null(X))
    #if it is a left+right formula
    if (is(Y, "formula") || length(Y)==3){ 
      X <- Y[c(1,3)]
      Y <- Y[c(1,2)] 
      if( !( (length(attr(terms(X, data=data), "term.labels"))==0) & 
               (length(attr(terms(Y, data=data), "term.labels"))==0)  )) {
        dup <- attr(terms(Y, data=data), "term.labels") %in% attr(terms(X, data=data), "term.labels")
        if (any(dup)) 
          Y <- formula(terms(Y,data=data)[!dup])
      }
    } else X <- ~1 
  
  if (missing(Z)) Z=NULL
  
  # remove terms from X that are also in Z
  if (is(Z, "formula") && is(X, "formula") && 
        identical(environment(Z), environment(X))) {
    if( !( (length(attr(terms(X, data=data), "term.labels"))==0) & (length(attr(terms(Z, data=data), "term.labels"))==0)  )) {
      dup <- attr(terms(X, data=data), "term.labels") %in% attr(terms(Z, data=data), "term.labels")
      if (all(dup)) stop("all covariates in X also in Z")
      if (any(dup)) 
        X <- formula(terms(X,data=data)[!dup])
    }
  }
  #browser()  
  # evaluate Y, which may be one of the colnames of data
  if(is(Y,"formula")) Y <- model.frame(Y, data, drop.unused.levels = TRUE,na.action=na.pass)
  Y <- eval(Y, data, parent.frame(sys.nframe()))
  
  n <- nrow(Y)
  
  # get Z and X
  X <- .getAlternative(X, data= if(is.null(data)) data.frame(Y) else data, n,dummyfy=dummyfy,forceFactor=forceFactor)
  attrXassign=attributes(X)$assign
  attrXfactors=attributes(X)$factors
  #browser()
  if(!is.null(Z)){
    Z <- .getNull(Z, data, n)
    offset <- Z$offset   
    Z <- Z$Z
    Z <- Z[,apply(Z,2,function(x) !all(x==0)),drop=FALSE]
    attrXassign=attrXassign[setdiff(colnames(X),colnames(Z))]
    X <- X[, setdiff(colnames(X),colnames(Z)),drop=FALSE]
  }  #else 
  #    if((ncol(X)>0) && rotationTest) {
  # 	 Z=matrix(rep(1,n))
  # 	 attrXassign=attrXassign[!.getIntercept(X)]
  # 	 X <- X[, !.getIntercept(X),drop=FALSE]
  # 	 }
  if(!is.null(Strata)){
    Strata <- .getStrata(Strata, data, n)
  } else Strata=NULL
  
  # # Adjust input due to levels argument
  # if ((!is.null(levels)) && is.factor(Y)) {
  # if (!all(levels %in% levels(Y)))
  # stop("argument \"levels\" does not match levels(Y)")
  # if (length(levels) > 1) {
  # select <- Y %in% levels
  # Y <- factor(Y[select], levels=levels)
  # X <- X[select,,drop=FALSE]
  # Z$Z <- Z$Z[select,, drop=FALSE]
  # if (!is.null(Z$offset)) 
  # Z$offset <- Z$offset[select]
  # if (length(levels) == 2)
  # model <- "logistic"
  # } else {
  # Y <- factor(Y == levels)
  # levels(Y) <- c("other", levels)
  # model <- "logistic"
  # }
  # }
  
  
  # # conservatively impute missing values in X
  # all.na <- apply(is.na(X), 2, all)
  # some.na <- apply(is.na(X), 2, any) & !all.na  
  # if (is.null(Z) || ncol(Z) == 0) {
  # X[is.na(X)] <- 0
  # } else {
  # if(any(some.na))
  # X[,some.na] <- apply(X[,some.na, drop=FALSE], 2, function(cov) {
  # fit <- lm(cov ~ 0 + Z, x = TRUE)
  # coefs <- coef(fit)
  # coefs[is.na(coefs)] <- 0
  # cov[is.na(cov)] <- drop(Z %*% coefs)[is.na(cov)]
  # cov
  # })
  # X[,all.na] <- 0 
  # }
  
  
  # if (is(Y, "formula")) {
  # name.Y <-  as.character(eval(Y)[[2]]) #serve sta roba?
  # Y <- eval(attr(terms(Y, data=data), "variables"), data, environment(Y))[[attr(terms(Y, data=data), "response")]]
  # } else {
  # name.Y <- deparse(call$Y) #serve sta roba?
  # }
  
  # keep NAs
  old.na.action <- options()$na.action  
  options(na.action="na.pass")
  #browser() 
  if(dummyfy$Y) {
    if(Yordinal) Y=as.data.frame(lapply(Y,factor,ordered=TRUE))
    Y =  .makeContrasts(~.,data=data.frame(Y),excludeRefCat=FALSE,excludeIntercept=TRUE,forceFactor=forceFactor)
  }
  #browser()
  # restore default
  options(na.action = old.na.action)
  attributes(X)$assign=attrXassign
  attributes(X)$factors=attrXfactors
  options(ref.cat =oldRefCat)
  if((is.null(Z) && !any(.getIntercept(X))) & rotationTest )   Z=matrix(1,nrow(X))
  return(list(Y=Y,X=X,Z=Z,Strata=Strata,intercept=FALSE))
}

##########################
# set the contrast for factors
#######################
###TODO : sistemare la funzione .makeContrasts!!!!!! mi pare che non funzioni bene excludeRefCat , 
#non funziona neppure con interazioni & excludeRefCat (non esclude l'interazione di riferimento!)
.makeContrasts <- function(formu, data=data,excludeRefCat=options()$exclude.ref.cat,
                           excludeIntercept=FALSE,forceFactor=FALSE){ 
  #excludeRefCat is used only for NOT ordered factors
  # make appropriate contrasts
  opt.ref.cat=options()$ref.cat
  options(ref.cat="first")
  mframe <- model.frame(formu, data=data,na.action = na.pass)
  if(forceFactor){
    toForce <- names(mframe)[!sapply(mframe, is.factor)]
    for(i in toForce)    mframe[,i]=factor(mframe[,i])
    attributes(attributes(mframe)$terms)$dataClasses[toForce]="factor"
  }
  if(length(mframe)>0){
    factors <- names(mframe)[sapply(mframe, is.factor)]
    contrs <- lapply(factors, function(fac) {
      levs <- levels(mframe[[fac]])
      k <- length(levs)
      if (is.ordered(mframe[[fac]])) {
        contr <- array(0, c(k, k-1), list(levs, paste("[",levs[-k], "<", levs[-1],"]", sep="")))
        contr[ lower.tri(contr)] <- 1
      } else {
        contr <- diag(k)
        rownames(contr) <- levs
        colnames(contr) <- paste(".",levs,".",sep="")
        
        if(excludeRefCat) contr <- .leaveRefCat(contr)
        
      }
      contr
    })
    names(contrs) <- factors
  }
  # make the design matrix
  formu <- terms(formu, data=data,na.action = na.pass)
  # if (length(attr(formu, "term.labels")) == 0)
  # stop("empty formu")
  if(excludeIntercept) attr(formu, "intercept") <- 0 else attr(formu, "intercept") <- 1
  # inutile :
  #	attributes(attributes(mframe)$terms)$dataClasses[attributes(attributes(mframe)$terms)$dataClasses=="ordered"]="factor"
  # ords=rep(FALSE,dim(data)[2])
  # for(i in 1:length(ords)) ords[i]=is.ordered(data[,i])
  # for(i in which(ords)) data[,i]=factor(data[,i],ordered=FALSE)
  
  formu <- model.matrix(formu, contrasts.arg=contrs, data=data,na.action = na.pass)
  if(exists("factors")) attributes(formu)$factors=factors
  #	if(!all(colnames(formu) == "(Intercept)" ) ) { #if only the intercept is present
  #	formu <- formu[,colnames(formu) != "(Intercept)",drop=FALSE]    # ugly, but I've found no other way
  #   }
  options(ref.cat=opt.ref.cat)
  formu
}

########################
.leaveRefCat <- function(D){
  if(is.null(options()$ref.cat)) return(D) else
    if(options()$ref.cat=="first") out=1 else
      if(options()$ref.cat=="last") out=ncol(D) else
        out=options()$ref.cat
  D=D[,-out,drop=FALSE]
  D
}
############################
# Get the X design matrix
############################
.getAlternative <- function(X, data, n,dummyfy=list(Y=TRUE,X=TRUE),forceFactor=FALSE) {
  # coerce X into a matrix
  if (is.data.frame(X) || is.vector(X)) {
    if (all(sapply(X, function(x) is.numeric(x)| is.logical(x)))) {
      X <- as.matrix(X)
    } else {
      stop("argument \"X\" could not be coerced into a matrix")
    }
  }
  
  if (is(X, "formula")) {
    # keep NAs
    old.na.action <- options()$na.action  
    options(na.action="na.pass")
    if(is.list(dummyfy) && dummyfy$X) {  X=.makeContrasts(X,data=data,forceFactor=forceFactor) }
    # restore default
    options(na.action = old.na.action)
  }
  #check dimensions and names
  if (nrow(X) != n) {
    stop("the length of \"Y\" (",n, ") does not match the row count of \"X\" (", nrow(X), ")")
  }
  # if (is.null(colnames(X)))
  # stop("colnames missing in X design matrix")
  if(is.null(colnames(X))) colnames(X)=paste("X",1:ncol(X),sep="")
  X
}

############################
# Get the Z design matrix
############################
.getNull <- function(Z, data, n) {
  
  # coerce Z into a matrix and find the offset term
  offset <- NULL
  if (is.data.frame(Z) || is.vector(Z)) {
    if (all(sapply(Z, is.numeric))) {
      Z <- as.matrix(Z)
    } else {
      stop("argument \"Z\" could not be coerced into a matrix")
    }
  }
  if (is(Z, "formula")) {
    if (is.null(data)) {
      tnull <- terms(Z)
      # prevent problems for input ~1 or ~0:
      if (((attr(tnull, "Y") == 0)|| is.null(attr(tnull, "Y")))  && (length(attr(tnull, "term.labels")) == 0)
          && (length(attr(tnull, "offset")) == 0)) {
        if (attr(tnull, "intercept") == 1)
          tnull <- terms(numeric(n) ~ 1)
        else
          tnull <- terms(numeric(n) ~ 0)
      }
      offset <- model.offset(model.frame(tnull))
    } else {
      offset <- model.offset(model.frame(Z, data=data))
      tnull <- terms(Z, data=data)
    }
    data <- model.frame(tnull, data, drop.unused.levels = TRUE)
    Z <- model.matrix(tnull, data)
    
    # # suppress intercept if necessary (can this be done more elegantly?)
    # if (model == "cox") Z <- Z[,names(Z) != "(Intercept)"]
  }
  
  # check dimensions
  if (nrow(Z) != n) {
    stop("the length of \"Y\" (",n, ") does not match the row count of \"Z\" (", nrow(Z), ")")
  }
  list(Z = Z, offset = offset)
}

############################
# Get the Strata vector
############################
.getStrata <- function(Strata, data, n) {
  # coerce Strata into a matrix and find the offset term
  #stop("argument \"Strata\" could not be coerced into a matrix")
  
  if (is(Strata, "formula")) {
    if (is.null(data)) {
      tnull <- terms(Strata)
    } else {
      tnull <- terms(Strata, data=data)
    }
    if(attr(tnull, "intercept") == 1) attr(tnull, "intercept") = NULL
    Strata <- model.frame(tnull, data, drop.unused.levels = TRUE)
    if(ncol(Strata)>1) Strata=as.matrix(apply(Strata,1,paste,collapse="-") )
  } else if(is.vector(Strata)) Strata=as.matrix(Strata)
  
  # check dimensions
  if (nrow(Strata) != n) {
    stop("the length of \"Y\" (",n, ") does not match the row count of \"Strata\" (", nrow(Strata), ")")
  }
  Strata
}

##############################################
.getIntercept <- function(X) apply(X,2,function(x)length(unique(x))==1)

##############################################
orthoZ <- function(Y, X=NULL, Z=NULL, data=NULL,returnGamma=FALSE){ 
  data=.orthoZ (list(Y=Y, X=X, Z=Z),returnGamma=returnGamma)
}

#####
.orthoZ <- function(data,returnGamma=FALSE){  
  if(is.null(data$Z) || (length(data$Z)==0)) return(data)
  attrsYassign<-attributes(data$Y)$assign
  attrsXassign<-attributes(data$X)$assign
  ZZ= try(solve(t(data$Z) %*% data$Z),silent=TRUE)
  if(is(ZZ,"try-error")) {warning("Data can not be orthoganalized"); return(data)}
  IP0 <- diag(nrow(data$Z)) - data$Z %*% ZZ %*% t(data$Z)
  IP0 <- (IP0 + t(IP0))/2
  ei=eigen(IP0)
  if(any(is.complex(ei$values))) {warning("Data can not be orthoganalized"); return(data)}
  ei$vectors <- ei$vectors[,(ei$values > 1e-1)] #gli autovalori sono tutti 0 o 1
  data$Y <- t(ei$vectors)%*%data$Y
  data$X <- t(ei$vectors)%*%data$X
  colnames(data$Y)=.getYNames(data$Y)
  colnames(data$X)=.getXNames(data$X)
  attributes(data$Y)$assign<-attrsYassign
  attributes(data$X)$assign<-attrsXassign
  data$Z=NULL
  if(returnGamma) data$Gamma=ei$vectors
  data
}



################################# conservatively impute missing values in alternative

.fillX <- function(alternative,null){
  #WARNING: have to be a 0-centered matrix:
  alternative= scale(alternative,scale =FALSE) 
  null =scale (null,scale=FALSE)
  
  all.na <- apply(is.na(alternative), 2, all)
  some.na <- apply(is.na(alternative), 2, any) & !all.na
  if (missing(null) || is.null(null) || (ncol(null) == 0)) {
    alternative[is.na(alternative)] <- 0
  } else {
    alternative[,some.na] <- apply(alternative[,some.na, drop=FALSE], 2, function(cov) {
      fit <- lm(cov ~ 0 + null, x = TRUE)
      coefs <- coef(fit)
      coefs[is.na(coefs)] <- 0
      cov[is.na(cov)] <- drop(null %*% coefs)[is.na(cov)]
      cov
    })
    alternative[,all.na] <- 0 
  }
  alternative <- alternative+matrix(attr(alternative,"scaled:center"),byrow=TRUE,nrow=dim(alternative)[1],ncol=dim(alternative)[2])
}


#####################################################
.getSubsetWeights<-function(weights, subsets,colNames.permSpace){
  if(missing(subsets)) {subsets=NULL; many.subsets=FALSE}
  if(missing(weights)) {weights=NULL; many.weights=FALSE}
  
  #subsets and weights
  if (!is.null(subsets) && !is.list(subsets))
    subsets <- list(subsets)
  many.subsets <- !is.null(subsets)
  one.weight <- (!is.null(weights)) && (!is.list(weights)) && (length(weights)==length(colNames.permSpace)) && many.subsets
  many.weights <- (!is.null(weights)) && (!one.weight)
  if (many.weights && !is.list(weights))
    weights <- list(weights)
  
  # check weights and subsets lengths
  if (many.subsets && many.weights) {
    if (length(subsets) != length(weights))
      stop("lengths of \"subsets\" and \"weights\" do not match")
    if (!((!is.null(names(subsets))) && (!is.null(names(weights))) && (!all(names(subsets)==names(weights)))))
      #if (is.null(alias))
      # alias <- names(weights)
      #else
      warning("names of subsets and weights do not match")
  }
  
  
  # make sure subsets is a list of names, compatible with colnames(tail)
  if (many.subsets) {
    osl <- sapply(subsets, length)
    subsets <- lapply(subsets, function(sst) {
      if (!is.character(sst)) 
        colNames.permSpace[sst]
      else
        intersect(sst, colNames.permSpace)
    })
  }
  
  # make sure that weights is a named list, compatible with colnames(tail)
  if (many.weights) {
    names.weights <- names(weights)
    weights <- lapply(1:length(weights), function (i) {
      wt <- weights[[i]]
      if (!is.null(names(wt)))
        wt <- wt[names(wt) %in% colnames(tail)]
      else 
        if (many.subsets && length(wt) == length(subsets[[i]]))
          names(wt) <- subsets[[i]]
      else if (length(wt) == ncol(tail))
        names(wt) <- colnames(tail)
      wt
    })
    names(weights) <- names.weights
    if (any(sapply(lapply(weights, names), is.null)))
      stop("weights input is not compatible with variables input.")
  }
  
  # make subsets and weights compatible
  if (many.subsets && many.weights) {
    weights <- lapply(1:length(weights), function(i) {
      if (all(subsets[[i]] %in% names(weights[[i]])))
        weights[[i]][subsets[[i]]]
      else 
        stop("names of weights input incompatible with subsets input.")
    })
  }
  
  # trim zero weights
  if (many.weights) {
    if (any(unlist(weights)==0)) {
      weights <- lapply(weights, function(wt) wt[wt != 0])
      many.subsets <- FALSE   # to redo subsets
    }
  }
  
  # make subsets in case of short named weights
  if (many.weights && !many.subsets && any(sapply(weights, length) != length(colNames.permSpace))) {
    subsets <- lapply(weights, names)
    many.subsets <- TRUE
  }
  
  # check missing values in subsets
  if (many.subsets && any(sapply(subsets, function(x) any(is.na(x))))) {
    stop("missing values in \"subsets\"")
  }
  
  return(list(one.weight=one.weight,many.subsets=many.subsets, many.weights=many.weights, subsets=subsets, weights=weights))
}

########################################################### tail util

#make "hight" (depending on the tail) values of the statistics to be significative
.fitTail <- function(permT,tail){
  if (missing(tail)||is.null(tail)) {
    tail = rep(0, ncol(permT))
  }     else if (length(tail) != ncol(permT)) {
    attrs=attributes(tail)$center
    tail <- rep(tail,len = ncol(permT))
    if(!is.null(attrs)) attributes(tail)$center=rep(attrs,len = ncol(permT))
  }
  tail
}

.setTail <- function(permT, tail){
  tail=.fitTail(permT,tail)
  if (!is.null(attributes(tail)$center)) {
    #tail==0 and need to be centered:
    center=attributes(tail)$center
    inters = intersect(colnames(permT)[tail==0],names(center)) 
    if(length(inters)>0) {
      center = attributes(tail)$center[inters] 
      permT[,names(center)] <- scale(permT[,names(center)],center=center)
    }
  }
  permT[, tail < 0] <- -permT[, tail < 0]
  permT[, tail == 0] <- abs(permT[, tail == 0])
  permT[is.na(permT)] <- 0
  permT
}

.setTailOut <- function(permT=permT, tail=tail){
  dir=as.character(sign(c(1,-1)%*%.setTail(matrix(c(1,-1),2,ncol(permT)),tail)))
  dir[dir=="1"]=">"
  dir[dir=="-1"]="<"
  dir[dir=="0"]="><"
  dir
}


###########################################################

################### get results
.getOut <- function(type="flip",res=NULL, data=NULL, call=NULL, flipReturn=list(permT=TRUE,call.env=TRUE),
                    separatedX=TRUE,extraInfoPre=NULL,extraInfoPost=NULL,call.env=NULL,...){ 
  colnames(res$permT)=.getTNames(data$Y,data$X,permT=res$permT) 
  
  # test statistic and std dev
  stat=res$permT[1,]
  #### da rivedere
  #stDev=apply(res$permT,2,sd,na.rm=TRUE)
  #pseudoZ=.t2stdt(res$permT)
  p=t2p(res$permT,obs.only=TRUE,tail=res$tail)
  
  
  if(type=="flip"){
    # tails of the test
    dir=.setTailOut (permT=res$permT, tail=res$tail)		
    #build the results table
    TAB=data.frame(Stat=as.vector(stat),#sd.permT=as.vector(stDev), pseudoZ=as.vector(pseudoZ),
                   tail=as.vector(dir), p=as.vector(p),stringsAsFactors =FALSE)
  } else if(type=="npc"){
    #build the results table
    TAB=data.frame(Stat=as.vector(stat),#sd.permT=as.vector(stDev), pseudoZ=as.vector(pseudoZ), 
                   p=as.vector(p),stringsAsFactors =FALSE)
    colnames(TAB)[colnames(TAB)=="nvar"]="#Vars"
  }
  
  if((!is.null(res$extraInfoPre))) {TAB=cbind(data.frame(res$extraInfoPre,row.names = NULL,stringsAsFactors =FALSE),TAB)}
  if((!is.null(res$extraInfoPost))) {TAB=cbind(TAB,data.frame(res$extraInfoPost,row.names = NULL,stringsAsFactors =FALSE))}
  colnames(TAB)[colnames(TAB)=="p"]="p-value"
  rownames(TAB)=colnames(res$permT)
  
  out <- new("flip.object")  
  out @res = TAB
  out @permSpace=if( (!is.null(flipReturn$permSpace)&&flipReturn$permSpace)||
                       (!is.null(flipReturn$permID)&&flipReturn$permID)  ) 
    res$perms else res$perms[-which(names(res$perms)=="permID")]
  out @call = if(!is.null(call)) call
  out @permP=if(!is.null(flipReturn$permP))if(flipReturn$permP) t2p(res$permT, obs.only=FALSE,tail=res$tail)
  out @permT=if(!is.null(flipReturn$permT))if(flipReturn$permT) res$permT
  out @data = if(!is.null(flipReturn$data))if(flipReturn$data) data
  out @call.env = if(is.null(flipReturn$call.env) || flipReturn$call.env) call.env
  if(!is.null(res$tail)) 
    out @tail = as.matrix(res$tail)
  out	
}


.getTNames <- function(Y,X=NULL,permT=NULL,checkUnique=FALSE){
  if(!is.null(colnames(permT))) {
    colnames(permT)[colnames(permT)==""] = paste("V",1:ncol(permT),sep="")[colnames(permT)==""]
    if(checkUnique){
      temp=table(colnames(permT))
      for(v in names(temp)[temp>1]){
        substitute=which(colnames(permT)==v)
        colnames(permT)[substitute]=paste(v,sep=".",1:length(substitute))
      }
    }
    
    return(colnames(permT))
  } else	if(is.null(Y)){ 
    return(paste("V",1:ncol(permT),sep=""))
  } else {	
    if(!is.null(X)) if(ncol(X)==0) X=NULL
    
    colnames(Y) = .getYNames(Y)
    colnames(X) = .getXNames(X)
    TNames=paste(
      rep(colnames(Y),rep(max(ncol(X),1),ncol(Y))),if((!is.null(X))&&(ncol(X)>1)) paste("_|_",colnames(X),sep="") else "" ,sep="")
    
    return(TNames)
  }
}

.getNames <- function(Y,prefix=".") {
  if(!is.null(Y)) {if(!is.null(colnames(Y))) colnames(Y) else	paste(prefix,sep="",if(ncol(Y)>1)1:ncol(Y)) } else NULL
}

.getYNames <- function(Y) {
  if(!is.null(Y)) {if(!is.null(colnames(Y))) colnames(Y) else	paste("Y",sep="",if(ncol(Y)>1)1:ncol(Y)) } else NULL
}

.getXNames <- function(X) {
  if(!is.null(X)) {if(!is.null(colnames(X)))  colnames(X) else paste("X",sep="",if(ncol(X)>1)1:ncol(X)) } else NULL
}

.getTRowNames <- function(permT){
  c("Tobs", paste("T*",1:(nrow(permT)-1),sep=""))
}




###################################
###################################
#####tests utilities
###################################

#################################
.prod.perms <- function(data,perms,testType="permutation"){
  if(testType%in%c("permutation","symmetry","rotation"))  {
    if(is.null(perms$permID)){
       digitsK=trunc(log10(perms$B))+1
      envOrig<-environment(perms$rotFunct)
      environment(perms$rotFunct) <- sys.frame(sys.parent())
      obs=as.vector(t(data$X)%*%data$Y)
      permT=matrix(,perms$B,length(obs))
      permT[1,]=obs
      rm(obs)     
      for(i in 1:(perms$B-1)) {
#         if (i%%10==0) {
#           cat(rep("\b", 2*digitsK+10), i, " / ", perms$B, sep="")
#           flush.console()
#         }
        permT[i+1,]=as.vector(t(data$X)%*%perms$rotFunct())
      }
      cat(rep("\b", 2*digitsK+1));  flush.console()
      
      environment(perms$rotFunct) <- envOrig
      colnames(permT)=.getTNames(data$Y,data$X)
      rownames(permT)=.getTRowNames(permT)
    } else { #permutation test uding IDs
      m=ncol(data$Y)
      q=ncol(data$X)
      digitsK=trunc(log10(m))+1
      if(testType=='symmetry')  { #on.exit(browser())
        permT=matrix(,nrow(perms$permID)+1,m*q)
        for( i in 1:ncol(data$X)){
          XY=diag(data$X[,i])%*%data$Y
          permT[,(0:(m-1))*q+i]=rbind(rep(1,perms$n),perms$permID) %*% XY
      }
        permT = rbind(permT,-permT[nrow(permT):1,,drop=FALSE])
      } else  {
        permT=matrix(,perms$B,m*q)
        for( i in 1:ncol(data$Y)){
          permT[,(i-1)*q+(1:q)]=(matrix(data$Y[rbind(1:perms$n,perms$permID),i],ncol=perms$n,nrow=(perms$B)) %*% data$X)
        }
      }
    }
  } else {warning("test type not implemented (yet?)"); return(NULL)}
  cat(rep("\b", 2*digitsK+3));  flush.console()
  permT
}


.prod.perms.P <-function(data,perms,testType="permutation",P=P){
  if(testType%in%c("permutation","rotation"))  {
    if(is.null(perms$permID)){
      digitsK=trunc(log10(perms$B))+1
      envOrig<-environment(perms$rotFunct)
      environment(perms$rotFunct) <- sys.frame(sys.parent())
      obs=as.vector(apply((t(data$Y)%*%P)^2,1,sum))
      permT=matrix(,perms$B,length(obs))
      permT[1,]=obs
      rm(obs)
      for(i in 1:(perms$B-1)){ 
        if (i%%10==0) {
          cat(rep("\b", 2*digitsK+10), i, " / ", perms$B, sep="")
          flush.console()
        }
        # R is random matrix of independent standard-normal entries 
        # Z shall be a random matrix with the same mean and covariance structure as Y 
        permT[i+1,]=apply((t(perms$rotFunct())%*%P)^2,1,sum)
      }
      environment(perms$rotFunct) <- envOrig
      colnames(permT)=.getTNames(data$Y)
      rownames(permT)=.getTRowNames(permT)
    } else { #permutation test using IDs
      m=ncol(data$Y)
      digitsK=trunc(log10(m))+1
      permT <- matrix(,perms$B,m)
      for(i in 1:m){
        #       if (i%%10==0) {
        # 		    cat(rep("\b", 2*digitsK+5), i, " / ", m, sep="")
        # 		    flush.console()
        # 		  }
        permT[,i]=rowSums((matrix(data$Y[rbind(1:perms$n,perms$permID),i],ncol=perms$n,nrow=(perms$B))%*%P)^2)
      }
      colnames(permT)=.getTNames(data$Y)
      rownames(permT)=.getTRowNames(permT)
    }
  } else {warning("test type not implemented (yet?)"); return(NULL)}
  cat(rep("\b", 2*digitsK+3));  flush.console()
  permT
}



.prod2sum <- function(permT,data){
  N=nrow(data$Y)
  colnames(permT)=.getTNames(data$Y,data$X)
  rownames(permT)=.getTRowNames(permT)
  permT
}	

.prod2t <- function(permT,data){
  N=nrow(data$Y)
  if(data$intercept==TRUE){
    sYs=colSums(data$Y)
    sXs=colSums(data$X)
    data$X=scale(data$X,center=TRUE,scale=FALSE)
    data$Y=scale(data$Y,center=TRUE,scale=FALSE)
    permT=scale(permT,center=as.vector(sXs%*%t(sYs)/N),scale=FALSE)
  }
  permT=scale( permT,center=FALSE,scale= rep(sqrt(colSums(data$Y^2)),rep(ncol(data$X),ncol(data$Y)))*rep(sqrt(colSums(data$X^2)),ncol(data$Y)))
  permT=permT/sqrt(1-permT^2)*sqrt(N-ncol(data$X))
  
  colnames(permT)=.getTNames(data$Y,data$X)
  rownames(permT)=.getTRowNames(permT)
  permT
}				

.prod2F <- function(permT,data){
  #sumY2=colSums(data$Y)^2
  for(i in 1:ncol(permT)){ 
    # (permT[,i]-sumY2[i])/(sum(data$Y[,i]^2) - sumY2[i]-permT[,i])} #non funziona. allora centrare le Y e poi
    permT[,i]=permT[,i]/(sum(data$Y[,i]^2)-permT[,i])
  }
  permT=permT*(nrow(data$Y)-ncol(data$X))/ncol(data$X)
  if(is.null(ncol(permT))) permT=as.matrix(permT)
  permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi    		
  colnames(permT)=.getTNames(data$Y)
  rownames(permT)=.getTRowNames(permT)	
  permT
}


.get.eigenv.proj.mat <- function(data){
  P=data$X%*%solve(t(data$X)%*%data$X)%*%t(data$X)
  P = eigen((P + t(P))/2)
  P <- P$vectors[,(P$values > 1e-1),drop=FALSE] #gli autovalori sono tutti 0 o 1
  P
}
# 
# .getDummiesGroups <- function(X,excludeRefCat=excludeRefCat) {
# 	grps=apply(X,1,paste,collapse="")
# 	if(length(unique(grps))<=1) {
# 		warning("At least two groups are required.") 
# 		return(X)
# 	} else if(length(unique(grps))==2 ){
# 		return(X) 
# 	} else {
# #          dummy.data.frame <- function(x){
# #               colnam=unique(x)
# #               out=outer(x,colnam,"==")*1
# #               colnames(out)=colnam
# #               out
# #          }
# 		res=dummy.data.frame(apply(X,1,paste,collapse=""))
# 		if(excludeRefCat>0) res <- res[,-excludeRefCat,drop=FALSE]
# 	}
# }
