###{{{ parsedesign

##' @export
parsedesign <- function(coef,x,...) {
    if (!is.vector(coef)) coef <- stats::coef(coef)
    if (is.numeric(coef) && !is.null(names(coef))) coef <- names(coef)
    dots <- substitute(list(...))[-1]
    expr <- suppressWarnings(inherits(try(x,silent=TRUE),"try-error"))
    if (expr) {
        ee <- c(deparse(substitute(x)), sapply(dots, deparse))
    } else {
        ee <- c(deparse(x), sapply(dots, function(x) deparse(x)))
    }
    res <- c()
    for (e in ee) {
        e0 <- gsub(" ","",e)
        ff <- strsplit(e0,'\"')[[1]]
        Val <- rbind(rep(0,length(coef)))
        for (i in seq(length(ff)/2)) {
            val0 <- gsub("[*()]","",ff[2*(i-1)+1])            
            suppressWarnings(val <- as.numeric(val0))
            if (is.na(val)) {
                val <- switch(val0,"-"=-1,1)
            }
            par0 <- ff[2*i]
            par0int <- suppressWarnings(as.integer(par0))
            if (is.na(par0int)) par0int <- match(par0,coef)
            if (par0int<=length(Val)) Val[par0int] <- val
        }
        if (any(Val!=0)) res <- rbind(res,Val)
    }
    res
}

###}}}

###{{{ contr

##' @export
contr <- function(p,n,...) {
  if (length(p)==1) {
    B <- matrix(0,ncol=p*n,nrow=p*(n-1))
    pos <- 0
    for (i in seq_len(p)) {
      for (j in seq_len(n-1)) {
        pos <- pos+1
        B[pos,i] <- 1;  B[pos,j*p+i] <- -1
      }
    }
    return(B)
  }
  if (missing(n)) n <- max(p)
  B <- matrix(0,ncol=n,nrow=length(p)-1)
  B[,p[1]] <- 1
  B[cbind(seq(nrow(B)),p[-1])] <- -1
  B
}

###}}} contr

###{{{ substArg

substArg <- function(x,env,...) {
  if (!missing(env)) {
    a <- with(env,substitute(x))
#    a <- substitute(x,environment(env))
  } else {
    a <- substitute(x)
  }
  myclass <- tryCatch(class(eval(a)),error=function(e) NULL)
  if (is.null(myclass) || myclass=="name") {
#  if (is.null(myclass)) {
    res <- unlist(sapply(as.character(a),
                         function(z) {
                           trimmed <- gsub(" ","",z,fixed=TRUE)
                           val <- strsplit(trimmed,"+",fixed=TRUE)
                           if (val[1]=="") val <- NULL
                           val
                         })); attributes(res)$names <- NULL
    return(res)
  }
  return(eval(a))
}

## g <- function(zz,...) {
##   env=new.env(); assign("x",substitute(zz),env)
##   substArg(zz,env=env)
## }
## h <- function(x,...) {
##   env=new.env(); assign("x",substitute(x),env)
##   substArg(x,env=TRUE)
## }

###}}}

###{{{ procrandomslope

procrandomslope <- function(object,data=object$data,...) {
  Xfix <- FALSE
  xfix <- myfix <- list()
  xx <- object
  for (i in seq_len(object$ngroup)) {
    x0 <- object$lvm[[i]]
    data0 <- data[[i]]
    xfix0 <- colnames(data0)[(colnames(data0)%in%parlabels(x0,exo=TRUE))]
    xfix <- c(xfix, list(xfix0))
    if (length(xfix0)>0) { ## Yes, random slopes
      Xfix<-TRUE
    }
    xx$lvm[[i]] <- x0
  }
  if (Xfix) {
    for (k in seq_len(object$ngroup)) {
      x1 <- x0 <- object$lvm[[k]]
      data0 <- data[[k]]
      nrow <- length(vars(x0))
      xpos <- lapply(xfix[[k]],function(y) which(regfix(x0)$labels==y))
      colpos <- lapply(xpos, function(y) ceiling(y/nrow))
      rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
      myfix0 <- list(var=xfix[[k]], col=colpos, row=rowpos)
      myfix <- c(myfix, list(myfix0))
      for (i in seq_along(myfix0$var))
        for (j in seq_along(myfix0$col[[i]]))
          regfix(x0,
                 from=vars(x0)[myfix0$row[[i]][j]],to=vars(x0)[myfix0$col[[i]][j]]) <-
                   colMeans(data0[,myfix0$var[[i]],drop=FALSE],na.rm=TRUE)
      index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
      object$lvm[[k]] <- x0
      yvars <- endogenous(x0)
      #parkeep <- c(parkeep, parord[[k]][coef(x1,mean=TRUE)%in%coef(x0,mean=TRUE)])
    }
#    parkeep <- sort(unique(parkeep))
    object <- multigroup(object$lvm,data,fix=FALSE,exo.fix=FALSE)
  }
  return(list(model=object,fix=myfix))
}

###}}} procrandomslope

###{{{ kronprod

## ' Calculate matrix product with kronecker product
## '
## ' \deqn{(A\crossprod B) Y}
## ' @title Calculate matrix product with kronecker product
## ' @param A
## ' @param B
## ' @param Y
## ' @author Klaus K. Holst
kronprod <- function(A,B,Y) {
  rbind(apply(Y,2,function(x) B%*%matrix(x,nrow=ncol(B))%*%t(A)))
}

###}}} kronprod

###{{{ izero

izero <- function(i,n) { ## n-1 zeros and 1 at ith entry
  x <- rep(0,n); x[i] <- 1
  x
}

###}}}

###{{{ Debug

`Debug` <-
  function(msg, cond=lava.options()$debug) {
    if (cond)
      print(paste(msg, collapse=" "))
  }

###}}}

###{{{ categorical2dummy

categorical2dummy <- function(x,data,silent=TRUE,...) {
  x0 <- x
  X <- intersect(index(x)$exogenous,colnames(data))
  catX <- c()
  for (i in X) {
    if (!is.numeric(data[,i])) catX <- c(catX,i)
  }
  if (length(catX)==0) return(list(x=x,data=data))
  f <- as.formula(paste("~ 1+", paste(catX,collapse="+")))
  opt <- options(na.action="na.pass")
  M <- model.matrix(f,data)

  options(opt)
  Mnames <- colnames(M)
  Mpos <- attributes(M)$assign
  A <- index(x)$A
  F <- regfix(x)
  count <- 0
  for (i in catX) {
    count <- count+1
    mnames <- Mnames[Mpos==count]
    kill(x0) <- i
    Y <- colnames(A)[A[i,]==1]
    if (length(mnames)==1) {
      fix <- as.list(F$labels[i,])
      fixval <- F$values[i,]
      fix[which(!is.na(fixval))] <- fixval[na.omit(fixval)]
      regression(x0,to=Y,from=mnames,silent=silent) <- fix[Y]
    } else {
      x0 <- regression(x0,to=Y,from=mnames,silent=silent)
    }
  }
  index(x0) <- reindex(x0,zeroones=TRUE,deriv=TRUE)
  return(list(x=x0,data=cbind(data,M)))
}

###}}}

###{{{ procdata.lvm

`procdata.lvm` <-
  function(x,data,categorical=FALSE,
##           na.method=ifelse(any(is.na(data[,intersect(colnames(data),exogenous(x))])),"pairwise.complete.obs","complete.obs")
           na.method=ifelse(any(is.na(data[,intersect(colnames(data),manifest(x))])),"pairwise.complete.obs","complete.obs")
##           na.method=c("pairwise.complete.obs")
           ) {
    if (is.numeric(data) & !is.list(data)) {
      data <- rbind(data)
    }
     if (is.data.frame(data) | is.matrix(data)) {
      nn <- colnames(data)
      data <- as.data.frame(data); colnames(data) <- nn; rownames(data) <- NULL
      obs <- setdiff(intersect(vars(x), colnames(data)),latent(x))
      Debug(obs)
      mydata <- subset(data, select=obs)
      if (NROW(mydata)==0) stop("No observations")
      for (i in seq_len(ncol(mydata))) {
        if (inherits(mydata[,i],"Surv"))
          mydata[,i] <- mydata[,i][,1]
        if (is.character(mydata[,i]) | is.factor(mydata[,i]))
          mydata[,i] <- as.numeric(as.factor(mydata[,i]))-1
      }

##      mydata <- data[,obs]
##      if (any(is.na(mydata))) {
##        warning("Discovered missing data. Going for a complete-case analysis. For data missing at random see 'missingMLE'.\n", immediate.=TRUE)
##        mydata <- na.omit(mydata)
##      }
      S <- NULL
      n <- nrow(mydata)
      if (n==1) {
        S <- diag(nrow=ncol(mydata)); colnames(S) <- rownames(S) <- obs
      }
      if (na.method=="pairwise.complete.obs") {
        mu <- colMeans(mydata,na.rm=TRUE)
        if (is.null(S)) {
          S <- (n-1)/n*cov(mydata,use=na.method)
          S[is.na(S)] <- 1e-3
        }
      }
      if (na.method=="complete.obs") {
        mydata <- na.omit(mydata)
        n <- nrow(mydata)
        mu <- colMeans(mydata)
        if (is.null(S))
          S <- (n-1)/n*cov(mydata) ## MLE variance matrix of observed variables
      }
    }
    else
      if (is.list(data)) {
        if ("cov"%in%names(data)) data$S <- data$cov
        if ("var"%in%names(data)) data$S <- data$var
        if ("mean"%in%names(data)) data$mu <- data$mean
        n <- data$n
        S <- reorderdata.lvm(x,data$S)
        mu <- reorderdata.lvm(x,data$mu)
        ##      if (is.null(n)) stop("n was not specified");
      }
      else
        stop("Unexpected type of data!");
    if (nrow(S)!=ncol(S)) stop("Wrong type of data!");
    return(list(S=S,mu=mu,n=n))
  }

###}}}

###{{{ reorderdata.lvm

`reorderdata.lvm` <-
  function(x, data) {
    if (is.vector(data)) {
      nn <- names(data)
      ii <- na.omit(match(index(x)$manifest, nn))
      data[ii,drop=FALSE]
    } else {
      nn <- colnames(data)
      ii <- na.omit(match(index(x)$manifest, nn))
      data[ii,ii,drop=FALSE]
    }
  }

###}}}

###{{{ symmetrize

`symmetrize` <-
function(M, upper=TRUE) {
  if (length(M)==1) return(M)
  if (!is.matrix(M) | ncol(M)!=nrow(M)) stop("Only implemented for square matrices.")
  if (upper) {
    for (i in seq_len(ncol(M)-1))
      for (j in seq(i+1,nrow(M)))
        M[i,j] <- M[j,i]
    return(M)
  } else {
    for (i in seq_len(ncol(M)))
      for (j in seq_len(nrow(M)))
        if (M[i,j]==0)
          M[i,j] <- M[j,i]
        else
          M[j,i] <- M[i,j]
    return(M)
  }
}

###}}}

###{{{ Inverse/pseudo

##' @export
Inverse <- function(X,tol=lava.options()$itol,det=TRUE,names=!chol,chol=FALSE) {
    n <- NROW(X)
    if (n==1L) {
        res <- 1/X
        if (det) attributes(res)$det <- X
        if (chol) attributes(res)$chol <- X
        return(res)
    }
    if (chol) {
        L <- chol(X)
        res <- chol2inv(L)
        if (det) attributes(res)$det <- prod(diag(L)^2)
        if (chol) attributes(res)$chol <- X        
    } else {
        svdX <- svd(X)
        id0 <- numeric(n)
        idx <- which(svdX$d>tol)
        id0[idx] <- 1/svdX$d[idx]
        res <- with(svdX, v%*%diag(id0,nrow=length(id0))%*%t(u))
        if (det)
            attributes(res)$det <- prod(svdX$d[svdX$d>tol])
        attributes(res)$pseudo <- (length(idx)<n)
        attributes(res)$minSV <- min(svdX$d)
    }
    if (names && !is.null(colnames(X))) dimnames(res) <- list(colnames(X),colnames(X))
    return(res)
}

###}}}

###{{{ naiveGrad

naiveGrad <- function(f, x, h=1e-9) {
  nabla <- numeric(length(x))
  for (i in seq_along(x)) {
    xh <- x; xh[i] <- x[i]+h
    nabla[i] <- (f(xh)-f(x))/h
  }
  return(nabla)
}

###}}}

###{{{ CondMom

# conditional on Compl(idx)
CondMom <- function(mu,S,idx,X) {
  idxY <- idx

  idxX <- setdiff(seq_len(ncol(S)),idxY)
  SXX <- S[idxX,idxX,drop=FALSE];
  SYY <- S[idxY,idxY,drop=FALSE]
  SYX <- S[idxY,idxX,drop=FALSE]
  iSXX <- solve(SXX)
  condvar <- SYY-SYX%*%iSXX%*%t(SYX)
  if (missing(mu)) return(condvar)

  muY <- mu[,idxY,drop=FALSE]
  muX <- mu[,idxX,drop=FALSE]
  if (is.matrix(mu))
    Z <- t(X-muX)
  else
    Z <- apply(X,1,function(xx) xx-muX)
  SZ  <- t(SYX%*%iSXX%*%Z)
##  condmean <- matrix(
  if (is.matrix(mu))
    condmean <- SZ+muY
  else
    condmean <- t(apply(SZ,1,function(x) muY+x))
##  ,ncol=ncol(SZ),nrow=nrow(SZ))
  return(list(mean=condmean,var=condvar))
}

###}}} CondMom

###{{{ Depth-First/acc (accessible)

DFS <- function(M,v,explored=c()) {
  explored <- union(explored,v)
  incident <- M[v,]
  for (v1 in setdiff(which(incident==1),explored)) {
    explored <- DFS(M,v1,explored)
  }
  return(explored)
}

acc <- function(M,v) {
  if (is.character(v)) v <- which(colnames(M)==v)
  colnames(M)[setdiff(DFS(M,v),v)]
}

###}}} Depth-First/acc (accessible)

## Trace operator
tr <- function(x) sum(diag(x))

npar.lvm <- function(x) {
  return(index(x)$npar+ index(x)$npar.mean+index(x)$npar.ex)

}

as.numeric.list <- function(x,...) {
  res <- list()
  asnum <- as.numeric(x)
  lapply(x,function(y) ifelse(is.na(as.numeric(y)),y,as.numeric(y)))
}

edge2pair <- function(e) {
  sapply(e,function(x) strsplit(x,"~"))
}
numberdup <- function(xx) { ## Convert to numbered list
  dup.xx <- duplicated(xx)
  dups <- xx[dup.xx]
  xx.new <- numeric(length(xx))
  count <- 0
  for (i in seq_along(xx)) {
    if (!dup.xx[i]) {
      count <- count+1
      xx.new[i] <- count
    } else {
      xx.new[i] <- xx.new[match(xx[i],xx)[1]]
    }
  }
  return(xx.new)
}


##' @export
logit <- function(p) log(p/(1-p))

##' @export
expit <- function(z) 1/(1+exp(-z))

##' @export
tigol <- expit

extractvar <- function(f) {
    yy <- getoutcome(f)
    xx <- attributes(terms(f))$term.labels
    myvars <- all.vars(f)
    return(list(y=yy,x=xx,all=myvars))
}

##' @export
getoutcome <- function(formula,sep) {
  aa <- attributes(terms(formula))
  if (aa$response==0) {
    res <- NULL
  } else {
    res <- paste(deparse(formula[[2]]),collapse="")
  }
  if (!missing(sep)) {
      attributes(res)$x <- lapply(strsplit(aa$term.labels,"\\|")[[1]],
                                  function(x) as.formula(paste0("~",x)))
  } else {
      attributes(res)$x <- aa$term.labels
  }
  return(res)
}


##' @export
Specials <- function(f,spec,split2="+",...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x)
  res <- unlist(strsplit(st,"[()]"))[2]
  if (is.null(split2)) return(res)
  unlist(strsplit(res,"+",fixed=TRUE))
}


##' @export
decomp.specials <- function(x,pattern="[()]",pattern2=NULL, pattern.ignore=NULL, sep=",",reverse=FALSE,...) {
  st <- gsub(" |^\\(|)$","",x) # Remove white space and leading/trailing parantheses
  if (!is.null(pattern.ignore)) {
      if (grepl(pattern.ignore,st,...)) return(st)
  }
  if (!is.null(pattern)) {
    st <- rev(unlist(strsplit(st,pattern,...)))[1]
  }
  if (!is.null(pattern2)) {
    st <- (unlist(strsplit(st,pattern2,...)))
    if (reverse) st <- rev(st)
  }
  unlist(strsplit(st,sep,...))
}

Decomp.specials <- function(x,pattern="[()]") {
  st <- gsub(" ","",x)
  st <- gsub("\n","",st)
  mysplit <- rev(unlist(strsplit(st,pattern)))
  type <- mysplit[2]
  vars <- mysplit[1]
  res <- unlist(strsplit(vars,","))
  if (type=="s" | type=="seq") {
    return(paste0(res[1],seq(as.numeric(res[2]))))
  }
  unlist(strsplit(vars,","))

}

printline <- function(n=70) {
    cat(rep("_", n), "\n", sep="");
}
