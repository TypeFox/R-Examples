###{{{ print.fix

##' @export
print.fix <- function(x,exo=FALSE,...) {
  switch(attributes(x)$type,
        reg = cat("Regression parameters:\n"),
        cov = cat("Covariance parameters:\n"),
        mean = cat("Intercept parameters:\n"))
  M <- linconstrain(x,print=TRUE)
  invisible(x)
}

linconstrain <- function(x,print=TRUE,indent="  ",exo=FALSE,...) {
  idx <- seq_len(attributes(x)$nvar)
  idx0 <- setdiff(idx,attributes(x)$exo.idx)
  if (!exo & attributes(x)$type!="reg")
    idx <- idx0
  if (attributes(x)$type=="mean") {
    M <- rbind(unlist(x[idx]))
    rownames(M) <- ""
    M[is.na(M)] <- "*"
  } else {
    M <- x$rel[idx,idx,drop=FALSE]
    M[M==0] <- NA
    M[M==1] <- "*"
    M[which(!is.na(x$labels[idx,idx]))] <- x$labels[idx,idx][which(!is.na(x$labels[idx,idx]))]
    M[which(!is.na(x$values[idx,idx]))] <- x$values[idx,idx][which(!is.na(x$values[idx,idx]))]
    if (attributes(x)$type=="reg")
        M <- t(M[,idx0,drop=FALSE])
  }
  if (print) {
    M0 <- M
    if (NROW(M)>0)
      rownames(M0) <- paste(indent,rownames(M))
    print(M0,quote=FALSE,na.print="",...)
  }
  invisible(M)
}

###}}} print.fix

###{{{ intfix

##' @export
"intfix" <- function(object,...) UseMethod("intfix")
##' @export
"intfix<-" <- function(object,...,value) UseMethod("intfix<-")

##' Fix mean parameters in 'lvm'-object
##'
##' Define linear constraints on intercept parameters in a \code{lvm}-object.
##'
##'
##' The \code{intercept} function is used to specify linear constraints on the
##' intercept parameters of a latent variable model. As an example we look at
##' the multivariate regression model
##'
##' \deqn{ E(Y_1|X) = \alpha_1 + \beta_1 X} \deqn{ E(Y_2|X) = \alpha_2 + \beta_2
##' X}
##'
##' defined by the call
##'
##' \code{m <- lvm(c(y1,y2) ~ x)}
##'
##' To fix \eqn{\alpha_1=\alpha_2} we call
##'
##' \code{intercept(m) <- c(y1,y2) ~ f(mu)}
##'
##' Fixed parameters can be reset by fixing them to \code{NA}.  For instance to
##' free the parameter restriction of \eqn{Y_1} and at the same time fixing
##' \eqn{\alpha_2=2}, we call
##'
##' \code{intercept(m, ~y1+y2) <- list(NA,2)}
##'
##' Calling \code{intercept} with no additional arguments will return the
##' current intercept restrictions of the \code{lvm}-object.
##'
##' @aliases intercept intercept<- intercept.lvm intercept<-.lvm intfix intfix
##' intfix<- intfix.lvm intfix<-.lvm
##' @param object \code{lvm}-object
##' @param vars character vector of variable names
##' @param value Vector (or list) of parameter values or labels (numeric or
##' character) or a formula defining the linear constraints (see also the
##' \code{regression} or \code{covariance} methods).
##' @param \dots Additional arguments
##' @usage
##' \method{intercept}{lvm}(object, vars, ...) <- value
##' @return
##'
##' A \code{lvm}-object
##' @note
##'
##' Variables will be added to the model if not already present.
##' @author Klaus K. Holst
##' @seealso \code{\link{covariance<-}}, \code{\link{regression<-}},
##' \code{\link{constrain<-}}, \code{\link{parameter<-}},
##' \code{\link{latent<-}}, \code{\link{cancel<-}}, \code{\link{kill<-}}
##' @keywords models regression
##' @export
##' @examples
##'
##'
##' ## A multivariate model
##' m <- lvm(c(y1,y2) ~ f(x1,beta)+x2)
##' regression(m) <- y3 ~ f(x1,beta)
##' intercept(m) <- y1 ~ f(mu)
##' intercept(m, ~y2+y3) <- list(2,"mu")
##' intercept(m) ## Examine intercepts of model (NA translates to free/unique paramete##r)
##'
##'
"intercept" <- function(object,...) UseMethod("intercept")

##' @export
##' @export
intercept.lvm <- intfix.lvm <- function(object,value,...) {
    if (!missing(value)) {
        intercept(object,...) <- value
        return(object)
    }
    res <- object$mean; attr(res,"type") <- "mean"
    attr(res,"exo.idx") <- index(object)$exo.idx
    attr(res,"nvar") <- length(res)
    class(res) <- "fix"
    return(res)
}

##' @export
"intercept<-" <- function(object,...,value) UseMethod("intercept<-")

##' @export
##' @export
"intercept<-.lvm" <- "intfix<-.lvm" <- function(object, vars,...,value) {
  if (!missing(vars) && inherits(value,"formula")) value <- all.vars(value)
  if (inherits(value,"formula")) {
    lhs <- getoutcome(value)
    yy <- decomp.specials(lhs)
    if ((inherits(value[[3]],"logical") && is.na(value[[3]]))) {
      intfix(object,yy) <- NA
      return(object)
    }
    tt <- terms(value)
    xf <- attributes(terms(tt))$term.labels
    res <- lapply(xf,decomp.specials)[[1]]
    myvalue <- suppressWarnings(as.numeric.list(as.list(res)))
    myvalue <- lapply(myvalue, function(x) ifelse(x=="NA",NA,x))
    intfix(object,yy) <- myvalue
    object$parpos <- NULL
    return(object)
  }
  if (inherits(vars,"formula")) {
    vars <- all.vars(vars)
  }
  object$mean[vars] <- value
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  return(object)
}

###}}} intfix

###{{{ covfix

##' @export
"covfix" <- function(object,...) UseMethod("covfix")

##' @export
covfix.lvm <- function(object,...) {
  res <- list(rel=object$cov, labels=object$covpar, values=object$covfix); attr(res,"type") <- "cov"
  attr(res,"exo.idx") <- index(object)$exo.idx
  attr(res,"nvar") <- NROW(res$rel)
  class(res) <- "fix"
  return(res)
}


##' @export
"covfix<-" <- function(object,...,value) UseMethod("covfix<-")

##' @export
"covfix<-.lvm" <- function(object, var1, var2=var1, pairwise=FALSE, exo=FALSE, ..., value) {

  if (inherits(var1,"formula")) {
      var1 <- all.vars(var1)
  }
  if (inherits(var2,"formula")) {
      var2 <- all.vars(var2)
  }
  object <- addvar(object,c(var1,var2),reindex=FALSE,...)

  allvars <- c(var1,var2)
  xorg <- exogenous(object)
  exoset <- setdiff(xorg,allvars)

  if (!exo & length(exoset)<length(xorg)) {
    exogenous(object) <- exoset
  }

  if (inherits(value,"formula")) value <- all.vars(value)
  if (pairwise) {
    p <- 0
    K <- length(var1)*(length(var1)-1)/2
    if (length(value)==1)
      value <- rep(value,K)
    if (length(value)!=K) stop("Wrong number of parameters")
    for (i in seq_len(length(var1)-1)) {
      for (j in seq(i+1,length(var1))) {
        p <- p+1
        valp <- suppressWarnings(as.numeric(value[[p]]))
        if (is.na(value[[p]]) | value[[p]]=="NA") {
          object$covfix[var1[i],var1[j]] <- object$covpar[var1[i],var1[j]] <- NA
          object$covfix[var1[j],var1[i]] <- object$covpar[var1[j],var1[i]] <- NA
        }
        else {
          object$cov[var1[i],var1[j]] <-  object$cov[var1[j],var1[i]] <- 1
          if (is.numeric(value[[p]]) | !is.na(valp)) {
            object$covfix[var1[i],var1[j]] <- object$covfix[var1[j],var1[i]] <- valp
            object$covpar[var1[i],var1[j]] <- object$covpar[var1[j],var1[i]] <- NA
          } else {
            object$covpar[var1[i],var1[j]] <- object$covpar[var1[j],var1[i]] <- value[[p]]
            object$covfix[var1[i],var1[j]] <- object$covfix[var1[j],var1[i]] <- NA
          }
        }
      }
    }
    newindex <- reindex(object)
    object$parpos <- NULL
    index(object)[names(newindex)] <- newindex
    return(object)
  }


  if (is.null(var2)) {
    if (length(value)==1)
      value <- rep(value,length(var1))
    if (length(value)!=length(var1)) stop("Wrong number of parameters")
    for (i in seq_along(var1)) {
      vali <- suppressWarnings(as.numeric(value[[i]]))
      if (is.na(value[[i]]) | value[[i]]=="NA") {
        object$covfix[var1[i],var1[i]] <- object$covpar[var1[i],var1[i]] <- NA
      }
      else {
        if (is.numeric(value[[i]]) | !is.na(vali)) {
          object$covfix[var1[i],var1[i]] <- vali
          object$covpar[var1[i],var1[i]] <- NA
        } else {
          object$covfix[var1[i],var1[i]] <- NA
          object$covpar[var1[i],var1[i]] <- value[[i]]
        }
      }
    }
    newindex <- reindex(object)
    object$parpos <- NULL
    index(object)[names(newindex)] <- newindex
    return(object)
  }

  if (length(var1)==length(var2) & length(var1)==length(value)) {
    p <- 0
    for (i in seq_along(var1)) {
      p <- p+1
      valp <- suppressWarnings(as.numeric(value[[p]]))
      if (is.na(value[[p]]) | value[[p]]=="NA") {
        object$covfix[var1[i],var2[i]] <- object$covpar[var1[i],var2[i]] <- NA
        object$covfix[var2[i],var1[i]] <- object$covpar[var2[i],var1[i]] <- NA
      }
      else {
        object$cov[var1[i],var2[i]] <-  object$cov[var2[i],var1[i]] <- 1
        if (is.numeric(value[[p]]) | !is.na(valp)) {
          object$covfix[var1[i],var2[i]] <- object$covfix[var2[i],var1[i]] <- valp
          object$covpar[var1[i],var2[i]] <- object$covpar[var2[i],var1[i]] <- NA
        } else {
          object$covpar[var1[i],var2[i]] <- object$covpar[var2[i],var1[i]] <- value[[p]]
          object$covfix[var1[i],var2[i]] <- object$covfix[var2[i],var1[i]] <- NA
        }
      }
    }
    newindex <- reindex(object)
    object$parpos <- NULL
    index(object)[names(newindex)] <- newindex
    return(object)
  }


  K <- length(var1)*length(var2)
  if (length(value)==1)
    value <- rep(value,K)
  if (length(value)!=K) stop("Wrong number of parameters")

  p <- 0
  for (i in seq_along(var1)) {
    for (j in seq_along(var2)) {
      if (!pairwise | var1[i]!=var2[j]) {
        p <- p+1
        valp <- suppressWarnings(as.numeric(value[[p]]))
        if (is.na(value[[p]]) | value[[p]]=="NA") {
          object$covfix[var1[i],var2[j]] <- object$covpar[var1[i],var2[j]] <- NA
          object$covfix[var2[j],var1[i]] <- object$covpar[var2[j],var1[i]] <- NA
        }
        else {
          object$cov[var1[i],var2[j]] <-  object$cov[var2[j],var1[i]] <- 1
          if (is.numeric(value[[p]]) | !is.na(valp)) {
            object$covfix[var1[i],var2[j]] <- object$covfix[var2[j],var1[i]] <- valp
            object$covpar[var1[i],var2[j]] <- object$covpar[var2[j],var1[i]] <- NA
          } else {
            object$covpar[var1[i],var2[j]] <- object$covpar[var2[j],var1[i]] <- value[[p]]
            object$covfix[var1[i],var2[j]] <- object$covfix[var2[j],var1[i]] <- NA
          }
        }
      }
    }
  }
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  return(object)
}

###}}} covfix

###{{{ regfix

##' @export
"regfix" <- function(object,...) UseMethod("regfix")

##' @export
regfix.lvm <- function(object,...) {
  res <- list(rel=index(object)$M, labels=object$par, values=object$fix); attr(res,"type") <- "reg"
  attr(res,"exo.idx") <- index(object)$exo.idx
  attr(res,"nvar") <- NROW(res$rel)
  class(res) <- "fix"
  return(res)
}

##' @export
"regfix<-" <- function(object,...,value) UseMethod("regfix<-")

##' @export
"regfix<-.lvm" <- function(object, to, from, exo=lava.options()$exogenous, variance, y,x, ..., value) {
    if (!missing(y)) {
        if (inherits(y,"formula")) y <- all.vars(y)
        to <- y
    }
    if (!missing(x)) {
        if (inherits(x,"formula")) x <- all.vars(x)
        from <- x
    }
    if (is.null(to)) stop("variable list needed")
    
  if (inherits(to,"formula")) {
      val <- procformula(object,to,exo=exo)
      object <- val$object
      ys <- val$ys
      xs <- val$xs      
      if (!missing(variance))
          covariance(object,ys) <- variance      
      to <- ys; from <- xs 
  } else {
      object <- addvar(object,c(to,from),reindex=FALSE,...)
      newexo <- from
      notexo <- to
      curvar <- index(object)$var  
      if (exo) {
          oldexo <- exogenous(object)
          newexo <- setdiff(newexo,c(notexo,curvar))
          exogenous(object) <- union(newexo,setdiff(oldexo,notexo))
      }
  }
    
  if (inherits(value,"formula")) value <- all.vars(value)

  if (length(from)==length(to) & length(from)==length(value)) {
    for (i in seq_along(from)) {
      if (object$M[from[i],to[i]]==0) { ## Not adjacent! ##!isAdjacent(Graph(object), from[i], to[i])) {
        object <- regression(object, to=to[i], from=from[i])
      }
      vali <- suppressWarnings(as.numeric(value[[i]]))
      if (is.na(value[[i]]) | value[[i]]=="NA") {
        object$fix[from[i],to[i]] <- object$par[from[i],to[i]] <- NA
      }
      else {
        if (is.numeric(value[[i]]) | !is.na(vali)) {
          object$fix[from[i],to[i]] <- vali
          object$par[from[i],to[i]] <- NA
        } else {
          object$par[from[i],to[i]] <- value[[i]]
          object$fix[from[i],to[i]] <- NA
        }
      }
    }
    newindex <- reindex(object)
    object$parpos <- NULL
    index(object)[names(newindex)] <- newindex
    return(object)
  }

  for (i in from) {
    for (j in to) {
      if (object$M[i,j]==0) { ##!isAdjacent(Graph(object), i, j)) {
        object <- regression(object,to=j,from=i)
      }
    }
  }

  K <- length(from)*length(to)
  if (length(value)==1)
    value <- rep(value,K)
  if (length(value)!=K) stop("Wrong number of parameters")

  for (j in seq_along(to)) {
    for (i in seq_along(from)) {
      p <- (j-1)*length(from) + i
      valp <- suppressWarnings(as.numeric(value[[p]]))
      if (is.na(value[[p]]) | value[[p]]=="NA")
        object$fix[from[i],to[j]] <- object$par[from[i],to[j]] <- NA
      else {
        if (is.numeric(value[[p]]) | !is.na(valp)) {
        object$fix[from[i],to[j]] <- valp
        object$par[from[i],to[j]] <- NA
      } else {
        object$par[from[i],to[j]] <- value[[p]]
        object$fix[from[i],to[j]] <- NA
      }
      }
    }
  }
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  index(object) <- reindex(object)
  return(object)
}

###}}} regfix

###{{{ parfix

##' @export
"parfix<-" <- function(x,...,value) UseMethod("parfix<-")

##' @export
"parfix<-.lvm" <- function(x,idx,...,value) {
  parfix(x,idx,value,...)
}

##' @export
"parfix" <- function(x,...) UseMethod("parfix")


## m <- lvm(c(y[m:v]~b*x))
## constrain(m,b~a) <- base::identity

##' @export
parfix.lvm <- function(x,idx,value,fix=FALSE,...) {
  object <- Model(x)
  if (fix)
    object <- fixsome(object)
  if (length(idx)!=length(value))
    value <- rep(value,length.out=length(idx))
  value <- as.list(value)
  I <- index(object)
  V <- with(I, matrices2(Model(object), seq_len(npar.mean+npar+npar.ex)))
  V$A[I$M0!=1] <- 0; V$P[I$P0!=1] <- 0
  v.fix <- which(V$v%in%idx) ## Intercepts
  vval <- V$v[v.fix]
  v.ord <- match(vval,idx)
  Pval <- V$P[V$P%in%idx] ## Variance/covariance
  P.fix <- which(matrix(V$P%in%idx,nrow=nrow(V$P)),arr.ind=TRUE)
  P.ord <- match(Pval,idx)
  Aval <- V$A[which(V$A%in%idx)] ## Regression parameters
  A.fix <- which(matrix(V$A%in%idx,nrow=nrow(V$A)),arr.ind=TRUE)
  A.ord <- match(Aval,idx)
  e.fix <- which(V$e%in%idx)
  eval <- V$e[e.fix]
  e.ord <- match(eval,idx)
  for (i in seq_len(length(e.fix))) {
      object$exfix[[e.fix[i]]] <- value[[e.ord[i]]]
  }
  for (i in seq_len(length(v.fix))) {
      object$mean[[v.fix[i]]] <- value[[v.ord[i]]]
  }
  for (i in seq_len(nrow(A.fix))) {
      if (is.numeric(value[[ A.ord[i] ]])){
          object$fix[A.fix[i,1],A.fix[i,2]] <- value[[A.ord[i]]]
          object$par[A.fix[i,1],A.fix[i,2]] <- NA
      } else {
          object$par[A.fix[i,1],A.fix[i,2]] <- value[[A.ord[i]]]
          object$fix[A.fix[i,1],A.fix[i,2]] <- NA
      }
  }
  for (i in seq_len(nrow(P.fix))) {
      if (is.numeric(value[[ P.ord[i] ]])) {
          object$covfix[P.fix[i,1],P.fix[i,2]] <- value[[P.ord[i]]]
          object$covpar[P.fix[i,1],P.fix[i,2]] <- NA
      } else {
          object$covpar[P.fix[i,1],P.fix[i,2]] <- value[[P.ord[i]]]
          object$covfix[P.fix[i,1],P.fix[i,2]] <- NA
      }
  }
  newindex <- reindex(object)
  object$parpos <- NULL
  index(object)[names(newindex)] <- newindex
  attributes(object)$fixed <- list(v=v.fix,A=A.fix,P=P.fix,e=e.fix)
  return(object)
}

###}}} parfix
