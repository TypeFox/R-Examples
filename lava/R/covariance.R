
##' Add covariance structure to Latent Variable Model
##'
##' Define covariances between residual terms in a \code{lvm}-object.
##'
##' The \code{covariance} function is used to specify correlation structure
##' between residual terms of a latent variable model, using a formula syntax.
##'
##' For instance, a multivariate model with three response variables,
##'
##' \deqn{Y_1 = \mu_1 + \epsilon_1}
##'
##' \deqn{Y_2 = \mu_2 + \epsilon_2}
##'
##' \deqn{Y_3 = \mu_3 + \epsilon_3}
##'
##' can be specified as
##'
##' \code{m <- lvm(~y1+y2+y3)}
##'
##' Pr. default the two variables are assumed to be independent. To add a
##' covariance parameter \eqn{r = cov(\epsilon_1,\epsilon_2)}, we execute the
##' following code
##'
##' \code{covariance(m) <- y1 ~ f(y2,r)}
##'
##' The special function \code{f} and its second argument could be omitted thus
##' assigning an unique parameter the covariance between \code{y1} and
##' \code{y2}.
##'
##' Similarily the marginal variance of the two response variables can be fixed
##' to be identical (\eqn{var(Y_i)=v}) via
##'
##' \code{covariance(m) <- c(y1,y2,y3) ~ f(v)}
##'
##' To specify a completely unstructured covariance structure, we can call
##'
##' \code{covariance(m) <- ~y1+y2+y3}
##'
##' All the parameter values of the linear constraints can be given as the right
##' handside expression of the assigment function \code{covariance<-} if the
##' first (and possibly second) argument is defined as well. E.g:
##'
##' \code{covariance(m,y1~y1+y2) <- list("a1","b1")}
##'
##' \code{covariance(m,~y2+y3) <- list("a2",2)}
##'
##' Defines
##'
##' \deqn{var(\epsilon_1) = a1}
##'
##' \deqn{var(\epsilon_2) = a2}
##'
##' \deqn{var(\epsilon_3) = 2}
##'
##' \deqn{cov(\epsilon_1,\epsilon_2) = b1}
##'
##' Parameter constraints can be cleared by fixing the relevant parameters to
##' \code{NA} (see also the \code{regression} method).
##'
##' The function \code{covariance} (called without additional arguments) can be
##' used to inspect the covariance constraints of a \code{lvm}-object.
##'
#
##'
##' @aliases covariance covariance<- covariance.lvm covariance<-.lvm
##' covfix<- covfix covfix<-.lvm covfix.lvm
##' variance variance<- variance.lvm variance<-.lvm
##' @param object \code{lvm}-object
##' @param var1 Vector of variables names (or formula)
##' @param var2 Vector of variables names (or formula) defining pairwise
##' covariance between \code{var1} and \code{var2})
##' @param constrain Define non-linear parameter constraints to ensure positive definite structure
##' @param pairwise If TRUE and \code{var2} is omitted then pairwise correlation is added between all variables in \code{var1}
##' @param \dots Additional arguments to be passed to the low level functions
##' @param value List of parameter values or (if \code{var1} is unspecified)
##' @usage
##' \method{covariance}{lvm}(object, var1=NULL, var2=NULL, constrain=FALSE, pairwise=FALSE,...) <- value
##' @return A \code{lvm}-object
##' @author Klaus K. Holst
##' @seealso \code{\link{regression<-}}, \code{\link{intercept<-}},
##' \code{\link{constrain<-}} \code{\link{parameter<-}}, \code{\link{latent<-}},
##' \code{\link{cancel<-}}, \code{\link{kill<-}}
##' @keywords models regression
##' @export
##' @examples
##'
##' m <- lvm()
##' ### Define covariance between residuals terms of y1 and y2
##' covariance(m) <- y1~y2
##' covariance(m) <- c(y1,y2)~f(v) ## Same marginal variance
##' covariance(m) ## Examine covariance structure
##'
##'
`covariance` <- function(object,...) UseMethod("covariance")

##' @export
"variance<-" <- function(object,...,value) UseMethod("covariance<-")

##' @export
`variance` <- function(object,...) UseMethod("variance")

##' @export
"variance.lvm" <- function(object,...) covariance(object,...)

##' @export
"variance<-.lvm" <- function(object,...,value) {
    covariance(object,...) <- value
    return(object)
}

##' @export
"covariance<-" <- function(object,...,value) UseMethod("covariance<-")

##' @export
"covariance<-.lvm" <- function(object, var1=NULL, var2=NULL, constrain=FALSE, pairwise=FALSE, ..., value) {

  if (!is.null(var1)) {
    if (inherits(var1,"formula")) {
      lhs <- getoutcome(var1)
      xf <- attributes(terms(var1))$term.labels
      xx <- unlist(lapply(xf, function(x) x[1]))
      if (length(lhs)==0) {
        covfix(object,var1,var2,pairwise=pairwise,...) <- value
        object$parpos <- NULL
        return(object)
      }
      else {
        yy <- decomp.specials(lhs)
      }
    } else {
      yy <- var1; xx <- var2
    }
    covfix(object,var1=yy,var2=xx,pairwise=pairwise,...) <- value
    object$parpos <- NULL
    return(object)
  }

  if (is.list(value)) {
    for (v in value) {
        covariance(object,pairwise=pairwise,constrain=constrain,...) <- v
      }
      return(object)
  }

  if (inherits(value,"formula")) {
    lhs <- getoutcome(value)
    if (length(lhs)==0) {
      return(covariance(object,all.vars(value),constrain=constrain,pairwise=pairwise,...))
    }
    yy <- decomp.specials(lhs)

    tt <- terms(value, specials=c("f","v"))
    xf <- attributes(terms(tt))$term.labels
    res <- lapply(xf,decomp.specials)
    nx <- length(xf)
    if (nx==1) {
      if(is.null(attr(tt,"specials")$f) | length(res[[1]])<2) {
        if(is.null(attr(tt,"specials")$v) & is.null(attr(tt,"specials")$f))

          {
            for (i in yy)
              for (j in res[[1]])
                object <- covariance(object, c(i,j), pairwise=TRUE, constrain=constrain, ...)
          } else {
            covfix(object,var1=yy,var2=NULL) <- res[[1]]
          }
      } else {
        covfix(object,var1=yy,var2=res[[1]][1]) <- res[[1]][2]
      }
      object$parpos <- NULL
      return(object)
    }

    xx <- unlist(lapply(res, function(z) z[1]))
    for (y in yy)
      for (i in seq_along(xx)) {
        if (length(res[[i]])>1) {
          covfix(object, var1=y, var2=res[[i]][1]) <- res[[i]][2]
        } else if ((i+1)%in%attr(tt,"specials")$f | (i+1)%in%attr(tt,"specials")$v) {
          covfix(object, var1=y, var2=NULL) <- res[[i]]
        } else {
          object <- covariance(object,c(y,xx[i]),pairwise=TRUE,...)
        }
      }

    object$parpos <- NULL
    return(object)
  }
  else covariance(object,value,pairwise=pairwise,...)
}

##' @export
`covariance.lvm` <-
    function(object,var1=NULL,var2=NULL,exo=FALSE,pairwise=FALSE,constrain=FALSE,value,...) {

        if (!missing(value)) {
            covariance(object,var1=var1,var2,exo=exo,pariwise=pairwise,constrain=constrain,...) <- value
            return(object)
        }
        if (!is.null(var1)) {
            if (inherits(var1,"formula")) {
                covariance(object,constrain=constrain,
                 pairwise=pairwise,exo=exo,...) <- var1
                return(object)
            }
    allvars <- var1
    if (!missing(var2)) {
      if (inherits(var2,"formula"))
        var2 <- all.vars(var2)
      allvars <- c(allvars,var2)
    }
    if (constrain) {
      if (length(allvars)!=2) stop("Constraints only implemented for pairs")
      return(covarianceconst(object,allvars[1],allvars[2],...))
    }

    object <- addvar(object, allvars, silent=TRUE, reindex=FALSE)

    xorg <- exogenous(object)
    exoset <- setdiff(xorg,allvars)
    if (!exo & length(exoset)<length(xorg)) {
      exogenous(object) <- exoset
    }

    if (!missing(var2)) {
      for (i in seq_len(length(var1))) {
        c1 <- var1[i]
        for (j in seq_len(length(var2))) {
          c2 <- var2[j]
          object$cov[c1,c2] <- object$cov[c2,c1] <- 1
          object$parpos <- NULL
          index(object) <- reindex(object)
        }
      }
    }
    else {
      if (pairwise) {
        for (i in seq_len(length(var1))) {
          c1 <- var1[i]
            for (j in seq_len(length(var1))) {
              c2 <- var1[j]
              object$cov[c1,c2] <- object$cov[c2,c1] <- 1
              object$parpos <- NULL
              index(object) <- reindex(object)
            }
        }
      }
    }
    return(object)
  }
  else
    return(covfix(object))
}

covarianceconst <- function(object,var1,var2,cname=NA,rname=NA,vname=NA,v2name=vname,lname=NA,l2name=lname,...) {
  if (inherits(var1,"formula")) {
    var1 <- getoutcome(var1)
    var2 <- attributes(var1)$x
  }
  curpar <- parlabels(object)
  if (is.na(cname)) {
    cname <- object$covpar[var1,var2]
  }

  if (is.na(v2name)) {
    v2name <- object$covpar[var2,var2]
  }
  if (is.na(vname)) {
    vname <- object$covpar[var1,var1]
  }

  nvarname <- c("rname","cname","vname","v2name","lname","l2name")[is.na(c(rname,cname,vname,v2name,lname,l2name))]

  nprefix <- sapply(nvarname, function(x) substr(x,1,1))
  for (i in seq_len(length(nvarname))) {
    count <- 0
    repeat {
      count <- count+1
      curname <- paste0(nprefix[i],count)
      if (!(curname%in%curpar)) break;
    }
    curpar <- c(curname,curpar)
    assign(nvarname[i],curname)
  }
  covariance(object,c(var1,var2)) <- c(vname,v2name)
  ff <- function(x) exp(x)
  constrain(object,vname,lname) <- ff
  if (vname!=v2name)
    constrain(object,v2name,l2name) <- ff
  covariance(object,var1,var2) <- cname
  cpar <- unique(c(vname,v2name,rname))
  constrain(object,cname,cpar) <- function(x)
    prod((x[seq(length(cpar)-1)]))^(1/(length(cpar)-1))*tanh(x[length(cpar)])

  return(structure(object,rname=rname,cname=cname))
}
