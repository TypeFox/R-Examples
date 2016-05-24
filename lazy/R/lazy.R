# ---------------------------------------- -*- mode: r; mode: font-lock -*- #
# lazy.R                                 Lazy learning for local regression #
# ------------------------------------------------------------------------- #
                                                                             
# ========================================================================= #
# Lazy learning for local regression                                        #
# ------------------------------------------------------------------------- #
# Copyright (C) 1999, 2003 Mauro Birattari and Gianluca Bontempi            #
# ========================================================================= #
# This program is free software; you can redistribute it and/or modify it   #
# under the terms of the GNU General Public License as published by the     #
# Free Software Foundation; either version 2 of the License, or (at your    #
# option) any later version.                                                #
#                                                                           #
# This program is distributed in the hope that it will be useful, but       #
# WITHOUT ANY WARRANTY; without even the implied warranty of                #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         #
# General Public License for more details.                                  #
#                                                                           #
# You should have received a copy of the GNU General Public License along   #
# with this program; if not, write to the Free Software Foundation, Inc.,   #
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.                  #
# ========================================================================= #

# ========================================================================= #
#             Mauro Birattari                     Gianluca Bontempi         #
#                IRIDIA                     Departement d'Informatique      #
#    Universite' Libre de Bruxelles       Universite' Libre de Bruxelles    #
#             mbiro@ulb.ac.be                     gbonte@ulb.ac.be          #
# ========================================================================= #

# $Id: lazy.R,v 1.8 2003/12/09 15:32:18 mbiro Exp $ #


lazy<-function(formula, data=NULL, weights, subset, na.action,
                control = lazy.control(...),...){
  mt <- terms(formula, data = data)
  mf <- match.call(expand.dots=FALSE)
  mf$control <- mf$... <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  na.act <- attr(mf, "na.action")
  y <- model.response(mf, "numeric")
  w <- model.weights(mf)
  if(is.null(w))
    w <- rep(1, length(y))
  else
    warning("Parameter `weights' will be ignored")
  nmx <- as.character(attr(mt, "variables"))[-(1:2)]
  x <- mf[, nmx, drop=FALSE]
  if(any(sapply(x, is.factor))) stop("predictors must all be numeric")
  x <- as.matrix(x)
  D <- ncol(x)
  nmx <- colnames(x)
  names(nmx) <- nmx

  obj<-control
  obj$call <- match.call()
  obj$terms <- mt
  obj$xnames <- nmx
  obj$x <- x
  obj$y <- y
  obj$weights <- w
  if(!is.null(na.act)) obj$na.action <- na.act
  class(obj)<-"lazy"
  obj
}


print.lazy<-function(x, digits=max(3, getOption("digits")-3), ...)
{
    if(!is.null(cl <- x$call)) {
	cat("Call:\n")
	dput(cl)
    }
    cat("\nNumber of observations:",nrow(x$x),"\n")
    invisible(x)
}



lazy.control <- function(conIdPar=NULL, linIdPar=1, quaIdPar=NULL,
                         distance=c("manhattan","euclidean"),
                         metric=NULL, cmbPar=1, lambda=1e+06 )
  list(conIdPar=conIdPar,linIdPar=linIdPar,quaIdPar=quaIdPar,
       distance=match.arg(distance),metric=metric,
       cmbPar=cmbPar,lambda=lambda)


summary.lazy <- function(object,...){
  class(object) <- "summary.lazy"
  object
}

print.summary.lazy <- function(x, digits=max(3, getOption("digits")-3), ...)
{
    if(!is.null(cl <- x$call)){
	cat("Call:\n")
	dput(cl)
    }
    cat("\nNumber of Observations:",nrow(x$x),"\n")
    invisible(x)
}

predict.lazy <- function(object, newdata = NULL,t.out=FALSE,k.out=FALSE,
                         S.out=FALSE,T.out=FALSE,I.out=FALSE,...){
    if(!inherits(object, "lazy"))
	stop("First argument must be a lazy object")
    if(is.null(newdata)) return(object)

    vars <- as.character(attr(delete.response(terms(object)),
                              "variables"))[-1]
    nvars<-length(vars)
    varstring<-if (nvars==1) vars
               else if (nvars==2) paste(vars,collapse="and")
               else paste(paste(vars[-nvars],collapse=", "),
                          ", and ",vars[nvars],sep="")

    if(is.vector(newdata)&&(length(newdata)==ncol(object$x))){
      warning(paste("   Coercing ", deparse(substitute(newdata)),
                    " to a 1x",length(newdata),
                    " data frame\n",
                    "   with columns named ",varstring,
                    "\n",sep=""))
      newdata<-matrix(newdata,1)
    }else if(is.matrix(newdata)&&(ncol(newdata)==ncol(object$x))){
      warning(paste("   Coercing ", deparse(substitute(newdata)),
                    " to a ", nrow(newdata),
                    "x", ncol(newdata),
                    " data frame\n",
                    "   with columns named ",varstring,
                    "\n",sep=""))
    }else{
      if(length(vars) > 1 || NCOL(newdata) > 1) {
        if(any(!match(vars, colnames(newdata), FALSE)))
          stop("newdata does not contain the variables needed")
        newdata<-as.matrix(newdata[, vars, drop=FALSE])
      }else{
        newdata<-as.matrix(newdata)
      }
    }
    
    .Call("lazy",object$x,object$y,newdata,
          object$conIdPar,object$linIdPar,object$quaIdPar,
          object$distance,object$metric,object$cmbPar,object$lambda,
          t.out,k.out,S.out,T.out,I.out,PACKAGE="lazy")
}

