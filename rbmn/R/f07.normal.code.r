
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cor4var <- function(ma)
#TITLE returns the correlation matrix from the variance
#DESCRIPTION returns the correlation matrix from the variance
# preserving possible variable names
#DETAILS
# Zero variances are detected and accepted (all associated correlation
# coefficients are forced to be zero.>>
#KEYWORDS 
#INPUTS
#{ma}<< The variance matrix.>>
#[INPUTS]
#VALUE
# The correlation matrix
#EXAMPLE
# cor4var(rbmn0mn.04$gamma);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_12_05
#REVISED 12_07_31
#--------------------------------------------
{
  dd <- 1/sqrt(diag(ma));
  dd[diag(ma)<=0] <- 0;
  dd <- diag(dd,nrow=length(dd));
  res <- dd %*% ma %*% dd;
  dimnames(res) <- dimnames(ma);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
var2pre <- function(ma)
#TITLE returns the precision matrix from the variance
#DESCRIPTION returns the precision matrix from the variance
# preserving possible variable names
#DETAILS
# Non full rank matrices are accepted, a generalized inverse
# is returned and a warning is issued.
#KEYWORDS 
#INPUTS
#{ma}<< The variance matrix.>>
#[INPUTS]
#VALUE
# The precision matrix
#EXAMPLE
# var2pre(rbmn0mn.04$gamma);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_12_17
#REVISED 11_12_17
#--------------------------------------------
{
  # constant
  eps <- 10^-10;
  # a look to the rank
  sivadi <- eigen(ma,symmetric=TRUE);
  nbv <- sum(sivadi$values > eps);
  if (nbv < nrow(ma)) {
    # issuing the warning
    r.erreur(NULL,message=
           paste("The matrix was found singular (rank = ",
                 nbv," < ",nrow(ma),") then a generalized inverse was provided",
                 sep=""));
    # computing the generalized inverse
    res <- sivadi$vectors[,r.bc(nbv)] %*%
           diag(1/sivadi$values[r.bc(nbv)],nrow=nbv,ncol=nbv) %*%
           t(sivadi$vectors[,r.bc(nbv)]);
  } else {
    res <- solve(ma);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
condi4joint <- function(mn,par,pour,x2=rep(0,length(pour))) 
#TITLE computes some conditional distribution of
# a multinormal vector
#DESCRIPTION returns the expectation and
# variance of a sub-vector conditionned
# with another (non overlapping) sub-vector
# from an initial random vector described 
# by \code{mn}.
#DETAILS
# when no names are given to \code{mn$mu},
# \code{par} and \code{pour} are supposed
# containing indices and default sequential
# names are provided.
#KEYWORDS 
#INPUTS
#{mn}<< list defining the distribution of
#       the initial vector with \code{$mu},
#       its expectation, and \code{$gamma},
#       its variance matrix.>>
#{par}<< names (or indices) of the sub-vector
#        to give the distribution.>>
#{pour}<< names (or indices) of the conditionning
#         sub-vector (can be \code{NULL} when
#         for non conditionning.>>
#
#[INPUTS]
#{x2} <<values to consider for the conditioning
#       sub-vector. When \code{NULL} the general
#       form is supplied, not a /mn/ object.  >>
#VALUE
# A list:\cr
# when \code{x2} provides the values taken by the
#      conditioning part, it is a /mn/ object with its
#      two components: \code{$mu} for the expectation vector
#      and \code{$gamma} for the variance matrix.\cr
# when \code{x2} is \code{NULL} the list has got three
#      components: \code{$a} for the fixed part of the
#      expectation vector, \code{$b} for the regression
#      coefficients to be associated to the non precised
#      \code{x2} values, varying part of the expectation
#      and \code{$gamma} for the variance matrix.\cr
#EXAMPLE
# print8mn(condi4joint(rbmn0mn.04,c("1.1","2.2","1.2","2.1"),NULL));
# print8mn(condi4joint(rbmn0mn.04,c("1.1","2.2","1.2","2.1"),"C",0));
# print(condi4joint(rbmn0mn.04,c("1.1","2.2","1.2","2.1"),"C",NULL));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_12_05
#REVISED 12_07_31
#--------------------------------------------
{
  # lengths of the involved vectors
  n <- length(mn$mu); n1 <- length(par); n2 <- length(pour);
  if (n1+n2>n) { stop("mu.parpour");}
  # naming if necessary
  if (is.null(names(mn$mu))) {
    va <- paste("V",as.character(1:n),sep="");
    par <- va[par];
    pour <- va[pour];
    names(mn$mu) <- va;
    dimnames(mn$gamma) <- list(va,va);
  }
  if (length(intersect(par,pour)) > 0) { stop("parpour");}
  #  
  if (n1==0) {
    if (is.null(x2)) {
      res <- list(mu=numeric(n1),
                  rho=matrix(NA,n1,n2,dimnames=list(par,pour)),
                  gamma=matrix(NA,n1,n1)
                 );
    } else {
      res <- list(mu=numeric(n1),
                  gamma=matrix(NA,n1,n1)
                 );
    }
  } else {
    if (n2==0) {
      if (is.null(x2)) {
        res <- list(mu=mn$mu[par],
                    rho=matrix(NA,n1,n2,dimnames=list(par,pour)),
                    gamma=mn$gamma[par,par]);
      } else {
        res <- list(mu=mn$mu[par],gamma=mn$gamma[par,par]);
      }
    } else {
      mu1 <- mn$mu[par];
      mu2 <- mn$mu[pour];
      s11 <- mn$gamma[par,par,drop=FALSE];
      s12 <- mn$gamma[par,pour,drop=FALSE];
      s22 <- mn$gamma[pour,pour,drop=FALSE];
      ss12 <- s12 %*% solve(s22);
      si <- s11 - ss12 %*% t(s12);
      if (is.null(x2)) {
        ac <- mu1 - ss12 %*% matrix(mu2,n2,1);
        ac <- as.vector(ac);
        names(ac) <- par;
        res <- list(mu=ac,
                    rho=ss12,
                    gamma=si);
      } else {
        mu <- mu1 + ss12 %*% matrix(as.numeric(x2-mu2),length(pour),1);
        mun <- as.vector(mu);
        names(mun) <- par;
        res <- list(mu=mun,gamma=si);
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
mn4joint1condi <- function(lmar,lcon) 
#TITLE computes a joint distribution from a marginal and a conditional
#  one for multinormal distributions
#DESCRIPTION returns the expectation and
# variance of the multinormal normal distribution
# defined through a marginal subcomponent and
# a conditional distribution.
#DETAILS
# The conditional distribution is defined
# with a list having \code{$a} for the constant
# part of the expectation; \code{$b} for the
# regression coefficient part of the expectation;
# and \code{$S} for the residual variance matrix.
#KEYWORDS 
#INPUTS
#{lmar}<< list defining the distribution of
#       the marginal part with \code{$mu},
#       its expectation, and \code{$gamma},
#       its variance matrix (in fact a /mn/ object).>>
#{lcon}<< list defining the distribution of
#       the conditional part (see the \emph{Details} section).>>
#[INPUTS]
#VALUE
# A list:
#{\$mu}<<The expectation vector.>>
#{\$gamma}<<The joint variance matrix.>>
# that is a /mn/ object.
#EXAMPLE
# lcon <- list(a=c(D=2,E=4),
#              b=matrix(1:6,2,dimnames=list(LETTERS[4:5],
#                                           LETTERS[1:3])),
#              S=matrix(c(1,1,1,2),2));
#
# print8mn(mn4joint1condi(rbmn0mn.01,lcon));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE (after the paper publication: change the components names of
# \code{lcon} to be consistent with the output of function
# \code{condi4joint} when regression coefficients are outputted.
#AUTHOR J.-B. Denis
#CREATED 11_12_14
#REVISED 11_12_14
#--------------------------------------------
{
  # getting constants
  nv <- length(lmar$mu);
  nu <- length(lcon$a);
  # degenerate cases
  if (nu == 0) {
    return(lmar);
  }
  if (nv == 0) {
    return(list(mu=lcon$a,gamma=lcon$S));
  }
  # non matrix cases
  if (nv == 1) {
    if (!is.matrix(lmar$gamma)) {
      lmar$gamma <- matrix(lmar$gamma,nv,nv);
    }
    if (!is.matrix(lcon$b)) {
      lcon$b <- matrix(lcon$b,nu,nv);
    }
  }
  if (nu == 1) {
    if (!is.matrix(lcon$S)) {
      lcon$S <- matrix(lcon$S,nu,nu);
    }
    if (!is.matrix(lcon$b)) {
      lcon$b <- matrix(lcon$b,nu,nv);
    }
  }
  # general case
  if (is.null(names(lmar$mu))) {
    vam <- paste("M",as.character(1:nv),sep="");
    names(lmar$mu) <- vam;
    dimnames(lmar$gamma) <- list(vam,vam);
  } else {
    vam <- names(lmar$mu);
  }
  if (is.null(names(lcon$a))) {
    vac <- paste("C",as.character(1:nu),sep="");
    names(lcon$a) <- vac;
    dimnames(lcon$b) <- list(vac,vam);
    dimnames(lcon$S) <- list(vac,vac);
  } else {
    vac <- names(lcon$a);
  }
  # a limited check
  if (length(intersect(vac,vam)) > 0) { stop("Overlap!");}
  # the computation
  n <- nv+nu;
  mu <- lcon$a + lcon$b %*% lmar$mu;
  mu <- c(lmar$mu,mu);
  s22 <- lcon$S + lcon$b %*% lmar$gamma %*% t(lcon$b);
  s12 <- lcon$b %*% lmar$gamma;
  gamma <- rbind(t(s12),s22);
  gamma <- cbind(rbind(lmar$gamma,s12),
                 gamma);
  names(mu) <- c(vam,vac);
  dimnames(gamma) <- list(c(vam,vac),c(vam,vac));
  #
  res <- list(mu=mu,gamma=gamma);
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
simulate8mn <- function(mn,nbs,tol=1e-7) 
#TITLE simulates a multinormal vector
#DESCRIPTION returns a matrix of simulated
# values with the variable in columns and the 
# simulations in rows.
#DETAILS
# Just a call to the basic function \code{mvrnorm}.
# Names of the variables are taken from those of
# \code{mn$mu}, when these does not exist, standard
# ones are provided.
#KEYWORDS 
#INPUTS
#{mn}<< list defining the distribution of
#       the initial vector with \code{$mu},
#       its expectation, and \code{$gamma},
#       its variance matrix.>>
#{nbs}<< number of simulations to return.>>
#[INPUTS]
#{tol} << tolerance value to be transmitted
#         to \code{mvrnorm}.>>
#VALUE
# A matrix/data frame of size : \code{nbs x length(mn$mu)}
#EXAMPLE
# print(simulate8mn(rbmn0mn.01,12));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_04_27
#REVISED 13_02_06
#--------------------------------------------
{
  # number of variables and their names
  nbv <- length(mn$mu);
  # 
  if (is.null(names(mn$mu))) {
    va <- paste("V",as.character(r.bc(nbv)),sep="");
  } else {
    va <- names(mn$mu);
  }
  # number of simulations
  nbs <- round(max(0,nbs));
  # simulating
  if (nbv*nbs > 1) {
    res <- mvrnorm(nbs,mn$mu,mn$gamma,tol=tol);
    if (nbs == 1) {
      res <- matrix(res,nbs,nbv);
    }
  } else {
    res <- matrix(NA,nbs,nbv);
  }
  # adding the variable names
  dimnames(res) <- list(NULL,va);
  # returning
  as.data.frame(res);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
simulate8gmn <- function(loi,cova,nbs,tol=1e-7) 
#TITLE simulates a multinormal vector with varying expectation
#DESCRIPTION returns a matrix of simulated
# values with the variable in columns and the 
# simulations in rows.
#DETAILS
# Just a call to the function \code{simulate8mn},
# adding the terms to the expectation due to the regression...
#KEYWORDS 
#INPUTS
#{loi}<< list defining the distribution of
#       the initial vector with \code{$mu},
#       its expectation, \code{$gamma},
#       its variance matrix and \code{$rho} a
#       matrix of regression coefficients for
#       the covariables modifying the expectation.>>
#{cova}<< Values to give to the covariables.
#   Must be a matrix with \code{nbs} rows and \code{ncol(loi$rho)}
#   columns or a vector with \code{ncol(loi$rho)} values to be
#   used for all simulations (i.e to replace a matrix with
#   identical rows..>>
#{nbs}<< number of simulations to return.>>
#[INPUTS]
#{tol} << tolerance value to be transmitted
#         to \code{mvrnorm}.>>
#VALUE
# A matrix of size : \code{nbs x length(loi$mu)}
#EXAMPLE
# loi <- list(mu=c(D=2,E=4),
#             rho=matrix(1:6,2,dimnames=list(LETTERS[4:5],
#                                           LETTERS[1:3])),
#             gamma=matrix(c(1,1,1,2),2));
# cova <- matrix(runif(36),12,dimnames=list(NULL,LETTERS[1:3]));
# print(simulate8gmn(loi,cova,12));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_02_06
#REVISED 13_02_18
#--------------------------------------------
{
  # simulating the constant part
  res <- simulate8mn(loi,nbs,tol);
  # adding the regression part
  if (!is.null(loi$rho)) { if (ncol(loi$rho) > 0) {
    nco <- dimnames(loi$rho)[[2]];
    if (!is.null(nco)) {
      nco <- r.bc(ncol(loi$rho));
    }
    if (is.matrix(cova)) {
      if (ncol(cova) != length(nco)) {
        stop("'cova' and 'loi$rho' column lengths are not consistent");
      }
      nco2 <- dimnames(loi$rho)[[2]];
      if (!is.null(nco2)) {
        nco2 <- r.bc(ncol(loi$rho));
      }
      if (!all(nco==nco2)) {
        stop("'cova' and 'loi$rho' column names are not consistent");
      }
      if (nrow(cova) != nbs) {
        stop("'cova' row number is not 'nbs'");
      }
    } else {
      if (length(cova) != length(nco)) {
        stop("'cova' and 'loi$rho' lengths are not consistent");
      }
      cova <- matrix(rep(as.numeric(cova),each=nbs),nbs);
    }
    res <- res + cova %*% t(loi$rho);
  }}
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
print8mn <- function(mn,what="msC",ordering=NULL,
                     digits=3,printed=TRUE)
#TITLE standard print function for a /mn/ object.
#DESCRIPTION prints a /mn/ object completely
# or a part of it.
#DETAILS
#KEYWORDS 
#INPUTS
#{mn} <<\code{mn} object to be printed.>>
#[INPUTS]
#{what} <<a \code{character(1)}; when comprising "m" the 
#         expectations are printed, "s" the standard deviations
#         are printed, "C" the correlation matrix is printed,
#         "S" the variance matrix is printed,
#         "P" the precision matrix is printed,
#         "p" the normalized precision matrix is printed.>>
#{ordering} << Nodes are given following the indices of "ordering" if \code{numeric}
#         or the names if it is \code{character}. \code{NULL} means the
#         identity permutation. Repetitions or missing nodes are accepted.>>
#{digits} << when not null, the number of digits for rounding the parameter values.>>
#{printed} << \code{TRUE} to issue a printing, if not the prepared matrix
#           is returned.>>
#VALUE
# The \code{mn} is printed or a matrix having \code{nn x ?} is returned
#  binding which elements precised in the argument \code{what}.
#EXAMPLE
# print8mn(rbmn0mn.01);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_19
#REVISED 12_07_31
#--------------------------------------------
{
  # checking
  # < to be done >
  # number of nodes
  nn <- length(mn$mu);
  nam <- names(mn$mu);
  if (nn!=length(nam)) {
    r.erreur(mn$mu,message="'mn$mu' must be a named vector");
  }
  # getting the ordering for the nodes
  ordering <- aux3(nam,ordering);
  nnr <- length(ordering);
  namr <- nam[ordering];
  # initializing
  cnam <- character(0);
  res <- matrix(NA,nnr,0);
  # printing each asked option
  if (r.expr3present("m",what)) {
    cnam <- c(cnam,"mu");
    res <- cbind(res,mn$mu[ordering]);
  }
  if (r.expr3present("s",what)) {
    cnam <- c(cnam,"s.d.");
    res <- cbind(res,sqrt(diag(mn$gamma[ordering,ordering,drop=FALSE])));
  }
  if (r.expr3present("C",what)) { if (nnr > 1) {
    cnam <- c(cnam,paste("C.",namr,sep=""));
    res <- cbind(res,cor4var(mn$gamma[ordering,ordering,drop=FALSE]));
  }}
  if (r.expr3present("S",what)) {
    cnam <- c(cnam,paste("V.",namr,sep=""));
    res <- cbind(res,mn$gamma[ordering,ordering,drop=FALSE]);
  }
  if (r.expr3present("P",what)) {
    cnam <- c(cnam,paste("P.",namr,sep=""));
    res <- cbind(res,var2pre(mn$gamma[ordering,ordering,drop=FALSE]));
  }
  if (r.expr3present("p",what)) {
    cnam <- c(cnam,paste("PC.",namr,sep=""));
    res <- cbind(res,cor4var(var2pre(mn$gamma[ordering,ordering,drop=FALSE])));
  }
  # dimnaming
  dimnames(res)[[2]] <- cnam;
  # rounding
  if (!is.null(digits)) {
    res <- round(res,digits);
  }
  # returning
  if (printed) {
    print(res);
    invisible();
  } else {
    return(res);
  }
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
dev4mn <- function(Y,EY,VY)
#TITLE Computes the deviance for a sample of multinormal vector
#DESCRIPTION 
# From the \code{n} observed values of a vector of size \code{p} (Y),
# their expectations (EY) and the
# variance matrix (VY) supposed identical for all vectors,
# returns the deviance, i.e. \code{-2*log(p(Y))}.
#DETAILS
# When \code{EY} is a vector with length \code{ncol(Y)} this
# supposes that all observations have the same expectation.
#KEYWORDS 
#INPUTS
#{Y} <<Matrix \code{nxp} of the \code{n} observed values of length \code{p}.>>
#{EY}<<Expectation of \code{Y} (matrix \code{nxp} or vector \code{p}).>>
#{VY}<<Matrix of the variance of each row of \code{Y} (matrix \code{pxp}).>>
#[INPUTS]
#VALUE
# A scalar 
#EXAMPLE
# dev4mn(matrix(runif(3),1),t(rbmn0mn.01$mu),rbmn0mn.01$gamma);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT to be made consistent with a /mn/ object!
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_06_19
#REVISED 13_07_30
#--------------------------------------------
{
  # checking 1
  if (!(is.matrix(Y) | is.data.frame(Y))) {
    stop("'Y' must be a data.frame or a matrix");
  }
  if (length(EY) == ncol(Y)) {
    EY <- matrix(rep(EY,each=nrow(Y)),nrow(Y));
  }
  if (!is.matrix(EY)) { stop("'EY' must be a matrix");}
  if (!is.matrix(VY)) { stop("'VY' must be a matrix");}
  # size of the vector
  pp <- ncol(Y);
  # checking 2
  if (!all(dim(Y)==dim(EY))) {
    stop("'Y' and 'EY' don't have the same dimensions");
  }
  if (!all(dim(VY)==rep(pp,2))) {
    stop("'VY' doesn't have the consistent dimensions");
  }
  # computing the precision matrix
  S1 <- solve(VY);
  # residuals
  del <- Y-EY;
  # computing the quadratic part
  qua <- 0;
  for (ii in r.bc(nrow(Y))) {
    deli <- matrix(as.numeric(del[ii,,drop=FALSE]),nrow=1);
    qua <- qua + deli %*% S1 %*% t(deli);
  }
  # computing the "constant" part
  kon <- determinant(VY)$modulus + pp*log(2*pi);
  # returning
  nrow(Y)*kon + qua;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
