## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### partial derivatives of the CDF wrt arguments ###############################

setGeneric("dCdu", function(cop, u) standardGeneric("dCdu"))

## Warning: This function assumes symmetry in u
dCduExplicitCopula <- function(cop, u)
{
    p <- cop@dimension
    alpha <- cop@parameters
    mat <- matrix(NA_real_, nrow(u),p)
    algNm <- paste(class(cop)[1], "cdfDerWrtArg.algr", sep=".")
    if(exists(algNm)) {
        der.cdf.u <- get(algNm)[p]
        unames0 <- paste0("u",1:p)
        for (j in 1:p)
        {
            unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
            colnames(u) <- unames
            mat[,j] <- eval(der.cdf.u, data.frame(u))
        }
    } else warning("there is no formula for dCdu*() for this copula")
    return(mat)
}

## this function is used for Khoudraji's device
dCduIndepCopula <- function(cop, u) {
  p <- cop@dimension
  mat <- matrix(0, nrow(u), p)
  for (j in 1:p) {
    mat[,j] <- apply(u[, -j, drop=FALSE], 1, prod)
  }
  mat
}

setMethod("dCdu", signature("archmCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("plackettCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("evCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("gumbelCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("indepCopula"), dCduIndepCopula)

dCduEllipCopula <- function(cop, u)
{
    p <- cop@dimension
    sigma <- getSigma(cop)

    ## quantile transformation
    v <- switch(class(cop),
		"normalCopula" = qnorm(u),
		"tCopula" = {
		    df <- cop@df
		    qt(u, df=df)
		},
		stop("not implemented for class ", class(cop)))
    n <- nrow(u)
    mat <- matrix(0,n,p)

    for (j in 1:p)
    {
	s <- sigma[-j,-j] - sigma[-j,j,drop=FALSE] %*% sigma[j,-j,drop=FALSE]

	switch(class(cop),
               "normalCopula" =
           {
               if (p == 2) {
                   rho <- cop@parameters
                   mat[,j] <- pnorm(v[,-j], rho * v[,j], sqrt(1 - rho^2))
               }
               else
                   for (i in 1:n)
                       mat[i,j] <- pmvnorm(lower = rep(-Inf, p - 1), upper = v[i,-j],
                                           mean = v[i,j] * sigma[-j,j],
                                           sigma = drop(s))
           },
               "tCopula" =
           {
               if (p == 2) {
                   rho <- cop@parameters
                   mat[,j] <-  pt(sqrt((df+1)/(df+v[,j]^2)) / sqrt(1 - rho^2)
                                  * (v[,-j] - rho * v[,j]), df=df+1)
               }
               else {
                   if(df != as.integer(df))
                       stop("'df' is not integer; therefore, dCdu() cannot be computed yet")
                   for (i in 1:n)
                       mat[i,j] <- pmvt(lower = rep(-Inf, p - 1),
                                        upper = drop(sqrt((df+1)/(df+v[i,j]^2)) *
                                        (v[i,-j] - v[i,j] * sigma[-j,j])),
                                        sigma = s, df = df + 1)
               }

           })
    }
    mat
}

setMethod("dCdu", signature("ellipCopula"), dCduEllipCopula)


### Plackett formula for elliptical copulas ####################################

setGeneric("plackettFormulaDim2", function(cop, x) standardGeneric("plackettFormulaDim2"))

plackettFormulaDim2NormalCopula <- function(cop, x)
  {
    rho <- cop@parameters
    ir2 <- 1 - rho^2
    as.matrix(exp(-(x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
                  (2 * ir2)) /
              (2 * pi * sqrt(ir2)))
  }

setMethod("plackettFormulaDim2", signature("normalCopula"), plackettFormulaDim2NormalCopula)

plackettFormulaDim2TCopula <- function(cop, x)
  {
    rho <- cop@parameters
    ir2 <- 1 - rho^2
    df <- cop@df
    as.matrix((1 + (x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
               (df * ir2))^(-df / 2) / (2 * pi * sqrt(ir2)))
  }

setMethod("plackettFormulaDim2", signature("tCopula"), plackettFormulaDim2TCopula)

setGeneric("plackettFormula",  function(cop, p, rho, s, m, x, i, j) standardGeneric("plackettFormula"))

plackettFormulaNormalCopula <- function(cop, p, rho, s, m, x, i, j)
{
    exp(-(x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) /
               (2 * (1 - rho^2))) / (2 * pi * sqrt(1 - rho^2)) *
           (if (p == 3) pnorm(drop((x[-c(i,j)] - m %*% x[c(i,j)])/sqrt(s)))
           else pmvnorm(lower = rep(-Inf, p - 2),
                        upper = drop(x[-c(i,j)] - m %*% x[c(i,j)]),
                        sigma = s))
}

setMethod("plackettFormula", signature("normalCopula"), plackettFormulaNormalCopula)

plackettFormulaTCopula <- function(cop, p, rho, s, m, x, i, j)
{
    stopifnot(p >= 3)
    df <- cop@df
    if(df != as.integer(df) && p > 3)
	stop("'df' is not integer; therefore, plackettFormula() cannot be computed yet")
    term <- 1 + (x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) / (df * (1 - rho^2))
    ## return:
    term^(-df / 2) / (2 * pi * sqrt(1 - rho^2)) *
	if (p == 3) pt(drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term * s)), df= df)
	else pmvt(df = df, lower = rep(-Inf, p - 2),
		  upper = drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term)),
		  sigma = s)
}

setMethod("plackettFormula", signature("tCopula"), plackettFormulaTCopula)


### partial derivatives of the CDF wrt parameters ##############################

setGeneric("dCdtheta", function(cop, u) standardGeneric("dCdtheta"))


dCdthetaExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    algNm <- paste(class(cop)[1], "cdfDerWrtPar.algr", sep=".")
    if(exists(algNm)) {
        der.cdf.alpha <- get(algNm)[p]
        colnames(u) <- paste0("u", 1:p)
        as.matrix(eval(der.cdf.alpha, data.frame(u)))
    } else {
        warning("there is no formula for dCdtheta*() for this copula")
        matrix(NA_real_, nrow(u),p)
    }
  }

dCdthetaEvCopula <- function(cop, u) {
  alpha <- cop@parameters
  loguv <- log(u[,1], u[,2])
  w <- log(u[,2]) / loguv
  der.cdf.alpha <- pCopula(u, cop) * loguv * dAdtheta(cop, w)
  return(as.matrix(der.cdf.alpha))
}

setMethod("dCdtheta", signature("archmCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("plackettCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("evCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("gumbelCopula"), dCdthetaExplicitCopula)

dCdthetaEllipCopula <- function(cop, u)
  {
    p <- cop@dimension

    ## quantile transformation
    v <- switch(class(cop),
		"normalCopula" = qnorm(u),
		"tCopula" = qt(u, df = cop@df),
		stop("not implemented for class ", class(cop)))

    if (p == 2)
      plackettFormulaDim2(cop, v)
    else
      {
        n <- nrow(u)
        sigma <- getSigma(cop)

        if (cop@dispstr %in% c("ex","ar1")) ## exchangeable or ar1
          {
            rho <- cop@parameters
            r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
            m <- sigma[-c(1,2),c(1,2)] %*% r
            s <- sigma[-c(1,2),-c(1,2)] - sigma[-c(1,2),c(1,2)] %*% r %*% sigma[c(1,2),-c(1,2)]

            mat <- matrix(0,n,1)

            if (cop@dispstr == "ex") ## exchangeable
              for (k in 1:n)
                for (j in 1:(p-1))
                  for (i in (j+1):p)
                    mat[k,1] <- mat[k,1] + plackettFormula(cop, p, rho, s, m, v[k,], i, j)
            else ## ar1
              for (k in 1:n)
                for (j in 1:(p-1))
                  for (i in (j+1):p)
                    mat[k,1] <- mat[k,1] + (i - j) * rho ^ (i - j - 1) *
                      plackettFormula(cop, p, rho, s, m, v[k,], i, j)

            mat
          }
	else # unstructured or toeplitz or ...
          {
            mat <- matrix(0,n,p*(p-1)/2)
            l <- 1
            for (j in 1:(p-1))
              for (i in (j+1):p)
                {
                  rho <- sigma[i,j]
                  r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
                  m <- sigma[-c(i,j),c(i,j)] %*% r
                  s <- sigma[-c(i,j),-c(i,j)] - sigma[-c(i,j),c(i,j)] %*% r %*% sigma[c(i,j),-c(i,j)]

                  for (k in 1:n)
                    mat[k,l] <- plackettFormula(cop, p, rho, s, m, v[k,], i, j)
                  l <- l + 1
                }
            if (cop@dispstr == "un") ## unstructured
                mat
            else if (cop@dispstr == "toep") {
                coef <- matrix(0,p*(p-1)/2,p-1)
                for (k in 1:(p-1))
                  {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- P2p(m)
                  }
                mat %*% coef
            }
            else stop("Not implemented yet for the dispersion structure ", cop@dispstr)
          }
      }
  }

setMethod("dCdtheta", signature("ellipCopula"), dCdthetaEllipCopula)


### Partial derivatives of the PDF wrt arguments for ellipCopula: DIVIDED BY PDF

setGeneric("derPdfWrtArgs", function(cop, u) standardGeneric("derPdfWrtArgs"))

derPdfWrtArgsExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    algNm <- paste(class(cop)[1], "pdfDerWrtArg.algr", sep=".")
    mat <- matrix(NA_real_, nrow(u),p)
    if(exists(algNm)) {
        der.pdf.u <- get(algNm)[p]
        unames0 <- paste0("u", 1:p)
        for (j in 1:p)
        {
            unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
            colnames(u) <- unames
            mat[,j] <- eval(der.pdf.u, data.frame(u))
        }
    } else warning("there is no formula for derPdfWrtArgs*() for this copula")
    return(mat)
  }

setMethod("derPdfWrtArgs", signature("archmCopula"), derPdfWrtArgsExplicitCopula)
setMethod("derPdfWrtArgs", signature("plackettCopula"), derPdfWrtArgsExplicitCopula)
setMethod("derPdfWrtArgs", signature("evCopula"), derPdfWrtArgsExplicitCopula)
setMethod("derPdfWrtArgs", signature("gumbelCopula"), derPdfWrtArgsExplicitCopula)

derPdfWrtArgsNormalCopula <- function(cop, u)
  {
    v <- qnorm(u)
    return( ( - v %*% solve(getSigma(cop)) + v)/ dnorm(v) )
  }

setMethod("derPdfWrtArgs", signature("normalCopula"), derPdfWrtArgsNormalCopula)

derPdfWrtArgsTCopula <- function(cop, u)
  {

    df <- cop@df
    v <- qt(u,df=df)
    w <- dt(v,df=df)
    m <- v %*% solve(getSigma(cop))

    return( - (df + cop@dimension) * m / ((df + rowSums(m * v)) * w) +
           (df + 1) * v / ((df +  v^2) * w) )
  }

setMethod("derPdfWrtArgs", signature("tCopula"), derPdfWrtArgsTCopula)


## Partial derivatives of the PDF wrt parameters for ellipCopula: DIVIDED BY PDF

setGeneric("derPdfWrtParams", function(cop, u) standardGeneric("derPdfWrtParams"))

derPdfWrtParamsExplicitCopula <- function(cop, u)
{
    p <- cop@dimension
    alpha <- cop@parameters
    algNm <- paste(class(cop)[1], "pdfDerWrtPar.algr", sep=".")
    if(exists(algNm) && !is.null((der.pdf.alpha <- get(algNm)[p])[[1]])) {
        colnames(u) <- paste0("u", 1:p)
        as.matrix(eval(der.pdf.alpha, data.frame(u)))
    } else {
        warning("there is no formula for derPdfWrtParam*() for this copula")
        matrix(NA_real_, nrow(u),p)
    }
}

setMethod("derPdfWrtParams", signature("archmCopula"), derPdfWrtParamsExplicitCopula)
setMethod("derPdfWrtParams", signature("plackettCopula"), derPdfWrtParamsExplicitCopula)
setMethod("derPdfWrtParams", signature("evCopula"), derPdfWrtParamsExplicitCopula)

derPdfWrtParamsEllipCopula <- function(cop, u)
{
    p <- cop@dimension

    ## quantile transformation
    v <- switch(clc <- class(cop),
		"normalCopula" = qnorm(u),
		"tCopula" = {
		    df <- cop@df
		    qt(u, df=df)
		},
		## else:
		stop("not implemented"))
    if (p == 2)
    {
	rho <- cop@parameters
	ir2 <- 1 - rho^2
        sv2 <- rowSums(v^2) # == v[,1]^2 + v[,2]^2
	as.matrix(switch(clc,
			 "normalCopula" =
			 (rho * ir2 - rho * sv2 + (rho^2 + 1) * v[,1] * v[,2])/ir2^2,
			 "tCopula" =
			 (1 + df) * rho / -ir2 + (2 + df) * (df * rho + v[,1] * v[,2])
			 / (df * ir2 + sv2 - 2 * rho * v[,1] * v[,2])))

    } else { ##  p >= 3
        n <- nrow(u)
        sigma <- getSigma(cop)
        detsig <- det(sigma)
        invsig <- solve(sigma)

        if (cop@dispstr %in% c("ex","ar1")) ## exchangeable or ar1
        {
            rho <- cop@parameters

            dersig <- matrix(1,p,p)
            if (cop@dispstr == "ex") ## ex
                diag(dersig) <- 0
            else ## ar1
                for (i in 1:p)
                    for (j in 1:p) {
                        ij <- abs(i - j)
                        dersig[i,j] <- ij * rho^(ij - 1)
                    }

            ## MM:  sum(diag(A %*% B)) == sum(A * t(B)) .. and B=dersig is symmetric here
            ## derdetsig <- detsig * sum(diag(invsig %*% dersig))
            derdetsig <- detsig *  sum(invsig * dersig)
            derinvsig <- - invsig %*% dersig %*% invsig
            firstterm <- derdetsig/detsig

	    mat <-
		switch(clc,
		       "normalCopula" =
		       - (firstterm + rowSums((v %*% derinvsig) * v))/2,
		       "tCopula" =
		       - (firstterm + (df + p) * rowSums((v %*% derinvsig) * v)
			  / (df +  rowSums((v %*% invsig) * v)) ) / 2)
            as.matrix(mat)
        }
	else # unstructured or toeplitz or ...
	  {
            mat <- matrix(0,n,p*(p-1)/2)
            l <- 1
            for (j in 1:(p-1))
                for (i in (j+1):p)
                {
                    derdetsig <- 2 * det(sigma[-i,-j,drop=FALSE]) * (-1)^(i+j)
                    derinvsig <- - invsig[,i] %*% t(invsig[,j]) - invsig[,j] %*% t(invsig[,i])
                    firstterm <- derdetsig/detsig

		    mat[,l] <-
			switch(clc,
			       "normalCopula" =
			       - (firstterm + rowSums((v %*% derinvsig) * v))/2,
			       "tCopula" =
			       - (firstterm + (df + p) * rowSums((v %*% derinvsig) * v)
				  / (df +  rowSums((v %*% invsig) * v)) ) / 2)
                    l <- l + 1
                }
            if (cop@dispstr == "un") ## unstructured
                mat
            else if (cop@dispstr == "toep") { ## toeplitz: p-1 parameters
                coef <- matrix(0,p*(p-1)/2,p-1)
                for (k in 1:(p-1))
                {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- P2p(m)
                }
                mat %*% coef
            }
            else stop("Not implemented yet for the dispersion structure ", cop@dispstr)
        }
    } ## p >= 3
}

setMethod("derPdfWrtParams", signature("ellipCopula"), derPdfWrtParamsEllipCopula)


### dCopula wrapper for influence coefficients #################################

setGeneric("dcopwrap",  function(cop, u, ...) standardGeneric("dcopwrap"))

dcopwrapExplicitCopula <- function(cop, u) dCopula(u,cop)

setMethod("dcopwrap", signature("archmCopula"),	   dcopwrapExplicitCopula)
setMethod("dcopwrap", signature("evCopula"),	   dcopwrapExplicitCopula)
setMethod("dcopwrap", signature("plackettCopula"), dcopwrapExplicitCopula)

dcopwrapEllipCopula <- function(cop, u) rep.int(1, NROW(u))

setMethod("dcopwrap", signature("ellipCopula"), dcopwrapEllipCopula)
