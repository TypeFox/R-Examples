#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

## cop is a copula of the desired family whose parameters will be used
## as starting values in fitCopula
gof <- function(cop, x, N=1000, method="kendall")
  {
    x <- as.matrix(x)

    ## number of observations
    n <- dim(x)[1]
    ## number of variables
    p <- dim(x)[2]

    if (p < 2)
      stop("The data should be at least of dimension 2")

    if (n < 2)
      stop("There should be at least 2 observations")

    if (cop@dimension != p)
      stop("The copula and the data should be of the same dimension")

    METHODS <- c("likelihood", "kendall", "spearman")
    method <-  pmatch(method, METHODS)
    if(is.na(method))
      stop("Invalid estimation method")
    if(method == -1)
      stop("Ambiguous estimation method")

    if ((method == 2 || method == 3) && p != 2)
      stop("Estimation based on Kendall's tau or Spearman's rho requires p = 2")

    if ((method == 2 || method == 3) && length(cop@parameters) != 1)
      stop("Estimation based on Kendall's tau or Spearman's rho will work only for one-parameter copulas")


    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    if (method==1)
      fcop <- fitCopula(u,cop,cop@parameters)@copula
    else if (method==2)
      fcop <- fitCopulaKendallsTau(cop,u)
    else if (method==3)
      fcop <- fitCopulaSpearmansRho(cop,u)

    ## compute the test statistic
    s <- .C("cramer_vonMises",
            as.integer(n),
            as.integer(p),
            as.double(u),
            as.double(pCopula(u, fcop)),
            stat = double(1),
            PACKAGE="copula")$stat

    ## simulation of the null distribution
    s0 <- numeric(N)
    for (i in 1:N)
      {
        cat(paste("Iteration",i,"\n"))
        x0 <- rCopula(n, fcop)
        u0 <- apply(x0,2,rank)/(n+1)

         ## fit the copula
        if (method==1)
          fcop0 <- fitCopula(u0,cop,fcop@parameters)@copula
        else if (method==2)
          fcop0 <- fitCopulaKendallsTau(cop,u0)
        else if (method==3)
          fcop0 <- fitCopulaSpearmansRho(cop,u0)

        s0[i] <- .C("cramer_vonMises",
                    as.integer(n),
                    as.integer(p),
                    as.double(u0),
                    as.double(pCopula(u0, fcop0)),
                    stat = double(1),
                    PACKAGE="copula")$stat
      }

    return(list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1)))
  }

##############################################################################

gofMobius <- function(cop, x, maxcard=ncol(x), use.empcop=TRUE, m=2000, N=100)
  {
    x <- as.matrix(x)

    ## number of observations
    n <- dim(x)[1]
    ## number of variables
    p <- dim(x)[2]

    if (p < 2)
      stop("The data should be at least of dimension 2")

    if (n < 2)
      stop("There should be at least 2 observations")

    if (cop@dimension != p)
      stop("The copula and the data should be of the same dimension")

    if (!is.numeric(maxcard) || (maxcard <- as.integer(maxcard)) < 2
        || maxcard > p)
      stop(paste("maxcard should be an integer greater than 2 and smaller than",p))

    ## number of subsets
    sb <- binom.sum(p,maxcard)
    nsubsets <- sb - p - 1

    ## power set of {1,...,dim} in integer notation
    subsets <-  .C("k_power_set",
                   as.integer(p),
                   maxcard,
                   subsets = integer(sb),
                   PACKAGE="copula")$subsets

    ## power set in character vector: {}, {1}, {2}, ..., {1,2}, ..., {1,...,dim}
    subsets.char <-  .C("k_power_set_char",
                        as.integer(p),
                        as.integer(sb),
                        as.integer(subsets),
                        sc = character(sb),
                        PACKAGE="copula")$sc

    ## remove emptyset and singletons
    subsets <- subsets[(p+2):sb]
    subsets.char <- subsets.char[(p+2):sb]

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    fcop <- fitCopula(u,cop,cop@parameters)@copula

    ## compute the test statistic
    if (use.empcop == TRUE)
      s <- .C("Mobius_cramer_vonMises_approx",
              as.integer(n),
              as.integer(p),
              as.integer(nsubsets),
              as.integer(subsets),
              as.double(u),
              as.double(rCopula(m, fcop)),
              as.integer(m),
              stat = double(nsubsets+1),
              PACKAGE="copula")$stat
    else
      {
        marg <- vector("list",nsubsets)
        pcop <- matrix(0,n,nsubsets)
        for (j in 1:nsubsets)
          {
            marg[[j]] <- sub('\\{', 'c(', subsets.char[j])
            marg[[j]] <- eval(parse(text=sub('\\}', ')', marg[[j]])))
            ones <- matrix(1,n,p)
            ones[,marg[[j]]] <- u[,marg[[j]]]
            pcop[,j] <- pCopula(ones, fcop)
          }

        s <- .C("Mobius_cramer_vonMises",
                as.integer(n),
                as.integer(p),
                as.integer(nsubsets),
                as.integer(subsets),
                as.double(u),
                as.double(pcop),
                stat = double(nsubsets+1),
                PACKAGE="copula")$stat
      }

    ## simulation of the null distribution
    s0 <- matrix(0,N,nsubsets+1)
    for (i in 1:N)
      {
        cat(paste("Iteration",i,"\n"))
        x0 <- rCopula(n, fcop)
        u0 <- apply(x0,2,rank)/(n+1)

        ## fit the copula
        fcop0 <- fitCopula(u0,cop,cop@parameters)@copula

        ## compute the test statistic
        if (use.empcop == TRUE)
          s0[i,] <- .C("Mobius_cramer_vonMises_approx",
                       as.integer(n),
                       as.integer(p),
                       as.integer(nsubsets),
                       as.integer(subsets),
                       as.double(u0),
                       as.double(rCopula(m, fcop0)),
                       as.integer(m),
                       stat = double(nsubsets+1),
                       PACKAGE="copula")$stat
        else
          {
            pcop0 <- matrix(0,n,nsubsets)
            for (j in 1:nsubsets)
              {
                ones <- matrix(1,n,p)
                ones[,marg[[j]]] <- u0[,marg[[j]]]
                pcop0[,j] <- pCopula(ones, fcop0)
              }

            s0[i,] <- .C("Mobius_cramer_vonMises",
                        as.integer(n),
                        as.integer(p),
                        as.integer(nsubsets),
                        as.integer(subsets),
                        as.double(u0),
                        as.double(pcop0),
                        stat = double(nsubsets+1),
                        PACKAGE="copula")$stat
        }
      }

    pval <-.C("Mobius_pvalues",
              as.integer(nsubsets),
              as.integer(N),
              as.double(rbind(s0,s)),
              pval = double(nsubsets+3),
              PACKAGE="copula")$pval

    test <- list(s0=s0, subsets=subsets.char, statistics=s[1:nsubsets],
                 pvalues = pval[1:nsubsets], global.statistic = s[nsubsets+1],
                 global.statistic.pvalue = pval[nsubsets+1], fisher.pvalue=pval[nsubsets+2],
                 tippett.pvalue=pval[nsubsets+3])

    #class(test) <- "ec.gof.test"
    return(test)
  }

##############################################################################

derivatives <- function(cop)
  {
    p <- cop@dimension

    ## partial derivatives of the cdf wrt arguments
    cdf.u <- vector("list",p)
    for (j in 1:p)
      cdf.u[[j]] <- D(cop@exprdist$cdf,paste("u",j,sep=""))

    ## partial derivative of the cdf wrt parameter
    cdf.alpha <- D(cop@exprdist$cdf,"alpha")

    ## partial derivative of the pdf wrt parameter
    pdf.alpha <- D(cop@exprdist$pdf,"alpha")

    ## partial derivatives of the (pdf wrt parameter) wrt arguments
    pdf.u <- vector("list",p)
    for (j in 1:p)
      pdf.u[[j]] <- D(cop@exprdist$pdf,paste("u",j,sep=""))

    return(list(cdf.u=cdf.u,cdf.alpha=cdf.alpha,
                pdf.alpha=pdf.alpha,pdf.u=pdf.u))
  }

##############################################################################

fitCopulaKendallsTau <- function(cop,u)
{
    cop@parameters <- iTau(cop,cor(u[,1],u[,2],method="kendall"))
    return(cop)
}

##############################################################################

fitCopulaSpearmansRho <- function(cop,u)
{
    cop@parameters <- iRho(cop,cor(u[,1],u[,2],method="spearman"))
    return(cop)
}

##############################################################################

## additional influence terms

influ.add <- function(x0, y0, influ1, influ2)
{
  M <- nrow(y0)
  o1 <- order(y0[,1], decreasing=TRUE)
  o1b <- ecdf(y0[,1])(x0[,1]) * M
  o2 <- order(y0[,2], decreasing=TRUE)
  o2b <- ecdf(y0[,2])(x0[,2]) * M
  return(c(0,cumsum(influ1[o1]))[M + 1 - o1b] / M - mean(influ1 * y0[,1]) +
         c(0,cumsum(influ2[o2]))[M + 1 - o2b] / M - mean(influ2 * y0[,2]))
}

##############################################################################

## grid == 0 => from the H0 copula but indep of x0 (what we used until now)
## grid == 1 => from the H0 copula but with x0 = g
## grid == 2 => unif on [0,1]^2
## grid == 3 => from the observed data
## grid == 4 => CVM from observed, H0 dist from generated

gofMultCLT <- function(cop,x, N=1000, method= c("kendall", "likelihood", "spearman"),
                       m=nrow(x), der=NULL, betavar=FALSE, center=FALSE, M=2500, grid=0)
{
    ## data
    x <- as.matrix(x)

    ## number of observations
    n <- dim(x)[1]
    ## number of variables
    p <- dim(x)[2]

    if (p < 2)
      stop("The data should be at least of dimension 2")

    if (n < 2)
      stop("There should be at least 2 observations")

    if (cop@dimension != p)
      stop("The copula and the data should be of the same dimension")

    method <- match.arg(method)

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    cop <- switch(method,
                  "likelihood" =
                  fitCopula(u,cop, iTau(cop, cor(u[,1],u[,2], method="kendall")))@copula,
                  "kendall" = fitCopulaKendallsTau(cop,u),
                  "spearman" = fitCopulaSpearmansRho(cop,u),
                  stop("impossible 'method' -- should never happen please report"))

    ## estimate of theta
    alpha <-  cop@parameters

    ## grid points where to evaluate the process
    if (grid == 0 || grid == 1)
      g <- rCopula(n, cop) ## default
    else if (grid==2)
      g <- matrix(runif(n*p),n,p)
    else if (grid >= 3)
      g <- u
    pcop <- pCopula(g, cop)

    ## compute the test statistic
    s <- .C("cramer_vonMises_2",
            as.integer(p),
            as.double(u),
            as.integer(n),
            as.double(g),
            as.integer(n),
            as.double(pcop),
            stat = double(1),
            PACKAGE="copula")$stat

    ## generate realizations under H0
    if (grid==1)
      x0 <- g
    else
      x0 <- rCopula(m, cop) ## default

    if (grid==4) ## a la Genest et al.
      {
        g <- apply(x0,2,rank)/(m+1) #pseudo-obs
        pcop <- pCopula(g, cop)
      }

    ## matrix of partial derivatives at g
    ## n lines  by p+1 columns
    ## (wrt arguments (p cols) + wrt parameter (1 col) + NEW
    dermat <- matrix(0,n,p+1)

    ## Now, prepare derivatives and influence coefficients

    #####################################################
    ##
    ## Copulas having an explicit cdf
    ##
    #####################################################

    if (!(class(cop) %in% c("normalCopula","tCopula")))
      {
        ## partial derivatives
        if (is.null(der))
          der <- derivatives(cop)

        v <- character()
        for (j in 1:p)
          v <- c(v,paste("u",j,sep=""))

        for (j in 1:p)
          assign(v[j], g[,j])
        for (j in 1:p)
          dermat[,j] <- eval(der$cdf.u[[j]], list(u1 = u1, u2 = u2))
        dermat[,p+1] <- eval(der$cdf.alpha, list(u1 = u1, u2 = u2))

        ## influence coefficients for estimation
        if (method == "likelihood")
          {
            for (j in 1:p)
              assign(v[j], x0[,j])

            ## influence: first part
            influ <- eval(der$pdf.alpha, list(u1 = u1, u2 = u2)) / dCopula(x0, cop)

            ## influence: second part
            ## integrals computed from M realizations by Monte Carlo
            y0 <- rCopula(M, cop)
            dcopy0 <- dCopula(y0, cop)
            for (j in 1:p)
              assign(v[j], y0[,j])
            influ0 <- eval(der$pdf.alpha, list(u1 = u1, u2 = u2)) / dcopy0
            beta <- if (betavar) var(influ0) else mean(influ0^2)
            influ1 <- influ0 * eval(der$pdf.u[[1]], list(u1 = u1, u2 = u2)) / dcopy0
            influ2 <- influ0 * eval(der$pdf.u[[2]], list(u1 = u1, u2 = u2)) / dcopy0

            ## influence: final
            influ <- (influ - influ.add(x0, y0, influ1, influ2))/beta
          }
      }

    #####################################################
    ##
    ## Metaelliptical copulas
    ##
    #####################################################

    else if (class(cop) == "normalCopula") #WARNING p=2 only !!!
      {
        ## derivatives
        v1 <- qnorm(g[,1])
        v2 <- qnorm(g[,2])
        dermat[,1] <- pnorm(v2, alpha * v1, sqrt(1 - alpha^2))
        dermat[,2] <- pnorm(v1, alpha * v2, sqrt(1 - alpha^2))
        dermat[,3] <- exp(-(v1^2 + v2^2 - 2 * alpha * v1 * v2)
                          /(2 * (1 - alpha^2)))/(2 * pi * sqrt(1 - alpha^2))

        ## influence coefficients for estimation
        if (method == "likelihood")
          {
            v1 <- qnorm(x0[,1])
            v2 <- qnorm(x0[,2])

            ## influence: first part
            influ <- (alpha * (1 - alpha^2) - alpha * (v1^2 + v2^2) +
                      (alpha^2 + 1) * v1 * v2)/(1 - alpha^2)^2

            ## influence: second part
            ## integrals computed from M realizations by Monte Carlo
            y0 <- rCopula(M, cop)
            v1 <- qnorm(y0[,1])
            v2 <- qnorm(y0[,2])
            influ0 <- (alpha * (1 - alpha^2) - alpha * (v1^2 + v2^2) +
                       (alpha^2 + 1) * v1 * v2)/(1 - alpha^2)^2
            beta <- if (betavar) var(influ0) else mean(influ0^2)
            influ1 <- influ0 * sqrt(2 * pi) * alpha * (alpha * v1 - v2) *
              exp(v1^2 / 2) / (alpha^2 - 1)
            influ2 <- influ0 * sqrt(2 * pi) * alpha * (alpha * v2 - v1) *
              exp(v2^2 / 2) / (alpha^2 - 1)

            ## influence: final
            influ <- (influ - influ.add(x0, y0, influ1, influ2))/beta
          }
      }
    else if (class(cop) == "tCopula") #WARNING p=2 only !!!
      {
        df <- cop@df

        ## derivatives
        v1 <- qt(g[,1], df=df)
        v2 <- qt(g[,2], df=df)
        dermat[,1] <- pt((v2 - alpha * v1) *
                         sqrt((df+1)/(df+v1^2)) / sqrt(1 - alpha^2), df=df+1)
        dermat[,2] <- pt((v1 - alpha * v2) *
                         sqrt((df+1)/(df+v2^2)) / sqrt(1 - alpha^2), df=df+1)
        dermat[,3] <- (1 + (v1^2 + v2^2 - 2 * alpha * v1 * v2) /
                       (df * (1 - alpha^2)))^(-df / 2) / (2 * pi * sqrt(1 - alpha^2))

        ## influence coefficients for estimation
        if (method== "likelihood")
          {
            v1 <- qt(x0[,1], df=df)
            v2 <- qt(x0[,2], df=df)

            ## influence: first part
            influ <- (1 + df) * alpha / (alpha^2 - 1) +
              (2 + df) * (df * alpha + v1 * v2) /
                (df * (1 - alpha^2) + v1^2 + v2^2 - 2 * alpha * v1 * v2)

            ## influence: second part
            ## integrals computed from M realizations by Monte Carlo
            y0 <- rCopula(M, cop)
            v1 <- qt(y0[,1], df=df)
            v2 <- qt(y0[,2], df=df)

            influ0 <-  (1 + df) * alpha / (alpha^2 - 1) +
              (2 + df) * (df * alpha + v1 * v2) /
                (df * (1 - alpha^2) + v1^2 + v2^2 - 2 * alpha * v1 * v2)
            beta <- if (betavar) var(influ0) else mean(influ0^2)

            influ1 <- - influ0 * sqrt(pi) * ((df+v1^2)/df)^((df-1)/2) *
              (df * (1 + (1 + df) * alpha^2) * v1 + v1^3 - df * alpha *
               (2 + df - v1^2) * v2 - (1 + df) * v1 * v2^2) * gamma(df/2) /
                 (sqrt(df) * (df * (1 - alpha^2) + v1^2 + v2^2 - 2 * alpha * v1 * v2) *
                  gamma((1+df)/2))
            influ2 <- - influ0 * sqrt(pi) * ((df+v2^2)/df)^((df-1)/2) *
              (df * (1 + (1 + df) * alpha^2) * v2 + v2^3 - df * alpha *
               (2 + df - v2^2) * v1 - (1 + df) * v2 * v1^2) * gamma(df/2) /
                 (sqrt(df) * (df * (1 - alpha^2) + v1^2 + v2^2 - 2 * alpha * v1 * v2) *
                  gamma((1+df)/2))

            ## influence: final
            influ <- (influ - influ.add(x0, y0, influ1, influ2))/beta
          }
      }
    else stop("H0 copula ",class(cop)," not implemented")

    #####################################################
    ##
    ## Influence coefficients for estimation by Kendall's tau
    ##
    #####################################################

    if (method == "kendall")
      {
        tauCtheta <- tau(cop)
        influ <- 4 * (2 * pCopula(x0, cop) - x0[,1] - x0[,2] + (1 - tauCtheta)/2)

        influ <-
            switch(class(cop),
                   "claytonCopula" = influ * (alpha+2)^2 / 2,
                   "frankCopula" = influ /
                   (4/alpha^2 + 4/(alpha * (exp(alpha) - 1)) - 8/alpha^2 * debye1(alpha)),
                   "gumbelCopula" =  influ * alpha^2,
                   "plackettCopula" =  ## added by JY; testing only
                   influ / dTauPlackettCopula(cop),
                   "normalCopula" = ,
                   "tCopula" = influ * pi * sqrt(1 - alpha^2) / 2,
                   ## otherwise:
                   stop("H0 copula (",class(cop),") not implemented"))
      }

    #####################################################
    ##
    ## Influence coefficients for estimation by Spearman's rho
    ##
    #####################################################

    if (method == "spearman")
      {
        rhoCtheta <- rho(cop)

        ## integrals computed from M realizations by Monte Carlo
        y0 <- rCopula(M, cop)
        influ <- 12 * (x0[,1] * x0[,2] + influ.add(x0, y0, y0[,2],y0[,1])) -
          3 - rhoCtheta

        influ <-
            switch(class(cop),
                   "claytonCopula" = ## added by JY; testing only
                   influ / dRhoClaytonCopula(cop),
                   "gumbelCopula" = ## added by JY; testing only
                   influ / dRhoGumbelCopula(cop),
                   "frankCopula" =
                   influ / (12 / (alpha * (exp(alpha) - 1)) -
                            36 / alpha^2 * debye2(alpha) +
                            24 / alpha^2 * debye1(alpha)),
                   "plackettCopula" =
                   influ * (alpha - 1)^3 / (2 * (2 - 2 * alpha
                                                 + (1 + alpha) * log(alpha))),
                   "normalCopula" = ,
                   "tCopula" = influ * pi * sqrt(4 - alpha^2) / 6,
                   ## otherwise:
                   stop("H0 copula (",class(cop),") not implemented"))
      }

    #####################################################
    ##
    ## Simulate under H0
    ##
    #####################################################

    if (center)
      influ <- scale(influ,scale=FALSE)

    s0 <- .C("simulate",
             as.integer(p),
             as.double(x0),
             as.integer(m),
             as.double(g),
             as.integer(n),
             as.double(pcop),
             as.double(dermat),
             as.double(influ),
             as.integer(N),
             s0 = double(N),
             PACKAGE="copula")$s0

    list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1), alpha=alpha)
}
