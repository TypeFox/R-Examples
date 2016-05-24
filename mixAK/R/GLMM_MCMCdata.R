##
##  PURPOSE:   Data manipulation for GLMM_MCMC function
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    04/08/2009
##
##  FUNCTIONS:  GLMM_MCMCdata
##
## ================================================================================================

## *************************************************************
## GLMM_MCMCdata
## *************************************************************
##
GLMM_MCMCdata <- function(y, dist, id, time, x, z, random.intercept)
{
##### The resulting object has the following
##### variables:
#####             y, id, time, x, z, random.intercept
#####                               adjusted values of original arguments  
#####  
#####             R                 number of response variables
#####             Rc                number of continuous response variables
#####             Rd                number of discrete response variables  
#####             ndist             numerical counterpart of 'dist' argument
#####             xempty            logical vector of length R
#####             zempty            logical vector of length R
#####             p                 numeric vector of length R  
#####             q                 numeric vector of length R
#####             fixed.intercept   logical vector of length R
#####             CrandomIntcpt     numerical counterpart of random.intercept
#####             CfixedIntcpt      numerical counterpart of fixed.intercept
#####             dimb              dimension of random effects (random intercepts included)
#####             LTb               length of lower triangle of matrix dimb x dimb
#####             naamLTb           names (indices) for elements of a lower triangle of matrix dimb x dimb  
#####             lalpha             length of alpha vector (fixed intercepts included)
#####             p_fi              numeric vector of length R
#####             q_ri              numeric vector of length R  
##### ------------------------------------------------------------------------------------------------------------------
#####
##### Other (potentially useful) variables created here:
#####             neworder          new ordering of rows of response and design matrices to conform to ordering of 'id'
#####
##### ------------------------------------------------------------------------------------------------------------------
    
    
  ##### Response
  ##### --------------------------------------------
  if (missing(y)) stop("y (response) must be given")
  if (is.null(dim(y))) y <- data.frame(y1=y)
  if (is.matrix(y)){
    y <- as.data.frame(y)
  }  
  R <- ncol(y)

  
  ##### Type of response
  ##### --------------------------------------------
  TAB.RESP <- c("gaussian", "binomial(logit)", "poisson(log)")
  if (length(dist) == 1) dist <- rep(dist, R)
  if (length(dist) != R) stop("dist argument has incorrect length")
  ndist <- pmatch(dist, table = c(TAB.RESP, "binomial (logit)", "poisson (log)"), nomatch = 0, duplicates.ok = TRUE)
  ndist[ndist == 4] <- 2  
  ndist[ndist == 5] <- 3
  if (any(!ndist)) stop(paste("Incorrect dist argument, possible dist values are: ", paste(paste("\"", TAB.RESP, "\"", sep=""), collapse=", "), ".", sep=""))
  ndist <- ndist - 1   ## 0 = gaussian etc.

  ### Make dist a full string (added on 20130311).
  ### Check whether binomial/poisson responses are correct (added here on 20140901),
  ### originally provided by the following code which causes error when there were two
  ### or more columns of the same type:
  ###
  ###     if (any(ndist == 1)) if (any(!(y[!is.na(y[, dist=="binomial(logit)"]), dist=="binomial(logit)"] %in% c(0, 1)))) stop("binomial response must be 0/1")
  ###     if (any(ndist == 2)) if (any(!(y[!is.na(y[, dist=="poisson(log)"]), dist=="poisson(log)"] >= 0))) stop("poisson response must be non-negative")
  ### (problem of indexing operator in the data.frame when we give more than one column index)
  ###
  dist <- as.character(rep(NA, length(ndist)))
  for (i in 1:length(ndist)){
    switch (as.character(ndist[i]),
            "0" = {
                dist[i] <- "gaussian"
             },
            "1" = {
                dist[i] <- "binomial(logit)"
                if (any(!(y[!is.na(y[, i]), i] %in% c(0, 1)))) stop("binomial response must be 0/1")                           
             },
            "2" = {
                dist[i] <- "poisson(log)"
                if (any(!(y[!is.na(y[, i]), i] >= 0))) stop("poisson response must be non-negative")                               
             }
            )
  }  
  
  Rc <- sum(ndist %in% c(0))
  Rd <- sum(ndist %in% c(1, 2))
  if (Rc & Rd){   ## check that continuous response comes first
    if (any(ndist[1:Rc] > 0)) stop("continuous response must preceed discrete response in the y matrix")
  }  
  
  
  ##### Cluster and time indicator
  ##### --------------------------------------------
  if (missing(id)) id <- 1:nrow(y)
  if (length(id) != nrow(y)) stop(paste("id must be of length ", nrow(y), sep=""))
  if (missing(time)){
    TAB.ID <- table(id)    
    time <- 1:TAB.ID[1]
    if (length(TAB.ID) > 1) for (i in 2:length(TAB.ID)) time <- c(time, 1:TAB.ID[i])
  }
  if (length(time) != length(id)) stop("time and id variables must be of the same length")  
  neworder <- order(id)  
  id <- id[neworder]
  time <- time[neworder]
  y <- y[neworder,]
  if (R == 1) y <- data.frame(y1=y)
  
  ##### Random and fixed intercept
  ##### -------------------------------------------
  if (missing(random.intercept)) random.intercept <- rep(FALSE, R)
  if (length(random.intercept) == 1) random.intercept <- rep(random.intercept, R)
  if (length(random.intercept) != R) stop(paste("random intercept must be of length ", R, sep=""))
  CrandomIntcpt <- as.numeric(random.intercept)

  fixed.intercept <- !random.intercept
  CfixedIntcpt <- as.numeric(fixed.intercept)  

  
  ##### Design matrices for fixed effects
  ##### --------------------------------------------
  if (missing(x)){
    x <- list(y1="empty")
    if (R >= 2) for (s in 2:R) x[[s]] <- "empty"
    names(x) <- paste(x, 1:R, sep="")
  }  
  if (is.matrix(x)) if (R >= 2) stop("x must be a list when there are two or more response variables")
  else              x <- list(x1=x)
  if (length(x) != R) stop(paste("x must be of length ", R, sep=""))

  xempty <- rep(FALSE, R)
  p <- rep(0, R)
  for (s in 1:R){
    if (!is.matrix(x[[s]]) & !is.data.frame(x[[s]])){
      if (is.character(x[[s]][1])){
        if (x[[s]][1] != "empty") stop("x[[", s, "]] is character not being equal to ", dQuote("empty"), ". The only possible character value is ", dQuote("empty"), ".", sep = "")
        xempty[s] <- TRUE
        p[s] <- 0
        next
      }
      x[[s]] <- data.frame(x1=x[[s]])
    }
    if (is.matrix(x[[s]])){
      x[[s]] <- as.data.frame(x[[s]])
    }
    if (nrow(x[[s]]) != nrow(y)) stop(paste("x[[", s, "]] must have ", nrow(y), " rows", sep=""))
    p[s] <- ncol(x[[s]])
    x[[s]] <- x[[s]][neworder,]
    if (p[s] == 1) x[[s]] <- data.frame(x1=x[[s]])
  }
    
  ##### Design matrices for random effects
  ##### -------------------------------------------
  if (missing(z)){
    z <- list(y1="empty")
    if (R >= 2) for (s in 2:R) z[[s]] <- "empty"
    names(z) <- paste(z, 1:R, sep="")
  }  
  if (is.matrix(z)) if (R >= 2) stop("z must be a list when there are two or more response variables")
  else              z <- list(z1=z)
  if (length(z) != R) stop(paste("z must be of length ", R, sep=""))
  
  zempty <- rep(FALSE, R)
  q <- rep(0, R)
  for (s in 1:R){
    if (!is.matrix(z[[s]]) & !is.data.frame(z[[s]])){
      if (is.character(z[[s]][1])){
        if (z[[s]][1] != "empty") stop("z[[", s, "]] is character not being equal to ", dQuote("empty"), ". The only possible character value is ", dQuote("empty"), ".", sep = "")
        zempty[s] <- TRUE
        q[s] <- 0
        next
      }  
      z[[s]] <- data.frame(z1=z[[s]])
    }
    if (is.matrix(z[[s]])){
      z[[s]] <- as.data.frame(z[[s]])
    }
    if (nrow(z[[s]]) != nrow(y)) stop(paste("z[[", s, "]] must have ", nrow(y), " rows", sep=""))
    q[s] <- ncol(z[[s]])
    z[[s]] <- z[[s]][neworder,]
    if (q[s] == 1) z[[s]] <- data.frame(z1=z[[s]])
  }

  ##### Dimension of random effects and length of alpha  vector
  ##### ------------------------------------------------------
  dimb <- sum(q) + sum(CrandomIntcpt)
  lalpha <- sum(p) + sum(CfixedIntcpt)
  LTb <- (dimb * (dimb + 1))/2
  if (dimb){
    Imat <- diag(dimb)
    rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
    colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE)] 
    naamLTb <- paste(".", rowsI, ".", colsI, sep="")
  }else{
    naamLTb <- "0"
  }    

  p_fi <- p + fixed.intercept
  q_ri <- q + random.intercept  

  RET <- list(y                = y,
              dist             = dist,
              id               = id,
              time             = time,
              x                = x,
              z                = z,
              random.intercept = random.intercept,
              R                = R,
              Rc               = Rc,
              Rd               = Rd,
              ndist            = ndist,
              xempty           = xempty,
              zempty           = zempty,
              p                = p,
              q                = q,
              fixed.intercept  = fixed.intercept,
              CrandomIntcpt    = CrandomIntcpt,
              CfixedIntcpt     = CfixedIntcpt,
              dimb             = dimb,
              LTb              = LTb,
              naamLTb          = naamLTb,
              lalpha           = lalpha,
              p_fi             = p_fi,
              q_ri             = q_ri)
  if (R > 1) RET$name.response <- names(x) else RET$name.response <- "Response"

  return(RET)  
}
