#######################################################################

plackettFormulaDim2 <- function(rho, df, x)
  {
    as.matrix((1 + (x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
               (df * (1 - rho^2)))^(-df / 2) / (2 * pi * sqrt(1 - rho^2)))
  }

plackettFormula <- function(dim, df, rho, s, m, x, i, j)
  {
    term <- 1 + (x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) / (df * (1 - rho^2))
    term^(-df / 2) / (2 * pi * sqrt(1 - rho^2)) *
      if (dim == 3)
        pt(drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term * s)), df= df)
      else
        pmvt(df = df, upper = drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term)),
             sigma = s)
  }

##################################################################################

gradCdfT <- function(param, x, df)
{
  dim <- NCOL(x)
  n <- NROW(x)
  mat <- matrix(NA,n,dim)
  mat2 <- matrix(NA,n,dim*(dim-1)/2)

  myargs <- param2args(param, "un", dim)
  sigma <- myargs$corr      
  
  ## standardize
  invdisp <-  diag(1/sqrt(param[(dim+1):(2*dim)]))
  x <- (x - matrix(param[1:dim],n,dim,byrow=TRUE)) %*% invdisp
  
  ## common to derivatives wrt means and with dispersions
  for (j in 1:dim)
    if (dim == 2)
      {
        rho <- param[5]
        mat[,j] <-  dt(x[,j], df=df) * pt(sqrt((df+1)/(df+x[,j]^2)) / sqrt(1 - rho^2)
                                * (x[,-j] - rho * x[,j]), df=df+1)
      }
    else
      {
        
        s <- sigma[-j,-j] - sigma[-j,j,drop=FALSE] %*% sigma[j,-j,drop=FALSE]
        for (i in 1:n)
          mat[i,j] <- dt(x[i,j], df=df) * pmvt(upper = drop(sqrt((df+1)/(df+x[i,j]^2)) *
                                                 (x[i,-j] - x[i,j] * sigma[-j,j])),
                                   sigma = s, df = df + 1)
      }
  
  ## derivatives wrt correlations
  if (dim == 2)
    {
      rho <- param[5]
      mat2 <- plackettFormulaDim2(rho, df, x)
    }
  else
    {
      l <- 1
      for (j in 1:(dim-1))
        for (i in (j+1):dim)
          {
            rho <- sigma[i,j]
            r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
            m <- sigma[-c(i,j),c(i,j)] %*% r
            s <- sigma[-c(i,j),-c(i,j)] - sigma[-c(i,j),c(i,j)] %*% r %*% sigma[c(i,j),-c(i,j)]
            
            for (k in 1:n) 
              mat2[k,l] <- plackettFormula(dim, df, rho, s, m, x[k,], i, j)
            l <- l + 1
          } 
    }
  cbind(- mat %*% invdisp, - (x %*% invdisp^2) * mat/2, mat2)
}

##################################################################################

gradLogPdfT <- function(param, x, df)
{
  dim <- NCOL(x)
  n <- NROW(x)
  mat2 <- matrix(NA,n,dim*(dim-1)/2)

  myargs <- param2args(param, "un", dim)
  sigma <- myargs$corr      
  
  ## standardize
  disp <- sqrt(param[(dim+1):(2*dim)])
  invdisp <-  diag(1/disp)
  x <- (x - matrix(param[1:dim],n,dim,byrow=TRUE)) %*% invdisp
  
  ## common to derivatives wrt means and with dispersions
  invsig <- solve(sigma)
  m <- x %*% invsig
  mat <- - (df + dim) * m / (df + rowSums(m * x))


  if (dim == 2)
    {
      rho <- param[5]
      mat2 <- as.matrix((1 + df) * rho / (rho^2 - 1) + (2 + df) * (df * rho + x[,1] * x[,2])
                        / (df * (1 - rho^2) + x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]))
    }
  else 
    {    
      detsig <- det(sigma)
      invsig <- solve(sigma)

      l <- 1
      for (j in 1:(dim-1))
        for (i in (j+1):dim)
          {
            derdetsig <- 2 * det(sigma[-i,-j,drop=FALSE]) * (-1)^(i+j)
            derinvsig <- - invsig[,i] %*% t(invsig[,j]) - invsig[,j] %*% t(invsig[,i])
            firstterm <- derdetsig/detsig
            
            mat2[,l] <- - (firstterm + (df + dim) * rowSums((x %*% derinvsig) * x)
                           / (df +  rowSums((x %*% invsig) * x)) ) / 2
            l <- l + 1
          }
    }
  cbind(- mat %*% invdisp, - matrix(1/(2*disp^2),n,dim,byrow=TRUE) - (x %*% invdisp^2) * mat/2, mat2)
}

####################################################
