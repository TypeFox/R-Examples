#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       bayesHistogram.design.R             ####
####                                                 ####
#### FUNCTIONS:  bayesHistogram.design               ####
#########################################################

### ======================================
### bayesHistogram.design
### ======================================
## Subfunction for 'bayesHistogram' function
##  - extract response matrices
##
## 26/10/2004
##
bayesHistogram.design <- function(y1, y2)
{
  if (missing(y1)) stop("Response y1 must be given")
  else if (missing(y2)) dim <- 1
       else             dim <- 2

  ## Extract response for each dimension
  extract.y <- function(y.y){
    if (!inherits(y.y, "Surv")) stop("Response must be a survival object. ")
    type <- attr(y.y, "type")
    if (type == 'counting') stop ("Invalid survival type ('counting' is not implemented). ")
    n.n <- nrow(y.y)
    n.y.y <- ncol(y.y)
    y.y.left <- y.y[,1]
    if (n.y.y <= 2){
      y.y.right <- rep(1, n.n)
      status.y.y <- y.y[,2]
    }
    else{
      y.y.right <- y.y[,2]
      status.y.y <- y.y[,3]
    }
    return(list(left=y.y.left, right=y.y.right, status=status.y.y, n=n.n))
  }    
  resp1 <- extract.y(y1)
  if (dim >= 2){
    resp2 <- extract.y(y2)
    if (resp1$n != resp2$n) stop ("Bivariate response of different length supplied")

    status <- rbind(resp1$status, resp2$status)
    y.left <- rbind(resp1$left, resp2$left)
    y.right <- rbind(resp1$right, resp2$right)    
  }
  else{
    status <- matrix(resp1$status, nrow=1)
    y.left <- matrix(resp1$left, nrow=1)
    y.right <- matrix(resp1$right, nrow=1)    
  }    

  return(list(y.left=y.left, y.right=y.right, status=status, dim=dim))
}
