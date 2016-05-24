f.anal.power <- function(coords, h, info)
  {
    ## created Oct 10th
    ## K <- function(d,params) -d^(2*params[1])/params[1]
    ## analytic results for the power model
    ## see write up and notes in binder for more details.

    ## OCT 13th 2003
    ## include an option for the mode of computing the 2F1
    ## values include NI, SUM and NR for
    ## numerical integration, summation and numerical recipes
    
    a <- info$rowwidth
    b <- info$colwidth
    ax <- coords[1]
    bx <- coords[2]
    i <- coords[3]
    j <- coords[4]
    ## 5
    ## 6
    ## 7 - which rows to evaluate - evalFactor
    ## 8 - lengths
    ## 9 - indicator for which rows to evaluate analytic results
    evaluate <- coords[9]
    ## 10 - lower limit for numerical integration
    ## 11 - upper limit for numerical integration

    ret <- 0
    if(evaluate) {
      h <- 2*h
      
      ## most common situation first
      if(i != 1 && j != 1) {
        ## different version, which takes care of possible problems like 0/0, or log(0)
        alphas <- c(ax, ax-a, ax, ax+a)
        betas <- c(bx,bx-b,bx,bx+b)
        asbs <- cbind(as.vector(matrix(alphas,4,4)),as.vector(matrix(1:4,4,4)),as.vector(t(matrix(betas,4,4))),as.vector(t(matrix(1:4,4,4))))
        ret <- 0
        
        ## function is defined below
        ret <- sum(apply(asbs,1,sepRowsepCol.power.res,h=h))
        ret <- ret/(4.*a*a*b*b)
      } else {
        
        if(i!=1 && j==1) {
          ## simplified same row result using mathematica
          ## problems when ax=a with log(ax-a) dealt with below
          alphas <- c(ax, ax-a, ax, ax+a)
          ret <- 0
          
          ## function is defined below
          ret <- sum(apply(cbind(alphas,1:4),1,sameRow.power.res,b,h))
          ret <- ret/(2.*a*a*b*b)
        } else {
          if(i==1 & j!=1) {
            ## SAME COLUMN
            alphas <- c(bx, bx-b, bx, bx+b)
            ret <- 0
            
            ## function is defined below
            ret <- sum(apply(cbind(alphas,1:4),1,sameRow.power.res,b=a,h=h))
            ret <- ret/(2.*a*a*b*b)
          } else {
            if(i==1 && j==1) {
              ## simplified diagaonal result using mathematica
              asqr <- a*a
              bsqr <- b*b
              dsqr <- asqr+bsqr
              ret <- -(2/(h+2)+2/(h+4))*(dsqr)^((h+4)/2)
              ret <- ret + (2/(h+2)-4/(h+3)+2/(h+4))*(a^(h+4)+b^(h+4))
              ret <- ret + 4/3*(dsqr)^(h/2)*(bsqr*bsqr*Hypergeometric2F1(-h/2,1,2.5,bsqr/(dsqr))+asqr*asqr*Hypergeometric2F1(-h/2,1,2.5,asqr/(dsqr)))
              ret <- ret + 4*asqr*bsqr/((h+2)*sqrt(dsqr))*(a^(1+h)*Hypergeometric2F1((3+h)/2,0.5,1.5,bsqr/(dsqr)) + b^(1+h)*Hypergeometric2F1((3+h)/2,0.5,1.5,asqr/(dsqr)))
              ret <- ret /(asqr*bsqr)
            }
          }
        }
      }
      ret <- -ret/(h/2)
    }
    ret
  }

sepRowsepCol.power.res <- function(abVal,h)
  {
    aVal <- abVal[1]
    ii <- abVal[2]
    bVal <- abVal[3]
    jj <- abVal[4]
    asqr <- aVal*aVal
    bsqr <- bVal*bVal
    dsqr <- asqr+bsqr
    sign <- (-1)^(ii+jj+1)
    ab3h <- c(aVal,bVal)^(3+h)
    dsqrh2 <- dsqr^(h/2)
    
    numSteps <- 2-((aVal==0) + (bVal==0))
    ret <- 0
    ret <- ret + sign*(2/(h+2)+2/(h+4))*(dsqr*dsqr*dsqrh2)
    if(is.nan(ret)) ret <- 0
    if(numSteps) {
      ret <- ret - sign*4/3*(dsqrh2)*(asqr*asqr*Hypergeometric2F1(-h/2,1,2.5,asqr/(dsqr))+bsqr*bsqr*Hypergeometric2F1(-h/2,1,2.5,bsqr/(dsqr)))
      ret <- ret - sign*4/((h+2)*sqrt(dsqr))*(ab3h[1]*bsqr*Hypergeometric2F1((3+h)/2,0.5,1.5,bsqr/(dsqr)) + ab3h[2]*asqr*Hypergeometric2F1((3+h)/2,0.5,1.5,asqr/(dsqr)))
    }
    ret
  }

sameRow.power.res <- function(abVal,b,h)
  {
    aVal <- abVal[1]
    ii <- abVal[2]
    asqr <- aVal*aVal
    bsqr <- b*b
    dsqr <- bsqr + asqr
    ret <- 0
    sign <- (-1)^(ii)
    ab1h <- c(aVal,b)^(1+h)
    dsqrh2 <- (dsqr)^(h/2)
    
    ret <- ret - sign*(2/(2+h)+2/(4+h))*(dsqr*dsqr*dsqrh2)
    ret <- ret + sign*(2/(2+h)-4/(3+h)+2/(4+h))*aVal^(4+h)
    
    ret <- ret + sign*4/3*(dsqrh2)*(bsqr*bsqr*Hypergeometric2F1(-h/2,1,2.5,bsqr/(dsqr))+asqr*asqr*Hypergeometric2F1(-h/2,1,2.5,asqr/(dsqr)))
    
    if(aVal) ret <- ret + sign*4*asqr*bsqr/((h+2)*sqrt(dsqr))*(ab1h[1]*Hypergeometric2F1((3+h)/2,0.5,1.5,bsqr/(dsqr)) + ab1h[2]*Hypergeometric2F1((3+h)/2,0.5,1.5,asqr/(dsqr)))            
    ret
  }
