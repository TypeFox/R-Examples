f.anal.ldt <- function(coords,info)
  {
    ## Copied over from ldt-analytic.r on Dec 10th
    ## trying to incorporate apply into my code
    
    ## created May 12th
    ## analytic results for the ldt model
    ## see write up and notes in binder for more details.

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
      if(i !=1 && j !=1)
        {
          ## most common situation:  different row, different column
          ## simplified different row, column result
          ## this is bad due to 0/0 and log(0) see below
          ## also its old, note the double for loops
          if(FALSE) {
            as <- c(-ax,ax-a,ax,ax+a)
            bs <- c(-bx,bx-b,bx,bx+b)
            ret <- 100*a^2*b^2
            for(ii in 1:4)
              for(jj in 1:4)
                {
                  ret <- ret + (-1)^(ii+jj+1)*(8*as[ii]*bs[jj]*(bs[jj]^2*atan(as[ii]/bs[jj])+as[ii]^2*atan(bs[jj]/as[ii]))-(as[ii]^4-6*as[ii]^2*bs[jj]^2+bs[jj]^4)*log(as[ii]^2+bs[jj]^2))
                } 
            ret <- ret/(48*a^2*b^2)
          }
          
          ## different version, takes care of 0/0, or log(0)
          ## doesn't have for loops
          ret <- 100*a^2*b^2
          as <- c(-ax,ax-a,ax,ax+a)
          bs <- c(-bx,bx-b,bx,bx+b)
          asbs <- cbind(as.vector(matrix(as,4,4)),as.vector(matrix(1:4,4,4)),as.vector(t(matrix(bs,4,4))),as.vector(t(matrix(1:4,4,4))))
                  
          ret <- 100*a^2*b^2
          ret <- ret + sum(apply(asbs,1,sepRowsepCol.ldt.res))
          ret <- ret/(48*a^2*b^2)        
        } else if(i == 1 && j !=1)
          {
            ## same column, different row
            ax <- bx
            bx <- 0
            a <- b
            b <- info$rowwidth
            ret <- (50*a^2*b^2-8*b^3*(ax-a)*atan((ax-a)/b)+
                    16*ax*b^3*atan(ax/b)-8*b^3*(ax+a)*atan((ax+a)/b)+
                    16*ax^3*b*atan(b/ax)-8*b*(ax-a)^3*atan(b/(ax-a))-
                    8*b*(ax+a)^3*atan(b/(ax+a))+4*ax^4*log(ax)-
                    2*(ax+a)^4*log(ax+a)+
                    (((ax-a)^2-b^2)^2-4*b^2*(ax-a)^2)*log((ax-a)^2+b^2)-
                    2*((ax^2-b^2)^2-4*b^2*ax^2)*log(ax^2+b^2)+
                    (((ax+a)^2-b^2)^2-4*b^2*(ax+a)^2)*log((ax+a)^2+b^2))
            if(ax>a) ret <- ret - 2*(ax-a)^4*log(ax-a)
            ret <- ret/(24*a^2*b^2)
          } else if(i!=1 && j==1)
            {
              ## simplified same row result using mathematica
              ## problems when ax=a with log(ax-a) dealt with below
              ret <- (50*a^2*b^2-8*b^3*(ax-a)*atan((ax-a)/b)+
                      16*ax*b^3*atan(ax/b)-8*b^3*(ax+a)*atan((ax+a)/b)+
                      16*ax^3*b*atan(b/ax)-8*b*(ax-a)^3*atan(b/(ax-a))-
                      8*b*(ax+a)^3*atan(b/(ax+a))+4*ax^4*log(ax)-
                      2*(ax+a)^4*log(ax+a)+
                      (((ax-a)^2-b^2)^2-4*b^2*(ax-a)^2)*log((ax-a)^2+b^2)-
                      2*((ax^2-b^2)^2-4*b^2*ax^2)*log(ax^2+b^2)+
                      (((ax+a)^2-b^2)^2-4*b^2*(ax+a)^2)*log((ax+a)^2+b^2))
              if(ax>a) ret <- ret - 2*(ax-a)^4*log(ax-a)
              ret <- ret/(24*a^2*b^2)
            } else {
              ## diagonal result, the most rare
              ret <- (25*a^2*b^2-8*a*b^3*atan(a/b)-8*a^3*b*atan(b/a)-2*a^4*log(a)-2*b^4*log(b)+(a^4-6*a^2*b^2+b^4)*log(a^2+b^2))/(12*a^2*b^2)
            }
    }
    ret
  }
   
sepRowsepCol.ldt.res <- function(abVal)
  {
    aVal <- abVal[1]
    ii <- abVal[2]
    bVal <- abVal[3]
    jj <- abVal[4]
    
    numSteps <- 2-((aVal==0) + (bVal==0))
    ret <- 0
    if(numSteps>=1) ret <- ret + (-1)^(ii+jj+2)*(aVal^4-6*aVal^2*bVal^2+bVal^4)*log(aVal^2+bVal^2)
    if(numSteps>1) ret <- ret + (-1)^(ii+jj+1)*(8*aVal*bVal*(bVal^2*atan(aVal/bVal)+aVal^2*atan(bVal/aVal)))
    ret
  }
