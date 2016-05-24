# polyfit computes and fits the OLS regression with polynomial terms

polyfitreg <- function(indnr, xv, yv, ch, sel, zv, vv) 
{  
  if (indnr == 2)
  {    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    
    # defining in- and output for the regression
    input <- cbind(rep(1, length(xv)), invx, invy, xv, yv, invx*invy, xv*invy, yv*invx,
                   xv*yv, xv^2, invx^2, yv^2, invy^2, xv^3, yv^3, invx^3, invy^3)
    inputs <- input[, sel]
    output <- ch  
    
    # computing the Ordinary Least Square Regression with no intercept
    res <- lm(output ~ inputs + 0, na.action=na.exclude) 
    B <- res$coef 
    
    # computing Error Sum of Squares (err_weights)
    if (length(sel) >= 2)
    {
      testoutput <- rowSums(inputs %*% diag(B))
    }  
    else
    {
      testoutput <- inputs * B
    }
    SEtest <- ((ch - testoutput)^2)    
    err_weights <- mean(SEtest)    
    
    # Log Likelihood and Rsquare
    LL <- sum(-SEtest/(2*var(ch)))
    print(c("LL:", LL), quote=F)
    Rsqr <- 1-mean(SEtest)/var(ch)
    print(c("Rsqr:", Rsqr), quote=F)
    
    return(list(B, err_weights))
  }
  
  ############################# indnr == 2 ends, indnr == 3 begins #################################
  
  if (indnr == 3)
  {
    # computing the inverse of x, y, z
    invx <- xv^(-1)
    invy <- yv^(-1)
    invz <- zv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    idx3 = which(is.infinite(invz)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    invz[idx3] <- NA
    
    # defining in- and output for the regression
    input <- cbind(rep(1, length(xv)), invx, invy, invz, xv, yv, zv, invx*invy, invy*invz, invx*invz, 
                   xv*yv, yv*zv, xv*zv, xv*invy, yv*invx, xv*invz, zv*invx, yv*invz, zv*invy, xv*invy*invz,
                   yv*invx*invz, zv*invx*invy, xv*yv*invz, yv*zv*invx, xv*zv*invy, xv*yv*zv, invx*invy*invz, 
                   xv^2, invx^2, yv^2, invy^2, zv^2, invz^2, xv^3, yv^3, zv^3, invx^3, invy^3, invz^3)
    inputs <- input[, sel]
    output <- ch    
    
    # computing the Ordinary Least Square Regression with no intercept
    res <- lm(output ~ inputs + 0, na.action=na.exclude) 
    B <- res$coef 
    # computing Error Sum of Squares (err_weights)
    if (length(sel) >= 2)
    {
      testoutput <- rowSums(inputs %*% diag(B))
    }  
    else
    {
      testoutput <- inputs * B
    }
    SEtest <- ((ch - testoutput)^2)    
    err_weights <- mean(SEtest)    
    
    # Log Likelihood and Rsquare
    LL <- sum(-SEtest/(2*var(ch)))
    print(c("LL:", LL), quote=F)
    Rsqr <- 1-mean(SEtest)/var(ch)
    print(c("Rsqr:", Rsqr), quote=F)
    
    return(list(B, err_weights))
  }
  
  ############################# indnr == 3 ends, indnr == 4 begins #################################
  
  if (indnr == 4)
  {
    # computing the inverse of x, y, z, v
    invx <- xv^(-1)
    invy <- yv^(-1)
    invz <- zv^(-1)
    invv <- vv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    idx3 = which(is.infinite(invz))
    idx4 = which(is.infinite(invv))
    invx[idx1] <- NA
    invy[idx2] <- NA
    invz[idx3] <- NA
    invv[idx4] <- NA
    
    # defining in- and output for the regression
    input <- cbind(rep(1, length(xv)), invx, invy, invz, invv, xv, yv, zv, vv, invx*invy, invy*invz, invx*invz, 
                   invx*invv, invy*invv, invz*invv, xv*yv, yv*zv, xv*zv, xv*vv, yv*vv, zv*vv, xv*invy, yv*invx, 
                   xv*invz, zv*invx, yv*invz, zv*invy, xv*invv, vv*invx, yv*invv, vv*invy, zv*invv, vv*invz,
                   xv*invy*invz, yv*invx*invz, zv*invx*invy, vv*invx*invy, vv*invx*invz, vv*invy*invz, xv*invy*invv, 
                   xv*invz*invv, yv*invx*invv, yv*invz*invv, zv*invx*invv, zv*invy*invv, xv*yv*invz, yv*zv*invx, 
                   zv*xv*invy, xv*yv*invv, yv*zv*invv, zv*xv*invv, xv*vv*invz, yv*vv*invz, yv*vv*invx, zv*vv*invx, 
                   vv*xv*invy, vv*zv*invy, xv*yv*zv, xv*yv*vv, xv*vv*zv, vv*yv*zv, invx*invy*invz, invx*invy*invv, 
                   invx*invv*invz, invv*invy*invz, xv*invv*invy*invz, yv*invx*invv*invz, zv*invx*invy*invv, 
                   vv*invx*invy*invz, xv*yv*zv*invv, xv*yv*vv*invz, xv*vv*zv*invy, vv*yv*zv*invx, xv*yv*invv*invz,
                   xv*zv*invv*invy, xv*vv*invy*invz, yv*zv*invv*invx, yv*vv*invz*invx, zv*vv*invx*invy, 
                   invx*invy*invz*invv, xv*yv*zv*vv, xv^2, invx^2, yv^2, invy^2, zv^2, invz^2, vv^2, invv^2,
                   xv^3, yv^3, zv^3, vv^3, invx^3, invy^3, invz^3, invv^3)
    inputs <- input[, sel]
    output <- ch    
    
    # computing the Ordinary Least Square Regression with no intercept
    res <- lm(output ~ inputs + 0, na.action=na.exclude) 
    B <- res$coef 
    
    # computing Error Sum of Squares (err_weights)
    if (length(sel) >= 2)
    {
      testoutput <- rowSums(inputs %*% diag(B))
    }  
    else
    {
      testoutput <- inputs * B
    }
    SEtest <- ((ch - testoutput)^2)    
    err_weights <- mean(SEtest)    
    
    # Log Likelihood and Rsquare
    LL <- sum(-SEtest/(2*var(ch)))
    print(c("LL:", LL), quote=F)
    Rsqr <- 1-mean(SEtest)/var(ch)
    print(c("Rsqr:", Rsqr), quote=F)
    
    return(list(B, err_weights))
  } 
}