# Bootstrap model stabilizer for multinomial logit models
# Version:            0.2
# Date:        2015-05-17
# Author: F.M., ctb: T.S.
# Note:   Needs nnet's multinom
# Further infos, references and credits:
#  See for nnet: Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth
#                Edition. New York: Springer.
# License: GPL-2 | GPL-3

BB.mod.stab.mlog <- function(data, BB.data, s.model, maxit.multi=3,...)
  {   
    options(warn=-1)
    regmod <- multinom(s.model,data=BB.data,trace=F,
                       na.action=na.exclude, ...=list(maxit=maxit.multi))
    options(warn=0)
    c.regmod <- multinom(s.model,data=data,trace=F,
                         na.action=na.exclude, ...=list(maxit=maxit.multi))
    
    misfactlevels <- !is.element(c.regmod$lev,regmod$lev)
    if (any(misfactlevels)==T) {
      pos.mis.lev <- which(misfactlevels)
      help.wts <- regmod$wts
      regmod$n <- c.regmod$n 
      regmod$nunits <- c.regmod$nunits
      regmod$nconn <- c.regmod$nconn
      regmod$conn <- c.regmod$conn
      multi1 <- c.regmod$n[1]+1
      regmod$wts <- c.regmod$wts
      fill.in.wts <- ((pos.mis.lev-1)*multi1+1):(multi1*pos.mis.lev)
      regmod$wts[fill.in.wts] <- -10e9
      regmod$wts[-fill.in.wts] <- help.wts
      regmod$lev <- c.regmod$lev
      regmod$lab <- c.regmod$lab
      regmod$edf <- c.regmod$edf
      
      help.fit <- regmod$fitted.values
      help.res <- regmod$residuals
      
      regmod$fitted.values <- c.regmod$fitted.values
      regmod$fitted.values[,pos.mis.lev] <- 0
      regmod$fitted.values[,-pos.mis.lev] <- help.fit
      
      regmod$residuals <- c.regmod$residuals
      regmod$residuals[,pos.mis.lev] <- 0
      regmod$residuals[,-pos.mis.lev] <- help.res
    }
    x <- list(model=regmod, c.model=c.regmod, misfactlevels=misfactlevels)
    return(x)
  }
