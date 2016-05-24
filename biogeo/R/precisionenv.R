precisionenv <-
function (dat, rst, x = "x", y = "y") 
{
  cn <- names(dat) #Column names of data
  f1 <- match(x, cn)
  f2 <- match(y, cn)
  dx <- abs(coord2numeric(dat[, f1]))
  dy <- abs(coord2numeric(dat[, f2]))
  xm <- dx - floor(dx)
  ym <- dy - floor(dy)
  xm <- round(xm, digits = 7)
  ym <- round(ym, digits = 7)  
  rs <- c(10,15,20,30,60) #Resolutions at which precision will be checked
  rr<-res(rst)
  resrst <- 1/rr[1]/0.6
  rs <- rs[rs>resrst]
  if(length(rs)>0){
    cn <- paste("p", rs, "m", sep = "") #Names of resulotions
    ers <- {
    } 
    ersc <- {
    }
    for (j in 1:length(rs)) {
      res <- rs[j]
      mp <- (seq(0, (60 - res), res) + (res/2))/60
      mc <- (seq(0, (60 - res), res))/60
      mp <- round(mp, digits = 7)
      mc <- round(mc, digits = 7)
      er <- rep(2, length(xm))
      erc <- er
      for (i in 1:length(xm)) {
        xp <- (is.element(xm[i], mp)) * 1
        yp <- (is.element(ym[i], mp)) * 1
        xc <- (is.element(xm[i], mc)) * 1
        yc <- (is.element(ym[i], mc)) * 1
        er[i] <- ifelse(yp == 1 & xp == 1, 1, 0)
        erc[i] <- ifelse(yc == 1 & xc == 1, 1, 0)
      }
      ers <- cbind(ers, er)
      ersc <- cbind(ersc, erc)
    }
    ers3 <- as.data.frame(ers + ersc)
    names(ers3) <- cn
    envpreci <- apply(ers3, MARGIN=1, FUN=function(x) any(x==1)*1)
#     dat$Exclude[envpreci==1] <- 1
#     dat$Reason <- as.character(dat$Reason)
#     dat$Reason[envpreci==1] <- 'Precision'
    z <- cbind(dat, ers3, envpreci)
    return(z)
  } else {
    stop('Raster resolution is too coarse for this function')
  }
}
