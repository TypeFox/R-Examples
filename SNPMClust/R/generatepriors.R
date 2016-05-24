generatepriors <-
function(x,y,calls,priorpoints=length(x)*.2,xm1=NA,xm2=NA,xm3=NA,
    ym1=NA,ym2=NA,ym3=NA,ranseed=ranseed) {
  set.seed(ranseed)
  x <- x[is.element(calls,c("AA","AB","BB"))]  
  y <- y[is.element(calls,c("AA","AB","BB"))]
  calls <- calls[is.element(calls,c("AA","AB","BB"))]
  if (is.na(xm1)) xmean1 <- mean(x[calls=="AA"]) else xmean1 <- xm1
  if (is.na(ym1)) ymean1 <- mean(y[calls=="AA"]) else ymean1 <- ym1
  if (is.na(xm2)) xmean2 <- mean(x[calls=="AB"]) else xmean2 <- xm2
  if (is.na(ym2)) ymean2 <- mean(y[calls=="AB"]) else ymean2 <- ym2
  if (is.na(xm3)) xmean3 <- mean(x[calls=="BB"]) else xmean3 <- xm3
  if (is.na(ym3)) ymean3 <- mean(y[calls=="BB"]) else ymean3 <- ym3
  xmeanvect <- x
  ymeanvect <- y
  xmeanvect[calls=="AA"]<-xmean1
  xmeanvect[calls=="AB"]<-xmean2
  xmeanvect[calls=="BB"]<-xmean3
  ymeanvect[calls=="AA"]<-ymean1
  ymeanvect[calls=="AB"]<-ymean2
  ymeanvect[calls=="BB"]<-ymean3
  xsd <- sqrt(sum((x-xmeanvect)^2)/length(x))
  ysd <- sqrt(sum((y-ymeanvect)^2)/length(y))
  if (length(table(calls))==3) {
    if (table(calls)["AA"] <= table(calls)["BB"]) {
      xpseudo <- c(rnorm(priorpoints,mean=xmean1,sd=xsd),rnorm(priorpoints,mean=xmean2,sd=xsd))
      ypseudo <- c(rnorm(priorpoints,mean=ymean1,sd=ysd),rnorm(priorpoints,mean=ymean2,sd=ysd))
    }
    if (table(calls)["BB"] < table(calls)["AA"]) {
      xpseudo <- c(rnorm(priorpoints,mean=xmean3,sd=xsd),rnorm(priorpoints,mean=xmean2,sd=xsd))
      ypseudo <- c(rnorm(priorpoints,mean=ymean3,sd=ysd),rnorm(priorpoints,mean=ymean2,sd=ysd))
    } 
  }
  if (length(table(calls))==2) {
    if (names(which(table(calls)==min(table(calls))))=="AA") {
      xpseudo <- rnorm(priorpoints,mean=xmean1,sd=xsd)
      ypseudo <- rnorm(priorpoints,mean=ymean1,sd=ysd)
    }
    if (names(which(table(calls)==min(table(calls))))=="AB") {
      xpseudo <- rnorm(priorpoints,mean=xmean2,sd=xsd)
      ypseudo <- rnorm(priorpoints,mean=ymean2,sd=ysd)
    }
    if (names(which(table(calls)==min(table(calls))))=="BB") {
      xpseudo <- rnorm(priorpoints,mean=xmean3,sd=xsd)
      ypseudo <- rnorm(priorpoints,mean=ymean3,sd=ysd)
    }
  }
  if (length(table(calls))<2) {
    xpseudo<-NA
    ypseudo<-NA
  } 	 
  if (!is.na(xpseudo) && !is.na(ypseudo)) {return(cbind(xpseudo,ypseudo))} else return(NA)
}
