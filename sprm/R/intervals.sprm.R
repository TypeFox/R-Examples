intervals.sprm <- function(x,data, optcomp, n, df, sigma2, covb, alpha){
  ynames <- names(x$YMeans)
  used.vars <- x$used.vars[[2*optcomp]]
  if(missing(data)){ 
    yp <- predict(x)
    y0 <- x$inputs$y0
    Xn <- as.matrix(x$inputs$Xs[,used.vars])
  } else { 
    yp <- predict(x,data)
    yindex <- which(colnames(data)==ynames)
    y0 <- data[,yindex]
    Xnames <- names(x$XMeans) #
    Xindex <- which(colnames(data)%in%Xnames) #
    if (length(Xindex)!=length(Xnames)){
      stop("Column names of data don't match variable names in the model.")
    }
    Xn <- scale(data[,Xindex], center=x$XMeans, scale=x$Xscales)[,used.vars] # attr(yp,"Xn.scaled")[,] #
    if (length(rownames(data))==0){
      rownames(Xn) <- 1:dim(Xn)[1]
    } else{
      rownames(Xn) <- rownames(data) 
    }
  }
  mT <- apply(x$scores,2,median)
  sT <- apply(x$scores,2,qn)
  Tn <- Xn%*%x$R[used.vars,]
  dt <- scale(Tn,center=mT,scale=sT)
  wtn <- sqrt(apply(dt^2, 1, sum))
  wtn <- wtn/median(wtn)
  fun <- x$inputs$fun 
  constants <- x$inputs$constants
  probp1 <- constants[1]
  probct <- qchisq(probp1,optcomp)
  if(fun=="Fair"){
    wte <- 1/(1 + abs(wtn/(probct*2)))
  } 
  if(fun=="Huber") {
    wte <- wtn
    wte[which(wtn <= probct)] <- 1
    wte[which(wtn > probct)] <- probct/abs(wtn[which(wtn > probct)]) 
    wte[which(wtn<wte)] <- wtn[which(wtn<wte)] 
  }
  if(fun=="Hampel") {
    hampelp2 <- constants[2]
    hampelp3 <- constants[3]
    hampelb <- qchisq(hampelp2, optcomp)
    hampelr <- qchisq(hampelp3, optcomp)
    wte <- wtn
    wte[which(wtn <= probct)] <- 1 
    wte[which(wtn > probct & wtn <= hampelb)] <- probct/abs(wtn[which(wtn > probct & wtn <= hampelb)])
    wte[which(wtn > hampelb & wtn <= hampelr)] <- probct*(hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)]))/(hampelr -hampelb)*1/abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
    wte[which(wtn > hampelr)] <- 0
  }
  outs <- wte
  outs[which(outs==0)] <- 2
  outs[which(outs>0 & outs < 1)] <- 1
  outs[which(outs==1)] <- 0
  
  sigma2 <- sqrt(sigma2+x$Yscales^2)
  covb <- covb*sigma2
  var.yp <- (sigma2*(1+1/n) + Xn%*%covb%*%t(Xn))/wte
  sd.yp <- sqrt(diag(var.yp))
  ulim <- yp + qt(1-alpha/2,df)*sd.yp 
  llim <- yp - qt(1-alpha/2,df)*sd.yp
  
  return(data.frame(y_original=y0,y_predicted=yp,nam=rownames(Xn),ulim=ulim,llim=llim, outs=outs))
}  