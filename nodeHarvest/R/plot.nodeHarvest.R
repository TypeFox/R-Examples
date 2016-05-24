plot.nodeHarvest <-
function(x,  XTEST=NULL, highlight=NULL, varnames= NULL, yoffset=0.12, labels="all",  cexfaclab=1, ...){

  Z <- x[["nodes"]]
  weight <- sapply(x[["nodes"]],attr,"weight")
  if(is.null(varnames)) varnames <- x[["varnames"]]
  cexfac <- 10
  colgrey <- rgb(0,0,0,0.33) 
  colline <- rgb(0.3,0.3,0.3,0.4)
  colred <- rgb(1,.1,0.1,0.6)

  attri <- attributes(Z)
  ord <- order(sapply(Z,attr,"n"))
  Z <- Z[ord]
  weight <- weight[ord]
  attributes(Z) <- attri
  
  
  isgrey <- rep(TRUE,length(Z))
  if(!is.null(highlight) ){
    if(highlight < 0 ){
      highlZ <- which(sapply(Z,attr,"leaf"))
    }else{
      IsignTEST <- abs(sign(getI(Z,XTEST[highlight,,drop=FALSE])$I))
      highlZ <- which( apply(IsignTEST,2,sum)!=0 )
    }
    if(length(highlZ)>0)  isgrey[ highlZ] <- FALSE
  }
    

  allmean <- sapply(Z,attr,"mean")
  alln <-  rank(sapply(Z,attr,"n") + seq(0.1,0,length=length(weight)) ,ties.method="first" )
  rang <- range(allmean)
  ylim <- c( min(alln)*0.9 , max(alln)*1.1 ) 
  xtmp <- allmean 
  xlimall <- mean(xtmp) + 1.1 * (1+(cexfaclab-1)/2) * range(xtmp-mean(xtmp))
  yvec <- alln
  cexx <- (cexfac*sqrt(weight/max(weight)))
  
  plot(xtmp , yvec, cex=cexx, ylim=ylim, type="n", axes=FALSE,xlab="RESPONSE",ylab="SAMPLES",xlim=xlimall)
  
  atvec <- round(seq(1,length(alln),length=6))
  axis(2,at=atvec,labels= sapply(Z,attr,"n")[atvec]  )
  

  for (k in 1:length(weight)){
    anc <- attr(Z[[k]],"ancestors")
    if(length(anc)>=1){
      choose <- anc[1]
      co <- c(xtmp[k] , yvec[k])
      coc <- c(xtmp[choose],yvec[choose])
      lines( (c(coc[1], co[1])), c(coc[2],co[2]),col=colline,lwd=1)
    }
  }



  points( xtmp , yvec ,cex=cexx,col="white", pch=20)
  

  if(all(isgrey)){
    points( xtmp , yvec ,cex=cexx,col=colgrey, pch=20)
  }else{
    indgrey <- which(isgrey)
    indred <- which(!isgrey)
    points( xtmp[indgrey] , yvec[indgrey] ,cex=cexx[indgrey],col=colgrey, pch=20)
    points( xtmp[indred] , yvec[indred] ,cex=cexx[indred],col=colred, pch=20)
  }

  axis(1)

  if(labels!="none"){
    drawlab <- if(labels=="all") 1:length(weight) else which(!isgrey)
    for (ky in drawlab){
      drawtext(Z,ky,offset= 0.00*abs(diff(xlimall)),varnames=varnames ,yoffset=1*yoffset,cex=0.55 * cexfaclab)    
    }
  }
    
   
}

