drawtext <-
function(Z,k,varnames= as.character(1:100) ,offset=0.1,yoffset=0.1,cex=0.75,alln=NULL,plot=TRUE){

 
  largersmaller <- Z[[k]]
  allmean <- sapply(Z,attr,"mean")
  if(is.null(alln)) alln <- sapply(Z,attr,"n")
  alln <-  rank(sapply(Z,attr,"n") + seq(0.1,0,length=length(Z)) ,ties.method="first" )

  
  offy <-  yoffset * (1:nrow(largersmaller))
  offy <- offy * max(alln) /5
  offy <- offy-mean(offy)

  if(!plot) rl <- list()
  for (kk in 1:nrow(largersmaller)){
    xval <- allmean[k] +  offset
    yval <-  alln[k] + offy[kk]
    var <- varnames[largersmaller[kk,"variable"]]
    larg <- signif(largersmaller[kk,"lower"],2)
    smal <- signif(largersmaller[kk,"upper"],2)
    nfact <- length(levelsu <- attr(Z,"levelvec")[[largersmaller[kk,"variable"]]])
    if(nfact<1){
      textplot <-  paste(if(!is.inf(larg)) paste( larg," <= ",sep="") else "", var, if(!is.inf(smal)) paste(" <= ", smal,sep="") else ""   ,sep="")
    }else{
      lev <- levelsu[dectobin(larg,nl=nfact)>0.5]
      if(length(lev)==1){
        textplot <- paste( var, " = ", lev,sep="")
      }else{
       if(length(lev)<=ceiling(nfact*0.55)){
         textplot <- paste( var, " in {", paste(lev,collapse=","),"}",sep="")
       }else{
         textplot <- paste( var, " not in {", paste(levelsu[dectobin(larg,nl=nfact)<0.5],collapse=","),"}",sep="")

       }
      }
    }
    if(plot){
      if(abs(larg)==Inf & abs(smal)==Inf){
        text(xval,yval, "root", cex=cex)
      }else{
        text(xval,yval,textplot,cex=cex)
      }
    }else{
      rl[[kk]] <- textplot
    }
  }
  if(!plot) return(rl)
  
}

