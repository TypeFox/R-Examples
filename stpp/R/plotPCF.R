plotPCF <- function(PCF,n=15,type="contour",legend=TRUE,which=NULL,main=NULL,...)
{
devl=dev.list()
if(is.null(devl))
{
old.par <- par(no.readonly = TRUE)    
on.exit(par(old.par))
}
else
{
dev.off(devl[length(devl)])
dev.new()
old.par <- par(no.readonly = TRUE)    
on.exit(par(old.par))
}


correc=c("none","isotropic","border","modified.border","translate")
correc2=PCF$correction
id <- match(correc2, correc, nomatch = NA)

if ((is.null(which) && length(id)>1) || any(is.na(match(which, correc, nomatch = NA))))
 {
    mess <- paste("Please specify the argument \'which\', among:", paste(dQuote(correc2), 
            collapse = ", "))
    stop(mess, call. = FALSE)
 }

if (is.matrix(PCF$pcf))
 { 
	if (is.null(which)) which=correc2
	else {
	if (!(is.null(which)) && which!=correc2) 
	{	
    	mess <- paste("Argument \'which\' should be", paste(dQuote(correc2), 
            collapse = ", "))
      stop(mess, call. = FALSE) 
	} }
 }

if (!is.matrix(PCF$pcf))
 {
  id <- match(which, correc2, nomatch = NA)
	if(is.na(id)) 
	{
    	mess <- paste("Please specify the argument \'which\', among:", paste(dQuote(correc2), 
            collapse = ", "))
    	stop(mess, call. = FALSE)
	}
	else PCF$pcf=PCF$pcf[[id]]
 }


 k=PCF$pcf
 K=PCF

if(!is.null(main)) {titl=main; subtitl=""}
else {
 titl <- expression(hat(g)* group("(",list(u,v),")") )
 subtitl <- paste("edge correction method: ",which,sep="")
     }


typeplot=c("contour","image","persp")
id <- match(type, typeplot, nomatch = NA)
  if (any(nbg <- is.na(id))) {
        mess <- paste("unrecognised plot type:", paste(dQuote(type[nbg]), 
            collapse = ", "))
        stop(mess, call. = FALSE)
    }
if((length(id)!=1)||is.na(id)) stop("Please specify one type among \"contour\", \"image\" and \"persp\" ")
typeplot=rep(0,3)
typeplot[id]=1

  colo <- colorRampPalette(c("red", "white", "blue"))
  M <- max(abs(range(k)))
  M <- pretty(c(-M,M),n=n)

  if (max(M)>15) M <- c(-max(M),-10, pretty(c(-5,5),n=n), 10, max(M))
  else if (max(M)>10 & max(M)<15) M <- c(-15,-10, pretty(c(-5,5),n=n), 10, 15)
  
  n <- length(M)
  COL <- colo(n)
  if (typeplot[3]==1)
    {
      mask <- matrix(0,ncol=length(K$times),nrow=length(K$dist))
      for(i in 1:length(K$dist)){ for(j in 1:length(K$times)){mask[i,j] <- COL[findInterval(x=k[i,j],vec=M)]}}
      COL <- mask[1:(length(K$dist)-1),1:(length(K$times)-1)]
      
      if(isTRUE(legend))
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=1,mar=c(0,0,3,0))
          par(fig=c(0,0.825,0,1))
          persp(x=K$dist, y=K$times, z=k, xlab="u",ylab="v", zlab="",expand=1, col=COL, ...)
          title(titl,cex.main=1.5,sub=subtitl,outer=TRUE,line=-1)
          par(fig=c(0.825,1,0,1))
          mini <- findInterval(x=min(k,na.rm=TRUE),vec=M)
          maxi <- findInterval(x=max(k,na.rm=TRUE),vec=M)
          legend("right",fill=colo(n)[maxi:mini],legend=M[maxi:mini],horiz=F,bty="n")
        }
      else
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=1)
          persp(x=K$dist, y=K$times, z=k, xlab="u",ylab="v", zlab="", expand=1, col=COL, ...)
          title(titl,cex.main=1.5,sub=subtitl)
        }
    }
   if(typeplot[1]==1)
    {
      if(isTRUE(legend))
        {
          par(cex.lab=1.5,cex.axis=1.5,font=2,plt=c(0,1,0,1),lwd=1,mar=c(0.5,0.5,2.5,0.5),las=1)
          par(fig=c(0.1,0.825,0.1,1))
          contour(K$dist, K$times, k, labcex=1.5,levels=M,drawlabels=F,col=colo(n),zlim=range(M),axes=F)
          box(lwd=2)
          at <- axTicks(1)
          axis(1,at=at[1:length(at)],labels=at[1:length(at)])
          at <- axTicks(2)
          axis(2,at=at[1:length(at)],labels=at[1:length(at)])
          title(titl,cex.main=1.5,sub=subtitl,outer=TRUE,line=-1)
          par(fig=c(0,1,0.1,1))
          mini <- findInterval(x=min(k,na.rm=TRUE),vec=M)
          maxi <- findInterval(x=max(k,na.rm=TRUE),vec=M)
          legend("right",fill=colo(n)[maxi:mini],legend=M[maxi:mini],horiz=F,bty="n")
        }
      else
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=2,las=1)
          contour(K$dist, K$times, k, labcex=1.5,levels=M,drawlabels=T,col=colo(n),zlim=range(M),axes=F)
          box(lwd=2)
          at <- axTicks(1)
          axis(1,at=at[1:length(at)],labels=at[1:length(at)])
          at <- axTicks(2)
          axis(2,at=at[1:length(at)],labels=at[1:length(at)])
          title(titl,cex.main=1.5,sub=subtitl)
        }
    }
  if(typeplot[2]==1)
    {
      if(isTRUE(legend))
        {
          par(cex.lab=1.5,cex.axis=1.5,font=2,lwd=1,plt=c(0,1,0,1),mar=c(0.5,0.5,2.5,0.5),las=1)
          par(fig=c(0.1,0.825,0.1,1))
          image(K$dist, K$times, k, col=colo(n),zlim=range(M),axes=F,xlab="",ylab="")
          box(lwd=2)
          at <- axTicks(1)
          axis(1,at=at[1:length(at)],labels=at[1:length(at)])
          at <- axTicks(2)
          axis(2,at=at[1:length(at)],labels=at[1:length(at)])
          title(titl,cex.main=1.5,sub=subtitl,outer=TRUE,line=-1)
          par(fig=c(0,1,0.1,1))
          mini <- findInterval(x=min(k,na.rm=TRUE),vec=M)
          maxi <- findInterval(x=max(k,na.rm=TRUE),vec=M)
          legend("right",fill=colo(n)[maxi:mini],legend=M[maxi:mini],horiz=F,bty="n")
        }
      else
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=2,las=1)
          image(K$dist, K$times, k, col=colo(n),zlim=range(M),axes=F,xlab="",ylab="")
          box(lwd=2)
          at <- axTicks(1)
          axis(1,at=at[1:length(at)],labels=at[1:length(at)])
          at <- axTicks(2)
          axis(2,at=at[1:length(at)],labels=at[1:length(at)])
          title(titl,cex.main=1.5,sub=subtitl)
        }
    }

par(old.par)
}
  
