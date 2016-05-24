plotK <- function(K,n=15,L=FALSE,type="contour",legend=TRUE,which=NULL,main=NULL,...)
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
correc2=K$correction
id <- match(correc2, correc, nomatch = NA)


if ((is.null(which) && length(id)>1) || any(is.na(match(which, correc, nomatch = NA))))
 {
    mess <- paste("Please specify the argument \'which\', among:", paste(dQuote(correc2), 
            collapse = ", "))
    stop(mess, call. = FALSE)
 }
 
if (isTRUE(K$infectious)) which="isotropic"


if (is.matrix(K$Khat))
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

if (!is.matrix(K$Khat))
 {
	id <- match(which, correc2, nomatch = NA)
	if(is.na(id)) 
	{
    	mess <- paste("Please specify the argument \'which\', among:", paste(dQuote(correc2), 
            collapse = ", "))
    	stop(mess, call. = FALSE)
	}
	else K$Khat=K$Khat[[id]]
 }



  if(!is.null(main)) {
   titl=main; subtitl=""
    if (isTRUE(L)) k <- K$Khat-K$Ktheo
    else k <- K$Khat 
   }
  else {
  if (isTRUE(L))
    {
      k <- K$Khat-K$Ktheo
      subtitl <- paste("edge correction method: ",which,sep="") 
      if (isTRUE(K$infectious))
        titl <- expression(hat(K)[ST] * group("(",list(u,v),")") - pi*u^2*v)
	else
        titl <- expression(hat(K)[ST] * group("(",list(u,v),")") - 2*pi*u^2*v) 
    }
  else
    {
      k <- K$Khat
      titl = expression(hat(K)[ST] * group("(",list(u,v),")"))
      subtitl <- paste("edge correction method: ",which,sep="") 
    }
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
  
