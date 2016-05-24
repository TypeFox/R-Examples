magerr <-
function(x, y, xlo, ylo, xhi=xlo, yhi=ylo, corxy, length=0.02, col='black',fill=FALSE,...){
if(length(col)==1){col=rep(col,length(x))}
if(!missing(corxy)){
  errbarsel=which(xlo==0 | ylo==0)
}else{
  errbarsel=1:length(x)
}
if(length(errbarsel)>0){
  if(!missing(xlo)){xlodraw=x-abs(xlo);doxlo=TRUE}else{xlodraw=0;doxlo=FALSE}
  if(!missing(xlo) | !missing(xhi)){xhidraw=x+abs(xhi);doxhi=TRUE}else{xhidraw=0;doxhi=FALSE}
  if(!missing(ylo)){ylodraw=y-abs(ylo);doylo=TRUE}else{ylodraw=0;doylo=FALSE}
  if(!missing(ylo) | !missing(yhi)){yhidraw=y+abs(yhi);doyhi=TRUE}else{yhidraw=0;doyhi=FALSE}
  
  if(doxlo & par()$xlog){
  	sel=which(xlodraw<=0)
  	xlodraw[sel]=1e-300
  }
  if(doxhi & par()$xlog){
    sel=which(xhidraw>0 & x<0)
    x[sel]=1e-300
  }
  if(doylo & par()$ylog){
  	sel=which(ylodraw<=0)
  	ylodraw[sel]=1e-300
  }
  if(doyhi & par()$ylog){
    sel=which(yhidraw>0 & y<0)
    y[sel]=1e-300
  }
  xarrow=x[errbarsel]
  yarrow=y[errbarsel]
  xlodraw=xlodraw[errbarsel]
  xhidraw=xhidraw[errbarsel]
  ylodraw=ylodraw[errbarsel]
  yhidraw=yhidraw[errbarsel]
  colarrow=col[errbarsel]
  if(doxlo & any(is.finite(xlodraw))){
    finalsel=xarrow>xlodraw
    arrows(xarrow[finalsel],yarrow[finalsel],xlodraw[finalsel],yarrow[finalsel],angle=90,length=length,col=colarrow[finalsel],...)
  }
  if(doxhi & any(is.finite(xhidraw))){
    finalsel=xarrow<xhidraw
    arrows(xarrow[finalsel],yarrow[finalsel],xhidraw[finalsel],yarrow[finalsel],angle=90,length=length,col=colarrow[finalsel],...)
  }
  if(doylo & any(is.finite(ylodraw))){
    finalsel=yarrow>ylodraw
    arrows(xarrow[finalsel],yarrow[finalsel],xarrow[finalsel],ylodraw[finalsel],angle=90,length=length,col=colarrow[finalsel],...)
  }
  if(doyhi & any(is.finite(yhidraw))){
    finalsel=yarrow<yhidraw
    arrows(xarrow[finalsel],yarrow[finalsel],xarrow[finalsel],yhidraw[finalsel],angle=90,length=length,col=colarrow[finalsel],...)
  }
  
}
if(!missing(corxy)){
  if(missing(xlo) | missing(ylo)){stop('For error ellipses xlo and ylo must be specified explicitly.')}
  if(any(is.finite(xlo)==FALSE) | any(is.finite(ylo)==FALSE)){stop('For error ellipses all xlo and ylo values must have real values.')}
  if(any(xlo!=xhi) | any(ylo!=yhi)){stop('xlo/ylo must equal xhi/yhi (i.e. the errors must be symmetric)')}
  n = length(x)
  for(i in 1:n){
    if(xlo[i]>0 & ylo[i]>0){
      Cov = matrix(c(xlo[i]^2,xlo[i]*ylo[i]*corxy[i],xlo[i]*ylo[i]*corxy[i],ylo[i]^2),2)
      E = eigen(Cov)
      a = sqrt(E$values[1])
      b = sqrt(E$values[2])
      angle = atan2(E$vector[2,1],E$vector[1,1])*180/pi
      if(fill){draw.ellipse(x[i],y[i],a=a,b=b,angle=angle,border=col[i],col=col[i],...)}else{draw.ellipse(x[i],y[i],a=a,b=b,angle=angle,border=col[i],col=NULL,...)}
    }
    if(xlo[i]>0 & ylo[i]==0){
      
    }
  }
}
}
