plot.bayesm.nmix=function(x,names,burnin=trunc(.1*nrow(probdraw)),Grid,bi.sel,nstd=2,marg=TRUE,
                          Data,ngrid=50,ndraw=200,...){
#
# S3 method to plot normal mixture marginal and bivariate densities
#     nmixlist is a list of 3 components, nmixlist[[1]]: array of mix comp prob draws,
#     mmixlist[[2]] is not used, nmixlist[[3]] is list of draws of components
#     P. Rossi 2/07
#     P. Rossi 3/07 fixed problem with dropping dimensions on probdraw (if ncomp=1)
#     P. Rossi 2/08 added marg flag to plot marginals
#     P. Rossi 3/08 added Data argument to paint histograms on the marginal plots
#
  nmixlist=x
  if(mode(nmixlist) != "list") stop(" Argument must be a list \n")
  probdraw=nmixlist[[1]]; compdraw=nmixlist[[3]]
  if(!is.matrix(probdraw)) stop(" First element of list (probdraw) must be a matrix \n")
  if(mode(compdraw) != "list") stop(" Third element of list (compdraw) must be a list \n")
  op=par(no.readonly=TRUE)
  on.exit(par(op))
  R=nrow(probdraw)
  if(R < 100) {cat(" fewer than 100 draws submitted \n"); return(invisible())}
  datad=length(compdraw[[1]][[1]]$mu)
  OneDimData=(datad==1)
  if(missing(bi.sel)) bi.sel=list(c(1,2))  # default to the first pair of variables
  ind=as.integer(seq(from=(burnin+1),to=R,length.out=max(ndraw,trunc(.05*R))))
  if(missing(names)) {names=as.character(1:datad)}
  if(!missing(Data)){
     if(!is.matrix(Data)) stop("Data argument must be a matrix \n")
     if(ncol(Data)!= datad) stop("Data matrix is of wrong dimension \n")      
  }
  if(mode(bi.sel) != "list") stop("bi.sel must be as list, e.g. bi.sel=list(c(1,2),c(3,4)) \n")

  if(missing(Grid)){
     Grid=matrix(0,nrow=ngrid,ncol=datad)
     if(!missing(Data))
	{for(i in 1:datad) Grid[,i]=c(seq(from=range(Data[,i])[1],to=range(Data[,i])[2],length=ngrid))}
     else
        {
	     out=momMix(probdraw[ind,,drop=FALSE],compdraw[ind])
         mu=out$mu
         sd=out$sd
         for(i in 1:datad ) Grid[,i]=c(seq(from=(mu[i]-nstd*sd[i]),
             to=(mu[i]+nstd*sd[i]),length=ngrid))
        }
  }
  #
  #  plot posterior mean of marginal densities
  #
  if(marg){
   mden=eMixMargDen(Grid,probdraw[ind,,drop=FALSE],compdraw[ind])
   nx=datad
   if(nx==1) par(mfrow=c(1,1)) 
   if(nx==2) par(mfrow=c(2,1))
   if(nx==3) par(mfrow=c(3,1))
   if(nx==4) par(mfrow=c(2,2))
   if(nx>=5) par(mfrow=c(3,2))

   for(index in 1:nx){
        if(index == 2) par(ask=dev.interactive())
        plot(range(Grid[,index]),c(0,1.1*max(mden[,index])),type="n",xlab="",ylab="density")
        title(names[index])
        if(!missing(Data)){
           deltax=(range(Grid[,index])[2]-range(Grid[,index])[1])/nrow(Grid)
           hist(Data[,index],xlim=range(Grid[,index]),
                freq=FALSE,col="yellow",breaks=max(20,.1*nrow(Data)),add=TRUE)
           lines(Grid[,index],mden[,index]/(sum(mden[,index])*deltax),col="red",lwd=2)}
        else
           {lines(Grid[,index],mden[,index],col="black",lwd=2)
           polygon(c(Grid[1,index],Grid[,index],Grid[nrow(Grid),index]),c(0,mden[,index],0),col="magenta")}
   }
  }
  #
  # now plot bivariates in list bi.sel
  #
  if(!OneDimData){
  par(ask=dev.interactive())
  nsel=length(bi.sel)
  den=array(0,dim=c(ngrid,ngrid,nsel))
  lstxixj=NULL
  for(sel in 1:nsel){
      i=bi.sel[[sel]][1]
      j=bi.sel[[sel]][2]
      xi=Grid[,i]
      xj=Grid[,j]
      lstxixj[[sel]]=list(xi,xj)
      for(elt in ind){
         den[,,sel]=den[,,sel]+mixDenBi(i,j,xi,xj,probdraw[elt,,drop=FALSE],compdraw[[elt]])
      }
      den[,,sel]=den[,,sel]/sum(den[,,sel])
   }     
  nx=nsel
  par(mfrow=c(1,1))

  for(index in 1:nx){
        xi=unlist(lstxixj[[index]][1])
        xj=unlist(lstxixj[[index]][2])
        xlabtxt=names[bi.sel[[index]][1]]
        ylabtxt=names[bi.sel[[index]][2]]
        image(xi,xj,den[,,index],col=terrain.colors(100),xlab=xlabtxt,ylab=ylabtxt)
        contour(xi,xj,den[,,index],add=TRUE,drawlabels=FALSE)
  }
  }

  invisible()
}

