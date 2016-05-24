"plot.spa" <-
function(x,xscale=NULL,g=NULL,vars=1:2,graph.grid=TRUE,plt.response,lcol=8,llwd=1,...){
  if(class(x)!="spa"){
    stop("Error:  `x' is not of type spa")
  }
  type=x$type
  y<-x$model$y
  lattice=FALSE
  control=x$control
  if(length(vars)>2){
    warning("To many varaibles set with vars ... only using first 2")
    vars<-vars[1:2]
  }  
  if(is.null(g)){
    stop("Error: must input a graph for the plot function")
  }
  if( x$control$dissimilar)
    g<-apply(is.finite(g),1,as.numeric)
  n<-x$model$dims[1]
  m<-x$model$dims[2]
  ##if(missing(xlim) & missing(ylim))
    plot(xscale[,vars],type="n",...)
  if(graph.grid){
    for(i in 1:n){
      ind<-which(g[,i]>0)
      piv1<-xscale[i,vars[1]]
      piv2<-xscale[i,vars[2]] 
      for(j in 1:length(ind)){  
        lines(c(piv1,xscale[ind[j],vars[1]]),c(piv2,xscale[ind[j],vars[2]]),col=lcol,lwd=llwd)
      }
    }
  }
  L=x$model$L
  y=as.numeric(x$model$y)
  points(xscale,cex=0.7,pch=16)
  if(nlevels(as.factor(y))==2 & plt.response){
    y=as.numeric(as.factor(y))
    points(xscale[c(L),vars],col=gray(c(0.2,0.8))[y],cex=2,pch=16)
    points(xscale[c(L),vars],col=1,cex=2,pch=1)
  }
  box()
}

