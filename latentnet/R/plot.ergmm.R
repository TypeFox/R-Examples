plot3d.ergmm<-function(x,...){
  plot.ergmm(x,rgl=TRUE,...)
}

plot.ergmm <- function(x, ..., vertex.cex=1, vertex.sides=16*ceiling(sqrt(vertex.cex)),
                       what="mkl",
                       main = NULL, xlab=NULL, ylab=NULL, zlab=NULL, xlim=NULL,ylim=NULL, zlim=NULL,
                       object.scale=formals(plot.network.default)[["object.scale"]],
                       pad=formals(plot.network.default)[["pad"]],
                       cluster.col=c("red","green","blue","cyan","magenta","orange","yellow","purple"),
                       vertex.col=NULL, print.formula=TRUE,
                       edge.col=8,
                       Z.ref=NULL,
                       Z.K.ref=NULL,
                       zoom.on=NULL,
                       pie = FALSE,
                       labels=FALSE,
                       rand.eff=NULL,
                       rand.eff.cap=NULL,
                       plot.means=TRUE,plot.vars=TRUE,
                       suppress.axes=FALSE,
                       jitter1D=1,curve1D=TRUE,
                       use.rgl=FALSE,
                       vertex.3d.cex=1/20,
                       edge.plot3d=TRUE,
                       suppress.center=FALSE,
                       density.par=list()){

  ## For convenience...
  Yg<-x[["model"]][["Yg"]]
  distances<-NULL
  n<-network.size(Yg)
  d<-x[["model"]][["d"]]
  if(d==0) stop("Plotting non-latent-space models is not available.")
  G<-x[["model"]][["G"]]
  if(G<1) pie<-FALSE

  if(use.rgl){
    if(!requireNamespace("rgl", quietly=TRUE)) stop("3D plots with use.rgl=TRUE option require the 'rgl' package.")
    if(pie) stop("3D plots cannot make pie charts.")
  }
  
  ## Set default axis labels.
  if(d==1){
    if (is.null(xlab)) 
      xlab <- ""
    if (is.null(ylab)) 
      ylab <- ""
  }else if(d==2){    
    if (is.null(xlab)) 
      xlab <- expression(Z[1])
    if (is.null(ylab)) 
      ylab <- expression(Z[2])
  }else if(d==3 && use.rgl){
    if (is.null(xlab)) 
      xlab <- expression(Z[1])
    if (is.null(ylab))
      ylab <- expression(Z[2])
    if (is.null(zlab)) 
      zlab <- expression(Z[3])
  }else if(d>3 && use.rgl){
    if (is.null(xlab)) 
      xlab <- "First principal component of Z"
    if (is.null(ylab)) 
      ylab <- "Second principal component of Z"
    if (is.null(zlab)) 
      zlab <- "Third principal component of Z"
  }else if(d>2){
    if (is.null(xlab)) 
      xlab <- "First principal component of Z"
    if (is.null(ylab)) 
      ylab <- "Second principal component of Z"
  }

  ## Find the requested plotting coordinates.
  ## Some "requests" require a substantially different code path, unfortunately.
  if(class(what)=="list"){
    summ<-what
    Z.pos <- summ[["Z"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    Z.pZK<-summ[["Z.pZK"]]
    if (is.null(main)) 
      main <- paste(deparse(substitute(what))," Latent Positions of ", 
                    deparse(substitute(x)),sep="")

  }else if(what=="start" || what=="burnin.start"){
    summ<-x[["start"]]
    Z.pos <- summ[["Z"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Initial Latent Positions of ", 
                    deparse(substitute(x)),sep="")

  }else if(what=="sampling.start"){
    summ<-x[["sampling.start"]]
    Z.pos <- summ[["Z"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Post-Burn-In Latent Positions of ", 
                    deparse(substitute(x)),sep="")

  }else if(what=="mle"){
    summ<-summary(x,point.est=c("mle"),se=FALSE, bic.eff.obs=NULL)
    Z.pos <- summ[["mle"]][["Z"]]
    summ<-summ[["mle"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    plot.means<-plot.vars<-FALSE
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Multistage MLEs of Latent Positions of", 
                    deparse(substitute(x)))

  }else if(what=="pmean"){
    summ<-summary(x,point.est=c("pmean"), bic.eff.obs=NULL)
    Z.pos <- summ[["pmean"]][["Z"]]
    summ<-summ[["pmean"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    Z.pZK<-summ[["Z.pZK"]]
    if (is.null(main)) 
      main <- paste("Posterior Mean Positions of", 
                    deparse(substitute(x)))

  }else if(what=="mkl"){
    summ<-summary(x,point.est=c("pmean","mkl"), bic.eff.obs=NULL)
    Z.pos <- summ[["mkl"]][["Z"]]
    if(!is.null(x[["mkl"]][["mbc"]])){
      Z.mean<-summ[["mkl"]][["mbc"]][["Z.mean"]]
      Z.var<-summ[["mkl"]][["mbc"]][["Z.var"]]
    }else{
      if(!is.null(summ[["pmean"]][["Z.mean"]])) Z.mean<-summ[["pmean"]][["Z.mean"]]
      else plot.means<-FALSE
      Z.var<-summ[["pmean"]][["Z.var"]]
    }
    Z.K<-summ[["pmean"]][["Z.K"]]
    Z.pZK<-summ[["pmean"]][["Z.pZK"]]
    summ<-summ[["mkl"]]
    if (is.null(main)) 
      main <- paste("MKL Latent Positions of", 
                    deparse(substitute(x)))

  }else if(what=="pmode"){
    summ<-summary(x,point.est=c("pmode"), bic.eff.obs=NULL)
    Z.pos <- summ[["pmode"]][["Z"]]
    summ<-summ[["pmode"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Posterior Mode Latent Positions of", 
                    deparse(substitute(x)))

  }else if(what=="cloud"){
    summ<-summary(x,point.est=c("pmean","mkl"), bic.eff.obs=NULL)
    Z.pos <- summ[["mkl"]][["Z"]]
    if(d!=2) stop("Cloud plots are only available for 2D latent space models.")
    
    if(!is.null(x[["mkl"]][["mbc"]])){
      Z.mean<-summ[["mkl"]][["mbc"]][["Z.mean"]]
      Z.var<-summ[["mkl"]][["mbc"]][["Z.var"]]
    }else{
      if(!is.null(summ[["pmean"]][["Z.mean"]])) Z.mean<-summ[["pmean"]][["Z.mean"]]
      else plot.means<-FALSE
      Z.var<-summ[["pmean"]][["Z.var"]]
    }
    Z.K<-summ[["pmean"]][["Z.K"]]
    Z.pZK<-summ[["pmean"]][["Z.pZK"]]
    summ<-summ[["mkl"]]
    if (is.null(main)) 
      main <- paste("MKL Latent Positions of", 
                    deparse(substitute(x)))
    plot(matrix(c(x[["sample"]][["Z"]]),ncol=2),pch=".")
    points(Z.pos,col=cluster.col[Z.K])
    points(Z.mean,col=cluster.col)
    
    return(invisible(NULL))

  }else if(what=="density"){
    if(is.null(density.par[["totaldens"]])) density.par[["totaldens"]] <- TRUE
    if(is.null(density.par[["subdens"]])) density.par[["subdens"]] <- TRUE
    if(is.null(density.par[["mfrow"]])){
      wanted<-density.par[["totaldens"]]+density.par[["subdens"]]*G
      density.par[["mfrow"]]<-rep(min(ceiling(sqrt(wanted)),4),2)
    }
    
    summ<-summary(x,point.est=c("pmean","mkl"), bic.eff.obs=NULL)
    Z.pos <- summ[["mkl"]][["Z"]]
    if(d!=2) stop("Density plots are only available for 2D latent space models.")

    if(!requireNamespace("KernSmooth",quietly=TRUE)){
      stop("The 'density' option requires the 'KernSmooth' package.")
    }

    old.par<-par(mfrow=density.par[["mfrow"]],mar=c(2.5,2.5,1,1))

    Z.all<-matrix(c(aperm(x[["sample"]][["Z"]],c(2,1,3))),ncol=2)

    if(density.par[["totaldens"]]){
      plot(Z.all,type='n',xlab=xlab,ylab=ylab,...)
      title(main=paste("Posterior density of",deparse(substitute(x))), cex.main=0.7, ...)
      Z.bkde <- KernSmooth::bkde2D(Z.all,0.2,c(201,201))
      image(Z.bkde[["x1"]],Z.bkde[["x2"]],Z.bkde[["fhat"]],col=grey(seq(1,0,length=255)),add=TRUE)
      box()
    }
    
    if(G>1 && density.par[["subdens"]]){
      Z.K.all <- c(t(x[["sample"]][["Z.K"]]))
      for(i in 1:G){
        plot(Z.all,main=paste("Class",i),type="n",...)
        Z.bkde <- KernSmooth::bkde2D(Z.all[Z.K.all==i,],0.2,c(101,101))
        col<-c(col2rgb(cluster.col[i])/255)
        image(Z.bkde[["x1"]],Z.bkde[["x2"]],Z.bkde[["fhat"]],add=TRUE,
              col=rgb(seq(1,col[1],length=255),
                seq(1,col[2],length=255),
                seq(1,col[3],length=255)))
        contour(Z.bkde[["x1"]],Z.bkde[["x2"]],Z.bkde[["fhat"]],add=TRUE, nlevels=4,
                drawlabels=FALSE,
                col="white")
        box()
      }
    }
    
    par(old.par)

    return(invisible(NULL))
    
  }else if(is.numeric(what) && round(what)==what){
    summ<-x[["sample"]][[what]]
    Z.pos <- summ[["Z"]]
    Z.mean<-summ[["Z.mean"]]
    Z.var<-summ[["Z.var"]]
    Z.K<-summ[["Z.K"]]
    if(pie) stop("Cannot make pie charts with the specified parameter type.")
    if (is.null(main)) 
      main <- paste("Iteration #",what," Latent Positions of ", 
                    deparse(substitute(x)),sep="")
  }else stop("Invalid latent space position estimate type.")

  if(!is.null(Z.ref)){
    R<-.procr(Z.pos,Z.ref,scale=FALSE,reflect=TRUE)
    Z.pos<-Z.pos%*%R
    if(G) Z.mean<-Z.mean%*%R
  }
  if(!is.null(Z.K.ref)){
    perm<-which.perm.nearest(Z.K.ref,Z.K)
    Z.K<-order(perm)[Z.K]
    Z.pZK<-Z.pZK[,perm]
    Z.mean<-Z.mean[perm,]
    Z.var<-Z.var[perm]
  }
  
  ## Transform coordinates for dimensions other than 2 (or 3, when using rgl).
  if(d==1){    
    Z.pos<-coords.1D(Z.pos,curve1D,jitter1D)
    if(G) Z.mean<-coords.1D(Z.mean,curve1D,jitter1D)
    if(curve1D){
      distances<-as.matrix(dist(Z.pos))
      distances<-distances/max(distances)
    }
  }else if(d>3 && use.rgl){
    ## Plot the first three principal components.
    prc<-prcomp(Z.pos)
    Z.pos<-predict(prc,Z.pos)[,1:3]
    if(G) Z.mean<-predict(prc,Z.mean)[,1:3]
  }else if(d>2 && !use.rgl){ ## I.e. high latent dimension.
    ## Plot the first two principal components.
    prc<-prcomp(Z.pos)
    Z.pos<-predict(prc,Z.pos)[,1:2]
    if(G) Z.mean<-predict(prc,Z.mean)[,1:2]
  }
  
  ## Set default vertex color.
  if(is.null(vertex.col)){
    if(is.latent.cluster(x) && !pie)
      vertex.col <- cluster.col[Z.K]
    else vertex.col<-cluster.col[1]
  }
  else if(length(vertex.col)==1 && is.character(vertex.col)){
    trycol <- as.numeric(as.factor(unlist(get.vertex.attribute(Yg,vertex.col))))
    if(!all(is.na(trycol))){
      vertex.col <- cluster.col[trycol]
    }
  }
  vertex.col<-rep(vertex.col,length.out=n)

  ## Set vertex sizes to correspond to random effect values.
  if(!is.null(rand.eff) && (rand.eff[1]=="total" || x[["model"]][rand.eff[1]][[1]])){
    if(rand.eff=="total")
      rand.eff.mul<-exp((summ["sender"][[1]]+summ["receiver"][[1]]))
    else      
      rand.eff.mul<-exp(summ[rand.eff][[1]])
    if(!is.null(rand.eff.cap)) rand.eff.mul<-pmax(pmin(rand.eff.mul,exp(rand.eff.cap)),exp(-rand.eff.cap))
    rand.eff.mul<-sqrt(rand.eff.mul)
    vertex.cex<-vertex.cex*rand.eff.mul
  }
  vertex.cex<-rep(vertex.cex,length.out=n)

  ## Find the bounds of the plotting region.
  xylim<-ergmm.plotting.region(if(!is.null(zoom.on)) Z.pos[zoom.on,] else Z.pos,
                                  if(plot.means && is.null(zoom.on)) Z.mean,if(plot.vars && is.null(zoom.on)) Z.var,!suppress.center && is.null(zoom.on),pad)
  if(is.null(xlim)) xlim<-xylim[["xlim"]] else xylim[["xlim"]]<-xlim
  if(is.null(ylim)) ylim<-xylim[["ylim"]] else xylim[["ylim"]]<-ylim
  if(is.null(zlim)) zlim<-xylim[["zlim"]] else xylim[["zlim"]]<-zlim

  if(!use.rgl){
    ## Go forth, and plot the network.
    old.warn<-options()[["warn"]]
    options(warn=-1)
    
    plot.network(Yg,coord=Z.pos,
                 main=main,
                 vertex.col=if(is.null(vertex.col)) 0 else vertex.col,
                 vertex.cex=vertex.cex,
                 xlab=xlab,
                 ylab=ylab,
                 xlim=xlim,
                 ylim=ylim,
                 object.scale=object.scale,
                 pad=pad,
                 suppress.axes=suppress.axes,
                 vertex.sides=16*sqrt(vertex.cex),
                 jitter=FALSE,
                 usecurve=curve1D,
                 edge.curve=distances,
                 displaylabels=labels&&!(x[["model"]][["d"]]==1 && curve1D==TRUE),
                 tick=!(x[["model"]][["d"]]==1 && curve1D==TRUE),
                 edge.col=edge.col,
                 ...)
  
    options(warn=old.warn)
  }else{
    vertex.radii <- vertex.cex*vertex.3d.cex
    rgl::plot3d(Z.pos,type="s",col= vertex.col,radius=vertex.radii,xlab=xlab,ylab=ylab,zlab=zlab,xlim=xlim,ylim=ylim,zlim=zlim,alpha=1,main=main)
    if(labels){
      Z.pos.r <- sqrt(rowSums(Z.pos^2))
      rgl::text3d(Z.pos*(Z.pos.r+vertex.radii*2)/Z.pos.r,texts=Yg %v% "vertex.names")
    }
    if(edge.plot3d){
      el<-as.matrix(Yg,matrix.type="edgelist")
      if(is.directed(Yg) && requireNamespace("heplots",quietly=TRUE)){
        medlen <- median(apply(el, 1, function(e){
          Z1 <- Z.pos[e[1],]
          Z2 <- Z.pos[e[2],]
          dZ <- Z2-Z1
          r1 <- vertex.radii[e[1]]
          r2 <- vertex.radii[e[2]]
          len <- sqrt(sum(dZ^2))
          len - r1 - r2
        }))

        apply(el, 1, function(e){
          Z1 <- Z.pos[e[1],]
          Z2 <- Z.pos[e[2],]
          dZ <- Z2-Z1
          r1 <- vertex.radii[e[1]]
          r2 <- vertex.radii[e[2]]
          len <- sqrt(sum(dZ^2))
          Z1edge <- Z1+r1/len*dZ
          Z2edge <- Z2-r2/len*dZ
          heplots::arrow3d(Z1edge,Z2edge,barblen=medlen*.1,n=4)
        })
      }else{
        if(is.directed(Yg)) message("3D plotting of directed edges with arrows requires package heplots, which is not installed. Using line segments instead.")
        rgl::segments3d(Z.pos[c(t(el)),])
      }
    }
  }

  if(!use.rgl){
    ## For 1D plots, only plot horizontal axis ticks.
    if(x[["model"]][["d"]]==1 && curve1D==TRUE) {
      axis(1)
    }
    
    ## Add the model formula as a subtitle.
    if(print.formula){
      fmla <- x[["model"]][["formula"]]   
      xformula <- paste(fmla[2],fmla[1],fmla[-c(1:2)],collapse=" ")
      if(!is.null(x[["model"]][["response"]])) xformula<-paste(xformula,"   (",x[["model"]][["response"]],")",sep="")
      title(main = xformula, line = 1, cex.main = 0.7)
  }
    
    ## Plot pie charts.
    if(pie){
      piesize<-rep(ergmm.plotting.vertex.radius(vertex.cex,xylim,object.scale),length=n)
      pie.order<-order(piesize,decreasing=TRUE)
      for(i in 1:n){
        ergmm.drawpie(Z.pos[pie.order[i],],piesize[pie.order[i]],Z.pZK[pie.order[i],],n=50,cols=cluster.col)
      }
    }

    ## Mark the events with "bullets" (small circles).
    if(is.bipartite(Yg)){
      bip<-Yg %n% "bipartite"
      points(Z.pos[-(1:bip),],pch=20,cex=vertex.cex[-(1:bip)],col=1)
    }
  }
  
  ## Mark the center.
  if(!suppress.center){
    if(!use.rgl)
      points(cbind(0,0),pch="+")
    else
      THra3d(0,0,0,color="black",size=1*vertex.3d.cex)
  }
  Z.mean<-if(G>0)Z.mean else cbind(0,0,if(use.rgl)0)

  ## Mark the cluster means.
  if(plot.means){
    if(!use.rgl)
      points(Z.mean,pch="+",col=cluster.col)
    else
      THra3d(Z.mean,color=cluster.col,size=1*vertex.3d.cex)
  }

  ## Plot the cluster standard deviations.
  if(plot.vars){
    if(!use.rgl)
      symbols(Z.mean,circles=sqrt(Z.var),fg=cluster.col,inches=FALSE,add=TRUE,asp=1)
    else
      for(c in seq_len(G)){
        rgl::spheres3d(Z.mean[c,],radius=sqrt(Z.var[c]),color=cluster.col[c],alpha=0.2)
      }
  }
  
  invisible(Z.pos)
}

coords.1D<-function(Z,curve1D,jitter1D){
  if(curve1D){
    Z.pos<-cbind(Z,rep(0,length(Z)))
  }else{
    jitter<-rnorm(length(Z),0,jitter1D*sd(Z)/sqrt(length(Z)))
    Z.pos<-cbind(Z,Z)+cbind(jitter,-jitter)
  }
  Z.pos
}

ergmm.plotting.region<-function(Z.pos,Z.mean,Z.var,include.center,pad){
  d3=(dim(Z.pos)[2]>=3)
  xylim<-t(apply(rbind(Z.pos+pad,Z.pos-pad,
                       Z.mean+pad,Z.mean-pad,
                       if(!is.null(Z.mean) && !is.null(Z.var)) Z.mean+sqrt(Z.var)+pad,
                       if(!is.null(Z.mean) && !is.null(Z.var)) Z.mean-sqrt(Z.var)-pad,
                       if(include.center) 0+pad,
                       if(include.center) 0-pad,
                       if(is.null(Z.mean) && length(Z.var)==1) matrix(sqrt(Z.var)+pad,ncol=2,byrow=FALSE),
                       if(is.null(Z.mean) && length(Z.var)==1) matrix(-sqrt(Z.var)-pad,ncol=2,byrow=FALSE)
                       ),2,range))
  xlim<-xylim[1,]
  ylim<-xylim[2,]
  if(d3)zlim<-xylim[3,]

  ### Shamelessly stolen from Carter's plot.network.default code (sans z).
  xrng<-diff(xlim)          #Force scale to be symmetric
  yrng<-diff(ylim)
  if(d3)zrng<-diff(zlim)
  
  xctr<-(xlim[2]+xlim[1])/2 #Get center of plotting region
  yctr<-(ylim[2]+ylim[1])/2
  if(d3)zctr<-(zlim[2]+zlim[1])/2
  
  if(!d3){
    if(xrng<yrng)
      xlim<-c(xctr-yrng/2,xctr+yrng/2)
    else
      ylim<-c(yctr-xrng/2,yctr+xrng/2)
  }
  ### End stolen part
  
  list(xlim=xlim,ylim=ylim,zlim=if(d3)zlim)
}

ergmm.plotting.vertex.radius<-function(vertex.cex,xylim,object.scale){
  baserad<-min(diff(xylim[["xlim"]]),diff(xylim[["ylim"]]))*object.scale
  vertex.cex*baserad
}

THra3d<-function(x,y=NULL,z=NULL,size=1,color="white",alpha=1,...){
  xyz<-xyz.coords(x=x,y=y,z=z)
  n<-length(xyz[["x"]])
  size<-rep(size,length.out=n)
  color<-rep(color,length.out=n)
  alpha<-rep(alpha,length.out=n)

  for(i in 1:n){
    with(xyz,THron3d(x[i],y[i],z[i],size=size[i],color=color[i],alpha=alpha[i]))
  }
}

THron3d<-function(x,y=NULL,z=NULL,size=1,...){
  xyz<-xyz.coords(x=x,y=y,z=z)
  xyz<-with(xyz,c(x,y,z))/size
  centroid<-c(1/2,1/(2*sqrt(2)),sqrt(3)/4)
  
  v.left<-c(0,0,0)+xyz-centroid
  v.right<-c(1,0,0)+xyz-centroid
  v.far<-c(1/2,0,sqrt(3)/2)+xyz-centroid
  v.top<-c(1/2,1/sqrt(2),sqrt(3)/4)+xyz-centroid
  rgl::triangles3d(rbind(v.left,v.right,v.top,
                    v.left,v.top,v.far,
                    v.left,v.right,v.far,
                    v.top,v.right,v.far)*size,y=NULL,z=NULL,...)
}
