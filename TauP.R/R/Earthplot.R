Earthplot <-
function(model,prop='vp',image.col=heat.colors(500),n=200,add=FALSE,...)
  {
    if(!is.null(model$rp)){
      rp=model$rp
    }else{
      rp=6371
    }
    if(!add){
      par(mar = c(1.1,1.1,4.1,1.1))
      plot(0,type='n',xlim = 1.15 * c(-rp,rp),ann=FALSE,axes=FALSE,asp=1)
    }
    
    if(is.character(prop)){
      x = seq(from = -rp, to = rp, length.out = n)
      y = seq(from = -rp, to = rp, length.out = n)
      L = meshgrid(x,y)
      H = rp - sqrt(L[[1]]^2 + L[[2]]^2)
      options(warn = -1)
      V = approx(model$z, model[[prop]], H)$y
      options(warn = 0)
      image(x, y, matrix(V,n,n), col = image.col, add = TRUE)
    }
    PolarPlot(0:360,rp,type='l',method=lines,degrees=TRUE,geographical=TRUE,col='black',...)
    polaraxis(rp)
    if(!is.null(model$conr)){
      PolarPlot(0:360,rp-model$conr,method=lines,degrees=TRUE,geographical=TRUE,lwd=.25,col='black',...)
    }
    if(!is.null(model$moho)){
      PolarPlot(0:360,rp-model$moho,method=lines,degrees=TRUE,geographical=TRUE,lwd=.25,col='black',...)
    }
    if(!is.null(model$d410)){
      PolarPlot(0:360,rp-model$d410,method=lines,degrees=TRUE,geographical=TRUE,col='black',...)
    }
    if(!is.null(model$d520)){
      PolarPlot(0:360,rp-model$d520,method=lines,degrees=TRUE,geographical=TRUE,col='black',...)
    }
    if(!is.null(model$d660)){
      PolarPlot(0:360,rp-model$d660,method=lines,degrees=TRUE,geographical=TRUE,col='black',...)
    }
    if(!is.null(model$cmb)){
      PolarPlot(0:360,rp-model$cmb,method=lines,degrees=TRUE,geographical=TRUE,col='black',...)
    }
    if(!is.null(model$icb)){
      PolarPlot(0:360,rp-model$icb,method=lines,degrees=TRUE,geographical=TRUE,col='black',...)
    }
  }

