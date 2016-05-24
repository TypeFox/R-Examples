vis.scam <- function(x,view=NULL,cond=list(),n.grid=30,too.far=0,col=NA,color="heat",
           contour.col=NULL,se=-1,type="link",plot.type="persp",zlim=NULL,nCol=50,...)
## hacked version of vis.gam() from mgcv package, (c) Simon N. Wood 23/2/03 
## predict.gam is simply replaced by predict.scam
# takes a scam object and plots 2D views of it, supply ticktype="detailed" to get proper axis anotation
{ fac.seq<-function(fac,n.grid)
  # generates a sequence of factor variables of length n.grid
  { fn<-length(levels(fac));gn<-n.grid;
    if (fn>gn) mf<-factor(levels(fac))[1:gn]
    else
    { ln<-floor(gn/fn) # length of runs               
      mf<-rep(levels(fac)[fn],gn)
      mf[1:(ln*fn)]<-rep(levels(fac),rep(ln,fn))
      mf<-factor(mf,levels=levels(fac))
    }
    mf
  }
  # end of local functions

  dnm <- names(list(...))

  ## basic issues in the following are that not all objects will have a useful `data'
  ## component, but they all have a `model' frame. Furthermore, `predict.gam' recognises
  ## when a model frame has been supplied

  v.names  <- names(x$var.summary) ## names of all variables

  ## Note that in what follows matrices in the parametric part of the model
  ## require special handling. Matrices arguments to smooths are different
  ## as they follow the summation convention. 
  if (is.null(view)) # get default view if none supplied
  { ## need to find first terms that can be plotted against
    k <- 0;view <- rep("",2) 
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) ok <- FALSE else
      if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]]))<=1) ok <- FALSE 
      } else {
        if (length(unique(x$var.summary[[i]]))==1) ok <- FALSE
      }
      if (ok) {
        k <- k + 1;view[k] <- v.names[i]
      }
      if (k==2) break;
    }
    if (k<2) stop("Model does not seem to have enough terms to do anything useful")
  } else { 
    if (sum(view%in%v.names)!=2) stop(
        paste(c("view variables must be one of",v.names),collapse=", "))
    for (i in 1:2) 
    if  (!inherits(x$var.summary[[view[i]]],c("numeric","factor"))) 
    stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }

  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]]))<=1) ok <- FALSE
  } else {
    if (length(unique(x$var.summary[[view[i]]]))<=1) ok <- FALSE 
  }
  if (!ok) 
  stop(paste("View variables must contain more than one value. view = c(",view[1],",",view[2],").",sep=""))

 
  # Make dataframe....
  if (is.factor(x$var.summary[[view[1]]]))
  m1<-fac.seq(x$var.summary[[view[1]]],n.grid)
  else { r1<-range(x$var.summary[[view[1]]]);m1<-seq(r1[1],r1[2],length=n.grid)}
  if (is.factor(x$var.summary[[view[2]]]))
  m2<-fac.seq(x$var.summary[[view[2]]],n.grid)
  else { r2<-range(x$var.summary[[view[2]]]);m2<-seq(r2[1],r2[2],length=n.grid)}
  v1<-rep(m1,n.grid);v2<-rep(m2,rep(n.grid,n.grid))
  
  newd <- data.frame(matrix(0,n.grid*n.grid,0)) ## creating prediction data frame full of conditioning values
  for (i in 1:length(x$var.summary)) { 
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) { 
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) ma <- ma[2] ## extract median
    }
    if (is.matrix(x$var.summary[[i]])) newd[[i]] <- matrix(ma,n.grid*n.grid,ncol(x$var.summary[[i]]),byrow=TRUE)
    else newd[[i]]<-rep(ma,n.grid*n.grid)
  }
  names(newd) <- v.names
  #row.names <- attr(newd,"row.names")
  #attributes(newd) <- attributes(x$model) # done so that handling of offsets etc. works
  #attr(newd,"row.names") <- row.names
  newd[[view[1]]]<-v1
  newd[[view[2]]]<-v2
  # call predict.gam to get predictions.....
  if (type=="link") zlab<-paste("linear predictor")
  else if (type=="response") zlab<-type
  else stop("type must be \"link\" or \"response\"")
  ## turn newd into a model frame, so that names and averages are valid
  #attributes(newd)<-attributes(x$model)
  #attr(newd,"row.names")<-as.character(1:(n.grid*n.grid))
  fv <- predict.scam(x,newdata=newd,se.fit=TRUE,type=type)
  z <- fv$fit # store NA free copy now
  if (too.far>0) # exclude predictions too far from data
  { ex.tf <- exclude.too.far(v1,v2,x$model[,view[1]],x$model[,view[2]],dist=too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf]<-NA
  }
  # produce a continuous scale in place of any factors
  if (is.factor(m1)) 
  { m1<-as.numeric(m1);m1<-seq(min(m1)-0.5,max(m1)+0.5,length=n.grid) }
  if (is.factor(m2)) 
  { m2<-as.numeric(m2);m2<-seq(min(m1)-0.5,max(m2)+0.5,length=n.grid) }
  if (se<=0)
  { old.warn<-options(warn=-1)
    av<-matrix(c(0.5,0.5,rep(0,n.grid-1)),n.grid,n.grid-1)
    options(old.warn)
    # z is without any exclusion of gridpoints, so that averaging works nicely
    max.z <- max(z,na.rm=TRUE)
    z[is.na(z)] <- max.z*10000 # make sure NA's don't mess it up
    z<-matrix(z,n.grid,n.grid) # convert to matrix
    surf.col<-t(av)%*%z%*%av   # average over tiles  
    surf.col[surf.col>max.z*2] <- NA # restore NA's
    # use only non-NA data to set colour limits
    if (!is.null(zlim))
    { if (length(zlim)!=2||zlim[1]>=zlim[2]) stop("Something wrong with zlim")
      min.z<-zlim[1]
      max.z<-zlim[2]
    } else
    { min.z<-min(fv$fit,na.rm=TRUE)
      max.z<-max(fv$fit,na.rm=TRUE)
    }
    surf.col<-surf.col-min.z
    surf.col<-surf.col/(max.z-min.z)  
    surf.col<-round(surf.col*nCol)
    con.col <-1
    if (color=="heat") { pal<-heat.colors(nCol);con.col<-3;}
    else if (color=="topo") { pal<-topo.colors(nCol);con.col<-2;}
    else if (color=="cm") { pal<-cm.colors(nCol);con.col<-1;}
    else if (color=="terrain") { pal<-terrain.colors(nCol);con.col<-2;}
    else if (color=="gray"||color=="bw") {pal <- gray(seq(0.1,0.9,length=nCol));con.col<-1}
    else stop("color scheme not recognised")
    if (is.null(contour.col)) contour.col<-con.col   # default colour scheme
    surf.col[surf.col<1]<-1;surf.col[surf.col>nCol]<-nCol # otherwise NA tiles can get e.g. -ve index
    if (is.na(col)) col<-pal[as.array(surf.col)]
    z<-matrix(fv$fit,n.grid,n.grid)
    if (plot.type=="contour")
    { stub <- paste(ifelse("xlab" %in% dnm, "" , ",xlab=view[1]"),
                    ifelse("ylab" %in% dnm, "" , ",ylab=view[2]"),
                    ifelse("main" %in% dnm, "" , ",main=zlab"),",...)",sep="")
      if (color!="bw")
      { txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)",stub,sep="") # assemble image() call
        eval(parse(text=txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)",
               ifelse("add" %in% dnm, "" , ",add=TRUE"),",...)" , sep="") # assemble contour() call
         eval(parse(text=txt))       
      } else
      { txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)",stub,sep="")  # assemble contour() call
        eval(parse(text=txt))
      }
    } else
    { stub <- paste(ifelse("xlab" %in% dnm, "" , ",xlab=view[1]"),
                    ifelse("ylab" %in% dnm, "" , ",ylab=view[2]"),
                    ifelse("main" %in% dnm, "" , ",zlab=zlab"),",...)",sep="")
      if (color=="bw")
      { op <- par(bg="white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ",stub,sep="") # assemble persp() call
        eval(parse(text=txt))
        par(op)
      } else
      { txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)",stub,sep="")  # assemble persp() call
        eval(parse(text=txt))
      }
    }
  } else # add standard error surfaces
  { if (color=="bw"||color=="gray") 
    { subs <- paste("grey are +/-",se,"s.e.") 
      lo.col <- "gray" ## ignore codetools claims about this
      hi.col <- "gray" ## ignore codetools 
    } else
    { subs<-paste("red/green are +/-",se,"s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim))
    { if (length(zlim)!=2||zlim[1]>=zlim[2]) stop("Something wrong with zlim")
      min.z<-zlim[1]
      max.z<-zlim[2]
    } else
    { z.max<-max(fv$fit+fv$se.fit*se,na.rm=TRUE)
      z.min<-min(fv$fit-fv$se.fit*se,na.rm=TRUE)
    }
    zlim<-c(z.min,z.max)
    z<-fv$fit-fv$se.fit*se;z<-matrix(z,n.grid,n.grid)
    if (plot.type=="contour") warning("sorry no option for contouring with errors: try plot.scam")

    stub <-  paste(ifelse("xlab" %in% dnm, "" , ",xlab=view[1]"),
                   ifelse("ylab" %in% dnm, "" , ",ylab=view[2]"),
                   ifelse("zlab" %in% dnm, "" , ",zlab=zlab"),
                   ifelse("sub" %in% dnm, "" , ",sub=subs"),
                   ",...)",sep="")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=lo.col"),
                 stub,sep="") # assemble persp() call
    eval(parse(text=txt))

    par(new=TRUE) # don't clean device
    z<-fv$fit;z<-matrix(z,n.grid,n.grid)

    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=\"black\""),
                 stub,sep="")
    eval(parse(text=txt))

    par(new=TRUE) # don't clean device
    z<-fv$fit+se*fv$se.fit;z<-matrix(z,n.grid,n.grid)
    
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim",
                 ifelse("border" %in% dnm, "" ,",border=hi.col"),
                 stub,sep="")
    eval(parse(text=txt))

  }
}
