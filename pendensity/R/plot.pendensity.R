#plotting the estimated density/densities, using the estimated density. Here, it's called 'obj'.
plot.pendensity <- function(x,plot.val=1,val=NULL,latt=FALSE,kernel=FALSE,confi=TRUE,main=NULL,sub=NULL,xlab=NULL,ylab=NULL,
               plot.base=FALSE,lwd=NULL,legend.txt=NULL,plot.dens=TRUE,...) {
  obj <- x
  if(plot.val==1) {
    weight <- obj$results$ck
    q <- obj$splines$q
    x.factor <- obj$values$covariate$x.factor
    help.env <- new.env()
    if(!is.null(x.factor)) {
      len.x.fac <- length(x.factor[1,])
      all.x <- 0
      all.x2 <- len.x.fac
    }
    max.bsp <-0
    max.kern <-0
    base <- obj$splines$base
    base.den <- obj$splines$base.den
    MeanW <- obj$splines$MeanW
    knots.spline <- obj$splines$knots.val$all
    knots.val <- obj$splines$knots.val
    Stand.abw <- obj$splines$Stand.abw
    help.degree <- obj$splines$help.degree
    K <- obj$splines$K
    N <- obj$splines$N
    Z <- obj$values$Z
    x <- obj$values$x
    m <- obj$splines$m
    Dm <- obj$splines$Dm
    h <- obj$splines$h
    sort <- obj$values$sort
    my.AIC <- obj$results$AIC$my.AIC
    lambda0 <- obj$results$lambda0
    var.par <- obj$results$variance.par
    eps <- 1e-4
    levels <- obj$values$covariate$levels
    how.combi <- obj$values$covariate$how.combi
    how.levels <- obj$values$covariate$how.levels
    
    lev1 <- c()
    for(i in 1:how.levels) lev1 <- c(lev1,as.numeric(as.vector(levels[[i]])))
    lev <- lev1[-1]
    if(base=="bspline") {
      h.help <- abs(knots.spline[1]-knots.spline[2])
      m <- length(knots.spline)
    }
    if(base=="gaussian") {
      h.help <- abs(MeanW[1]-MeanW[2])
      m <- length(MeanW)
    }

    max.bw.all <- c()
    
    if (is.null(obj$values$x)) {
      if(!is.null(val)) {
        y.val <- val
        cond <- (y.val>=min(knots.val$val) & y.val <=max(knots.val$val))
        if(any(cond==FALSE)) stop("calculation of density is only possible in the range of the response",call.=FALSE)
      }
      y <- obj$values$y
      if(!sort) {
        o <- order(y)
        base.den <- base.den[,o]
        y <- y[o]
      }
      y.r <- range(knots.val$val)
      y.help <- seq(y.r[1],y.r[2],length=200)
      list <- c()
      
      if(base=="bspline") {
        if(q>2) K.help <- K-q+2 else K.help <- 0
        help.base.den <- my.bspline(h,q,knots.val,y.help,K.help,plot.bsp=FALSE)$base.den
        if(!is.null(val)) {
          help.base.y.val <- my.bspline(h,q,knots.val,y.val,K.help,plot.bsp=FALSE)$base.den
          assign("base1.y.val",help.base.y.val,help.env)
          assign("y.help1.val",y.val,help.env)
        }
        assign(paste("base1"),help.base.den,envir=help.env)
        assign("y.help1",y.help,envir=help.env)
      
      }
      if(base=="gaussian") {
        for(i in 1:m) {
          if(MeanW[i] >=y.r[1] & MeanW[i]<=y.r[2]) list <- c(list,i)
        }
        nn <- matrix(1:length(y.help))
        help.base.den <- apply(nn,1,function(i,MeanW,Stand.abw,y.help) dnorm(y.help[i],MeanW,Stand.abw),MeanW,Stand.abw,y.help)[list,]
        assign(paste("y.help1",sep=""),y.help,envir=help.env)
        assign(paste("base1"),help.base.den,envir=help.env)
        if(!is.null(val)) {
          nn <- matrix(1:length(y.val))
          help.base.y.val <- apply(nn,1,function(i,MeanW,Stand.abw,y.val) dnorm(y.val[i],MeanW,Stand.abw),MeanW,Stand.abw,y.val)[list,]
          assign("y.help1.val",y.val,help.env)
          assign("base1.y.val",help.base.y.val,help.env)
        }
      }   
      sd.cal <- new.env()
      sd.cal1 <- variance.val(help.base.den,var.par,weight,K,x,Z,x.factor) #confidence intervals
      assign("sd.cal1",sd.cal1,sd.cal)
      if(!is.null(val)) {
        sd.cal1.y.val <- variance.val(help.base.y.val,var.par,weight,K,x,Z,x.factor) #confidence intervals
        assign("sd.cal1.y.val",sd.cal1.y.val,sd.cal)
      }
      
      assign("weight1",weight,envir=help.env)
      assign("y.later1",y,envir=help.env)
      assign("len.x.fac1",1,envir=help.env)
      assign("b.w1",bw <- get("weight1",envir=help.env)%*%get("base1",envir=help.env),help.env)
      if(!is.null(val)) {
        assign("b.w1.y.val",den.val <- get("weight1",envir=help.env)%*%get("base1.y.val",envir=help.env),help.env)
      }
      
      max.bw <- max(bw)
      cut <- 0.0005*max.bw

      y.help <- get("y.help1",help.env)
      ind <- which(bw<cut)
      assign("ind1",ind,help.env)

      if(length(ind)!=0){
        assign("b.w1",bw[-ind],envir=help.env)
        assign("y.help1",get("y.help1",help.env)[-ind],help.env)
      }
      max.bsp <- max.bw
      if(kernel) assign("kern.den1",density(get("y.later1",envir=help.env),kernel="gaussian",bw="ucv"),envir=help.env)
      if(kernel) max.kern <- max(get("kern.den1",envir=help.env)$y) else max.kern <- 0
      all.x2 <- 1
      list.len <- 1
    }
    
    if (!is.null(obj$values$x)) {
      y <- obj$values$y
      help.knots <- c()
      list <- c()
      list.len <- length(x.factor[,1])
      for(i in 1:list.len) {
        name <- paste("weight",i,sep="")
        obj2 <- weight[i,]
        assign(name,obj2,envir=help.env)
      }
      if(!sort) {
        o <- order(y)
        y.help <- y[o]
        base.den <- base.den[,o]
        Z <- Z[o,]
      }
      else y.help <- y
      
      x.factor.len <- length(x.factor)
      for(j in 1:list.len) {
        set <- c()
        y <- c()
        for (i in 1:length(obj$values$y)) {
          if (all.equal(as.vector(Z[i,]),as.vector(x.factor[j,]))==TRUE) {
            y <- c(y,y.help[i])
            set <- c(set,i)
          }
        }
        assign(paste("y.later",j,sep=""),y,envir=help.env)
      }

      if(!is.null(val)) {
        y.val <- val
        assign("y.val",y.val,help.env)
        cond <- (y.val>=min(knots.val$val) & y.val <=max(knots.val$val))
        if(any(cond==FALSE)) stop("calculation of density is only possible in the range of the response",call.=FALSE)
        den.val <- list()
        sd.up.y.val <- list()
        sd.down.y.val <- list()
      }

      for(j in 1:list.len) { #which covariate?
        knot.left <- 1
        knot.right <- K-help.degree
   
        y.help <- seq(knots.val$val[knot.left],knots.val$val[knot.right],length=200)
        assign(paste("y.help",j,sep=""),y.help,envir=help.env)
        if(base=="bspline") {
          if(q>2) K.help <- K-q+2 else K.help <- 0
          base.val<- my.bspline(h,q,knots.val,y.help,K.help,plot.bsp=FALSE)
          help.base.den <- base.val$base.den
          assign(paste("base",j,sep=""),help.base.den,envir=help.env)
          if(!is.null(val)) {
            help.base.y.val<- my.bspline(h,q,knots.val,y.val,K.help,plot.bsp=FALSE)$base.den
            assign(paste("base.y.val",j,sep=""),help.base.y.val,envir=help.env)
          }
        }
 
        if(base=="gaussian") {
          for(i in 1:m) if(MeanW[i] >=y.r[1] & MeanW[i]<=y.r[2]) list <- c(list,i)
          nn <- matrix(1:length(y.help))
          help.base.den <- apply(nn,1,function(i,k,MeanW,Stand.abw,y.help) dnorm(y.help[i],MeanW,Stand.abw),MeanW,Stand.abw,y.help)[list,]
          assign(paste("base",j,sep=""),help.base.den,envir=help.env)
          if(!is.null(val)) {
            nn <- matrix(1:length(y.val))
            help.base.y.val<- apply(nn,1,function(i,k,MeanW,Stand.abw,y.help) dnorm(y.val,MeanW,Stand.abw),MeanW,Stand.abw,y.val)[list,]
            assign(paste("base.y.val",j,sep=""),help.base.y.val,envir=help.env) 
          }
          
        }
        
        bw <- get(paste("weight",j,sep=""),envir=help.env)%*%get(paste("base",j,sep=""),envir=help.env)
        if(!is.null(val)) {
          den.val[[j]] <- c(get(paste("weight",j,sep=""),envir=help.env)%*%get(paste("base.y.val",j,sep=""),envir=help.env))
          assign(paste("b.w.y.val",j,sep=""),den.val[[j]],envir=help.env)
        }
                
        max.bw.all <- c(max.bw.all,max(bw))
        assign(paste("b.w",j,sep=""),bw,envir=help.env)
        
        if(plot.base) max.bsp <- max(max(get(paste("b.w",j,sep=""),envir=help.env)),max.bsp)
        if(kernel) assign(paste("kern.den",j,sep=""),density(get(paste("y.later",j,sep=""),envir=help.env),kernel="epanechnikov"),envir=help.env)
        if(kernel) max.kern <- max(max(get(paste("kern.den",j,sep=""),envir=help.env)$y),max.kern) else max.kern <- 0
        
      }

      if(is.null(val)) y.val <- NULL else y.val <- val
      sd.cal <- variance.val(help.env,var.par,weight,K,x,list.len,Z,x.factor,y.val=y.val)

      max.bw <- max(max.bw.all)
      cut <- 0.0005*max.bw

      max.x.all <- c()
      min.x.all <- c()

      for(i in 1:list.len) {
        bw <- get(paste("b.w",i,sep=""),help.env)
        y.help <- get(paste("y.help",i,sep=""),help.env)
        ind <- which(bw<cut)
        if(length(ind)!=0) {
          assign(paste("b.w",i,sep=""),bw[-ind],help.env)
          y.help <- y.help[-ind]
        }
        assign(paste("y.help",i,sep=""),y.help,help.env)
        assign(paste("ind",i,sep=""),ind,help.env)

        max.x.all <- c(max(y.help),max.x.all)
        min.x.all <- c(min(y.help),min.x.all)
      }

      max.x <- max(max.x.all)
      min.x <- min(min.x.all)
    }

    AIC <- round(my.AIC,2)
    
    if(is.null(lambda0)) lambda0 <- 0
    lam <- round(lambda0,2)
    help.lam <- substitute(lambda ==s,list(s=lam))

    if(is.null(obj$values$x)) {
      ind <- get("ind1",help.env)
      if(length(ind)!=0) {
        assign("conf.plus1",get("b.w1",help.env)+2*get("sd.cal1",sd.cal)[-ind],help.env)
        assign("conf.minus1",get("b.w1",help.env)-2*get("sd.cal1",sd.cal)[-ind],help.env)
        if(!is.null(val)) {
          sd.up.y.val <- c(get("b.w1.y.val",help.env)+2*get("sd.cal1.y.val",sd.cal),help.env)
          sd.down.y.val <- c(get("b.w1.y.val",help.env)-2*get("sd.cal1.y.val",sd.cal),help.env)
        }
      }
      else {
        assign("conf.plus1",get("b.w1",help.env)+2*get("sd.cal1",sd.cal),help.env)
        assign("conf.minus1",get("b.w1",help.env)-2*get("sd.cal1",sd.cal),help.env)
        if(!is.null(val)) {
          sd.up.y.val <- c(get("b.w1.y.val",help.env)+2*get("sd.cal1.y.val",sd.cal))
          sd.down.y.val <- c(get("b.w1.y.val",help.env)-2*get("sd.cal1.y.val",sd.cal))
        }
      }
      list.len <- 1
    }
    else {
      list.len <- length(x.factor[,1])
      c1 <- c()
      for(i in 1:list.len) {#which grouping of covariates?
        ind <- get(paste("ind",i,sep=""),help.env)
        if(length(ind)!=0) {
          assign(paste("conf.plus",i,sep=""),conf.plus <- get(paste("b.w",i,sep=""),help.env)+2*get(paste("sd.cal",i,sep=""),sd.cal)[-ind])
          c1 <- c(c1,conf.plus)
          assign(paste("conf.minus",i,sep=""),get(paste("b.w",i,sep=""),help.env)-2*get(paste("sd.cal",i,sep=""),sd.cal)[-ind])
          if(!is.null(val)) {
            assign(paste("conf.plus.y.val",i,sep=""),sd.up.y.val[[i]] <- c(get(paste("b.w.y.val",i,sep=""),help.env)+2*get(paste("sd.cal.y.val",i,sep=""),sd.cal)))
            c1 <- c(c1,conf.plus)
            assign(paste("conf.minus.y.val",i,sep=""),sd.down.y.val[[i]] <- c(get(paste("b.w.y.val",i,sep=""),help.env)-2*get(paste("sd.cal.y.val",i,sep=""),sd.cal)))
          }
        }
        else {
          assign(paste("conf.plus",i,sep=""),conf.plus <- get(paste("b.w",i,sep=""),help.env)+2*get(paste("sd.cal",i,sep=""),sd.cal))
          c1 <- c(c1,conf.plus)
          assign(paste("conf.minus",i,sep=""),get(paste("b.w",i,sep=""),help.env)-2*get(paste("sd.cal",i,sep=""),sd.cal))
          if(!is.null(val)) {
            assign(paste("conf.plus.y.val",i,sep=""),sd.up.y.val[[i]] <- c(get(paste("b.w.y.val",i,sep=""),help.env)+2*get(paste("sd.cal.y.val",i,sep=""),sd.cal)))
            c1 <- c(c1,conf.plus)
            assign(paste("conf.minus.y.val",i,sep=""),sd.down.y.val[[i]] <- c(get(paste("b.w.y.val",i,sep=""),help.env)-2*get(paste("sd.cal.y.val",i,sep=""),sd.cal)))
          }
        } 
      }
    }

    if(is.null(main)) main.title <- substitute("K= "*a*", AIC= "*b*", "*c*"="*d, list(a=K,b=AIC,c=parse(text="lambda")[[1]],d=lam))
    if(is.null(sub) & base=="bspline") sub.title <- paste("base is", base, " order is", m) else sub.title <- sub
    if(is.null(sub) & base=="gaussian") sub.title <- paste("base is", base) else sub.title <- sub
    if(is.null(xlab)) xlab.title <- "y" else xlab.title <- xlab
    if(is.null(ylab)) ylab.title <- "density" else ylab.title <- ylab
    if(is.null(lwd)) lwd.value <- 3 else lwd.value <- lwd
    if(!is.null(legend.txt) & length(legend.txt)!=list.len) stop("Input length not equal to length of groupings")
    if(!is.null(legend.txt) & kernel) legend.txt <- c(legend.txt,paste("kernel density estimation of ", legend.txt,sep=""))

    if(is.null(obj$values$x)) x.range <- range(get("y.help1",help.env))
    else  x.range <- c(min.x,max.x)
    x.range.add <- 0.05*diff(x.range)
    x.lim <- c(x.range[1]-x.range.add,x.range[2]+x.range.add)

    if(plot.val==1 & latt==FALSE & plot.dens) {
      y.temp <- obj$values$y
      if(is.null(obj$values$x)) plot(y.temp, dnorm(y.temp,mean=0,sd=1),xlab=xlab.title,ylab=ylab.title,xlim=x.lim,ylim=c(0,max(max(get("conf.plus1",help.env)),max.kern,max.bsp)),type="n", main=main.title, sub=sub.title,cex.axis=1.2,cex.lab=1.5,cex.main=1.8)
      
      if(!is.null(obj$values$x)) plot(y.temp, dnorm(y.temp,mean=0,sd=1),xlab=xlab.title,ylab=ylab.title,xlim=x.lim,ylim=c(0,max(c1,max.kern,max.bsp)),type="n", main=main.title, sub=sub.title,cex.axis=1.2,cex.lab=1.5,cex.main=1.8)
      
      if(is.null(obj$values$x)) {
        lines(get(paste("y.help1",sep=""),envir=help.env),get(paste("b.w1",sep=""),envir=help.env),type="l",lwd=lwd.value)
        if(confi)lines(get("y.help1",help.env),get("conf.plus1",help.env),type="l",lwd=lwd.value-1,lty=3)
        if(confi)lines(get("y.help1",help.env),get("conf.minus1",help.env),type="l",lwd=lwd.value-1,lty=3)
        if(kernel) lines(get(paste("kern.den1",sep=""),envir=help.env),col=3,lwd=lwd.value,lty=2)
        if(!is.null(legend.txt)&!kernel) legend("topright",legend.txt,lty=1,lwd=lwd.value,cex=1.2)
        if(!is.null(legend.txt)&kernel) legend("topright",legend.txt,lty=c(1,2),col=c(1,3),lwd=lwd.value,cex=1.2)
        if(plot.base) {
          help <- get(paste("base",i,sep=""),envir=help.env)
          help.w <- get(paste("weight",i,sep=""),envir=help.env)
          help.l <- length(help[,1])
          for(j in 1:help.l) lines(get(paste("y.help",i,sep=""),envir=help.env),help[j,]*help.w[j],col=1)
        }
      }
      else {
        for(i in 1:list.len) {#which grouping of covariates?
          if(kernel) lines(get(paste("kern.den",i,sep=""),envir=help.env),col=2,lwd=lwd.value,lty=1+i-1)
          conf.plus <- get(paste("conf.plus",i,sep=""),help.env)
          conf.minus <- get(paste("conf.minus",i,sep=""),help.env)
          lines(get(paste("y.help",i,sep=""),envir=help.env),t(get(paste("b.w",i,sep=""),envir=help.env)),type="l",col=1,lwd=lwd.value,lty=1+i-1)
          if(confi) lines(get(paste("y.help",i,sep=""),help.env),conf.plus,type="l",col=1,lty=1+i-1,lwd=lwd.value-1)
          if(confi) lines(get(paste("y.help",i,sep=""),help.env),conf.minus,type="l",col=1,lty=1+i-1,lwd=lwd.value-1)
          if(!is.null(legend.txt)&!kernel) legend("topright",legend.txt,lty=seq(1,list.len),lwd=lwd.value,cex=1.2)
          if(!is.null(legend.txt)&kernel) legend("topright",legend.txt,lty=c(kronecker(seq(1,1+list.len-1),matrix(1,list.len,1))),col=c(kronecker(matrix(1,1,list.len),seq(1,1+list.len-1))),lwd=lwd.value)
          if(plot.base) {
            help <- get(paste("base",i,sep=""),envir=help.env)
            help.w <- get(paste("weight",i,sep=""),envir=help.env)
            help.l <- length(help[,1])
            for(j in 1:help.l) lines(get(paste("y.help",i,sep=""),envir=help.env),help[j,]*help.w[j],col=1)
          }
        }
      }
    }
    
    if(plot.val==1 & latt==TRUE &plot.dens) {
      main.title <- as.expression(main.title)
      y <- c()
      fy <- c()
      fac <- c()
      confi.plus <- c()
      confi.minus <- c()
      for(i in 1:list.len) {
        y <- c(y,get(paste("y.help",i,sep=""),envir=help.env))
        fy <- c(fy,get(paste("b.w",i,sep=""),envir=help.env))
        if(confi) {
          y <- c(y,rep(get(paste("y.help",i,sep=""),help.env),2))
          fy <- c(fy,get(paste("conf.plus",i,sep=""),help.env))
          fy <- c(fy,get(paste("conf.minus",i,sep=""),help.env))
        }
        if(!is.null(obj$values$x)) {
          fac <- c(fac,rep(paste(i,".x=",paste(x.factor[i,],collapse=",",sep=""),sep=""),length(get(paste("b.w",i,sep=""),envir=help.env))))
          if(confi) fac <- c(fac,rep(paste(i,".confi+ x=",paste(x.factor[i,],collapse=",",sep=""),sep=""),length(get(paste("conf.plus",i,sep=""),help.env))),rep(paste(i,".confi- x=",paste(x.factor[i,],collapse=",",sep=""),sep=""),length(get(paste("conf.plus",i,sep=""),help.env))))
        }
        if(is.null(obj$values$x)) {
          fac <- c(fac,rep("y~1",length(get(paste("b.w",i,sep=""),envir=help.env))))
          if(confi) fac <- c(fac,rep("confi+ y~1",length(get(paste("conf.plus",i,sep=""),help.env))),rep("confi- y~1",length(get(paste("conf.plus",i,sep=""),help.env))))
        }
      }
      graph.sets <-list(superpose.line=list(col=sort(rep(1:list.len,3))),superpose.symbol = list(col = sort(rep(1:list.len,3))))
      datafr <- data.frame(y,fy,fac)
      help55 <- xyplot(fy~y,groups=fac,type="l",auto.key=list(space="right",title="grouping",sort=FALSE),main=main.title,sub=sub.title,data=datafr,xlab=xlab.title,ylab=ylab.title,par.settings=graph.sets,lty=rep(c(2,2,1),list.len))
      print(help55)
    }
    if(!is.null(val)) {
      if(!is.null(obj$values$x)) {
        name <- "baseline"
        col.names <- colnames(obj$values$covariate$Z)[-1]
        rownames(x.factor) <- c()
        x.factor <- x.factor[,-1]
        if(is.matrix(x.factor)) {
          c.x.factor <- 1:length(x.factor[1,])
          r.x.factor <- 1:length(x.factor[,1])
        }
        if(is.vector(x.factor)) {
          c.x.factor <- 1
          r.x.factor <- length(x.factor)
          x.factor <- matrix(x.factor,r.x.factor,c.x.factor)
          r.x.factor <- 1:length(x.factor)
        }
        for(i in r.x.factor) {
          name.help <- c()
          for(j in c.x.factor) {
            if(x.factor[i,j]==1) name.help <- c(name.help,col.names[j])
          }
          if(length(name.help)>1) name.help <- paste(name.help,collapse="&")
          name <- c(name,name.help)
        }
        name <- noquote(name)
        names(den.val) <- name
        names(sd.up.y.val) <- name
        names(sd.down.y.val) <- name
      }
      return(list(y=y.val,fy=den.val,sd.down.y.val=sd.down.y.val,sd.up.y.val=sd.up.y.val))
    }
  }
  
if(plot.val==2) {
  sort <- obj$values$sort
  y.list <- list()
  sum.list <- list()
  x.factor <- obj$values$covariate$x.factor

  if(is.null(obj$values$x)) {
    N <- obj$splines$N
    p <- obj$splines$p
    y <- obj$values$y
    base.den2 <- obj$splines$base.den2
    y.order <- order(obj$values$y)
    K <- length(base.den2[,1])-1
    weight <- c(obj$results$ck)
    row.help <- rep(0,length(base.den2[1,]))
    sum <- c(0)
    for(k in 1:K) {
      sum <- weight[k]*colSums(base.den2[(k:(K+1)),]) +sum
    }
    #y<-sort(y)
    y.list[[1]] <- y
    sum.list[[1]] <- sum
    if(!latt&plot.dens) plot(sort(y),sum,xlab="y",ylab="F(y)",main="Distribution of f(y)")
    if(latt&plot.dens) {
      datafr <- data.frame(y,sum)
      obj <- xyplot(sum~sort(y),data=datafr,ylab="F(y)",xlab="y",main="Distribution of f(y)")
      print(obj)
    }
    if(!is.null(val)) {
      r.y<-range(obj$splines$knots.val$val)
      val<-sort(val)
      if(any(val<r.y[1])|any(val>r.y[2])) {
          print("Any value of val outside of the range of the observed values, extended with 0 at the left and 1 at the right of the interval of observed values")
          ind<-which(val>=r.y[1]&val<=r.y[2])
      }
      else ind<-c(1:length(val))
      base.den2<- my.bspline(h=obj$splines$h,q=obj$splines$q,knots.val=obj$splines$knots.val,y=val[ind],K=obj$splines$K-1,plot.bsp=FALSE)$base.den2
      sum.add <- c(0)
      for(k in 1:(obj$splines$K-1)) {
        sum.add <- obj$results$ck[k]*colSums(base.den2[(k:(obj$splines$K+1)),]) +sum.add
      }
      if(any(val<r.y[1])) sum.add<-c(rep(0,length(which(val<r.y[1]))),sum.add)
      if(any(val>r.y[2])) sum.add<-c(sum.add,rep(1,length(which(val>r.y[2]))))
      y.list<-list(val)
      sum.list<-list(sum.add)
      if(!latt&plot.dens) points(val,sum.add,col=2)
      if(latt&plot.dens) {
        datafr <- data.frame(val,sum.add)
        obj <- xyplot(sum.add~val,data=datafr,ylab="F(val)",xlab="y",main="Distribution of f(val)")
        print(obj)
      }
    }
    return(list(y=y.list,sum=sum.list))
  }
  if(!is.null(obj$values$x)) {
    base.den2 <- obj$splines$base.den2
    K <- length(base.den2[,1])-1
    help.env <- new.env()
    base <- obj$splines$base   
    knots.spline <- obj$splines$knots.spline$all
    beta <- obj$results$beta.val
    x.factor <- obj$values$covariate$x.factor
    N <- obj$splines$N
    y <- obj$values$y
    x <- obj$values$x
    K <- obj$splines$K
    len.x.fac <- length(x.factor[,1])
    all.x <- 0
    all.x2 <- 1
    Z <- obj$values$covariate$Z
    weight <- obj$results$ck

    for(i in 1:len.x.fac) {
      name <- paste("weight",i,sep="")
      assign(name,weight[i,],envir=help.env)
    }
    y.help <- c()
    y.call <- c()
    if(!sort) {
      o <- order(y)
      y.help <- y[o]
      base.den2 <- base.den2[,o]
    }
    else {
      y.help <- y
    }
    x.factor.len <- length(x.factor[1,])
    for(j in 1:len.x.fac) {
      set <- c()
      y <- c()
      for (i in 1:length(obj$values$y)) {
          if (all.equal(as.vector(Z[i,]),as.vector(x.factor[j,]))==TRUE) {
          y <- c(y,y.help[i])
          set <- c(set,i)
         }
      }
      assign(paste("y.later",j,sep=""),y,envir=help.env)
      assign(paste("y.list",j,sep=""),set,envir=help.env)
    }

    if(!latt) par(mfrow=c(2,ceiling(len.x.fac/2)))
    else{
      y.lattice <- c()
      sum.lattice <- c()
      factor <- c()
    }
    for(j in 1:len.x.fac) { #which covariate?
      y.r <- range(get(paste("y.later",j,sep=""),envir=help.env))
      y.list <- get(paste("y.list",j,sep=""),envir=help.env)
      if(base=="bspline") {
        row.help <- rep(0,length(base.den2[1,]))
        base.den2 <- rbind(base.den2,row.help)
        sum <- c(0)
        weight <- get(paste("weight",j,sep=""),help.env)
        for(k in 1:K) {
          sum <- weight[k]*colSums(base.den2[(k:(K+1)),y.list]) +sum
        }
        assign(paste("distr",j,sep=""),sum,envir=help.env)
        if(latt) {
          y.lattice <- c(y.lattice,get(paste("y.later",j,sep=""),envir=help.env))
          sum.lattice <- c(sum.lattice,sum)
          factor <- c(factor,rep(paste("Distribution-No",j,sep=""),length(sum)))
          y.call[[j]] <- get(paste("y.later",j,sep=""),envir=help.env)
          sum.list[[j]] <- sum
        }
        else {
          y.call[[j]] <- get(paste("y.later",j,sep=""),envir=help.env)
          sum.list[[j]] <- sum
          if(plot.dens) plot(get(paste("y.later",j,sep=""),envir=help.env),sum,xlab="y",ylab=paste("F(y|x",j,")",sep=""),main=paste("Distribution of f(y|x",j,")",sep=""))
        }
      }  
    }
    if(latt&plot.dens) {
      datafr <- data.frame(y.lattice,sum.lattice,factor)
      print(xyplot(sum.lattice~y.lattice|factor,data=datafr,xlab="y",ylab="F(y)",autokey=TRUE))
    }
  return(list(y=y.call,sum=sum.list))
  }
}
  if(plot.val==3) {
    help.env <- distr.func.help(obj)
    func <- distr.func(yi=NULL,obj,help.env)
    len.b <- length(obj$splines$base.den[,1])-obj$splines$q+1
    x.factor <- get("x.factor",envir=func)
    if(!is.null(obj$values$x)) x.factor.len <- length(x.factor[1,]) else x.factor.len <- 1
    knots.val <- obj$splines$knots.val
    all.x2 <- get("allx",envir=func)
    par(mfrow=c(ceiling(sqrt(all.x2)),1))
    y <- obj$values$y
    Z <- obj$values$covariate$Z
    x <- obj$values$x
    eps <- 1e-10
    y.help <- c()
    y.lattice <- c()
    sum.lattice <- c()
    factor <- c()
    eps <- 1e-10
    for(i in 1:all.x2) {
      w <- c()
      tt <- c()
      if(!is.null(x)) com.h <- x.factor[i,]

      if(!is.null(x)) {
        for(j in 1:length(Z[,1])) {
          if(identical(as.vector(Z[j,]),as.vector(com.h))==TRUE) y.help <- c(y.help,y[j])
        }
      }
      else y.help <- y
      
      min.y <- min(y.help)
      max.y <- max(y.help)
      val.min <- c()
      val.max <- c()
      for(k in 1:(length(knots.val$val)-1)) {
        if(knots.val$val[k]-eps <= min.y & min.y <= knots.val$val[k+1]+eps) val.min <- k
        if(knots.val$val[k]-eps <= max.y & max.y <= knots.val$val[k+1]+eps) val.max <- k
      }
      
      for(j in val.min:val.max) {
        funcy <- get(paste("distr.func",i,".",j,sep=""),envir=func)
        eval(parse(text=funcy))
        if(j==val.min) xi <- seq(min.y,knots.val$val[j+1],length=100)
        if(j!=val.min | j!=val.max) xi <- seq(knots.val$val[j],knots.val$val[j+1],length=100)
        if(j==val.max) xi <- seq(knots.val$val[j],max.y,length=100)
        ti <- obj(xi)
        w <- c(w,xi)
        tt <- c(tt,ti)
      }
      y.lattice <- c(y.lattice,w)
      sum.lattice <- c(sum.lattice,tt)
      factor <- c(factor,rep(paste("Distribution-No",i,sep=""),length(tt)))
      assign(paste("x",i,sep=""),w,envir=func)
      assign(paste("F(x)",i,sep=""),tt,envir=func)
      if(!latt&plot.dens) plot(w,tt,xlab="y",ylab="F(y)",main=paste("Distribution function of f(y|x",i,")",sep=""))
    }
    datafr <- data.frame(y.lattice,sum.lattice,factor)
    if(latt&plot.dens)    print(xyplot(sum.lattice~y.lattice|factor,data=datafr,xlab="y",ylab="F(y)",autokey=TRUE)) 
    return(datafr)  
  }
}

