plot.paircopula <- function(x,val=NULL,marg=TRUE,plot.view=TRUE,int=FALSE,main.txt=NULL,
                    sub.txt=NULL,contour=FALSE,cuts=20,cex=1,cex.axes=1,
                    xlab=NULL,ylab=NULL,zlab=NULL,xlim=NULL,ylim=NULL,zlim=NULL,show.observ=FALSE,margin.normal=FALSE,...) {
  if(!class(x)=="paircopula") stop("obj has to be of class paircopula")
  env <- list()
  p <- get("p",x)
  q <- 0
  d <- get("d",x)
  ddb <- get("ddb",x)
  Index.basis.D <- get("Index.basis.D",x)
  ck <- get("ck.val",x)
  D <- get("D",x)
  alpha <- 0
  base <- get("base",x)
  index.b <- matrix(0:get("dd",x))
  #if(int & base=="Bernstein") stop("Distribution function is not supported for Bernstein polynomials")

  if(is.null(val))
    {
      if(p!=2) stop("geht nicht!")
      else {
        x.grid <- seq(0,1,length=51)
        grid <- expand.grid(y1=x.grid, y2=x.grid)  
        if(margin.normal) {
           x.grid.unif <- seq(-5,5,length=51)
           x.grid <- pnorm(x.grid.unif)
           grid <- expand.grid(y1=x.grid, y2=x.grid) 
        }
        tilde.Psi.d <- array(NA, dim=c(dim(grid)[1],get("ddb",x),p))
        for (j in 1:p)
          {
            if(base=="Bernstein") {
              if(int) tilde.Psi.d[,,j] <- int.bernstein(x,Y=grid[,j])
              else tilde.Psi.d[,,j] <- apply(index.b,1,bernstein,x=grid[,j],n=get("dd",x))
            }
            if(base=="B-spline") {
              if(int) tilde.Psi.d[,,j] <- int.bspline2(x,Y=grid[,j])
              else tilde.Psi.d[,,j] <- my.bspline(y=grid[,j],K=get("K",x)+get("q",x)-1,q=get("q",x),margin.normal=margin.normal)$base.den
            }  
          }
        tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
        
        for (j in 2:p)
          {        
            tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
          }
        grid[["plot"]] <- tilde.PSI.d.D%*%ck
        if(margin.normal) {
           plot<-grid$plot
           grid <- expand.grid(y1=x.grid.unif, y2=x.grid.unif)
           grid[["plot"]]<-plot*dnorm(grid$y1)*dnorm(grid$y2)
        }
        if(is.null(zlim)) zlim <- c(0,round(max(grid[["plot"]]),2))
        if(is.null(ylim)) ylim <- range(grid$y2)
        if(is.null(xlim)) xlim <- range(grid$y1)
        lam1 <- get("lambda",x)[1]               
        if(is.null(main.txt)) {
          main.txt <- substitute("K="*a*","*c*"="*d,list(a=D,c=parse(text="lambda")[[1]],d=lam1))
          main.txt <- as.expression(main.txt)
        }
        k <- dim(x$liste)[1]
        log.like <- round(get("log.like",x),3)
        pen.log.like <- round(get("pen.log.like",x),3)
        if(is.null(sub.txt)) sub.txt <- paste("log like=",log.like,", pen. log like= ",pen.log.like,", cAIC=",round(get("cAIC",x),3),sep="")
        if(int) z.txt <- "distribution" else z.txt <- "density"
        hh <- c("y1","y2")
        values <- as.formula(paste("plot~",paste(hh,collapse="*"),sep=""))
        if(!contour) obj1 <- wireframe(values,data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(cex=cex.axes,label=xlab),
                                       ylab=list(cex=cex.axes,label=ylab),scales=list(arrows=FALSE,col="black",font=3,x=list(cex=cex),y=list(cex=cex),z=list(cex=cex)),
                                       main=main.txt,shade=TRUE,zlim=zlim, par.settings = list(axis.line = list(col = "transparent")),
                                       par.box = c(col = "black"))
        if(contour) obj1 <- contourplot(values, data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),
                                 xlab=list(label=xlab,cex=cex.axes),ylab=list(label=ylab,cex=cex.axes),
                                 scales=list(arrows=FALSE,col="black",font=3,cex=cex,x=list(lim=xlim),y=list(lim=ylim)),
                                 main=main.txt,shade=TRUE,cuts=cuts)
        if(contour&show.observ) {
                                x1 <- get("Y",x)[,1]
                                x2 <- get("Y",x)[,2]
                                grid2 <- data.frame(x1=x1,x2=x2)
                                obj1 <- obj1+layer(panel.points(x=x1,y=x2,subscripts=TRUE,under=TRUE,pch=1),data=grid2)
        }
        if(marg) {
          T.marg <- get("T.marg",x)
          base <- get("tilde.Psi.knots.d",x)
          
          xx <- rep(seq(0,1,length=get("ddb",x)),p)
          density <- c()
          fac <- c()
          
          for(j in 1:p)
            {
              density <- c(density,round(base%*%(T.marg[,,j]%*%ck),5))
              fac <- c(fac,rep(j,get("ddb",x)))
            }
          datafr <- data.frame(xx,density,fac)
          graph.sets <-list(superpose.line=list(col=c(1:p),superpose.symbol = list(col = c(1:p))))
          
          obj2 <- xyplot(density~xx|fac,type="l",auto.key=list(space="right",title="marginal densities",sort=FALSE),par.settings=graph.sets)
          if(plot.view) {
            print(obj1,position=c(0,0.35,1,1),more=TRUE)
            print(obj2,position=c(0,0,1,0.35))
          }
          else return(list(density=obj1,marg.density=obj2))
        }
        if(plot.view) {
          print(obj1)
          #check <- any(grid$plot<0)
          #if(check) print(grid[grid$plot<0,])
        }
        else return(list(density=obj1,grid=grid))
      }
    }
  else
    {
      if(!is.matrix(val)) {
        if(is.data.frame(val)) val <- as.matrix(val) else stop("val has to be a data.frame or a matrix")
      }
        tilde.Psi.d <- array(NA, dim=c((length(val)/p),get("dd",x)+1,p))
        val <- matrix(val,(length(val)/p),p)
        for (j in 1:p)
          {
            if(base=="Bernstein") {
              if(int) tilde.Psi.d[,,j] <- int.bernstein(x,Y=val[,j])
              else tilde.Psi.d[,,j] <- apply(index.b,1,bernstein,x=val[,j],n=get("dd",x))
            }
            if(base=="B-spline") {
              if(int) tilde.Psi.d[,,j] <- int.bspline2(x,Y=val[,j])
              else tilde.Psi.d[,,j] <- my.bspline(y=val[,j],K=get("K",x)+get("q",x)-1,q=get("q",x))$base.den
            }  
          }
        tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
   
        for (j in 2:p)
          {  
            tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
          }
        datafr <- data.frame(val,tilde.PSI.d.D%*%ck)
        colnames(datafr)[p+1] <- "fit"
        return(datafr)
      }
}
