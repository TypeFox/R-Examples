 plot.pencopula <- function(x,val=NULL,marg=TRUE,plot=TRUE,int=FALSE,main.txt=NULL,
                    sub.txt=NULL,contour=FALSE,cond=NULL,cuts=20,cex=1,cex.axes=1,cex.contour=1,
                    xlab=NULL,ylab=NULL,zlab=NULL,biv.margin=NULL,show.observ=FALSE,...) {
  env <- list()
  p <- get("p",x)
  q <- get("q",x)
  d <- get("d",x)
  pen.order <- get("pen.order",x)
  ddb <- get("ddb",x)
  Index.basis.D <- get("Index.basis.D",x)
  ck <- get("ck.val",x)
  D <- get("D",x)
  alpha <- get("alpha",x)
  base <- get("base",x)
  DD <- get("DD",x)
  #if(int & base=="Bernstein") stop("Distribution function is not supported for Bernstein polynomials")
  if(!is.null(cond)) cond.true <- TRUE  else cond.true <-  FALSE

  if(is.null(val))
    {
      #browser()
      if(p!=2 & (!cond.true & is.null(biv.margin))) stop("geht nicht!")
      else {
        if(p==2) {
	    x.grid <- seq(0,1,length=21)
            grid <- expand.grid(y1=x.grid, y2=x.grid)
            show.ob <- c(1,2)
        }  
        if(cond.true & is.null(biv.margin)) {
          if((length(cond))!=p) stop("change in the input for cond")
          list.full <- which(cond==-1)
          list.cond <- (1:p)[-list.full]
          name <- c()
          for(j in 1:p) {
            name <- c(paste("y",j,sep=""))
            if(j%in%list.full) {
              env[[noquote(name)]] <- get("Y",x)[,j]
            }
            else env[[noquote(name)]] <- rep(cond[j],get("n",x))
          }

          env.extend <- list()
          for(j in 1:p) {
            name <- c(paste("y",j,sep=""))
            env.extend[[noquote(name)]] <- with(env, do.breaks(range(env[[j]]),30))
          }
          grid <- unique(expand.grid(env.extend)) 
        }
        if(is.null(biv.margin)) {
          tilde.Psi.d <-  array(NA, dim=c(dim(grid)[1],get("ddb",x),p))
          for (j in 1:p)
            {
              if(int) tilde.Psi.d[,,j] <-  hierarch.bs(grid[,j], d = d, plot.bsp = FALSE,typ=3,penden.env=x,int=TRUE)$int.B.tilde else tilde.Psi.d[,,j] <-  hierarch.bs(grid[,j], d = d, plot.bsp = FALSE,typ=3,penden.env=x,int=FALSE)$B.tilde
            }
          tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
          for (j in 2:p)
            {        
              tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
            }                                        #browser()
          grid[["plot"]] <- tilde.PSI.d.D%*%ck
        }
        if(!is.null(biv.margin)) {
   
          env.extend <- list()
          for(j in 1:p) {
            name <- c(paste("y",j,sep=""))
            env[[noquote(name)]] <- seq(0,1,length=21)
          }
 
          grid <- unique(expand.grid(env)) 
          tilde.Psi.d <-  array(NA, dim=c(dim(grid)[1],get("ddb",x),p))
          for (j in biv.margin)
            {
              if(int) tilde.Psi.d[,,j] <-  hierarch.bs(grid[,j], d = d, plot.bsp = FALSE,typ=3,penden.env=x,int=TRUE)$int.B.tilde else tilde.Psi.d[,,j] <-  hierarch.bs(grid[,j], d = d, plot.bsp = FALSE,typ=3,penden.env=x,int=FALSE)$B.tilde
            }
          k <- seq(1,p)[-biv.margin]

          for(j in k) tilde.Psi.d[,,j] <- matrix(rep(1,dim(grid)[1]*get("ddb",x),nrow=dim(grid)[1]))

          tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
          for(j in 2:p)  tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
            
        #browser()
          grid[["plot"]] <- tilde.PSI.d.D%*%ck
          show.ob <- biv.margin
        }

        if(!cond.true) {
          lam1 <- get("lambda",x)[1]
          lam2 <- get("lambda",x)[2]
        }
                
        if(is.null(main.txt) & !cond.true & base=="B-Spline") {
          main.txt <- substitute("d="*a*", D="*b*", "*c*"="*d*", "*e*"="*f*", pen.order="*g*", q="*h, list(a=d,b=D,c=parse(text="lambda[1]")[[1]],d=lam1, e=parse(text="lambda[2]")[[1]],f=lam2,g=pen.order,h=q))
          main.txt <- as.expression(main.txt)
        }
        if(is.null(main.txt) & cond.true & base=="B-Spline") {
          main.txt <- substitute("d="*a*", D="*b*", pen.order="*g, list(a=d,b=D,g=pen.order))
          main.txt <- as.expression(main.txt)
        }
        if(is.null(main.txt) & base=="Bernstein") {
          main.txt <- substitute("Marg. Coeff.="*a*", Number all Coeff.="*b,list(a=ddb,b=DD))
          main.txt <- as.expression(main.txt)
        }
           
        k <- dim(x$liste)[1]
        log.like <- round(get("log.like",x),3)
        pen.log.like <- round(get("pen.log.like",x),3)
        if(is.null(sub.txt)) sub.txt <- paste("log like=",log.like,", pen. log like= ",pen.log.like,", AIC=",round(get("AIC",x),3),", alpha=",alpha,sep="")
        if(int) z.txt <- "distribution" else z.txt <- "density"

        if(cond.true) {
          values <- c()
          for(j in 1:length(list.full)) values <- c(values,paste("y",list.full[j],sep=""))
          values <- as.formula(paste("plot~",paste(values,collapse="*"),sep=""))
        }
        else {
          hh <- c("y1","y2")
          values <- as.formula(paste("plot~",paste(hh,collapse="*"),sep=""))
        }
        if(!is.null(biv.margin)) {
          values <- c()
          for(j in biv.margin) values <- c(values,paste("y",j,sep=""))
          values <- as.formula(paste("plot~",paste(values,collapse="*"),sep=""))
        }
     
        if(!contour) obj1 <- wireframe(values,data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(cex=cex.axes,label=xlab),
                                       ylab=list(cex=cex.axes,label=ylab),scales=list(arrows=FALSE,col="black",font=3,x=list(cex=cex),y=list(cex=cex),z=list(cex=cex)),
                                       main=main.txt,shade=TRUE,zlim=c(0,max(grid)), par.settings = list(axis.line = list(col = "transparent")),
                                       par.box = c(col = "black"))
        if(contour&!show.observ) obj1 <- contourplot(values, data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(label=xlab,cex=cex.axes),
                                 ylab=list(label=ylab,cex=cex.axes),scales=list(arrows=FALSE,col="black",font=3,cex=cex),
                                 zlim=c(0,max(grid)),main=main.txt,shade=TRUE,cuts=cuts,subscripts=TRUE)

        if(contour&show.observ) {
          yy <- NA
          rm(yy)
          grid['xx'] <- c(get("Y",x)[,biv.margin[1]],rep(NA,dim(grid)[1]-length(get("Y",x)[,biv.margin[1]])))
          grid['yy'] <- c(get("Y",x)[,biv.margin[2]],rep(NA,dim(grid)[1]-length(get("Y",x)[,biv.margin[2]])))
 
          obj1 <- contourplot(values, data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),
                               xlab=list(label=xlab,cex=cex.axes), ylab=list(label=ylab,cex=cex.axes),scales=list(arrows=FALSE,col="black",font=3,cex=cex),
                               panel = function(contour,labels,cex,...)  panel.contourplot(contour=TRUE,labels=list(labels=TRUE,cex=cex.contour),...),
                              zlim=c(0,max(grid)),main=main.txt,shade=TRUE,cuts=cuts,subscripts=TRUE)+layer(panel.points(x=xx,y=yy,subscripts=TRUE,under=TRUE,pch=1),data=grid)
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
          if(plot) {
            print(obj1,position=c(0,0.35,1,1),more=TRUE)
            print(obj2,position=c(0,0,1,0.35))
          }
          else return(list(density=obj1,marg.density=obj2))
        }
        else if(plot) {
          print(obj1)
          check <- any(grid$plot<0)
          print(check)
          if(check) print(grid[grid$plot<0,])
        }
        else return(list(density=obj1,grid=grid))
      }
    }
  else
    {
        tilde.Psi.d <-  array(NA, dim=c((length(val)/p),(2**d)+q,p))
        if(!is.matrix(val)) val <- matrix(val,(length(val)/p),p)

        for (j in 1:p)
          {
            if(int) tilde.Psi.d[,,j] <-  hierarch.bs(val[,j], d = d, plot.bsp = FALSE,typ=3,penden.env=x,int=TRUE)$int.B.tilde
            else tilde.Psi.d[,,j] <-  hierarch.bs(val[,j], d = d, plot.bsp = FALSE,penden.env=x,typ=3,int=FALSE)$B.tilde
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
