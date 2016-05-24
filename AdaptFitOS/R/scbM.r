# Function to calculate the simultaneous confidence bands using the bias correction, resulted from the mixed model represenation of penalized splines.
scbM=function(object,select=NULL,drv=0,level=0.95,div=1000,calc.stdev=TRUE,bayes=FALSE)
{
   if (class(object)!="asp") stop("Only asp objects created by asp2() are supported.")
    fit=object
    x= object$info$pen$x 
    y=object$info$y
    n <- length(y)

    if (is.null(select)) nsmooth=1:length(object$info$pen$degree) else nsmooth=select
    if (!is.null(object$sigmax$fitted)) hetero=T else hetero=F
    fit.mat=list()
    fit.mat$model=data.frame(y=y,x=x)
    names(fit.mat$model)[2:ncol(fit.mat$model)]=object$info$pen$name
    sigma2 <- object$fit$sigma^2
    Pen<-sigma2/diag(object$aux$random.var)

    for (j in 1:length(object$info$pen$degree)){
      if (object$info$pen$basis=="os"){
        m=  object$info$pen$degree[[j]]
        k=length(object$knots[[1]][[j]])
        smooth=list(m=m,bs.dim=k,knots= object$knots[[1]][[j]],term=object$info$pen$name[j],p.order=m)
        class(smooth)= "ospline.smooth"
        xx= as.matrix(x)[,j]
        data.xx = data.frame(x = xx)
        names(data.xx) <- object$info$pen$name[j]
        fit.mat$smooth[[j]]= smooth.construct.os.smooth.spec(smooth,data.xx,knots=data.frame(object$knots[[1]][[j]]))
      } else {
        m= c(object$info$pen$degree[[j]],object$info$pen$degree[[j]]+1)
        k=length(object$knots[[1]][[j]])
        fit.mat$smooth[[j]]=list(m=m,bs.dim=k,knots= object$knots[[1]][[j]],term=object$info$pen$name[j])
        if (object$info$pen$basis=="trunc.poly") class(fit.mat$smooth[[j]])="tlspline.smooth"
        else class(fit.mat$smooth[[j]])= object$info$pen$basis
      }
    }

    crit=Stdev.fit= ucb=lcb=fitted= k0=cov.coef=WjXj=list()
    for (j in nsmooth){
      m=sum(sapply(fit.mat$smooth, FUN=function(x) return(x$m[1])))
      df <- n-m-1

     ##################
        Cnj=Znj=numeric()
        for (i in 1:ncol(object$x[[2]])){
          if (i==j) {
              if (i==1) firstZ=1 else firstZ=ncol(Znj)+1
              data.grid <- data.frame(x=x[,i])
              names(data.grid) <- names(fit.mat$model)[1+i]
              CZj= Predict.matrix.lme(fit.mat$smooth[[i]],data.grid)
              lastZ=firstZ+ncol(CZj$Z)-1
          }
          else{
              data.grid = data.frame(x = x[,i])
              names(data.grid) <- names(fit.mat$model)[1+i]
              CZ.temp <- Predict.matrix.lme(fit.mat$smooth[[i]],data.grid)
              if (fit.mat$smooth[[i]]$m[1]==1) CZ.temp$C=matrix(0,n,0)
              if (class(fit.mat$smooth[[j]])=="ospline.smooth") Cnj=cbind(Cnj,CZ.temp$C)
              else Cnj=cbind(Cnj,CZ.temp$C[,-1])
              Znj=cbind(Znj,CZ.temp$Z)
          }
        }


        if (class(fit.mat$smooth[[j]])=="ospline.smooth") Cj= CZj$C[,drop=F]  else Cj= CZj$C[,-1,drop=F]
        if (fit.mat$smooth[[j]]$m[1]==1) Cj=matrix(0,n,0)
        Zj= CZj$Z
        Xj=cbind(Cj,Zj)

        Cnj=cbind(matrix(rep(1,n)),object$info$lin$x,Cnj)
        Xnj=cbind(Cnj,Znj)

        Lambdaj  =  diag(c(rep(0,ncol(Cj)),Pen[(firstZ:lastZ)]))

        if (!hetero){
          if (ncol(object$x[[2]])>1|!is.null(object$info$lin$x)) {
            SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj)+diag(c(rep(0,ncol(Cnj)),Pen[-((firstZ:lastZ))]))),Xnj)%*%Xj)
            WjXj[[j]]=Xj-SnjXj
          }
          else  WjXj[[j]]= Xj
        }
        else {
          if (ncol(object$x[[2]])>1|!is.null(object$info$lin$x)) SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj,((object$sigmax$fitted^2)^(-1)*Xnj))+diag(c(rep(0,ncol(Cnj)),Pen[-((firstZ:lastZ))]))),((object$sigmax$fitted^2)^(-1)*Xnj))%*%Xj)
          else   SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj,((object$sigmax$fitted^2)^(-1)*Xnj))),((object$sigmax$fitted^2)^(-1)*Xnj))%*%Xj)
          WjXj[[j]]= (object$sigmax$fitted^(-2))*(Xj-SnjXj )
        }

        cov.coef[[j]]= solve(crossprod(Xj,WjXj[[j]]) +Lambdaj )
        cc.ev <- eigen(cov.coef[[j]],symmetric=F)
        cov.coef12 <- cc.ev$vectors%*%diag(sqrt(cc.ev$values))%*%t(cc.ev$vectors)

        integ <- function(xx,diffe) {
          data.grid <- data.frame(x=xx)
          names(data.grid) <- names(fit.mat$model)[1+j]
          if (drv==0) {
            CZj.grid <- Predict.matrix.lme(fit.mat$smooth[[j]],data.grid)
            if (class(fit.mat$smooth[[j]])=="ospline.smooth") Cj.grid= CZj.grid$C[-diffe,,drop=F]  else Cj.grid= CZj.grid$C[-diffe,-1,drop=F]
            Zj.grid= CZj.grid$Z[-diffe,]
            if (fit.mat$smooth[[j]]$m[1]==1) Cj.grid=matrix(0,nrow(Zj.grid),0)
            Xj.grid=  cbind(Cj.grid,Zj.grid)
          }
          else {
            if (class(fit.mat$smooth[[j]])=="ospline.smooth") {
              CZj.grid <- Predict.matrix.lme(fit.mat$smooth[[j]],data.grid,drv=drv,center=F)
              Cj.grid= CZj.grid$C[-diffe,,drop=F]
              Zj.grid= CZj.grid$Z[-diffe,]
              Xj.grid=  cbind(Cj.grid,Zj.grid)
            }
            else Xj.grid=drvbasis(xx[-diffe],degree=fit.mat$smooth[[j]]$m[1],knots=fit.mat$smooth[[j]]$knots,drv=drv,basis=class(fit.mat$smooth[[j]]))[,-1]
          }
          SX <-  tcrossprod(cov.coef12,(Xj.grid))
          SX.norm = sqrt(apply(SX^2, 2, sum))
          SX/SX.norm
        }



      sx=seq(min(x[,j]),max(x[,j]),length=div)
      k0[[j]]=sum(sqrt(apply((integ(sx,div)-integ(sx,1))^2,2,sum)))

      # finding the critical value
        crit[[j]] <- .C("scritval", k0 = as.numeric(c(k0[[j]],1)), d = as.integer(1), cov = as.numeric(level),
                   m = as.integer(2), rdf = as.numeric(df), x = numeric(1), k = as.integer(1),PACKAGE="AdaptFitOS")$x

      cat("Critical value for f(",names(fit.mat$model)[1+j],"): ", crit[[j]], "\n",sep="")

      if (calc.stdev){
if (drv>0) warning("Don't use SCBs for full sample in the case of derivative estimation (calc.stdev does currently not work!)")
        Sj= Xj%*%tcrossprod(cov.coef[[j]],(WjXj[[j]])) #t(Xj)%*%t(W[[j]])
        fitted[[j]]=Sj%*%y
        
        if (bayes==F){
            if (!hetero) SS <-tcrossprod(Sj)
            else SS <- (Sj)%*%((object$sigmax$fitted^2)*t(Sj))
            #     SS=Xj%*%tcrossprod(cov.coef[[j]],(WjXj[[j]]))%*%Xj
            Stdev.fit[[j]] <- sqrt(diag(SS)*(sigma2))
        }
        else{
            S= Xj%*%tcrossprod(cov.coef[[j]],Xj) #t(Xj)%*%t(W[[j]])
            Stdev.fit[[j]] <- sqrt(diag(S)*(sigma2))
        }

        ucb[[j]] <- fitted[[j]]+crit[[j]]*Stdev.fit[[j]]
        lcb[[j]] <- fitted[[j]]-crit[[j]]*Stdev.fit[[j]]
      }
      else Stdev.fit[[j]] =fitted[[j]]=ucb[[j]] <-lcb[[j]] <-NA
    }

  scbm <- list(crit=crit,aspobject=fit,fit=fit.mat,cov.coef=cov.coef,WjXj=WjXj,Stdev=Stdev.fit,sigma2=sigma2,df=df,k0=k0,drv=drv,select=nsmooth,hetero=hetero,sigmax.fitted=object$sigmax$fitted,fitted=fitted,lcb=lcb,ucb=ucb)
  class(scbm) <- "scbm"
  scbm
}
