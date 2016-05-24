`scbTest` <-
function(object,level=0.95,div=1000,select=NULL,drv=0,calc.stdev=TRUE)
{
   if (class(object)!="asp") stop("Only asp objects created by asp2() are supported.")
   if (object$info$pen$basis!="os") stop("Only B-spline basis functions supported.")
    fit=object
    x= object$info$pen$x #object$x[[2]]
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
        m=  object$info$pen$degree[[j]]#rep((object$info$pen$degree[j]),2)
        k=length(object$knots[[1]][[j]])
        smooth=list(m=m,bs.dim=k,knots= object$knots[[1]][[j]],term=object$info$pen$name[j],p.order=m)
        class(smooth)= "ospline.smooth"
        xx= as.matrix(x)[,j]
        data.xx = data.frame(x = xx)
        names(data.xx) <- object$info$pen$name[j]
        fit.mat$smooth[[j]]= smooth.construct.os.smooth.spec(smooth,data.xx,knots=data.frame(object$knots[[1]][[j]]))
    }

    crit=Stdev.fit= ucb=lcb=fitted= k0=cov.coef=WjXj=tstat=pval=list()
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
              Cnj=cbind(Cnj,CZ.temp$C)
              Znj=cbind(Znj,CZ.temp$Z)
          }
        }


        Cj= CZj$C[,drop=F]
        if (fit.mat$smooth[[j]]$m[1]==1) Cj=matrix(0,n,0)
        Zj= CZj$Z
        Cnj=cbind(matrix(rep(1,nrow(Cj))),object$info$lin$x,Cj,Cnj)
        Xnj=cbind(Cnj,Znj)

        Lambdaj  =  diag(c(Pen[(firstZ:lastZ)]))

        if (!hetero){
          if (ncol(object$x[[2]])>1|!is.null(object$info$lin$x)) {
            SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj)+diag(c(rep(0,ncol(Cnj)),Pen[-((firstZ:lastZ))]))),Xnj)%*%Zj)
            WjXj[[j]]=Zj-SnjXj
          }
          else{
            SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj)),Xnj)%*%Zj)
            WjXj[[j]]= Zj-SnjXj
          }
        }
        else {
          if (ncol(object$x[[2]])>1|!is.null(object$info$lin$x)) SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj,((object$sigmax$fitted^2)^(-1)*Xnj))+diag(c(rep(0,ncol(Cnj)),Pen[-((firstZ:lastZ))]))),((object$sigmax$fitted^2)^(-1)*Xnj))%*%Zj)
          else   SnjXj=Xnj%*%( tcrossprod(solve(crossprod(Xnj,((object$sigmax$fitted^2)^(-1)*Xnj))),((object$sigmax$fitted^2)^(-1)*Xnj))%*%Zj)
          WjXj[[j]]= (object$sigmax$fitted^(-2))*(Zj-SnjXj )
        }

        cov.coef[[j]]= solve(crossprod(Zj,WjXj[[j]]) +Lambdaj )
        cc.ev <- eigen(cov.coef[[j]],symmetric=F)
        cov.coef12 <- cc.ev$vectors%*%diag(sqrt(cc.ev$values))%*%t(cc.ev$vectors)

        integ <- function(xx,diffe) {
          data.grid <- data.frame(x=xx)
          names(data.grid) <- names(fit.mat$model)[1+j]
          if (drv==0) {
            CZj.grid <- Predict.matrix.lme(fit.mat$smooth[[j]],data.grid)
            Zj.grid= CZj.grid$Z[-diffe,]
          }
          else {
              CZj.grid <- Predict.matrix.lme(fit.mat$smooth[[j]],data.grid,drv=drv,center=F)
              Zj.grid= CZj.grid$Z[-diffe,]
           }
          SX <-  tcrossprod(cov.coef12,(Zj.grid))
          SX.norm = sqrt(apply(SX^2, 2, sum))
          SX/SX.norm
        }



      sx=seq(min(x[,j]),max(x[,j]),length=div)
      k0[[j]]=sum(sqrt(apply((integ(sx,div)-integ(sx,1))^2,2,sum)))

      # finding the critical value
        crit[[j]] <- .C("scritval", k0 = as.numeric(c(k0[[j]],1)), d = as.integer(1), cov = as.numeric(level),
                   m = as.integer(2), rdf = as.numeric(df), x = numeric(1), k = as.integer(1),PACKAGE="AdaptFitOS")$x


      if (calc.stdev){
if (drv>0) warning("Don't use SCBs for full sample in the case of derivative estimation (calc.stdev does currently not work!)")
        Sj= Zj%*%tcrossprod(cov.coef[[j]],(WjXj[[j]]))#t(Xj)%*%t(W[[j]])
        if (!hetero) SS <-tcrossprod(Sj)
        else SS <- (Sj)%*%((object$sigmax$fitted^2)*t(Sj))
        #     SS=Xj%*%tcrossprod(cov.coef[[j]],(WjXj[[j]]))%*%Xj
        fitted[[j]]=Sj%*%y
        Stdev.fit[[j]] <- sqrt(diag(SS)*(sigma2))
        ucb[[j]] <- fitted[[j]]+crit[[j]]*Stdev.fit[[j]]
        lcb[[j]] <- fitted[[j]]-crit[[j]]*Stdev.fit[[j]]
      }
      else Stdev.fit[[j]] =fitted[[j]]=ucb[[j]] <-lcb[[j]] <-NA

tstat[[j]]= max(abs((fitted[[j]]))*(Stdev.fit[[j]] )^(-1) )
pval[[j]] <- .C("stailp",crit=as.numeric(tstat[[j]]), k0 = as.numeric(c(k0[[j]],1)), d = as.integer(1), m = as.integer(2), rdf = as.numeric(df), x = numeric(1), k = as.integer(1),PACKAGE="AdaptFitOS")$x

  }
  scbt <- list(crit=crit, tstat=tstat,pval=pval, aspobject=fit,fit=fit.mat,cov.coef=cov.coef,WjXj=WjXj,Stdev=Stdev.fit,sigma2=sigma2,df=df,k0=k0,drv=drv,select=nsmooth,hetero=hetero,sigmax.fitted=object$sigmax$fitted,fitted=fitted,lcb=lcb,ucb=ucb)
  class(scbt) <- "scbTest"
  scbt
}
