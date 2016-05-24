##Model selection,specification by a step-wise like AIC procedure in a "forward" direction
model.selection.gwr<-function(DeVar=NULL,InDeVars=NULL, data=list(),bw=NULL,approach="CV",
                     adaptive=F,kernel="bisquare",dMat=NULL,p=2, theta=0, longlat=F)
{
   if(is.null(DeVar)||!is.character(DeVar)||is.null(InDeVars)||!is.character(InDeVars))
    stop("Input are not correct, please recheck!")
   ##Data points
  spdf<-data
  if (!is.null(data))
  {
    if (is(data, "Spatial"))
      {
        p4s <- proj4string(data)
        dp.locat<-coordinates(data)
        data <- as(data, "data.frame")
      }
    else
      {
        if (!is(data, "data.frame"))
           stop("Given regression data must be data.frame or Spatial*DataFrame")
      }
  }
  else stop("No regression data frame is avaiable!")
  #################################################################
   vars.df<-names(data)
   dp.n<-nrow(data)
   var.n<-length(InDeVars)
   InDeVars.Sub<-InDeVars
   model.list<-list()
   GWR.df<-c()
   
   if (is.null(dMat))
    {
      dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
    }
  #vars.idxs<-c()###Record indices used for each model
   varsindx.list<-list()
   level.vars<-c()
   tag<-1
   adapt<-NULL
   for (i in 1:var.n)
  {
    AICcs<-c()
    for (j in 1:(var.n-i+1))
    {
      vars.j<-c(level.vars,InDeVars.Sub[j])
      fml<-Generate.formula(DeVar,vars.j)
	  cat("Now calbrating the model: \n", fml,"\n")
	  matL<-extract.mat(fml, data)
	  y<-matL[[1]]
	  x<-matL[[2]]
      if (is.null(bw))
      {
        part1<-paste("bandwidth<-bw.gwr(",fml,sep="")
        part2<-"data=spdf,kernel=kernel,approach=approach,dMat=dMat)"
        #part2<-paste(paste(part2.1,part2.2,sep=","),")",sep="")
        #part1<-paste("bw<-bw.sel(",fml)
        #part2<-paste()
        expression<-paste(part1,part2,sep=",")
        print(expression)
        eval(parse(text=expression))
      }
      else
      {
        if (adaptive)
        {
          stopifnot(is.numeric(bw))
          stopifnot((bw >= 0))
        }
        else
        {
          stopifnot(is.numeric(bw))
          stopifnot((bw > min(dMat)))
        }
        bandwidth<-bw
       } 

        ##############Calibrate the GWR model
        S<-matrix(nrow=dp.n,ncol=dp.n)
        betas <-matrix(nrow=dp.n, ncol=var.n)
          for (i in 1:dp.n)
          {
            dist.vi<-dMat[,i]
            W.i<-gw.weight(dist.vi,bandwidth,kernel,adaptive)
            gw.resi<-gw.reg(x,y,W.i,hatmatrix=T,i)
            #betas[i,]<-gw.resi[[1]] ######See function by IG
            S[i,]<-gw.resi[[2]]
            #Ci<-gw.resi[[3]]
          }
	      
        
        v1 <- sum(diag(S))
        v2<-0
         for (i in 1:dp.n)
	        v2 <-v2+ sum(S[,i]^2)
        Q <- t(diag(dp.n)-S)%*%(diag(dp.n)-S)
  	    rss <- c(t(y)%*%Q%*%y)
        sigma2 <- rss/dp.n
        AIC <- dp.n*log(sigma2) + dp.n*log(2*pi) +dp.n+v1
        AICc <- dp.n*log(sigma2) + dp.n*log(2*pi) + dp.n *((dp.n + v1) / (dp.n - 2 - v1))
        model.list[[tag]]<-list(fml, vars.j)
        GWR.df<-rbind(GWR.df, c(bw, AIC, AICc, rss))
        AICcs<-c(AICcs, AICc)
        tag<-tag+1
    }
    idx<-which.min(AICcs)[1]
    level.vars<-c(level.vars,InDeVars.Sub[idx])
    InDeVars.Sub<-InDeVars.Sub[-idx]
  }
  res<-list(model.list,GWR.df)
  res
}

#######Extract model matrix
extract.mat<-function(formula, data=list())
{
	this.call <- match.call()
	if (!is.null(data))
    {
       if (is(data, "Spatial"))
	   {
	     data <- as(data, "data.frame")
	   }
       else
       {
         if (!is(data, "data.frame"))
           stop("Given regression data must be data.frame or Spatial*DataFrame")
       }
    }
    else stop("No regression data frame is avaiable!")
	mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    mat.list<-list(y,x)
	mat.list
}

model.view.gwr<-function(DeVar, InDeVars, model.list)
{
  n<-length(InDeVars)
  if (n>10)
  {
    cex<-10/n
  }
  else
  {
    cex<-1
  }

  #InDeVars<-sort(InDeVars)
  numModels<-length(model.list)
  alpha<-2*pi/numModels
  cols<-rainbow(n)
  pchs<-rep(c(8,9,10,15,16,17,18,23,24),length.out=n)
  plot(x=0,y=0,xlim=c(-3*n/4, n+6),ylim=c(-n-1, n+1), cex=2, axes=F, pch=22,xlab="",ylab="",main="View of GWR model selection with different variables")
  for (i in 1:numModels)
  {
    vars<-model.list[[i]][[2]]
    nvar<-length(vars)
    p1<-c(0,0)
    for (j in 1:nvar)
    {
      radius<-sqrt(n)*sqrt(j)
      var.idx<-which(InDeVars==vars[j])
      coord<-c(radius*cos((i-1)*alpha),radius*sin((i-1)*alpha))
      lines(x=c(p1[1], coord[1]),y=c(p1[2], coord[2]), col="grey",lwd=cex)
      points(x=coord[1], y=coord[2], col=cols[var.idx],pch=pchs[var.idx],cex=(cex*i/numModels+0.3))
      p1<-coord
    }
    text(x=(radius+0.5)*cos((i-1)*alpha),y=(radius+0.5)*sin((i-1)*alpha), as.character(i), cex=cex*0.6)
  }
  legend(x=n+2, y=n/2, col=c("black",cols),pch=c(22,pchs), c(DeVar, InDeVars),box.col="white")
}

model.sort.gwr<-function(Sorting.list , numVars, ruler.vector)
{
  n<-length(Sorting.list)
  numMoldes<-length(ruler.vector)
  indxs<-c()
  tag<-0
  for (i in numVars:1)
  {
    tmpV<-ruler.vector[(tag+1):(tag+i)]
    indx<-sort(tmpV, decreasing=T, index=T)$ix
    indxs<-c(indxs, indx+tag)
    tag<-tag+i
  }
  res<-list()
  for (i in 1:n)
    {
      list.i<-Sorting.list[[i]]
      if (is.list(list.i))
      {
         tmp.L<-list()
         for (j in 1:numMoldes) tmp.L[[j]]<-list.i[[indxs[j]]]
         res[[i]]<-tmp.L
      }
      else
      {
        tmp.V<-c()
        for (j in 1:numMoldes) tmp.V<-rbind(tmp.V,list.i[indxs[j],])
        res[[i]]<-tmp.V
      }
    }
    res
}

####Generate formula based on the given depedent and indepednent variables
Generate.formula<-function(DeVar,InDeVars)
{
  fml<-paste(paste(DeVar, "~", sep=""), InDeVars[1],sep="")
  var.n<-length(InDeVars)
  if (var.n>1)
  {
    for (i in 2:var.n)
      fml<-paste(fml, InDeVars[i], sep="+")
  }
  fml
}