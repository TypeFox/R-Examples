spatialPredict.copula <- function(object,...){
  if(!is.null(object$params$nmax)){
    searchradius<-min(object$params$nmax,dim(object$observations)[1],50)
  } else {
    searchradius<-min(dim(object$observations)[1],50)
  }
  estimation = object$copulaParams
  
  predictions<-bayesCopula(object,estimation,search=searchradius,object$outputWhat, 
                       testMean = object$params$testMean)
  names(predictions)[1:2] = c("mean","variance")
  for (i in 1:length(predictions)) {
    predictions[[i]][is.na(predictions[[i]])] = 99999999
  }
  object$predictions = SpatialPointsDataFrame(coordinates(object$predictionLocations),
     data=as.data.frame(predictions), match.ID = FALSE)
  if (gridded(object$predictionLocations)) gridded(object$predictions) = TRUE
  if (!is.na(proj4string(object$predictionLocations))) 
     proj4string(object$predictions) = proj4string(object$predictionLocations) 
  object
}

bayesCopula <- function(obj,estimates,search=10,calc=list(mean=TRUE,variance=TRUE),testMean = FALSE){
  if (is.null(testMean)) testMean = FALSE
  if (is.null(calc$mean)) calc$mean=FALSE
  if (is.null(calc$variance)) calc$variance=FALSE
  if (testMean==TRUE && calc$mean==FALSE) testMean=1

  data<-NULL
  data$x<-coordinates(obj$observations)[,1]
  data$y<-coordinates(obj$observations)[,2]
  depVar = as.character(obj$formulaString[[2]])
  data$z<-obj$observations@data[[depVar]]
  field<-obj$predictionLocations

  manyquantiles=length(calc[names(calc)=="quantile"])
  manyexcprob=length(calc[names(calc)=="excprob"])
  calcnew<-NULL
  if(manyquantiles>0){
    for(i in 1:manyquantiles){
      calcnew$quantiles[i]<-calc[names(calc)=="quantile"][[i]]
    }
  }
  if(testMean==1 && sum(calcnew$quantiles==0.5)>0) testMean<--1

  if(manyexcprob>0){
    for(i in 1:manyexcprob){
      calcnew$excprob[i]<-calc[names(calc)=="excprob"][[i]]
    }
  }
  calcnew$mean<-calc$mean
  calcnew$variance<-calc$variance

  prediction<-NULL
  tempfield<-NULL
  tempfield$x<-coordinates(field)[,1]
  tempfield$y<-coordinates(field)[,2]
  tempfield$F<-model.matrix(delete.response(terms(obj$formulaString)),obj$predictionLocations)
  prediction<-predictioncopula(tempfield,data,estimates,search,calcnew,obj$params$debug.level, obj$params$nclus)

  if(manyquantiles>0){
    for(i in 1:manyquantiles){
      name=paste("quantile",calcnew$quantile[i],sep="")
      prediction[[name]]<-prediction$quantiles[,i]
    }
  }
  prediction$quantiles<-NULL
  if(manyexcprob>0){
    for(i in 1:manyexcprob){
      name=paste("excprob",calcnew$excprob[i],sep="")
      prediction[[name]]<-prediction$excprob[,i]
    }
  }
  prediction$excprob<-NULL

#Test if the mean has reasonable values. If not, calculate the median.
  if(testMean==1){
    index<-!is.na(tempfield$x) & !is.na(tempfield$y)
    med=median(data$z)
    ma=max(data$z)
    mi=min(data$z)
    if((sum(is.na(prediction$mean[index]))>0 || max(prediction$mean[index])>ma+2*(ma-med) 
               || min(prediction$mean[index])<mi-2*(med-mi))){
      warning("Problem in bayescopula. Estimated mean values are nonsensical. Calculating median instead.", 
                    call. = FALSE, immediate. = TRUE)
 	    prediction$mean<-predictioncopula(tempfield,data,estimates,search,list(mean=FALSE,variance=FALSE,quantiles=0.5),obj$params$debug.level)$quantiles
 	    prediction$quantile0.5<-prediction$mean
    }
  } else if (testMean==-1){
  		warning("Problem in bayescopula. Estimated mean values are nonsensical. Calculating median instead.", call. = FALSE, immediate. = TRUE)
  		prediction$mean<-prediction$quantile0.5
  }
  prediction
}



`predictioncopula` <- function(locations,data,estimates,search=10,calc,debug.level, nclus){
  margin<-estimates$margin
  distribfunction<-get(paste("p",margin$name,sep=""),mode="function")
  quantilefunction<-get(paste("q",margin$name,sep=""),mode="function")
  densityfunction<-get(paste("d",margin$name,sep=""),mode="function")

  expected<-estimates$trend$F %*% estimates$trend$params
  if(length(margin$params)==1) newdata<-distribfunction(data$z,expected,margin$params[1]) else 
    newdata<-distribfunction(data$z,expected,margin$params[1],margin$params[2])

  copula<-estimates$copula

  if(copula$method=="chisq"){
    newdata<-qchisq(newdata,1,copula$params)
    ijmatrix=matrix(0,nrow=2^(search+1),ncol=search+1)
    for(j in 1:(search+1)) ijmatrix[,j]=kronecker(t(0:(2^(search-j+2)-1)%%2),matrix(1,nrow=2^(j-1),ncol=1))
    multeps=(-1)^ijmatrix
  } else if(copula$method=="norm") {
    newdata=qnorm(newdata,0,1)
    multeps = NULL
  }

  len<-length(data$x)

  anisotropy<-estimates$anisotropy
  if(!is.null(anisotropy$params)){
    xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),
         -anisotropy$params[2]*sin(anisotropy$params[1]),
         anisotropy$params[2]*cos(anisotropy$params[1])),
         ncol=2,byrow=TRUE) %*% rbind(t(data$x),t(data$y))
    data$x=xy[1,]
    data$y=xy[2,]
    xy<-matrix(c(cos(anisotropy$params[1]),sin(anisotropy$params[1]),
         -anisotropy$params[2]*sin(anisotropy$params[1]),
         anisotropy$params[2]*cos(anisotropy$params[1])),
         ncol=2,byrow=TRUE) %*% rbind(t(locations$x),t(locations$y))
    locations$x=xy[1,]
    locations$y=xy[2,]
  }

  h<-as.matrix(dist(cbind(data$x,data$y)))

  if(calc$mean==T) mean1<-rep(NA,length(locations$x)) else  mean1<-NULL
  if(is.null(calc$excprob))  exceedanceprob<-NULL else   exceedanceprob<-rep(NA,length(calc$excprob))

  numquantiles<-length(calc$quantiles)
  if(numquantiles==0)  quantiles<-NULL else   quantiles<-rep(NA,numquantiles)

  if(calc$variance==FALSE)  variance<-NULL else variance<-rep(NA,length(locations$x))

  correlation<-estimates$correlation

  
#  clusterExport(cl,list("data","debug.level","distribfunction", "quantilefunction","densityfunction, newdata,
#        copula, search, correlation, h, estimates, calc, margin, numquantiles, exceedanceprob,
#        quantiles, multeps))
#  splt = sample(1:nclus, nrow(coordinates(locations)), replace = TRUE)
#  newdlst = lapply(as.list(1:nclus), function(w) locations[splt == w,])

#  res <- do.call("rbind", parLapply(cl, newdlst, function(lst) 
#    intamap:::pfunc(data,locations[i,],debug.level,distribfunction, quantilefunction,densityfunction, newdata,
#        copula, search, correlation, h, estimates, calc, margin, numquantiles, exceedanceprob,
#        quantiles, multeps)))


  if (is.null(nclus)) nclus = 1
  loc = as.data.frame(locations)
  names(loc) = names(locations)
  if (nclus > 1 & dim(loc)[1] > 3) {
    if (!suppressMessages(suppressWarnings(require(doParallel))))
  	  stop("nclus is > 1, but package doParallel is not available")

    locs = vector("list", dim(loc)[1])
    for (i in 1:length(locs)) locs[[i]] = loc[i,]
    cl <- makeCluster(nclus)
    res <- clusterApply(cl, locs, fun = pfunc, data = data, debug.level = debug.level,
      distribfunction = distribfunction, quantilefunction = quantilefunction, 
        densityfunction = densityfunction, newdata = newdata,
          copula = copula, search = search, correlation = correlation, h = h, 
          estimates = estimates, calc = calc, margin = margin, numquantiles = numquantiles, 
          exceedanceprob = exceedanceprob, quantiles = quantiles, multeps = multeps)
    res = matrix(unlist(res), ncol = length(res[[1]]), byrow = TRUE)

 #   if (FALSE) {
#      clusterEvalQ(cl, intamap:::pfunc)
#      registerDoParallel(cl, nclus)
#      i = 1 # just to avoid check warnings
#      res <- foreach(i = 1:length(locations$x), .combine = rbind) %dopar% {
#        intamap:::pfunc(data,loc[i,],debug.level,distribfunction, quantilefunction,densityfunction, newdata,
#            copula, search, correlation, h, estimates, calc, margin, numquantiles, exceedanceprob,
#          quantiles, multeps)
#     } # end of foreach loop
#    }
    stopCluster(cl)
  } else {
    for (i in 1:length(locations$x)) {
      res0 = pfunc(data,loc[i,],debug.level,distribfunction, quantilefunction,densityfunction, newdata,
          copula, search, correlation, h, estimates, calc, margin, numquantiles, exceedanceprob,
          quantiles, multeps)
      if (i ==1) res = res0 else res = rbind(res,res0)
    }
    if (length(locations$x) == 1) res = matrix(res, nrow = 1)
  }   

  numexcprob = length(calc$excprob)
  if (numexcprob > 0) excprob = as.matrix(res[,3:(2+numexcprob)]) else excprob = NULL
  if (numquantiles > 0) quantiles = as.matrix(res[,(3+numexcprob):(2+numexcprob+numquantiles)]) else quantiles = NULL

  list(mean = res[,1],variance = res[,2], excprob = excprob, quantiles = quantiles)
}



`meanintgauss` <- function(loc,sigma,newdata,len,margin,quantilefunction,expected){
  loc[loc==1]=NaN
  loc[loc==0]=NaN
  if(length(margin$params)==1){
    xy<-quantilefunction(loc,expected,margin$params[1])
  } else {
    xy<-quantilefunction(loc,expected,margin$params[1],margin$params[2])
  }
  dens<-xy*condcopuladensgauss(loc,sigma,newdata,len)
  fin<-!is.finite(dens)
  if(sum(fin)>0) dens[fin]=0
  dens
}

`varintgauss` <- function(loc,sigma,newdata,len,margin,quantilefunction,m,expected){
  loc[loc==1]=NaN
  loc[loc==0]=NaN
  if(length(margin$params)==1){
    xy<-quantilefunction(loc,expected,margin$params[1])
  } else {
    xy<-quantilefunction(loc,expected,margin$params[1],margin$params[2])
  }
  dens<-(xy-m)^2*condcopuladensgauss(loc,sigma,newdata,len)
  fin<-!is.finite(dens)
  if(sum(fin)>0) dens[fin]=0
  dens
}

`condcopuladensgauss` <- function(loc,sigma,newdata,len){
  loc=qnorm(loc)
  div=dnorm(loc)
  sigma22inv=solve(sigma[2:len,2:len])
  x=sigma[1,2:len]%*%sigma22inv
  dnorm(loc,x%*%newdata,sqrt(1-x%*%sigma[2:len,1]))/div
}

`condcopuladenschi2` <- function(loc,lambda,sigma,multeps,newdata,search){
  distrib<-NULL
  loc<-qchisq(loc,1,lambda)
  div<-dchisq(loc,1,lambda)
  eps=multeps[1:2^search,1:search]*kronecker(matrix(1,nrow=2^search,ncol=1),t(sqrt(newdata)));
  c=sum(dmvnorm(eps,sqrt(lambda)*matrix(1,nrow=1,ncol=search),sigma[2:dim(sigma)[1],2:dim(sigma)[1]]))
  for(i in 1:length(loc)){
    eps=multeps*kronecker(matrix(1,nrow=2^(search+1),ncol=1),t(c(sqrt(loc[i]),sqrt(newdata))))
    distrib[i]=sum(dmvnorm(eps,sqrt(lambda)*matrix(1,nrow=1,ncol=search+1),sigma))/(2*sqrt(loc[i])*div[i]*c)
  }
  distrib
}

`meanintchi2` <- function(loc,lambda,sigma,multeps,newdata,len,margin,expected,quantilefunction){
  loc[loc==1]=NaN
  loc[loc==0]=NaN
  if(length(margin$params)==1){
    xy<-quantilefunction(loc,expected,margin$params[1])
  } else {
    xy<-quantilefunction(loc,expected,margin$params[1],margin$params[2])
  }
  dens<-xy*condcopuladenschi2(loc,lambda,sigma,multeps,newdata,len)
  fin<-!is.finite(dens)
  if(sum(fin)>0) dens[fin]=0
  dens
}

`findquantiles` <- function(loc,sigma,newdata,len,quant){
  if(loc<1 && loc>0){
    val<-quant-integrate(condcopuladensgauss,0,loc,sigma,newdata,len,subdivisions=100,rel.tol=.Machine$double.eps^0.25,abs.tol=1e-9,stop.on.error=FALSE)$value;
  } else {
    val<-1e10;
  }
  val
}

`exceedance` <- function(loc,sigma,newdata,len,margin,distribfunction,densityfunction,expected){
  loc[loc==1]=NaN
  loc[loc==0]=NaN
  if(length(margin$params)==1){
    dens<-densityfunction(loc,expected,margin$params[1])
    loc<-distribfunction(loc,expected,margin$params[1])
  } else {
    dens<-densityfunction(loc,expected,margin$params[1],margin$params[2])
    loc<-distribfunction(loc,expected,margin$params[1],margin$params[2])
  }
  ex<-condcopuladensgauss(loc,sigma,newdata,len)*dens
  ex[!is.finite(ex)]<-0
  ex
}




pfunc <- function(data,locations,debug.level,distribfunction, quantilefunction,densityfunction, newdata,
        copula, search, correlation, h, estimates, calc, margin, numquantiles, exceedanceprob,
        quantiles, multeps){
  if(!is.nan(locations$x)){
    numpoints <- search + 1
    distances = sqrt((data$x-locations$x)^2+(data$y-locations$y)^2)	#hnew[,i]
    sorted = sort(distances,index.return=TRUE)
    distances <- sorted$x[1:search]
    index <- sorted$ix[1:search]
    if(sum(distances == 0) > 0){
      index <- index[distances>0]
      distances <- distances[distances>0]
      numpoints <- search
    }
    sigma=covar(rbind(c(0,distances),cbind(distances,h[index,index])),
               correlation$params,correlation$model)
    expected<-locations$F %*% estimates$trend$params

#Predictive mean
    if(calc$mean==T || calc$variance==T){
      if(copula$method=="chisq"){
      integr=integrate(meanintchi2,0,1,copula$params,sigma,multeps,newdata[index],
          search,margin,expected,quantilefunction,subdivisions=100,
          rel.tol = .Machine$double.eps^0.25, abs.tol=1e-3,stop.on.error=FALSE)
      mean1<-ifelse(integr$message=="OK", integr$value, NaN)
      } else if(copula$method=="norm"){
        integr=integrate(meanintgauss,0,1,sigma,newdata[index],numpoints,
            margin,quantilefunction,expected,subdivisions=100,
            rel.tol = .Machine$double.eps^0.25, abs.tol=1e-3,stop.on.error=FALSE)
        mean1<-ifelse(integr$message=="OK", integr$value, NaN)
      }
    }

#Predictive quantiles
    if(numquantiles>0){
      if(copula$method=="norm"){
      	if(numquantiles<=9){
	      	test<-seq(from=0,to=1,len=numquantiles+2)
	      	test<-test[2:(numquantiles+1)]
      	} else {
      		test<-seq(from=0,to=1,len=11)
	      	test<-test[2:10]
	      }
        val<-NULL
        integr<-NULL
        for(zz in 1:length(test)){
	        integr[[zz]]<-integrate(condcopuladensgauss,0,test[zz],sigma,newdata[index],
             numpoints,subdivisions=100, rel.tol = .Machine$double.eps^0.25, abs.tol=5e-3,stop.on.error=FALSE)
        }
        for(j in 1:numquantiles){
          for(zz in 1:length(test)){
          	val[zz]<-if(integr[[zz]]$message=="OK"){calc$quantiles[j]-integr[[zz]]$value}else{Inf}
          }
          zz=which.min(val^2)
          q<-test[zz]

          if(length(test)==1){
          	if(val[zz]>=0){
	          	q1<-q
          		f1<-val[zz]
          		q2<-1
          		f2<-calc$quantiles[j]-1
	          } else {
          		q1<-0
	          	f1<-calc$quantiles[j]
	   	        q2<-q
    	      	f2<-val[zz]
  	        }
          } else {
            if(zz!=length(test) && zz!=1){
            	if(val[zz]>=0){
            		q1<-q
            		f1<-val[zz]
            		q2<-test[zz+1]
            		f2<-val[zz+1]
            	} else {
            		q1<-test[zz-1]
            		f1<-val[zz-1]
            		q2<-q
            		f2<-val[zz]
            	}
            } else {
            	if(zz==1){
            		if(val[zz]>=0){
            			q1<-q
            			f1<-val[zz]
            			q2<-test[zz+1]
            			f2<-val[zz+1]
            		} else {
            			q1<-0
            			f1<-calc$quantiles[j]
            			q2<-q
            			f2<-val[zz]
            		}
            	} else {
            		if(val[zz]>=0){
             			q1<-q
            			f1<-val[zz]
            			q2<-1
            			f2<-calc$quantiles[j]-1
             		} else {
            			q1<-test[zz-1]
             			f1<-val[zz-1]
            			q2<-q
            			f2<-val[zz]
            		}
            	}
            }
          }
          if(sign(f1)!=-sign(f2)){
          	q1<-0
          	q2<-1
          	f1<-calc$quantiles[j]
          	f2<-calc$quantiles[j]-1
          }
          quantiles[j]<-uniroot(findquantiles,c(q1,q2),f.lower=f1,f.upper=f2,tol=.Machine$double.eps^0.4,sigma=sigma,newdata=newdata[index],len=numpoints,quant=calc$quantiles[j])$root
        }
        if(length(margin$params)==1){
          quantiles<-quantilefunction(quantiles,expected,margin$params[1])
        } else {
          quantiles<-quantilefunction(quantiles,expected,margin$params[1],margin$params[2])
        }
      } #TODO: copula$method="chisq"
    }

#Exceedance above Thresholds
    if(!is.null(calc$excprob)){
      if(copula$method=="norm"){
        for(j in 1:length(calc$excprob)){
          integr<-integrate(exceedance,calc$excprob[j],Inf,sigma,newdata[index],numpoints,margin,distribfunction,densityfunction,expected,subdivisions=100,rel.tol=.Machine$double.eps^0.25,abs.tol=1e-9,stop.on.error=FALSE)
          exceedanceprob[j]<-if(integr$message=="OK"){integr$value}else{NaN}
        }
      } #TODO: copula$method="chisq"
    }

#Prediction Variance
    if(calc$variance==TRUE){
      if(copula$method=="norm"){
        integr=integrate(varintgauss,0,1,sigma,newdata[index],numpoints,margin,quantilefunction,mean1,expected,subdivisions=100, rel.tol = .Machine$double.eps^0.25, abs.tol=1e-5,stop.on.error=FALSE)
        variance<-if(integr$message=="OK"){integr$value}else{NaN}
      } #TODO: copula$method="chisq"
    }
  }
  return(c(mean1,variance,exceedanceprob,quantiles))
}
