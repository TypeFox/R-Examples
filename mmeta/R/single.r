######################################################################################
### Purpose: Create object "singletable"
### Input:   data(y1,n1,y2,n2), hyperparemeters(a1,a2,b1,b2,rho),
###          comaparative measure, model, method, significance level (alpha),
###          number of samples (nsam)
### Output:  S3 object "singletable", refer to the help file of singletable
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################
singletable <- function(y1=y1,n1=n1,y2=y2,n2=n2,measure=measure,model="Sarmanov",
                        method="exact",a1=0.5,b1=0.5,a2=0.5,b2=0.5,rho=0,alpha=0.05,
                        nsam=10000) {
  
  if (measure=="RD"& method=="exact") {
    print("only sampling based mehtod is available for RD")
    method <- "sampling"
  }

  if(length(y1)>=2 |length(n1)>=2 |length(y2)>=2 |length(n2)>=2 )
    stop ("only for single table analysis \n")

  if(y1<0|y2<0) stop("y1,y2 should be greater than 0")
  if(n1<=0|n2<=0) stop("n1,n2 should not be less or equal to 0")
   
  if ((trunc(y1)!=y1)|(trunc(n1)!=n1)|(trunc(y2)!=y2)|(trunc(n2)!=n2))
    stop("y1,n1,y2,n2 should be interger")
  
  if(y1>n1) stop("y1 should be less than n1")
  if(y2>n2) stop("y2 should be less than n2")

  if (any(y1==n1)) {
    index=which(y1==n1)
    for(i in 1:length(index)) {
      y1[index[i]]=y1[index[i]]-0.02
    }
  }

  if (any(y2==n2)) {##if y1=n2, the bugs will down
    index=which(y2==n2)
    for(i in 1:length(index)) {
      y2[index[i]]=y2[index[i]]-0.02
    }
  }	   
   
  # check the range of parameters
  if(model=="Sarmanov") {
    cc <- sqrt(a1*a2*b1*b2)/sqrt((a1+b1+1)*(a2+b2+1))
    upper.bound <- cc/max(a1*b2, a2*b1)
    lower.bound <- -cc/max(a1*a2, b1*b2)
    rho.range<-c(lower.bound,upper.bound)
    names(rho.range)<-c("lower.bound","upper.bound")
    cat("Range of corelation (rho)",fill=T)
    cat("Lower bound:", rho.range[1],fill=T)
    cat("Upper bound:",rho.range[2],fill=T)
    if (rho > upper.bound | rho < lower.bound) stop(paste("rho is out of bound: ",
    lower.bound, upper.bound))
  }
  
  #save parameter
  if(model=="Independent") {
    parameter <- c(a1,b1,a2,b2)
    names(parameter) <- c("a1","b1","a2","b2")
  }

  if(model=="Sarmanov") {
    parameter<-c(a1,b1,a2,b2,rho)
    names(parameter) <- c("a1","b1","a2","b2","rho")
  }
                 
  # save data
  dataset <- c(y1,n1,y2,n2)
  names(dataset) <- c("y1","n1","y2","n2")

  if (measure=="OR") measurename<-"Odds ratio"
  if (measure=="RR") measurename<-"Relative risk"
  if (measure=="RD") measurename<-"Risk difference"

  mysample <- list()
  temp <- sampling(a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,n1=n1,y1=y1,
                            n2=n2,y2=y2,measure=measure,model=model,nsam=nsam)
  upper=quantile(temp,prob=0.9995,na.rm=T)
  mysample[[1]]=temp[temp<=upper]  
  temp<- sampling(a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,n1=0,y1=0,n2=0,
                            y2=0,measure=measure,model=model,nsam=nsam)
  upper=quantile(temp,prob=0.9995,na.rm=T)
  mysample[[2]]=temp[temp<=upper]
  
  dens <- list()
  if (method=="sampling") {
  #emprial density
    if (measure!="RD") dens <- lapply(mysample,density,from=0,n=2048)
    if (measure=="RD") dens <- lapply(mysample,density,from=-1,to=1,n=2048)
  }

  if (method=="exact") {
    xmin.post <- quantile(mysample[[1]],prob=0.0025,na.rm=T)
    xmax.post <- quantile(mysample[[1]],prob=0.9975,na.rm=T)
    dens[[1]] <- dens.post(y1=y1, n1=n1, y2=y2, n2=n2, a1=a1, b1=b1, a2=a2, b2=b2,
                           rho=rho,grid.start=xmin.post,grid.end=xmax.post,grid.num=1000,
                           measure=measure,model=model)
    dens[[2]] <- dens.post(y1=0, n1=0, y2=0, n2=0, a1=a1, b1=b1, a2=a2, b2=b2,
                           rho=rho, grid.start=xmin.post, grid.end=xmax.post, grid.num=1000,
                           measure=measure,model=model)
  }
  names(dens) <- c("Posterior","Prior")
  studynames <- c("Posterior","Prior")
  result <- list(measure=measure,model=model,method=method,dataset=dataset,
                 parameter=parameter,alpha=alpha,density=dens,sample=mysample,
                 studynames=studynames,measurename=measurename)

  class(result) <- "singletable"
  invisible(result)
}     
