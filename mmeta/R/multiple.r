######################################################################################
### Purpose: Wrapper function to create object "multipletables"
### Input:   data.frame(y1,n1,y2,n2),other(measure, model, method,alpha, nsam)
### Output:  S3 object "multipletables", refer to the help file of multipletables
### Note:    Implement is function "multipletables_sar"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################

multipletables <- function(data=NULL, measure=NULL, model="Sarmanov",
                           method="sampling", nsam=10000, alpha=0.05) {
 
  if (is.null(measure)) stop("measure is missing")
  if (is.null(data)) stop("data is missing")
  y1 <- data$y1; y2 <- data$y2
  n1 <- data$n1; n2 <- data$n2
  studynames <- data$studynames

  if(length(y1)<2 |length(n1)<2 |length(y2)<2 |length(n2)<2 )
    stop ("only for multiple tables analysis \n")

  if (is.null(studynames)) {
    studynames<-seq(1:length(y1))
    studynames<-as.character(studynames)
  }
 
  if(!(length(y1)==length(y2) & length(y1)==length(n1) & length(y1)==length(n2)
       & length(y2)==length(n2) & length(y2)==length(n1) & length(n2)==length(n1)
       & length(n2)==length(studynames)))
    stop ("length of y1,n1,y2,n2 or studynames are non-conformable. \n")

  if(all(all(n1<=0),all(n2<=0)))
    stop("n1,n2 should be not be less or euqual to 0")
		
  if(all(all(y1<0),all(y2<0))) 
    stop("y1,n1,y2,n2 should be greater than 0")
		
  if (all(all(trunc(y1)!=y1),all(trunc(n1)!=n1),all(trunc(y2)!=y2),all(trunc(n2)!=n2)))
    stop("y1,n1,y2,n2 should be interger")
  if(all(y1>n1)) stop("y1 should be less than n1")
  if(all(y2>n2)) stop("y2 should be less than n2")
  if (any(y1==n1)){
    index <- which(y1==n1)
    for(i in 1:length(index)) {
      y1[index[i]] <- y1[index[i]]-0.02
    }
  }

  if (any(y2==n2)){##if y1=n2, the bugs will down
    index <- which(y2==n2)
    for(i in 1:length(index)){
      y2[index[i]] <- y2[index[i]]-0.02
    }
  }

  if (measure=="RD" & method=="exact"){
    cat("Only sampling based method available for RD \n")
    method <- "sampling"
  }

  out <- multipletables_sar(y1=y1,n1=n1,y2=y2,n2=n2,studynames=studynames,measure=measure,
                            model=model,method=method,nsam=nsam,alpha=alpha)
  invisible(out)
}

######################################################################################
### Purpose: Create object "multipletables"
### Input:   data(y1,n1,y2,n2),other(measure, model, method,alpha, nsam)
### Output:  S3 object "multipletables", refer to the help file of multipletables
### Note:    Wrapper is function "multipletables"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu 
### Data:    7/13/2012
######################################################################################
multipletables_sar <- function(y1=y1,n1=n1,y2=y2,n2=n2,studynames=studynames,measure=measure,
                               model=model,method=model,nsam=nsam,alpha=alpha) {
  nstudy <- length(y1)
  ##parameters estimates
  mypar <- MLE.function(y1=y1,n1=n1,y2=y2,n2=n2,model=model)
  MLE <- mypar$MLE
  a1 <- MLE[1]
  b1 <- MLE[2]
  a2 <- MLE[3]
  b2 <- MLE[4]
  rho <- MLE[5]  
  chi2 <- mypar$chi2; pvalue<-mypar$pvalue;
  hessian <- mypar$hessian
  cov.matrix <- inverse.matrix.func(-hessian)

  if (model=="Sarmanov") {
    colnames(hessian) <- colnames(cov.matrix)<-c("loga1","logb1","loga2","logb2","eta")
    rownames(hessian) <- rownames(cov.matrix)<-c("loga1","logb1","loga2","logb2","eta")
  }
  if (model=="Independent") {
    colnames(hessian) <- colnames(cov.matrix)<-c("loga1","logb1","loga2","logb2")
    rownames(hessian) <- rownames(cov.matrix)<-c("loga1","logb1","loga2","logb2")
  }
  #estimate overall measure
  myoverall <- overall(MLE,hessian,measure,model,alpha=alpha)                                 
  #study-specific sample
  mysample <- list()
  for(i in 1:nstudy){
    temp <- sampling(a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,n1=n1[i],
                              y1=y1[i],n2=n2[i],y2=y2[i], model=model,measure=measure,
                              nsam=nsam)
	upper=quantile(temp,prob=0.9995,na.rm=T)
    mysample[[i]] =temp[temp<=upper]  
	}						  
   priordsample <- sampling(a1=a1,b1=b1,a2=a2,b2=b2,rho=rho,n1=0,
                           y1=0,n2=0,y2=0, model=model,measure=measure,
                           nsam=nsam)
	upper=quantile(temp,prob=0.9995,na.rm=T)
     priordsample=temp[temp<=upper]  
				   
   dens <- list()	 
  ###compute the exact density, if method="exact"
  if (method=="exact") {
    range.overlap.data <- matrix(NA,ncol=2,nrow=nstudy)
    study.min <- sapply(mysample,quantile,probs=0.025,na.rm=TRUE)
    study.max<-sapply(mysample,quantile,probs=0.975,na.rm=TRUE)                                          
    xmin <- min(study.min)
    xmax <- max(study.max)
    for (i in 1:nstudy) {
      dens[[i]] <- dens.post(y1=y1[i],n1=n1[i],y2=y2[i],n2=n2[i],a1=a1,b1=b1,a2=a2,b2=b2,
                             rho=rho,grid.start=xmin,grid.end=xmax,grid.num=1000,
                             measure=measure,model=model)
    }
    priordens <- dens.post(y1=0,n1=0,y2=0,n2=0,a1=a1,b1=b1,a2=a2,b2=b2,
                           rho=rho,grid.start=xmin,grid.end=xmax,grid.num=1000,
                           measure=measure,model=model)
	}
	
    ###Compute the empirical density: if method="sampling"
    if(method=="sampling") {
      if (measure=="RD"){
        dens <- lapply(mysample,density,from=-1,to=1,n=3072)
        priordens <- density(priordsample,from=-1,to=1,n=3072)
      }
      if (measure=="OR" | measure=="RR") {
        dens <- lapply(mysample,density,from=0,n=3072)
        priordens <- density(priordsample,from=0,n=3072)
      }
    }
    ###create object "multiple"
    #data
    dataset <- matrix(NA,ncol=length(y1),nrow=4)
    dataset[1,] <- y1
    dataset[2,] <- n1
    dataset[3,] <- y2
    dataset[4,] <- n2
    rownames(dataset) <- c("y1","n1","y2","n2")

    if (measure=="OR") measurename <- "Odds ratio"
    if (measure=="RR") measurename <- "Relative risk"
    if (measure=="RD") measurename <- "Risk difference"

    out <- list(measure=measure,model=model,dataset=dataset,studynames=studynames,
                measurename=measurename,alpha=alpha,chi2=chi2,pvalue=pvalue,
                MLE=MLE,hessian=hessian,cov.matrix=cov.matrix,overall=myoverall,
                sample=mysample,density=dens,priordens=priordens,priordsample=priordsample)
    class(out) <- "multipletables"
    invisible(out)
}
