#' @name ResamplingVariance
#' @aliases ResamplingVariance
#' @title Resampling variance of randomized response models
#' 
#' @description To estimate the variance of the randomized response estimators using resampling methods.
#' 
#' @usage ResamplingVariance(output,pi,type=c("total","mean"),option=1,N=NULL,pij=NULL,str=NULL,
#' clu=NULL,srswr=FALSE)
#' @param output output of the qualitative or quantitative method depending on the variable of interest
#' @param pi vector of the first-order inclusion probabilities. By default it is NULL
#' @param type the estimator type: total or mean 
#' @param option method used to calculate the variance (1: Jackknife, 2: Escobar-Berger, 3: Campbell-Berger-Skinner). By default it is 1 
#' @param N size of the population
#' @param pij matrix of the second-order inclusion probabilities. This matrix is necessary for the Escobar-Berger and Campbell-Berger-Skinner options. By default it is NULL
#' @param str strata ID. This vector is necessary for the Jackknife option. By default it is NULL
#' @param clu cluster ID. This vector is necessary for the Jackknife option. By default it is NULL
#' @param srswr variable indicating whether sampling is with replacement. By default it is NULL
#' 
#' @details 
#' Functions to estimate the variance under stratified, cluster and unequal probability sampling by resampling methods (Wolter, 2007). 
#' The function ResamplingVariance allows us to choose from three models:
#' 
#' - The Jackknife method (Quenouille, 1949)
#' 
#' - The Escobar-Berger method (Escobar and Berger, 2013)
#' 
#' - The Campbell-Berger-Skinner method (Campbell, 1980; Berger and Skinner, 2005).
#' 
#' The Escobar-Berger and Campbell-Berger-Skinner methods are implemented using the functions defined in samplingVarEst package:
#' 
#' VE.EB.SYG.Total.Hajek, VE.EB.SYG.Mean.Hajek; 
#' 
#' VE.Jk.CBS.SYG.Total.Hajek, VE.Jk.CBS.SYG.Mean.Hajek 
#' 
#' (see López, E., Barrios, E., 2014, for a detailed description of its use).
#' 
#' Note: Both methods require the matrix of the second-order inclusion probabilities. When this matrix is not an input, the program will give a warning and, by default, a jackknife
#' method is used.
#' 
#' @return The resampling variance of the randomized response technique
#' 
#' @references Berger, Y.G., Skinner, C.J. (2005).
#' \emph{A jackknife variance estimator for unequal probability sampling.}
#'  Journal of the Royal Statistical Society B, 67, 79-89.
#' 
#' @references Campbell, C. (1980).
#' \emph{A different view of finite population estimation.}
#'  Proceedings of the Survey Research Methods Section of the American Statistical Association, 319-324.
#'   
#' @references Escobar, E.L., Berger, Y.G. (2013).
#' \emph{A new replicate variance estimator for unequal probability sampling without replacement.}
#'  Canadian Journal of Statistics 41, 3, 508-524.
#'   
#' @references López, E., Barrios, E. (2014).
#' \emph{samplingVarEst: Sampling Variance Estimation.}
#' R package version 0.9-9. Online http://cran.r-project.org/web/packages/survey/index.html 
#'
#' @references Quenouille, M.H. (1949).
#' \emph{Problems in Plane Sampling.}
#' The Annals of Mathematical Statistics 20, 355-375.
#' 
#' @references Wolter, K.M. (2007).
#' \emph{Introduction to Variance Estimation.}
#'  2nd Edition. Springer.
#' 
#' @seealso \code{\link{Warner}}
#' @seealso \code{\link{ChaudhuriChristofides}}
#' @seealso \code{\link{EichhornHayre}}
#' @seealso \code{\link{SoberanisCruz}}
#' @seealso \code{\link{Horvitz}}
#' 
#' @keywords Randomized_response Resampling Variance Jackknife Escobar_Berger Campbell_Berger_Skinner
#' @examples
#' N=417
#' data(ChaudhuriChristofidesData)
#' dat=with(ChaudhuriChristofidesData,data.frame(z,Pi))
#' mu=c(6,6) 
#' sigma=sqrt(c(10,10))
#' cl=0.95
#' data(ChaudhuriChristofidesDatapij)
#' out=ChaudhuriChristofides(dat$z,mu,sigma,dat$Pi,"mean",cl,pij=ChaudhuriChristofidesDatapij)
#' out
#' ResamplingVariance(out,dat$Pi,"mean",2,N,ChaudhuriChristofidesDatapij)
#' 
#' #Resampling with strata
#' data(EichhornHayreData)
#' dat=with(EichhornHayreData,data.frame(ST,z,Pi))
#' mu=1.111111
#' sigma=0.5414886
#' cl=0.95
#' out=EichhornHayre(dat$z,mu,sigma,dat$Pi,"mean",cl)
#' out
#' ResamplingVariance(out,dat$Pi,"mean",1,str=dat$ST)
#' 
#' #Resampling with cluster 
#' N=1500
#' data(SoberanisCruzData)
#' dat=with(SoberanisCruzData, data.frame(CL,z,Pi))
#' p=0.7
#' alpha=0.5
#' cl=0.90
#' out=SoberanisCruz(dat$z,p,alpha,dat$Pi,"total",cl)
#' out
#' ResamplingVariance(out,dat$Pi,"total",2,N,samplingVarEst::Pkl.Hajek.s(dat$Pi))
#' 
#' #Resampling with strata and cluster
#' N=1500
#' data(HorvitzDataStCl)
#' dat=with(HorvitzDataStCl, data.frame(ST,CL,z,Pi))
#' p=0.6
#' alpha=0.5
#' cl=0.95
#' out=Horvitz(dat$z,p,alpha,dat$Pi,"mean",cl,N)
#' out
#' ResamplingVariance(out,dat$Pi,"mean",3,N,samplingVarEst::Pkl.Hajek.s(dat$Pi))
#'
#' @export
ResamplingVariance=function(output,pi,type=c("total","mean"),option=1,N=NULL,pij=NULL,str=NULL,clu=NULL,srswr=FALSE){
  
  if(!is.list(output)){stop("output must be a list.")}
  if(any(is.na(output))){stop("There are missing values in output.")}
  
  if(!is.character(type)){stop("type must be a character")}
  if((type!="total")&(type!="mean")){stop("The value of type must be total or mean.")}                
  
  if((option<1)|(option>3)){stop("The value of option must be between 1 and 3.")}
  
  if(!is.logical(srswr)){srswr=as.logical(srswr)} 
  if((srswr!=TRUE)&(srswr!=FALSE)){stop("The value of srswr must be TRUE or FALSE")}
  
  if((is.null(pij))&(option!=1)){
    warning("If you do not introduce pij you must choose the jackknife method. The output given is done with the jackknife method without strata nor cluster.")
    option=1
  }
  
  if(is.null(N)){
    N=round(sum(1/pi))
  }

  
  if(option==1){
    strata=FALSE
    if(!is.null(str)){
      strata=TRUE
      if(!is.factor(str)){str=as.factor(str)}
    }
    cluster=FALSE
    if(!is.null(clu)){
      cluster=TRUE
      if(!is.factor(clu)){clu=as.factor(clu)}
    }

    if((strata==FALSE)&(cluster==FALSE)){
      Rve=JackknifeVariance(length(pi),output$TransformedVariable,pi)
    }
    
    if((strata==TRUE)&(cluster==FALSE)){
     nh=table(str)
     ve=vector()
     for(i in levels(str)){
       ve[i]=JackknifeVariance(nh[i],output$TransformedVariable[str==i],pi[str==i],strata)
     }
     Rve=sum(ve)
    }
   
    if((strata==FALSE)&(cluster==TRUE)){
     Rve=JackknifeVariance(length(levels(clu)),output$TransformedVariable,pi,cluster=cluster,clu=clu)
    }
    
    if((strata==TRUE)&(cluster==TRUE)){
     ve=vector()
     for(i in levels(str)){
       cl_sti=droplevels(clu[str==i])         
       ve[i]=JackknifeVariance(nlevels(cl_sti),output$TransformedVariable[str==i],pi[str==i],strata,cluster,cl_sti)
     }
     Rve=sum(ve)
    }  
    
    if(type=="mean"){
      if(length(N)!=1){stop("N must be a scalar.")}   
      if(N<0){stop("N must be a positive number.")}
      
      Rve=Rve/N^2
    }
  }
    
  if(option==2){
    Rve=EscobarBergerVariance(N,output$TransformedVariable,pi,pij,type)
  }

  if(option==3){
    Rve=CampbellBergerSkinnerVariance(N,output$TransformedVariable,pi,pij,type)
  }
  
  if(srswr==TRUE){
    Rve=Rve/(1-mean(pi))  
  }

  return(Rve)
}