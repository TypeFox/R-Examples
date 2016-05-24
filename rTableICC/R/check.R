check <-
function(p, theta, M, N=NULL, K=NULL, row.margins=NULL, col.margins=NULL, lambda=NULL, sampling, 
                           ICC=TRUE, structure="2x2xK"){

  if ((sampling!="Multinomial") & (sampling!="Product") & (sampling!="Poisson")) {
    stop("The parameter 'sampling' can take 'Product' for product multinomial sampling, 'Multinomial' for multinomial sampling, and 'Poisson' for Possion sampling schemes!")
  }
  
  if (sampling=="Poisson"){
    if ((is.finite(lambda)==FALSE) & (min(lambda)<0)){
      stop("Mean number of individuals must be either a positive scalar or a positive vector under Poisson samlping plan!")     
    }
  }

  if (ICC==TRUE){
    if ((is.na(theta)==TRUE) | (is.finite(theta)==FALSE) | (min(theta)<0)) {
      stop("Pre-determined intraclass correlations must be entered as a positive value!")
    }
    if (structure=="2x2xK"){
      if ((is.na(p)==TRUE) | (is.finite(p)==FALSE) | (length(dim(p))!=2) |  (ncol(p)!=4) | (min(p)<=0) | (max(p)>=1) | (sum(p)!=1)){
        stop("Parameter p must be entered as a Kx4 finite matrix with elements within the interval [0,1] and sum up to 1!")
      }
      if ((is.na(M)==TRUE) | (is.finite(M)==FALSE) | (min(M)<0)) {
        stop("Number of clusters under each center (M>0) must be entered!")
      }       
    }else if (structure=="RxC"){
      if ((length(M)!=1) | (is.finite(M)==FALSE) | (min(M)<=0)){
        stop("Total number of clusters (M>0) must be a scalar value for RxC tables!")
      }       
    }    
    
  }else if (ICC==FALSE){
    if (structure=="2x2xK"){
      if (sampling=="Multinomial"){ 
        if ((length(N)!=1) | (is.finite(N)==FALSE)){ 
          stop("Total number of observation should be entered as a scalar under multinomial samlping plan.")      
        }
        if (((is.na(p)==TRUE) | (is.finite(p)==FALSE) | (length(dim(p))!=3) | (dim(p)[1]!=2) | (dim(p)[2]!=2) | (min(p)<=0) | (max(p)>=1))){
          stop("Parameter p must be entered as a 2x2xK finite matrix with elements within the interval [0,1]!")
        }
      } 
      if (sampling=="Poisson"){
        if (length(dim(lambda))==3){
          if ((is.na(lambda)==TRUE) | (is.finite(lambda)==FALSE) |  (dim(lambda)[1]!=2) | (dim(lambda)[2]!=2) | (min(lambda)<=0)){
            stop("Parameter lambda must be entered as a 2x2xK finite matrix!")
          }
        }else if (is.null(dim(lambda))==TRUE){
          if ((is.na(lambda)==TRUE) | (is.finite(lambda)==FALSE) | (is.null(K)==TRUE)| (min(lambda)<=0)){
            stop("The scalar entered as the value of parameter lambda must be a finite scalar and number of centers K must be entered!")
          }
        } else {          
          stop("Parameter lambda must be either a 2x2xK finite matrix or a scalar!")
        }      
      }
      if (sampling=="Product"){
        if ((is.na(N)==TRUE) | (is.finite(N)==FALSE) | (min(N)<0)) {
          stop("Center margins must be entered as scalar positive values!")
        } else if (length(N)<1){
          stop("At least one center must be entered!")
        } 
        if (((is.na(p)==TRUE) | (is.finite(p)==FALSE) | (length(dim(p))!=3) | (dim(p)[1]!=2) | (dim(p)[2]!=2) | (min(p)<=0) | (max(p)>=1))){
          stop("Parameter p must be entered as a 2x2xK finite matrix with elements within the interval [0,1]!")
        }
      }
    } else if (structure=="RxC"){
      if (sampling=="Multinomial"){ 
        if ((length(N)!=1) | (is.finite(N)==FALSE)){ 
          stop("Total number of observation should be entered as a scalar under multinomial samlping plan.")      
        }
        if ((is.na(p)==TRUE) | (is.finite(p)==FALSE) | (length(dim(p))!=2) | (min(p)<=0) | (max(p)>=1)){
          stop("Parameter p must be entered as an RxC finite matrix with elements within the interval [0,1]!")
        }      
        if (min(dim(p))<=1){
          stop("Minimum colum or row length of p must be greater than 1!")
        }
      }
      if (sampling=="Product"){
        if ((is.null(row.margins)==FALSE) & (is.null(col.margins)==FALSE)){
          stop("Row or columun margins must be a scalar vector under product multinomial sampling plan!")
        } else if ((length(row.margins)>1) & (length(col.margins)>1)){
          stop("Number of either row.margins or columun margins must be greater than one. Both cannot be greater than one!")
        } else if ((length(row.margins)<2) & (length(col.margins)<2)){
          stop("At least number of fixed row or columun margins must be greater than one under product multinomial sampling!")
        }
        if ((is.na(p)==TRUE) | (is.finite(p)==FALSE) | (length(dim(p))!=2) | (min(p)<=0) | (max(p)>=1)){
          stop("Parameter p must be entered as an RxC finite matrix with elements within the interval [0,1]!")
        }  
      }
      if (sampling=="Poisson"){
        if ((is.na(lambda)==TRUE) | (is.finite(lambda)==FALSE) | (length(dim(lambda))!=2) | (min(lambda)<=0)){
          stop("Parameter lambda must be entered as an RxC finite matrix!")
        }  
      }      
    }
  }
}
