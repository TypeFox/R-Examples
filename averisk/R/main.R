#' Calculate average attributable fractions
#'
#' This function calculates average attributable fractions and confidence intervals for discrete riskfactors
#' @param f A formula object specifying a logistic regression model for the risk factors.  See Details for more information.
#' @param the.data A dataframe containing the disease indicator, riskfactors and possible confounders.
#' @param ref_cat A character vector indicating the reference values for each risk factor. Defaults to \code{NULL}. In the case that this argument is not specified, the default R assignment for the reference level is used for risk factors coded as \code{character} and \code{factor} variables (see \code{?factor}) whereas the minimum value is assigned as the reference level for risk factors coded as \code{numeric} vectors. Attributable fractions calculate the proportional change in disease prevalence that might be expected if the entire population had the reference value for each risk factor.
#' @param refy Response value that specifies controls.  Defaults to \code{0}.
#' @param cat_confounders A character vector indicating categorical variables that need to be adjusted for (excluding risk factors).  Defaults to \code{NULL}.
#' @param cont_confounders A character vector indicating categorical variables that need to be adjusted for (excluding risk factors).  Defaults to \code{NULL}.
#' @param prev A proportion specifying the percentage of the population that have the disease.  Defaults to \code{NA}.  \code{NA} is appropriate for prospective designs. Prevalence should be specified for case control designs.
#' @param allperm \code{TRUE} gives an exact calculation of the sample average fraction.  \code{FALSE} gives an approximate calculation.  Defaults to \code{TRUE}.
#' @param nsample_perm How many permutations are used when calculating approximate average attributable fractions? (only necessary to specify when \code{allperm=FALSE}).  If \code{approx_error} is specified, then \code{nsample_perm} is the number of sampled permutations for the construction of the standard error (if \code{ci=TRUE}), but is overrided in the construction of the point estimate.
#' @param approx_error Specifying this option will calculate average fractions using the number of permutations necessary to approximate the true sample average fraction point estimate to within approx_error.  Defaults to \code{NA}. 
#' @param ci Is a confidence interval required? Defaults to \code{FALSE}.  
#' @param conf_level The confidence level specified as a proportion.  i.e. \code{conf_level=0.95} would imply 95 percent confidence. Only necessary to specify when \code{ci=TRUE}
#' @param nsample_var The number of monte carlo iterates of the average fraction that are used when calculating the confidence interval.  Only necessary to specify when \code{ci=TRUE} 
#' @param correction_factor Whether an extra correction term is subtracted from the estimated standard error due to the monte carlo simulated AFs and the point estimate for the AF being based on different numbers of permutations.  Defaults to \code{TRUE}. (Only necessary to specify when \code{ci=TRUE} and \code{allperm=FALSE})
#' @param quantile_int confidence interval is calculated from empirical quantiles of the monte carlo simulated AFs.  Defaults to \code{FALSE}.   Only necessary to specify when \code{ci=TRUE}.  See Details for more information.
#' @param sep_est The point estimate of the AF in a separate calculation to the confidence interval.  Defaults to \code{TRUE}   (Only necessary to specify when \code{ci=TRUE})
#' @keywords getAF 
#' @details The model formula \code{f} is specified using traditional R notation.  For instance, in the situation where the binary response is 'y' and there are 3 risk factors: x1, x2 and x3, \code{f} would be specified as \code{y~x1+x2+x3}.  Interactions are not permitted in the model.  Confounders (either categorical or continuous) are added as separate arguments in character vectors.  When a confidence interval is requested, a symmetric interval around the point estimate AF is given by default.  In the case that \code{nsample_var} is large, a possibly assymetric interval may instead be requested using \code{quantile_int=TRUE}.  In this case, the interval is generated from percentiles of the Monte Carlo simulates of the AF.  Since estimating percentiles of a distribution is more difficult than estimating the overall variance, \code{nsample_var} should be increased if quantile based confidence intervals are desired, and doing so will significantly increase run-time.  Quantile based confidence intervals maybe superior to symmetric intervals when the distribution of the simulated AFs is skewed.  However, in our experience, the distribution of Monte Carlo simulates is usually relatively symmetric. 
#' @return If \code{ci=TRUE}, a 3 x (K+1) matrix, where K is the number of risk factors.  The first row represents the point estimate for the AF, and the second and third rows lower and upper confidnece bounds.  If \code{ci=FALSE}, a (K+1)-dimensional vector with the calculated point estimates for the AF is returned. 
#' @importFrom MASS mvrnorm
#' @importFrom stats as.formula glm model.frame predict qt quantile sd var weighted.mean
#' @importFrom utils flush.console
#' @author John Ferguson (john.ferguson@@nuigalway.ie)
#' @export 
#' @references Eide, Geir Egil and Olaf Gefeller.  Sequential and Average Attributable Fractions as Aids in the Selection of Preventative Strategies.  J. Clinical Epidemiology, Vol. 48, No. 5, pp. 645-655, 1995.
#' @examples
#' # the following example is from Eide and Gefeller, 1995
#' # simulate data
#' 
#' ex_probs <- c(.06732,.02976,.01570,.01787,.01445,.01008,.06986,.06553,.03,.05766,
#'              .09680,.04194,.02741,.02194,.02474,.01031,.12410,.09537,.08408,.09509) # P(E|D)
#' disease_probs <- c(.036,.0621,.0236,.0411,.0507,.0864,.1066,.1745,.1867,.2891,.0514,
#'              .0875,.0339,.0584,.0718,.1206,.1474,.2345,.2497,.3708) # P(D|E)
#' pe <- ex_probs/disease_probs ## marginal P(E)
#' pe <- pe/sum(pe)
#' nond_exposure_probs <- (1-disease_probs)*pe  # P(E|not D)
#' nond_exposure_probs <- nond_exposure_probs/sum(nond_exposure_probs)
#' ex_probs <- ex_probs/sum(ex_probs)
#' the.mat <- cbind(c(rep(0,10),rep(1,10)),rep(rep(1:5,each=2),2),rep(c(0,1),10))
#' ncase <- 500
#' ncontrol <- 500
#' casemat <- the.mat[sample(1:20,size=ncase,replace=TRUE,prob=ex_probs),]
#' case_rows <- cbind(rep(1,ncase),casemat)
#' controlmat <- the.mat[sample(1:20,size=ncase,replace=TRUE,prob=nond_exposure_probs),]
#' control_rows <- cbind(rep(0,ncontrol),controlmat)
#' the.d <- rbind(case_rows,control_rows)
#' colnames(the.d) <- c("y","urban.rural","smoking.category","occupational.exposure")
#' 
#' # Just get the estimate (no confidence interval)
#' getAF(y~urban.rural+smoking.category+occupational.exposure,the.d,prev=0.09)
#' 
#' ## other examples -  uncomment these lines to run them.
#' ## find the average fraction and associated monte-carlo calculated 99% confidence 
#' ##  No need for approximation here.  Assume population prevalence is 9 percent.
#' 
#' #getAF(y~urban.rural+smoking.category+occupational.exposure,the.d,prev=0.09,ci=TRUE,conf_level=0.99)
#'
#'
#' ## genetic simulation using more risk factors (disease prevalence = 0.01)
#'
#'#thevec <- dbinom(0:40, size= 40, prob=0.2, log = FALSE)
#'#bin_fun <- function(beta_0){
#' # sum(thevec*exp(beta_0+.1*(0:40))/(1+exp(beta_0+.1*(0:40))))-0.01
#'#}
#'#beta_0 <- uniroot(bin_fun,lower=-8,upper=5)$root
#'#total_risk <- (0.01-exp(beta_0)/(1+exp(beta_0)))/0.01
#'#risk_per_snp <- total_risk/20
#'#case_probabilities <- (thevec*exp(beta_0 + (0:40)*0.1)/(1+exp(beta_0 + (0:40)*0.1)))/0.01
#'#control_probabilities <- thevec*1/(1+exp(beta_0 + (0:40)*0.1))/0.99
#'#simdata_genetic <- function(ncase,ncontrol){ 
#'# numbersnps_case <- sample(0:40,ncase,prob=case_probabilities,replace=TRUE)
#'#  numbersnps_control <- sample(0:40,ncase,prob=control_probabilities,replace=TRUE)
#'#  case_rows <- cbind(rep(1,ncase),matrix(0,nrow=ncase,ncol=20))
#'#  control_rows <- cbind(rep(0,ncase),matrix(0,nrow=ncontrol,ncol=20))
#'#  for(i in 1:ncase){
#'#    if(numbersnps_case[i]>0){   
#'#      positions <- sample(1:40,numbersnps_case[i])  
#'#      positions <- ceiling(positions/2)
#'#     for(j in 1:length(positions)) case_rows[i,positions[j]+1] <- case_rows[i,positions[j]+1]+1
#'#    }
#'#  }
#'#  for(i in 1:ncontrol){
#'#    if(numbersnps_control[i]>0){   
#'#      positions <- sample(1:40,numbersnps_control[i])
#'#      positions <- ceiling(positions/2)
#'#      for(j in 1:length(positions)){
#'#          control_rows[i,positions[j]+1]<- control_rows[i,positions[j]+1]+1 
#'#      }                                                       
#'#    }
#'# }
#'#  return(rbind(case_rows,control_rows))
#'#}
#'#the.d <- simdata_genetic(ncase=250, ncontrol=250)
#'#colnames(the.d) <- c("y",paste("SNP",1:20,sep=""))
#'
#'## Here we just calculate the approximate average fraction
#'## from 50 permutations and no confidence interval.
#'## If CI desired add the argument ci=TRUE and nsample_var to the function. 
#'## 50 permuations is chosen for speed. In reality, 1000 maybe needed
#'
#'#thesnp <- paste0("SNP", 1:20,sep="")
#'#(fmla <- as.formula(paste("y ~ ", paste(thesnp, collapse= "+"))))
#'#getAF(fmla, the.d,prev=0.01, allperm=FALSE,nsample_perm=50,ci=FALSE)
#'
#'## Instead of specifying the number of permutations, 
#'## you can specify an estimated approximation error. 
#'## The approximation error will be within this bound with 95% confidence 
#'## approximation error of 0.01 specified for reasons of speed.
#'## In reality, you may want to use a smaller value for approx_error.
#' 
#'#getAF(fmla, the.d,prev=0.01, allperm=FALSE,approx_error=0.01,ci=FALSE)

getAF <- function(f, the.data, ref_cat = NULL,refy=0, cat_confounders=NULL, cont_confounders=NULL, prev=NA,allperm=TRUE, nsample_perm=1000, approx_error = NA, ci=FALSE,conf_level=0.99, nsample_var=100, correction_factor=TRUE, quantile_int=FALSE, sep_est=TRUE){
  y <- as.character(f)[2]
  cat_riskfactors <- unlist(strsplit(as.character(f)[3],split=" + ",fixed=TRUE))
  # make sure ref_cat is in same order as variables in dataframe
  positions <- numeric(length(cat_riskfactors))
  for(i in 1:length(positions)) positions[i] <- grep(paste("^",cat_riskfactors[i],"$",sep=""),colnames(the.data),perl=TRUE)
  for(i in 1:length(positions)){
    if(colnames(the.data)[i]!=y){
      # transform categorical variables to numeric codes
      if(is.character(the.data[,positions[i]])){
        if(is.null(ref_cat)) the.data[,positions[i]] <- as.numeric(factor(the.data[,positions[i]]))
        if(!is.null(ref_cat)){
          the.data[,positions[i]] <- as.numeric(factor(the.data[,positions[i]],levels=c(ref_cat[i],setdiff(unique(the.data[,positions[i]]),ref_cat[i]))))
          ref_cat[i] <- 1
        }
      }    
      if(is.factor(the.data[,positions[i]])){
        if(is.null(ref_cat)) the.data[,positions[i]] <- as.numeric(the.data[,positions[i]])
        if(!is.null(ref_cat)){
          cha <- levels(the.data[,positions[i]])[the.data[,positions[i]]]
          the.data[,positions[i]] <- as.numeric(factor(cha,levels=c(ref_cat[i],setdiff(unique(cha),ref_cat[i]))))
          ref_cat[i] <- 1
        }
      } 
    }
  }

  if(!is.null(cat_confounders)){
    positions <- numeric(length(cat_confounders)) 
    for(i in 1:length(positions)) positions[i] <- grep(paste("^",cat_confounders[i],"$",sep=""),colnames(the.data),perl=TRUE)
    for(i in 1:length(positions)){
      if(is.character(the.data[,positions[i]])) the.data[,positions[i]] <- as.numeric(factor(the.data[,positions[i]]))
      if(as.factor(the.data[,positions[i]])) the.data[,positions[i]] <- as.numeric(the.data[,positions[i]])
    }
    cat_confounders <- as.matrix(the.data[,colnames(the.data)%in%cat_confounders])
    
  }
  if(!is.null(cont_confounders)) cont_confounders <- as.matrix(the.data[,colnames(the.data)%in%cont_confounders])
  the.colnames <- colnames(the.data)[colnames(the.data)%in%cat_riskfactors]
  cat_riskfactors <- as.matrix(the.data[,colnames(the.data)%in%cat_riskfactors])
  colnames(cat_riskfactors) <- the.colnames
  y <- as.character(f)[2]
  y <- the.data[,colnames(the.data)%in%as.vector(y)]
  ##  remove nas
  cat_confounders_n <- cat_confounders
  cont_confounders_n <- cont_confounders
  cat_riskfactors_n <- cat_riskfactors
  y_n <- y
  if(is.null(cat_confounders_n)) cat_confounders_n <- matrix(rep(1,length(y)),ncol=1)
  if(is.null(cont_confounders_n)) cont_confounders_n <- matrix(rep(1,length(y)),ncol=1)
  if(is.null(cat_riskfactors_n)) cat_riskfactors_n <- matrix(rep(1,length(y)),ncol=1)
   if(is.null(y_n)) y_n <- matrix(rep(1,length(y)),ncol=1)
  bigmat <- cbind(cat_confounders_n,cont_confounders_n,cat_riskfactors_n,y_n)
  thenas <-  apply(bigmat,1,function(x){sum(is.na(x))>0})
  cat_confounders <- cat_confounders[!thenas,,drop=FALSE]
  cont_confounders <- cont_confounders[!thenas,,drop=FALSE]
  cat_riskfactors <- cat_riskfactors[!thenas,,drop=FALSE]
   y <- y[!thenas]
  nvar_cat <- ncol(cat_riskfactors)
  the.colnames <- colnames(cat_riskfactors)
  if(is.null(the.colnames)) the.colnames <- paste("riskfactor",1:ncol(cat_riskfactors),sep="")
  if(is.list(y)) y <- unlist(y)
  if(!is.na(refy)) y <- recode(y,refy)  
  d <- matrix(0,nrow=length(y),ncol=0)  
  if(!is.null(cont_confounders)){
    d <- cbind(d,cont_confounders)
    colnames(d) = paste('cont_adjust',1:ncol(cont_confounders),sep="")
    cont_confounders<- as.matrix(cont_confounders)
  }
  
  if(!is.null(cat_confounders)){
    new.mat <- matrix(0,nrow=length(y),ncol=0)  
    for(i in 1:ncol(cat_confounders)){
      stuff <- recode(cat_confounders[,i],min(cat_confounders[,i],na.rm=TRUE))
      colnames(stuff) <- paste('cat_adjust_',i,"_",1:ncol(stuff),sep="")
      new.mat <- cbind(new.mat,stuff)
    }
    cat_confounders <- new.mat
    d <- cbind(d,cat_confounders)
  }
  
  if(!is.null(cat_riskfactors)){
    new.mat <- matrix(0,nrow=length(y),ncol=0)  
    for(i in 1:ncol(cat_riskfactors)){
      if(is.null(ref_cat)) stuff <- recode(cat_riskfactors[,i],min(cat_riskfactors[,i],na.rm=TRUE)) ## minimum value assumed to be the reference category
      ###  need to ensure that risk factors are added to the model in the right way for this to be true...
      if(!is.null(ref_cat)) stuff <- recode(cat_riskfactors[,i],ref_cat[i]) ## minimum value assumed to be the reference category
      colnames(stuff) <- paste('o',i,"_",1:ncol(stuff),sep="")
      new.mat <- cbind(new.mat,stuff)
    }
    cat_riskfactors <- new.mat
    d <- cbind(d,cat_riskfactors)
  }
  
  d <- cbind(y, d)
  
  df_vec <- 0
  
  d <<- as.data.frame(d)
  xnam <- colnames(d)[2:ncol(d)]
  w <- rep(1,length(y))
  thenas <- apply(matrix(d[,2:ncol(d)],nrow=length(y)),1,function(x){sum(is.na(x))>0})
  if(!is.na(prev)){ 
    
    
    w[y==0] <- (sum(y==1 & !thenas)/sum(y==0 & !thenas))*((1-prev)/prev)
    w[y==1] <-  1
    
  }
  
  
  (fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))
  f <- glm(fmla,family="binomial",data=as.data.frame(d), weights=w)
  
  quantile_est <- quantile_int
  ###  don't call function when K=1
  if(nvar_cat  == 1){
    startcol <- grep("(c|o)[0-9]_",colnames(d),perl=TRUE)[1]
    obs.cases <- sum(y==1)
    d_new=model.frame(f,d)
    if(colnames(d_new)[ncol(d_new)]=="(weights)") d_new <- d_new[,1:(ncol(d_new)-1)]
    d_new[,startcol:ncol(d_new)] <- 0
    the_coeffs <- f$coefficients
    orig_probs <- predict(f,type="response")
    nlp <- as.matrix(cbind(1,d_new[,2:ncol(d_new)]))%*%matrix(the_coeffs,ncol=1)
    new_probs <- exp(nlp)/(1+exp(nlp))
    w <- w[!thenas]
    af <- sum(w*(orig_probs-new_probs))
    newm <- glm(d[,1]~as.matrix(d[,2:ncol(d)]),family='binomial')
    cov_beta <-  summary(newm)$cov.unscaled
    if(!ci) return(af/obs.cases)
    if(ci){
      afvec <- numeric(nsample_var)
      beta_mat <- mvrnorm(n = nsample_var, mu=the_coeffs, Sigma=cov_beta, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      for(i in 1:nsample_var){
        beta_new <- beta_mat[i,]
        nlp <- as.matrix(cbind(1,d_new[,2:ncol(d_new)]))%*%matrix(beta_new,ncol=1)
        new_probs <- exp(nlp)/(1+exp(nlp))
        afvec[i] <- sum(w*(orig_probs-new_probs))
      }
      the.df <- nsample_var - 1
      the.sd <- sd(afvec)
      critval <- qt(1-(1-conf_level)/2,the.df)
      if(quantile_est) return(rbind(af/obs.cases,quantile(afvec/obs.cases,(1-conf_level)/2),quantile(afvec/obs.cases,1-(1-conf_level)/2)))
      out <- rbind(af/obs.cases,(af-critval*the.sd)/obs.cases,(af + critval*the.sd)/obs.cases)
      colnames(out) <- the.colnames
      rownames(out) <- c("estimate",paste(100*conf_level,"% lower bound",sep=""),paste(100*conf_level,"% upper bound",sep="")) 
      return(out)
    }
  }
  
  if(!ci){
    if(!allperm) out <- affun(f,d,nsample_perm =nsample_perm, prev=prev, allperm=allperm, approx_error=approx_error)
    if(allperm) out <- AF_exact(f,d,prev=prev)
    names(out) <- c(the.colnames,"total")
    return(out)
  }
  
  the.num <- nsample_var ###  make sure that the chunks are at least 100
  out <- affun_ci(f,d,nsample_perm =nsample_perm,nsample_var=the.num, prev=prev,conf_level=conf_level, sep_est=sep_est, correction_factor=correction_factor,allperm=allperm, quantile_est=quantile_est, approx_error=approx_error)
  colnames(out) <- c(the.colnames,"total")
  rownames(out) <- c("estimate",paste(100*conf_level,"% lower bound",sep=""),paste(100*conf_level,"% upper bound",sep=""))
  return(out)
}