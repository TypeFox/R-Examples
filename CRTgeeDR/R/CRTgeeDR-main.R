#' Doubly Robust Inverse Probability Weighted Augmented GEE Estimator 
#' 
#' This function implements a GEE estimator. It implements classical GEE, IPW-GEE, augmented GEE and IPW-Augmented GEE (Doubly robust). \cr
#' 
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#' @param id a vector which identifies the clusters. The length of "id" should be the same as the number of observations. Data are assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which CRTgeeDR is called.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See family for details of family functions.)
#' @param corstr a character string specifying the correlation structure. The following are permitted: '"independence"', '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'
#' @param Mv for "m-dependent", the value for m
#' @param weights A vector of weights for each observation. If an observation has weight 0, it is excluded from the calculations of any parameters. Observations with a NA anywhere (even in variables not included in the model) will be assigned a weight of 0.
#' @param aug A list of vector (one for A=1 treated, one for A=0 control) for each observation representing E(Y|X,A=a).
#' @param pi.a A number, Probability of treatment attribution P(A=1)
#' @param corr.mat The correlation matrix for "fixed". Matrix should be symmetric with dimensions >= the maximum cluster size. If the correlation structure is "userdefined", then this is a matrix describing which correlations are the same.
#' @param init.beta an optional vector with the initial values of beta. If not specified, then the intercept will be set to InvLink(mean(response)). init.beta must be specified if not using an intercept.
#' @param init.alpha an optional scalar or vector giving the initial values for the correlation. If provided along with Mv>1 or unstructured correlation, then the user must ensure that the vector is of the appropriate length.
#' @param init.phi an optional initial overdispersion parameter. If not supplied, initialized to 1.
#' @param scale.fix if set to TRUE, then the scale parameter is fixed at the value of init.phi.
#' @param sandwich if set to TRUE, the sandwich variance is provided together with the naive estimator of variance.
#' @param maxit maximum number of iterations.
#' @param tol tolerance in calculation of coefficients.
#' @param print.log if set to TRUE, a report is printed.
#' @param typeweights a character string specifying the weights implementation. The following are permitted: "GENMOD" for $W^{1/2}V^{-1}W^{1/2}$, "WV" for $V^{-1}W$
#' @param nameTRT Name of the variable containing information for the treatment 
#' @param model.weights an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted for the propensity score. Must model the probability of being observed. 
#' @param model.augmentation.trt an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted for the ouctome model for the treated group (A=1). 
#' @param model.augmentation.ctrl an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted for the ouctome model for the control group (A=0). 
#' @param stepwise.weights if set to TRUE, a stepwise for the propensity score is performed during the fit of the augmentation model for the OM
#' @param stepwise.augmentation if set to TRUE, a stepwise for the augmentation model is performed during the fit of the augmentation model for the OM
#' @param nameMISS Name of the variable containing information for the Missing indicator 
#' @param nameY Name of the variable containing information for the outcome 
#' @param sandwich.nuisance if set to TRUE, the nuisance adjusted sandwich variance is provided.
#' @param fay.adjustment if set to TRUE, the small-sample nuisance adjusted sandwich variance with Fay's adjustement is provided.
#' @param fay.bound if set to 0.75 by default, bound value used for Fay's adjustement.
#'        
#' @export

#' @author Melanie Prague [based on R packages 'geeM' L. S. McDaniel, N. C. Henderson, and P. J. Rathouz.  Fast Pure R Implementation of GEE: Application of the Matrix Package. The R Journal, 5(1):181-188, June 2013.]
#' 
#'
#' @details 
#' The estimator is founds by solving:
#'  \deqn{ 0= \sum_{i=1}^M \Bigg[ \boldsymbol D_i^T \boldsymbol V_i^{-1} \boldsymbol  W_i(\boldsymbol X_i, A_i, \boldsymbol \eta_W) \left( \boldsymbol Y_i - \boldsymbol B(\boldsymbol X_i, A_i, \boldsymbol \eta_B) \right) }
#'  \deqn{ \qquad + \sum_{a=0,1} p^a(1-p)^{1-a} \boldsymbol D_i^T \boldsymbol V_i^{-1}  \Big( \boldsymbol B(\boldsymbol X_i,A_i=a, \boldsymbol \eta_B) -\boldsymbol \mu_i(\boldsymbol \beta,A_i=a)\Big) \Bigg]}
#'  where \eqn{\boldsymbol D_i=\frac{\partial \boldsymbol \mu_i(\boldsymbol \beta,A_i)}{\partial \boldsymbol \beta^T}} is the design matrix, \eqn{\boldsymbol V_i} is the covariance matrix equal to \eqn{\boldsymbol U_i^{1/2} \boldsymbol C(\boldsymbol \alpha)\boldsymbol U_i^{1/2}$ with $\boldsymbol U_i} a diagonal matrix with elements \eqn{{\rm var}(y_{ij})} and \eqn{\boldsymbol C(\boldsymbol \alpha)} is the working correlation structure with non-diagonal terms \eqn{\boldsymbol \alpha}.
#'  Parameters \eqn{\boldsymbol \alpha} are estimated using simple moment estimators from the Pearson residuals. 
#'  The matrix of weights \eqn{\boldsymbol W_i(\boldsymbol X_i, A_i, \boldsymbol \eta_W)=diag\left[R_{ij}/\pi_{ij}(\boldsymbol X_i, A_i, \boldsymbol \eta_W)\right]_{j=1,\dots,n_{i}}}, where \eqn{\pi_{ij}(\boldsymbol X_i, A_i, \boldsymbol \eta_W)=P(R_{ij}|\boldsymbol X_i, A_i)} is the Propensity score (PS).
#'  The function \eqn{\boldsymbol B(\boldsymbol X_i,A_i=a,\boldsymbol \eta_B)}, which is called the Outcome Model (OM), is a function linking \eqn{Y_{ij}} with \eqn{\boldsymbol X_i} and \eqn{A_i}.
#'  The \eqn{\boldsymbol \eta_B} are nuisance parameters that are estimated.
#'  The estimator is most efficient if the OM is equal to \eqn{E(\boldsymbol Y_i|\boldsymbol X_i,A_i=a)}
#'  The estimator denoted \eqn{\hat{\beta}_{aug}} is found by solving the estimating equation.
#'  Although  analytic solutions sometimes exist, coefficient estimates are generally obtained using an iterative procedure such as the Newton-Raphson method. 
#'  Automatic implementation is such that, \eqn{\hat{ \boldsymbol \eta}_W$ in $\boldsymbol W_i(\boldsymbol X_i, A_i, \hat{ \boldsymbol \eta}_W)} are obtained using a logistic regression and \eqn{\hat{ \boldsymbol \eta}_B$ in $\boldsymbol B(\boldsymbol X_i,A_i,\hat{ \boldsymbol \eta}_B)} are obtained using a linear regression. 
#'  
#'  
#'  The variance of \eqn{\hat{\boldsymbol \beta}_{aug}} is estimated by the sandwich variance estimator. 
#'  There are two external sources of variability that need to be accounted for: estimation of \eqn{\boldsymbol \eta_W} for the PS and of \eqn{\boldsymbol \eta_B} for the OM. 
#'  We denote \eqn{\boldsymbol \Omega=(\boldsymbol \beta, \boldsymbol \eta_W,\boldsymbol \eta_B)} the estimated parameters of interest and nuisance parameters. 
#'  We can stack estimating functions and score functions for \eqn{\boldsymbol \Omega}:
#'  \deqn{\small \boldsymbol U_i(\boldsymbol \Omega)= \left( \begin{array}{c} \boldsymbol \Phi_i(\boldsymbol Y_i,\boldsymbol X_i,A_i,\boldsymbol \beta, \boldsymbol \eta_W, \boldsymbol \eta_B) \\ \boldsymbol S^W_i(\boldsymbol X_i, A_i, \boldsymbol \eta_W)\\ \boldsymbol S^B_i(\boldsymbol X_i, A_i, \boldsymbol \eta_B)\\ \end{array} \right)}
#'  where \eqn{\boldsymbol S^W_i} and \eqn{\boldsymbol S^B_i} represent the score equations for patients in cluster \eqn{i} for the estimation of \eqn{\boldsymbol \eta_W} and \eqn{\boldsymbol \eta_B}  in the PS and the OM. 
#'  A standard Taylor expansion paired with Slutzky's theorem and the central limit theorem give the sandwich estimator adjusted for nuisance parameters estimation in the OM and PS:
#'  \deqn{Var(\boldsymbol \Omega)={{E\left[\frac{\partial   \boldsymbol U_i(\boldsymbol \Omega)}{\partial \boldsymbol \Omega}\right]}^{-1}}^{T} \underbrace{{E\left[ \boldsymbol U_i(\boldsymbol \Omega)\boldsymbol U_i^T(\boldsymbol \Omega) \right]}}_{\boldsymbol \Delta_{adj}} \underbrace{E\left[\frac{\partial   \boldsymbol U_i(\boldsymbol \Omega)}{\partial \boldsymbol \Omega}\right]^{-1} }_{\boldsymbol \Gamma^{-1}_{adj}}.}
#'  
#'  
#'  
#'  
#' @references Details regarding implementation can be found in 
#' \itemize{
#'       \item 'Augmented GEE for improving efficiency and validity of estimation in cluster randomized trials by leveraging cluster-and individual-level covariates' - 2012 - Stephens A., Tchetgen Tchetgen E. and De Gruttola V. : Stat Med 31(10) - 915-930.
#'       \item 'Accounting for interactions and complex inter-subject dependency for estimating treatment effect in cluster randomized trials with missing at random outcomes' - 2015 - Prague M., Wang R., Stephens A., Tchetgen Tchetgen E. and De Gruttola V. : in revision.
#'       \item 'Fast Pure R Implementation of GEE: Application of the Matrix Package' - 2013 - McDaniel, Lee S and Henderson, Nicholas C and Rathouz, Paul J : The R Journal 5(1) - 181-197.
#'       \item 'Small-Sample Adjustments for Wald-Type Tests Using Sandwich Estimators' - 2001 - Fay, Michael P and Graubard, Barry I : Biometrics 57(4) - 1198-1206.
#' }
#' 
#' 
#' @return An object of type 'CRTgeeDR' \cr 
#' @return $beta Final values for regressors estimates \cr 
#' \itemize{
#' \item $phi scale parameter estimate\cr 
#' \item $alpha Final values for association parameters in the working correlation structure when exchangeable\cr 
#' \item $coefnames Name of the regressors in the main regression \cr
#' \item $niter Number of iteration done by the algorithm before convergence
#' \item $converged convergence status
#' \item $var.naiv Variance of the estimates model based (naive)\cr 
#' \item $var Variance of the estimates sandwich\cr 
#' \item $var.nuisance Variance of the estimates nuisance adjusted sandwich\cr 
#' \item $var.fay Variance of the estimates nuisance adjusted sandwich with Fay correction for small samples
#' \item $call Call function
#' \item $corr Correlation structure used
#' \item $clusz Number of unit in each cluster
#' \item $FunList List of function associated with the family 
#' \item $X design matrix for the main regression
#' \item $offset Offset specified in the regression
#' \item $eta predicted values
#' \item $weights Weights vector used in the diagonal term for the IPW
#' \item $ps.model Summary of the regression fitted for the PS if computed internally
#' \item $om.model.trt Summary of the regression fitted for the OM for treated if computed internally
#' \item $om.model.ctrl Summary of the regression fitted for the OM for control if computed internally
#' }
#' 
#' @import MASS
#' @import Matrix
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom stats as.formula binomial fitted gaussian glm median model.frame model.matrix model.offset model.response na.pass pnorm predict quantile step terms weights
#' @importFrom methods as
#'  @examples
#'  
#'  data(data.sim)
#'  \dontrun{
#'  #### STANDARD GEE
#'  geeresults<-geeDREstimation(formula=OUTCOME~TRT,
#'                                id="CLUSTER" , data = data.sim,
#'                                family = "binomial", corstr = "independence")
#'  summary(geeresults)                            
#'  #### IPW GEE
#'  ipwresults<-geeDREstimation(formula=OUTCOME~TRT,
#'                                id="CLUSTER" , data = data.sim,
#'                                family = "binomial", corstr = "independence",
#'                                model.weights=I(MISSING==0)~TRT*AGE)
#'  summary(ipwresults)     
#'  #### AUGMENTED GEE
#'  augresults<-geeDREstimation(formula=OUTCOME~TRT,
#'                                id="CLUSTER" , data = data.sim,
#'                                family = "binomial", corstr = "independence",
#'                                model.augmentation.trt=OUTCOME~AGE,
#'                                model.augmentation.ctrl=OUTCOME~AGE, stepwise.augmentation=FALSE)
#'  summary(augresults)   
#'  }
#'  #### DOUBLY ROBUST
#'  drresults<-geeDREstimation(formula=OUTCOME~TRT,
#'                                id="CLUSTER" , data = data.sim,
#'                                family = "binomial", corstr = "independence",
#'                                model.weights=I(MISSING==0)~TRT*AGE,
#'                                model.augmentation.trt=OUTCOME~AGE,
#'                                model.augmentation.ctrl=OUTCOME~AGE, stepwise.augmentation=FALSE)
#'  summary(drresults)                            

geeDREstimation <- function(formula, id,data = parent.frame(), family = gaussian, corstr = "independence", Mv = 1, weights = NULL, aug=NULL,pi.a=1/2, corr.mat = NULL, init.beta=NULL, init.alpha=NULL, init.phi = 1, scale.fix=FALSE, sandwich=TRUE, maxit=20, tol=0.00001,print.log=FALSE,typeweights="VW",nameTRT="TRT",model.weights=NULL,model.augmentation.trt=NULL,model.augmentation.ctrl=NULL,stepwise.augmentation=FALSE,stepwise.weights=FALSE,nameMISS="MISSING",nameY="OUTCOME",sandwich.nuisance=FALSE,fay.adjustment=FALSE,fay.bound=0.75){
  if(print.log)print("********************************************************************************************")
  if(print.log){print("DESCRIPTION: Doubly Robust Inverse Probability Weighted Augmented GEE estimator")}
  if(print.log)print("********************************************************************************************")
  
  call <- match.call()
  
  ### Get the information from the family: link function for the outcome
  FunList <- getfam(family)
  LinkFun <- FunList$LinkFun
  VarFun <- FunList$VarFun
  InvLink <- FunList$InvLink
  InvLinkDeriv <- FunList$InvLinkDeriv
  
  ### Check that all the attributes are ok
  if(scale.fix & is.null(init.phi)){
    stop("If scale.fix=TRUE, then init.phi must be supplied")
  }
  
  if((!(sum(!(unique(data[,nameTRT])%in%c(0,1)))==0))&(!(is.null(aug)&is.null(model.augmentation.trt)&is.null(model.augmentation.ctrl)))){
    stop("Augmentation is requested whereas more than two level of treatment exist. Implementation not available yet.")
  }
  
  if(is.null(weights)&is.null(model.weights)){
    typeweights <- NULL
  }else{
    w.vec <- c("GENMOD", "VWR")
    w.match <- charmatch(typeweights, w.vec)
    if(is.na(w.match)){stop("Unsupported type of weights specification")}  
  }
  
  if(!is.null(model.weights))dat.ps <- data
  if(!is.null(model.augmentation.trt)){
    dat.om.trt <- data
    dat.om.trt<-cleandata(dat=dat.om.trt,type="OM model",nameY=nameY,cc=FALSE,formula=model.augmentation.trt)
    #dat.om.trt<-dat.om.trt[which(!is.na(dat.om.trt[,nameY])),]
  }
  if(!is.null(model.augmentation.ctrl)){
    dat.om.ctrl <- data
    dat.om.ctrl<-cleandata(dat=dat.om.ctrl,type="OM model",nameY=nameY,cc=FALSE,formula=model.augmentation.ctrl)
    #dat.om.ctrl<-dat.om.ctrl[which(!is.na(dat.om.ctrl[,nameY])),]
  }

  ## Compute the PS and the weights if needed
  propensity.score<-NULL
  if(!is.null(model.weights)){
    dat.ps<-cleandata(dat=dat.ps,type="PS model",nameY=nameMISS,cc=FALSE,formula=model.weights)
    if(print.log){print("------------------------------------------------------------> Information for PS")}
    if(print.log){print("PS is computed internally...")}
    warning("Warning: Propensity score is computed internally, make sure that the formula given in model.weights models the probability of being observed. Information given in attribute 'weights' will be ingnored.")
    if(stepwise.weights){
      propensity.score<-step(glm(model.weights,data=dat.ps,family=binomial(link = "logit")),trace=0)
    }else{
      propensity.score<-glm(model.weights,data=dat.ps,family=binomial(link = "logit"))
    }
    dat.ps$weights.inside<-1/(fitted(propensity.score))
    weights <- "weights.inside"
    if(print.log){
      print("Details for PS regression")
      print(summary(propensity.score))}
  }
  weightsname<-weights
  
  ### Get the data - either from the environement or from the attribute dataset
  dat <- model.frame(formula, data, na.action=na.pass) 
  if(!is.null(model.weights)){
    data$weights.inside<-dat.ps$weights.inside
  }
  nn <- dim(dat)[1]
  
  if(typeof(data) == "environment"){
    id <- id
    if(!is.null(weights)) weights <- weights
    dat$id <- id
  }else{
    if(length(id) == 1){
      subj.col <- which(colnames(data) == id)  
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      }else{
        id <- parent.frame()$id
      }
    }else if(is.null(id)){
      id <- 1:nn
    }
    
    if(!is.null(weights)){
      if(length(weights) == 1){
        weights.col <- which(colnames(data) == weights)  
        if(length(weights.col) > 0){
          weights <- data[,weights.col]
        }else{          
          weights <- parent.frame()$weights
        }
        dat$weights<-weights
      }else if(is.null(weights)){
        weights <- NULL
      }
    }
  }
  dat$id <- id 
  if(!is.null(weights)){
    dat$weights <- weights
  }else{
    dat$weights <- 1
  }

  ## Clean the dataset for missing data in the covariates
  dat<-cleandata(dat=dat,type="marginal model",nameY=nameY,cc=TRUE,formula=formula)
  weights<-dat$weights
  
  includedvec <- weights>0
  inclsplit <- split(includedvec, id)

  #Drop cluster with no observations
  dropid <- NULL
  allobs <- T
  if(any(!includedvec)){
    allobs <- F
    for(i in 1:length(unique(id))){
      if(all(!inclsplit[[i]])){
        dropid <- c(dropid, i)
      }    
    }
  }

  
  if(length(dropid)>0){
    dropind <- which(is.element(id, dropid))
    dat <- dat[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]    
    id <- id[-dropind]
  }
  ### Get the numbers of clusters and individuals
  nn <- dim(dat)[1]
  K <- length(unique(id))
  
  modterms <- terms(formula)
	X <- model.matrix(formula,dat)
  X.t <- X.c <- X
  if(nameTRT %in% colnames(dat)){
  X.t[,nameTRT]<-1.0000
  X.c[,nameTRT]<-0.0000
  }else{
    stop("User need to provide the name of the treatment variable in nameTRT. Default nameTRT='TRT' does not exist in the dataset.")   
  }
  Y <- as.matrix(model.response(dat))
  
  
  
  ## Compute the OM and the augmentation terms if needed
  B<-NULL
  om.t<-NULL
  om.c<-NULL
  
  if(!(is.null(model.augmentation.trt)|is.null(model.augmentation.ctrl))){
    if(print.log){print("------------------------------------------------------------> Information for OM")}
    if(print.log){print("OM are computed internally...")}
    warning("Warning: Outcome model for augmentation is computed internally. Information given in attribute 'aug' will be ingnored.")
  
  data.trt<-dat.om.trt[which((dat.om.trt[,nameTRT]==1)),]
  data.ctrl<-dat.om.ctrl[which((dat.om.ctrl[,nameTRT]==0)),]
  data.t<-dat.om.trt
  data.t[,nameTRT]<-1
  data.c<-dat.om.ctrl
  data.c[,nameTRT]<-0 
 
  if(stepwise.augmentation){
    om.t<-step(glm(model.augmentation.trt,data=data.trt,family=family),trace=0)
    om.c<-step(glm(model.augmentation.ctrl,data=data.ctrl,family=family),trace=0)
  }else{
    om.t<-glm(model.augmentation.trt,data=data.trt,family=family)
    om.c<-glm(model.augmentation.ctrl,data=data.ctrl,family=family)     
  }
  data$B1<-InvLink(predict(om.t,newdata=data.t))
  data$B0<-InvLink(predict(om.c,newdata=data.c))    
  aug<-c(ctrl="B0",trt="B1")
  if(print.log){
    print("Details for OM regression in treated")
    print(summary(om.t))
    print("Details for OM regression in control")
    print(summary(om.c))
  } 
  }

  if(!is.null(aug)){
    if(length(aug)!=2){
      stop("If augmentation is requested, then aug must be supplied (length=2)")
    }  
   
    ###  Alert the user if there are covariates in the main regression and augmentation is used -- Stop had been removed because theoritical result is still valid.
    if(formula!=as.formula(paste(paste(as.character(formula)[2],as.character(formula)[1],sep=""),nameTRT,sep=""))){
      warning("Warning: Augmentation approach is used with a marginal model including covariates.")
    }
    B.c<-data[,aug["ctrl"]]
    B.t<-data[,aug["trt"]]
    temp<-cbind(X,B.c,B.t)
    Bi<-ifelse(temp[,nameTRT]==1,temp[,"B.t"],temp[,"B.c"])
    B<-cbind(B.c,B.t,Bi)
  }
  
	### if no offset is given, then set to zero
  offset <- model.offset(dat)			
  p <- dim(X)[2]
	if(is.null(offset)){
 		off <- rep(0, nn)
	}else{
		off <- offset
	}
  
	# Is there an intercept column?
	interceptcol <- apply(X==1, 2, all)
	
	## Basic check to see if link and variance functions make any kind of sense
	linkOfMean <- LinkFun(mean(Y,na.rm=T))
	if( any(is.infinite(linkOfMean) | is.nan(linkOfMean)) ){
		stop("Infinite or NaN in the link of the mean of responses.  Make sure link function makes sense for these data.")
	}
	if( any(is.infinite( VarFun(mean(Y))) | is.nan( VarFun(mean(Y)))) ){
		stop("Infinite or NaN in the variance of the mean of responses.  Make sure variance function makes sense for these data.")
	}
	
	if(is.null(init.beta)){
		if(any(interceptcol)){
			#if there is an intercept and no initial beta, then use link of mean of response
			init.beta <- rep(0, dim(X)[2])
			init.beta[which(interceptcol)] <- linkOfMean
		}else{
			stop("Must supply an initial beta if not using an intercept.")
		}
	}
  

  # Number of included observations for each cluster
  includedlen <- rep(0, K)
  len <- rep(0,K)
  uniqueid <- unique(id)

  tmpwgt <- as.numeric(includedvec)
  idspl <-ifelse(tmpwgt==0, NA, id)
  includedlen <- as.numeric(summary(split(Y, idspl, drop=T))[,1])
  len <- as.numeric(summary(split(Y, id, drop=T))[,1])

  W <- Diagonal(x=weights)
  sqrtW <- sqrt(W)
  included <- Diagonal(x=(as.numeric(weights>0)))
  
	# Figure out the correlation structure
	cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
	cor.match <- charmatch(corstr, cor.vec)
	
  if(is.na(cor.match)){stop("Unsupported correlation structure")}  
    
	# Set the initial alpha value
	if(is.null(init.alpha)){
		alpha.new <- 0.2
		if(cor.match==4){
			# If corstr = "m-dep"
			alpha.new <- 0.2^(1:Mv)
		}else if(cor.match==5){
			# If corstr = "unstructured"
			alpha.new <- rep(0.2, sum(1:(max(len)-1)))
		}else if(cor.match==7){
			# If corstr = "userdefined"
			alpha.new <- rep(0.2, max(unique(as.vector(corr.mat))))
		}
	}else{
		alpha.new <- init.alpha	
	}
	#if no initial overdispersion parameter, start at 1
	if(is.null(init.phi)){
		phi <- 1
	}else{
    phi <- init.phi
	}
	
	beta <- init.beta
	

	
	#Set up matrix storage
	StdErr <- Diagonal(nn)
	dInvLinkdEta <- Diagonal(nn)
	Resid <- Diagonal(nn)


	# Initialize for each correlation structure
  if(print.log){print(paste("Initialize the correlation structure:",corstr,sep=" "))}
	if(cor.match == 1){
		# INDEPENDENCE
		R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
		BlockDiag <- getBlockDiag(len)$BDiag
	}else if(cor.match == 2){
		# AR-1
		tmp <- buildAlphaInvAR(len)
		# These are the vectors needed to update the inverse correlation
		a1<- tmp$a1
		a2 <- tmp$a2
		a3 <- tmp$a3
		a4 <- tmp$a4
		# row.vec and col.vec for the big block diagonal of correlation inverses
		# both are vectors of indices that facilitate in updating R.alpha.inv
		row.vec <- tmp$row.vec
		col.vec <- tmp$col.vec
		BlockDiag <- getBlockDiag(len)$BDiag

	}else if(cor.match == 3){
		# EXCHANGEABLE
		# Build a block diagonal correlation matrix for updating and sandwich calculation
		# this matrix is block diagonal with all ones.  Each block is of dimension cluster size.
		tmp <- getBlockDiag(len)
		BlockDiag <- tmp$BDiag
		
		# Create a vector of length number of observations with associated cluster size for each observation
		n.vec <- vector("numeric", nn)
		index <- c(cumsum(len) - len, nn)
		for(i in 1:K){
			n.vec[(index[i]+1) : index[i+1]] <-  rep(len[i], len[i])
		}
	}else if(cor.match == 4){
		# M-DEPENDENT, check that M is not too large
		if(Mv >= max(len)){
			stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
		}
		
		# Build block diagonal similar to in exchangeable case, also get row indices and column
		# indices for fast matrix updating later.		
		tmp <- getBlockDiag(len)
		BlockDiag <- tmp$BDiag
		row.vec <- tmp$row.vec
		col.vec <- tmp$col.vec
		R.alpha.inv <- NULL
	}else if(cor.match == 5){
		# UNSTRUCTURED
		if( max(len^2 - len)/2 > length(len)){
			stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
		}
		tmp <- getBlockDiag(len)
		BlockDiag <- tmp$BDiag
		row.vec <- tmp$row.vec
		col.vec <- tmp$col.vec	
	}else if(cor.match == 6){
		# FIXED
		# check if matrix meets some basic conditions
		corr.mat <- checkFixedMat(corr.mat, len)

		R.alpha.inv <- as(getAlphaInvFixed(corr.mat, len), "symmetricMatrix")/phi
		BlockDiag <- getBlockDiag(len)$BDiag
	}else if(cor.match == 7){
		# USERDEFINED
		corr.mat <- checkUserMat(corr.mat, len)

		# get the structure of the correlation matrix in a way that
		# I can use later on.
		tmp1 <- getUserStructure(corr.mat)
		corr.list <- tmp1$corr.list
		user.row <- tmp1$row.vec
		user.col <- tmp1$col.vec
		struct.vec <- tmp1$struct.vec
		
		# the same block diagonal trick.
		tmp2 <- getBlockDiag(len)
		BlockDiag <- tmp2$BDiag
		row.vec <- tmp2$row.vec
		col.vec <- tmp2$col.vec

	}else if(cor.match == 0){
		stop("Ambiguous Correlation Structure Specification")
	}else{
		stop("Unsupported Correlation Structure")	
	}

  stop <- F
	converged <- F	
	count <- 0
	beta.old <- beta
	unstable <- F
	phi.old <- phi
	
  
if(print.log){print("------------------------------------------------------------> Dataset description")}
if(print.log){print(paste("Number of CLUSTERS:",length(unique(id)),sep=" "))}
if(print.log){print(paste("Number of INDIVIDUAL:",nn,sep=" "))}	
if(print.log){print(paste("Number of Observations included:",sum(includedvec),sep=" "))}
if(print.log){print(paste("Variable for PS:",weightsname,sep=" "))}
if(print.log){print(paste("Variable for OUTCOME:",nameY,sep=" "))}
if(print.log){print(paste("Variable for MISSING:",nameMISS,sep=" "))}
if(print.log){print(paste("Variable for TRT:",nameTRT,sep=" "))}

# Main fisher scoring loop
if(print.log){print("------------------------------------------------------------> Estimations")}
if(max(diag(sqrtW))==1){
  if(print.log)print("NON-IPW Analysis")
}else{
  if(print.log)print("IPW Analysis")
}
if(is.null(B)){
  if(print.log)print("NON-AUGMENTED Analysis")
}else{
  if(print.log)print("AUGMENTED Analysis")
}

	while(!stop){		
		count <- count+1
		if(print.log){print(paste("Main loop for estimation:",count,sep=" "))}
    
		eta <- as.vector(X %*% beta) + off
				
		mu <- InvLink(eta)
		
		diag(StdErr) <- sqrt(1/VarFun(mu))
		
    	if(!scale.fix){
			  phi <- updatePhi(Y, mu, VarFun, p, StdErr, included, includedlen)
    	}
    	phi.new <- phi
	   	if(print.log){print(paste("Phi",phi.new,sep=" "))}
		
    ## Calculate alpha, R(alpha)^(-1) / phi
		if(cor.match == 2){
			# AR-1
			alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs)
			R.alpha.inv <- getAlphaInvAR(alpha.new, a1,a2,a3,a4, row.vec, col.vec)/phi
		}else if(cor.match == 3){
			#EXCHANGEABLE		  
			alpha.new <- updateAlphaEX(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, included, includedlen)
			R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
		}else if(cor.match == 4){
			# M-DEPENDENT
			if(Mv==1){
				alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, len, StdErr, p, included, includedlen, includedvec, allobs)
			}else{
				alpha.new <- updateAlphaMDEP(Y, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, Mv, included, includedlen, allobs)
				if(sum(len>Mv) <= p){
					unstable <- T
				}
			}
			if(any(alpha.new >= 1)){
				stop <- T
				warning("Some estimated correlation is greater than 1, stopping.")
			}
			R.alpha.inv <- getAlphaInvMDEP(alpha.new, len, row.vec, col.vec)/phi		
		}else if(cor.match == 5){
			# UNSTRUCTURED
			alpha.new <- updateAlphaUnstruc(Y, mu, VarFun, phi, id, len, StdErr, Resid,  p, BlockDiag, included, includedlen, allobs)
			# This has happened to me (greater than 1 correlation estimate)
			if(any(alpha.new >= 1)){
				stop <- T
				warning("Some estimated correlation is greater than 1, stopping.")
			}
			R.alpha.inv <- getAlphaInvUnstruc(alpha.new, len, row.vec, col.vec)/phi
		}else if(cor.match ==6){
			# FIXED CORRELATION, DON'T NEED TO RECOMPUTE
			R.alpha.inv <- R.alpha.inv*phi.old/phi
		}else if(cor.match == 7){
			# USER SPECIFIED
			alpha.new <- updateAlphaUser(Y, mu, phi, id, len, StdErr, Resid, p, BlockDiag, user.row, user.col, corr.list, included, includedlen, allobs)
			R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, user.row, user.col, row.vec, col.vec)/phi
		}else if(cor.match == 1){
			# INDEPENDENT
			R.alpha.inv <-  Diagonal(x = rep.int(1/phi, nn))
			alpha.new <- "independent"
		}
			
		beta.list <- updateBeta(Y=Y, X=X,X.t=X.t,X.c=X.c, B=B, beta=beta, off=off, 
                             InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun,
                             R.alpha.inv=R.alpha.inv, StdErr=StdErr, dInvLinkdEta=dInvLinkdEta,
                             tol=tol, sqrtW=sqrtW,W=W,included=included,typeweights=typeweights,pi.a=pi.a)	
		if(print.log){print(paste(cat("beta",unlist(beta.list$beta))," "))}
		if(print.log){print(paste("alpha",alpha.new,sep=" "))}   
    #print("LOOP")
    beta <- beta.list$beta
		phi.old <- phi
		if( max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol ){converged <- T; stop <- T}
		if(count >= maxit){stop <- T}		
		beta.old <- beta		
	}
	biggest <- which.max(len)[1]
	index <- cumsum(len[biggest])
	biggest.R.alpha.inv <- R.alpha.inv[(index+1):(index+len[biggest]) , (index+1):(index+len[biggest])]
	eta <- as.vector(X %*% beta) + off
  if(print.log){print("------------------------------------------------------------> Variance estimation")}

	sandvar.list <- list()
	sandvar.list$sandvar <- NULL
	if(sandwich){
		sandvar.list <- getSandwich(Y=Y, X=X,X.t=X.t,X.c=X.c, B=B, beta=beta, off=off, 
                                id=id, R.alpha.inv=R.alpha.inv, phi=phi,
                                InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun,
                                hessMat=beta.list$hess, StdErr=StdErr, dInvLinkdEta=dInvLinkdEta, 
                                BlockDiag=BlockDiag, sqrtW=sqrtW,W=W,included=included,typeweights=typeweights,pi.a=pi.a,
                                print.log=print.log)
	}else{
		sandvar.list <- list()
		sandvar.list$sandvar <- NULL
		if(print.log){print("No sandwich")}   
	}

  if(is.null(model.augmentation.trt)&is.null(model.augmentation.ctrl)&is.null(model.weights)){
    sandwich.nuisance<-FALSE
  }
	
#	sandvarnuis.list <- list()
#	sandvarnuis.list$sandadjvar <- NULL
  if(sandwich.nuisance){
     dat.nuis<-cleandata(dat=data,type="nuisance",nameY=nameY,cc=FALSE,formula=formula,print=FALSE)
     tryCatch({
     sandvarnuis.list <- getSandwichNuisance(Y=Y, X=X,X.t=X.t,X.c=X.c, B=B, beta=beta, off=off, 
                                            id=id, R.alpha.inv=R.alpha.inv, phi=phi,
                                            InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun,
                                            hessMat=beta.list$hess, StdErr=StdErr, dInvLinkdEta=dInvLinkdEta, 
                                            BlockDiag=BlockDiag, sqrtW=sqrtW,W=W,included=included,typeweights=typeweights,pi.a=pi.a,
                                            nameTRT=nameTRT,propensity.score=propensity.score,om.t=om.t,om.c=om.c,
                                            data=dat.nuis,nameY=nameY,nameMISS=nameMISS,print.log=print.log)
     }, error=function(e){
       cat("There was an error in the nuisance variance computation \n")
      })
  }else{
    sandvarnuis.list <- list()
    sandvarnuis.list$sandadjvar <- NULL
    if(print.log){print("No nuisance-adjusted sandwich")}   
  }
	
	sandvarfay.list <- list()
	sandvarfay.list$sandadjfay <- NULL
	if(fay.adjustment){
	 tryCatch({
	    sandvarfay.list <- getFay(formula=formula,id=id,family=family,data=data,corstr=corstr,b=fay.bound,beta=beta,alpha=alpha.new,scale=phi.new,Y=Y, X=X,hessMAT=beta.list$hess,
	                              X.t=X.t,X.c=X.c, B=B, off=off, 
	                              R.alpha.inv=R.alpha.inv, phi=phi,
	                              InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun,
	                              StdErr=StdErr, dInvLinkdEta=dInvLinkdEta, 
	                              BlockDiag=BlockDiag, sqrtW=sqrtW,W=W,included=included,typeweights=typeweights,pi.a=pi.a,
	                              nameTRT=nameTRT,propensity.score=propensity.score,om.t=om.t,om.c=om.c,
	                              nameY=nameY,nameMISS=nameMISS,print.log=print.log)
	   },
	  error=function(e){
	    cat("There was an error in the fay variance computation \n")} 
	 )
	}else{
	  sandvarfay.list <- list()
	  sandvarfay.list$sandadjfay <- NULL
	  if(print.log){print("No fay adjustment sandwich for small sample")}   
	}

	if(!converged){warning("Did not converge")}
	if(unstable){warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")}

	# Create object of class CRTgeeDR with information about the fit
  dat <- model.frame(formula, data, na.action=na.pass)
  X <- model.matrix(formula, dat)
  
  if(alpha.new == "independent"){alpha.new <- 0}
	results <- list()
	results$beta <- as.vector(beta)
	results$phi <- phi
	results$alpha <- alpha.new
	if(cor.match == 6){
		results$alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat,1)!=0)])
	}
	results$coefnames <- colnames(X)
	results$niter <- count
	results$converged <- converged
	results$var.naiv <- solve(beta.list$hess)##*phi.new  ## call model-based
	results$var <- sandvar.list$sandvar
  results$var.nuisance <- sandvarnuis.list$sandadjvar
  results$var.fay <- sandvarfay.list$sandadjfay
	results$call <- call
	results$corr <- cor.vec[cor.match]
	results$clusz <- len
	results$FunList <- FunList
	results$X <- X
	results$offset <- off
	results$eta <- eta
  results$dropped <- dropid
  results$weights <- weights
  results$ps.model<-propensity.score
  if(is.null(propensity.score)){
    results$used.weights<- NULL
  }else{
    results$used.weights<-diag(W)
  }
  results$om.model.trt<-om.t
  results$om.model.ctrl<-om.c
 	class(results) <- "CRTgeeDR"
  if(print.log){
    print("********************************************************************************************")
    print("RESULTS: Doubly Robust Inverse Probability Weighted Augmented GEE estimator")
    print("********************************************************************************************")
    print(results)
    print(summary(results))
  }
	return(results)
}


#' The data.sim Dataset.
#' 
#' HIV risk of infection after STI/HIV intervention in a cluster randomized trial.
#'
#' @details A dataset containing the HIV risk scores and presence of risky behaviors (yes/no) and other covarites of 10000 subjects among 100 communities.
#' The variables are as follows:
#'
#' \itemize{
#'   \item IDPAT subject id
#'   \item CLUSTER cluster id
#'   \item TRT treatment status, 1 is received STI/HIV intervention
#'   \item X1 A covariate following a N(0,1)
#'   \item JOB employement status
#'   \item MARRIED marital status
#'   \item AGE age
#'   \item HIV.KNOW Score for HIV knowlege
#'   \item RELIGION religiosity score
#'   \item OUTCOME Binary outcome - 1 if the subject is at high risk of HIV infection, 0 otherwise. NA if missing.
#'   \item MISSING 1 if the ouctome is missing - 0 otherwise.
#' }
#'
#' @format A data frame with 10000 rows and 8 variables
#' @name data.sim
NULL




#' Doubly Robust Inverse Probability Weighted Augmented GEE estimator 
#'
#' The CRTgeeDR package allows you to estimates parameters in a regression model (with possibly a link function). 
#' It allows treatment augmentation and IPW for missing data alone.  
#'
#' The only function you're likely to need from \pkg{CRTgeeDR} is
#' \code{\link{geeDREstimation}}. Otherwise refer to the help documentation.
#'
#' @docType package
#' @name CRTgeeDR
NULL





