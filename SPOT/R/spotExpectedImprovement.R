##############################################################
#' Lower Confidence Bound Infill Criterion
#'
#' Calculate lower confidence bound for all objectives and all observations.
#' Also used by\code{\link{spotInfillLcbHyperVolume}}.
#'
#' @param resy predicted objective values
#' @param resvar predicted variance
#' @param y the current Pareto front
#' @param ref unused
#'
#' @return returns the LCB for each row and column in resy.
#'
#' @export
##############################################################
spotInfillLcbMulti <- function(resy,resvar,y,ref=NULL){
	p=NULL
	for(i in 1:ncol(resy)){
		p=cbind(p,resy[,i]-((0.5^(1/length(y)))*resvar[,i]/2)) #todo multiply gain alpha to s
	}
	p #does not consider std or mse. no EI criterion. just an aggregation.
}

##############################################################
#' Hypervolume Lower Confidence Bound Infill Criterion
#'
#' This multi objective infill criterion is similar to the SMS-EGO infill criterion by Ponweiser (2008).
#' It aggregates the objective values for each point by calculating the hypervolume contribution.
#' As a first step the lower confidence bound is calculated, decreasing the predicted objective values by their predicted variance.
#' Unlike SMS-EGO, epsilon dominance is not employed here. Also, the penalties for dominated points are calculated differently: The hypervolume between the dominated
#' points and the current true Pareto front is used.
#'
#' @note An optimizer like pso will work signif. better than cmaes with this infill criterion.
#'
#' @param resy predicted objective values
#' @param resvar predicted variance
#' @param y the current Pareto front
#' @param ref reference point, if not given will be chosen as maximum of observed values plus one
#'
#' @return returns the contribution (or penalty) for each row in resy
#' 
#' @references W. Ponweiser, T. Wagner, D. Biermann, and M. Vincze. Multiobjective optimization
#' on a limited budget of evaluations using model-assisted -metric selection. In
#' PPSN, pages 784-794, 2008.
#'
#' @export
##############################################################
spotInfillLcbHyperVolume <- function(resy,resvar,y,ref=NULL){
	p=spotInfillLcbMulti(resy,resvar,y,ref)
	nn=nrow(p)
	y= matrix(unlist(y),ncol=length(y))
	p=as.matrix(p)
	colnames(p)=colnames(y)
	pfront=rbind(p,y) #todo better than y would be reevaluated y on model.
	if(is.null(ref)){ref=apply(pfront,2,max)+1}
	result <- NULL
	for(i in 1:nn){
		if(is_dominated(t(pfront))[i])
			result<-c(result,(dominated_hypervolume(t(y),p[i,]))) #penalty for dominated points
		else
			result<-c(result,(dominated_hypervolume(t(y),ref)-dominated_hypervolume(t(rbind(p[i,],y)),ref)))
	}
	result[1:nn] #does not consider std or mse. no EI criterion. just an aggregation.
}

##############################################################
#' MCO Infill Criterion
#'
#' Same as \code{\link{spotInfillLcbHyperVolume}}, but without using lower confidence bound.
#' Simply aggregates contributed hypervolumes based on predicted means, penalizing dominated points.
#'
#' @param resy predicted objective values
#' @param resvar not used
#' @param y the current Pareto front
#' @param ref reference point, if not given will be chosen as maximum of observed values plus one
#'
#' @return returns the contribution (or penalty) for each row in resy
#'
#' @export
##############################################################
spotInfillHyperVolume <- function(resy,resvar,y,ref=NULL){
	p=resy
	nn=nrow(p)
	y= matrix(unlist(y),ncol=length(y))
	p=as.matrix(p)
	colnames(p)=colnames(y)
	pfront=rbind(p,y) #todo better than y would be reevaluated y on model.
	if(is.null(ref)){ref=apply(pfront,2,max)+1}
	result <- NULL
	for(i in 1:nn){
		if(is_dominated(t(pfront))[i])
			result<-c(result,(dominated_hypervolume(t(y),p[i,]))) #penalty for dominated points
		else
			result<-c(result,(dominated_hypervolume(t(y),ref)-dominated_hypervolume(t(rbind(p[i,],y)),ref)))
	}
	result[1:nn] #does not consider std or mse. no EI criterion. just an aggregation.
}

##############################################################
#' Single objective lower confidence bound
#'
#' LCB=mean-sd 
#'
#' @param mean predicted mean values
#' @param sd predicted st. deviation
#' @param min the currently known optimum value
#'
#' @return Returns the difference of mean and sd
#' @export
##############################################################
spotInfillLcbSingle<- function(mean,sd,min){ #NegProbImp from Forrester
	mean-sd #SPOT minimizes, large values are better, so negate it
}

##############################################################
#' Probability of Improvement Infill Criterion
#'
#' This is the single objective infill criterion "Probability of Improvement" as introduced by Forrester (2008).
#'
#' @param mean predicted mean values
#' @param sd predicted st. deviation
#' @param min the currently known optimum value
#'
#' @return Returns the negative of the Probability of Improvement.
#' 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#'
#' @export
##############################################################
spotInfillProbImp<- function(mean,sd,min){ #NegProbImp from Forrester
	-pnorm((min-mean)/sd) #SPOT minimizes, large values are better, so negate it
}

##############################################################
#' Neg. Log. of Expected Improvement Infill Criterion
#'
#' This is the single objective infill criterion "Expected Improvement" as introduced by Forrester (2008).
#'
#' @param mean predicted mean values
#' @param sd predicted st. deviation
#' @param min the currently known optimum value
#'
#' @return Returns the negative log. of the Expected Improvement.
#' 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#'
#' @export
##############################################################
spotInfillExpImp <- function(mean,sd,min){ #NegLogExpImp from Forrester
	EITermOne=(min-mean)*pnorm((min-mean)/sd)
	EITermTwo=sd*(1/sqrt(2*pi))*exp(-(1/2)*((min-mean)^2/(sd^2)))
	-log10(EITermOne+EITermTwo+(.Machine$double.xmin)) #SPOT minimizes, large values are better, so negate it
}

##############################################################
#' Neg. SD Infill Criterion
#'
#' Infill criterion: negative standard deviation only.
#' This type of infill criterion is useful when pure exploration is desired.
#'
#' @param mean predicted mean values
#' @param sd predicted st. deviation
#' @param min the currently known optimum value
#'
#' @return Returns the of \code{sd}.
#'
#' @export
##############################################################
spotInfillSD <- function(mean,sd,min){ 
	-sd
}


##############################################################
#' spotExiPsi
#'
#' helper function for spotSExI2d
#'
#' @param a a
#' @param b b
#' @param m m
#' @param s s
#'
#' @keywords internal
##############################################################
spotExiPsi <- function(a,b,m,s){
  s*dnorm((b-m)/s) + (a-m)*pnorm((b-m)/s)
}

##############################################################
#' S-metric Expected Improvement SExI Infill Criterion
#'
#' This two-objective infill criterion is the Expected Improvement of the S-metric. That is,
#' it aggregates the predicted objective values by an exact calculation of the Expected Improvement EI
#' in hypervolume. As this gets more complex and time-consuming for higher dimensional objective spaces,
#' this is only implemented for the two-objective case. An approximation approach for higher dimensional problems
#' exists, but is not yet implemented in SPOT. 
#'
#' @param resy predicted objective values
#' @param resvar predicted variance
#' @param y the current Pareto front
#' @param ref reference point, if not given will be chosen as maximum of observed values plus one
#'
#' @return returns the EI for each row in resy
#' 
#' @references M. Emmerich, A.H. Deutz, J.W. Klinkenberg: The computation of the expected improvement in dominated hypervolume of Pareto front approximations , LIACS TR-4-2008, Leiden University, The Netherlands 
#'
#' @seealso \code{\link{spotSExI2d}}
#'
#' @export
##############################################################
spotInfillSExI2d <- function(resy,resvar,y,ref=NULL){
	EI=NULL
	y1=y[[1]]
	y2=y[[2]]
	if(is.null(ref)){ref=c(max(y1),max(y2))+1}
	Pfront=data.frame(y1=y1,y2=y2)	
	Pfront=Pfront[!is_dominated(t(Pfront)),]
	#% prediction of each objective
	pred1=resy[,1]
	pred2=resy[,2]
	#% MSE of each objective
	s1=resvar[,1]
	s2=resvar[,2]	
	for(i in 1:length(pred1)){
		EI=c(EI,-spotSExI2d(Pfront,ref,c(pred1[i],pred2[i]),sqrt(c(s1[i],s2[i])))) #sqrt is used to compute stddeviation, that means models are expected to return MSE or Variance.
	}
	EI
}

##############################################################
#' S-metric Expected Improvement SExI Infill Criterion
#'
#' This two-objective infill criterion is the Expected Improvement of the S-metric. That is,
#' it aggregates the predicted objective values by an exact calculation of the Expected Improvement EI
#' in hypervolume. As this gets more complex and time-consuming for higher dimensional objective spaces,
#' this is only implemented for the two-objective case. An approximation approach for higher dimensional problems
#' exists, but is not yet implemented in SPOT. 
#'
#' @param P Approximation set: provide f1,f2 coordinates of current Pareto front approximation 
#' @param r Reference point: used for computing the hypervolume
#' @param mu  Mean vector: mean value of predictive distribution (e.g. from Gaussian process), f1, f2
#' @param s Standard deviations of predictive distribution
#'
#' @return returns the EI for each row in resy
#'
#' @author (c) Michael Emmerich and Andre Deutz, LIACS, Leiden University, 2010 \cr
#'    \email{emmerich@@liacs.nl}, \email{deutz@@liacs.nl} \cr
#'    R port by Patrick Koch, Cologne University of Applied Sciences \cr
#'	   \email{patrick.koch@@fh-koeln.de}
#'
#' @references M. Emmerich, A.H. Deutz, J.W. Klinkenberg: The computation of the expected improvement in dominated hypervolume of Pareto front approximations , LIACS TR-4-2008, Leiden University, The Netherlands 
#'
#' @seealso \code{\link{spotInfillSExI2d}}
#'
#' @examples
#' print(spotSExI2d(data.frame(x1=c(0,1,2),x2=c(2,1,0)),c(3,3),c(0,0),c(0.1,0.1)))
#' ##should be approx. 3.08
#' print(spotSExI2d(data.frame(x1=c(1,2),x2=c(2,1)),c(11,11),c(10,10),c(4,4)))
#' ##should be approx. 0.0726
#'
#' @export
##############################################################
spotSExI2d <- function(P,r,mu,s){

  S = P[order(P[1]),]
  k = nrow(S)
  #contr = data.frame()
  contr = matrix(0,k+1,k+1) #avoid dynamic extension of contr - this costs time

  c2 = sort(S[,2])
  c1 = sort(S[,1])
 
  cL1=c(-Inf,c1)
  cL2=c(-Inf,c2)
  cU1=c(c1,r[1])
  cU2=c(c2,r[2])
  fMax1=c(r[1],rev(c1))
  fMax2=c(r[2],rev(c2))
  
  #Cumulative Gaussian over length for correction constant
  GaussCDF1 = pnorm((cU1-mu[1])/s[1]) - pnorm((cL1-mu[1])/s[1])
  #Cumulative Gaussian over length for correction constant
  GaussCDF2 = pnorm((cU2-mu[2])/s[2]) - pnorm((cL2-mu[2])/s[2])
  #Marginal integration over the length of a cell
  Psi1 = outer(fMax1, 0:k, function(fMax1,j) {spotExiPsi(fMax1,cU1[j+1],mu[1],s[1]) - spotExiPsi(fMax1,cL1[j+1],mu[1],s[1]) } )
  #Marginal integration over the height of a cell
  Psi2 = outer(0:k, fMax2, function(i,fMax2) {spotExiPsi(fMax2,cU2[i+1],mu[2],s[2]) - spotExiPsi(fMax2,cL2[i+1],mu[2],s[2]) } )

  for (i in 0:k){  # hight from below to above
    for (j in 0:(k-i)){ # first coordinate (length) of cell grid
        SM = S[ which(cU1[j+1] <= S[,1] & cU2[i+1] <= S[,2]) , ]
		sPlus=if(nrow(SM)==0){0    #MZ: the hypervolume function from emoa package does not allow empty matrices, but its faster
		}else{dominated_hypervolume(t(SM), c(fMax1[i+1],fMax2[j+1]))}
        #ExI Kontribution fuer die aktuelle Zelle
        contr[i+1,j+1]= Psi1[i+1,j+1]*Psi2[i+1,j+1]-sPlus*GaussCDF1[j+1]*GaussCDF2[i+1]
    }
  }
  sum(contr)
}
