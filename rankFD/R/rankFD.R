#' Rank-based tests for general factorial designs
#' 
#' The rankFD() function calculates the Wald-type statistic (WTS) and the ANOVA-type 
#' statistic (ATS) for general factorial designs.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'    contains the response variable and the right hand side contains the factor
#'    variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'    \code{formula}. The default option is \code{NULL}.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param CI.method Either "Logit" or "Normal", specifying the method used for
#'    calculation of the confidence intervals.
#' @param effect Should the weighted or unweighted effects be calculated?
#' @param hypothesis The null hypothesis to test, either "H0F" or "H0p".
#' @param Factor.Information Logical. If TRUE, descriptive statistics for 
#'    the different factor level combinations is printed.
 
#' @details The package provides the Wald-type as well as the ANOVA-type statistic 
#'    for rank-based factorial designs, i.e., even for ordinal data.  It is
#'    implemented for crossed designs and allows
#'    for an arbitrary number of factor combinations as well as different sample
#'    sizes.
#'   
#' @return An \code{rankFD} object containing the following components:
#' \item{Descriptive}{Some descriptive statistics of the data for all factor
#'    level combinations. Displayed are the number of individuals per factor
#'    level combination (size), the relative effect (pd), variance and 100*(1-alpha)\% confidence
#'    intervals.}
#'  \item{WTS}{The value of the WTS along with degrees of freedom of the central chi-square distribution and p-value}
#'  \item{ATS}{The value of the ATS, degrees of freedom of the central F distribution and the corresponding p-value.}
#' 
#' @examples
#' data(Coal)
#' model <- rankFD(Acidity ~ NaOH * Type, data = Coal, CI.method = "Normal",
#'  effect = "unweighted", hypothesis = "H0F")
#' 
#' data(Muco)
#' model.oneway <- rankFD(HalfTime ~ Disease, data = Muco, CI.method = "Logit",
#'  effect = "weighted", hypothesis = "H0p")
#' 
#' 
#' @references Friedrich, S., Konietschke, F., Pauly, M.(2015) GFD - An R-package
#' for the Analysis of General Factorial Designs. Submitted to Journal of Statistical Software.
#' 
#' Pauly, M., Brunner, E., Konietschke, F.(2015) Asymptotic Permutation Tests in General Factorial Designs. Journal of the Royal Statistical Society - Series B \bold{77}, 461--473.
#'
#' @importFrom lattice xyplot panel.superpose panel.arrows panel.points panel.xyplot
#' @importFrom stats formula model.frame pchisq pf qnorm terms var aggregate as.formula confint cov median pnorm pt qt quantile sd
#' @importFrom utils read.table
#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix
#' @importFrom coin wilcox_test pvalue statistic
#' @importFrom graphics abline axis box plot points polygon title
#' 
#' @export


rankFD <- function(formula, data,alpha=0.05, CI.method=c("Logit","Normal"),
                     effect=c("unweighted","weighted"),hypothesis=c("H0F","H0p"),
                    Factor.Information=FALSE){
  
  effect=match.arg(effect)
  hypothesis = match.arg(hypothesis)
  CI.method = match.arg(CI.method)

  #-------- Determine the model------#
  
  dat.Model0 <- model.frame(formula, data)
  
  #--------Numbers of factors---------#
  
  for(i in length(dat.Model0):2) dat.Model0[,i] <-as.factor(dat.Model0[,i])
  nf = ncol(dat.Model0)-1
  
  n.levels = c()
  names.levels=list()
  
  for (i in 2:(nf+1)){n.levels[i-1]= nlevels(dat.Model0[,i])
                      names.levels[[i-1]] = levels(dat.Model0[,i])
  }
  names(names.levels) <- names(dat.Model0[,2:(nf+1)])
  
  
  #-------------Hypotheses matrices--------#

  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  nr_hypo <- attr(terms(formula), "factors")
  fac_names <- colnames(nr_hypo)
  Hypotheses0 = HC(n.levels,"Hyp",perm_names,fac_names)
  Hypotheses=Hypotheses0[[1]]
  CI.Matrices = HC(n.levels,"CI",perm_names,fac_names)[[1]]
  
  
  n.hypotheses = length(Hypotheses)
  n.levels.prod=prod(n.levels)
  #Output.names <- attr(terms(formula), "term.labels")
  Output.names <- Hypotheses0[[2]]
  
  #-------Sort Data according to Factors------#
  for (i in length(dat.Model0):2) {
    dat.Model0 <- dat.Model0[order(dat.Model0[,i]),]}
  
  #-----------------Introduce Pseudo Factor -----#
  dat.Model0$Interaction = interaction(dat.Model0[,2:length(dat.Model0)],sep=":")
  dat.response <- dat.Model0[,1]
  
  #----------------Compute means and variances etc---------#
  n <- aggregate(formula,data=dat.Model0,length)
  for(i in (length(n)-1):1) {
    n <-n[order(n[,i]),]
  }
  colnames(n)[ncol(n)]<-"Size"
  
  #-------------Compute Inference Methods------#
  dat.Model0$INum <- as.factor(rep(1:n.levels.prod, n$Size))
  
  
  #--------------Compute the relative Effects----#
  H0pW <-Effects(dat.response, dat.Model0$INum,effect)
  n$pd <- c(H0pW$pd)
  n$Var <- c(diag(H0pW$VBF))
  
  
  CI <- Limits(c(H0pW$pd),H0pW$VBF,alpha,c(H0pW$N))
  n$L.Normal <- CI$Normal[,1]
  n$U.Normal <- CI$Normal[,2]
  n$L.Logit <- CI$Logit[,1]
  n$U.Logit <- CI$Logit[,2]
  
  WTS = matrix(0,n.hypotheses,3)
  ATS = matrix(0,n.hypotheses, 4)
  ATSp = matrix(0,n.hypotheses,4)
  Descriptive.Factors = list()
  Levels.Factors = list()
  
  
  for(i in 1:n.hypotheses){
    
    
    WTS[i,] = Wald(c(H0pW$pd),Hypotheses[[i]],H0pW$VH0F)
    ATS[i,] =ANOVATYP(c(H0pW$pd),Hypotheses[[i]],H0pW$VH0F,n$Size)
    ATSp[i,]= ANOVATYPH0P(c(H0pW$pd),Hypotheses[[i]],H0pW$VBF,n$Size,H0pW$dfATS)
    CILimits <-Limits(c(CI.Matrices[[i]]%*%c(H0pW$pd)),CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]]),alpha,H0pW$N)
    
    Descriptives <-data.frame(pd=CI.Matrices[[i]]%*%n$pd,
                              Var= c(diag(CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]]))),
                              L.Normal=CILimits$Normal[,1],
                              U.Normal = CILimits$Normal[,2],
                              L.Logit = CILimits$Logit[,1],
                              U.Logit = CILimits$Logit[,2])
    Output.namesi <-Output.names[i] 
    formula.act <- as.formula(paste(names(dat.Model0)[1], Output.namesi, sep=" ~ "))
    aha <- data.frame(aggregate(formula.act,data=dat.Model0,mean))
    
    for(ii in (length(aha)-1):1) {aha <-aha[order(aha[,ii]),]}
    Descriptive.Factors[[i]] <-data.frame(aha,Descriptives)
    Descriptive.Factors[[i]] <-Descriptive.Factors[[i]][,-length(aha)]
    if(length(grep(":", Output.names[i]))<1){
      pos <- which(names(dat.Model0)==Output.names[i])
      Levels.Factors[[i]] = data.frame(X=levels(dat.Model0[,pos]))
    }
    
    if(length(grep(":", Output.names[i]))>=1){
      facs.singles <- c(strsplit(Output.names[i], ":")[[1]])
      
      Levels.Factors[[i]] = data.frame(n[,facs.singles])
      
      
    }
  }
  
  names(Descriptive.Factors) <- Output.names
  
  rownames(WTS) <- Output.names
  rownames(ATS) <- Output.names
  rownames(ATSp)<- Output.names
  colnames(WTS) <- c("Statistic", "df", "p-Value")
  colnames(ATS) <- c("Statistic", "df1", "df2", "p-Value")
  colnames(ATSp) <- c("Statistic", "df1", "df2", "p-Value")

  
  if(hypothesis=="H0F"){
    result <- list(Descriptive=n, Wald.Type.Statistic = WTS, ANOVA.Type.Statistic=ATS)
  }
  
  if(hypothesis=="H0p"){
    result <- list(Descriptive=n,  ANOVA.Type.Statistic=ATSp)
  } 
  class(result) <- "rankFD"

  result$plotting <- list(nf = nf, fac_names = fac_names, n.hypotheses = n.hypotheses,
                          Descriptive.Factors = Descriptive.Factors, CI.method = CI.method)
  if (Factor.Information==TRUE){
    print(Descriptive.Factors)
  }
  return(result)
}

