# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the lm function from the Stats package the lme function of the nlme package
# and functions from the lmeSpline, reshape and gdata packages
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Moleculesral Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Moleculesral Public License for more details.
#
# You should have received a copy of the GNU Moleculesral Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#' Data-driven linear mixed effect model spline modelling
#' 
#' Function that models a linear or limear mixed model depending on the best fit. Alternatively, the function can return THE derivation information of the fitted models
#' for the fixed (original) times points and a chosen \code{basis}.
#' 
#' @import methods 
#' @importFrom nlme lme lmeControl pdIdent pdDiag
#' @importFrom parallel parLapply detectCores makeCluster clusterExport stopCluster
#' @importFrom gdata drop.levels
#' @importFrom lmeSplines smspline approx.Z
#' @importFrom reshape2 melt dcast
#' @importFrom stats lm predict.lm predict anova quantile na.exclude
#' @usage lmmSpline(data, time, sampleID, timePredict, deri, basis, knots, keepModels,numCores)
#' @param data \code{data.frame} or \code{matrix} containing the samples as rows and features as columns
#' @param time \code{numeric} vector containing the sample time point information.
#' @param sampleID \code{character}, \code{numeric} or \code{factor} vector containing information about the unique identity of each sample
#' @param timePredict \code{numeric} vector containing the time points to be predicted.  By default set to the original time points observed in the experiment.
#' @param deri \code{logical} value. If \code{TRUE} returns the predicted derivative information on the observed time points.By default set to \code{FALSE}.
#' @param basis \code{character} string. What type of basis to use, matching one of \code{"cubic"}, \code{"p-spline"} or \code{"cubic p-spline"}. The \code{"cubic"} basis (\code{default}) is the cubic smoothing spline as defined by Verbyla \emph{et al.} 1999, the \code{"p-spline"} is the truncated p-spline basis as defined by Durban \emph{et al.} 2005.
#' @param knots Alternatively an \code{integer}, the number of knots used for the \code{"p-spline"} or \code{"cubic p-spline"} basis calculation. Otherwise calculated as proposed by Ruppert 2002. Not used for the "cubic" smoothing spline basis as it used the inner design points.
#' @param keepModels alternative \code{logical} value if you want to keep the model output. Default value is FALSE
#' @param numCores Alternative \code{numeric} value indicating the number of CPU cores to be used. Default value is automatically estimated.
#' @details  
#' The first model (\code{modelsUsed}=0) assumes the response is a straight line not affected by individual variation. 
#' 
#' Let \eqn{y_{ij}(t_{ij})} be the expression of a feature for individual (or biological replicate) \eqn{i} at time \eqn{t_{ij}}, where \eqn{i=1,2,...,n}, \eqn{j=1,2,...,m_i}, \eqn{n} is the sample size and \eqn{m_i} is the number of observations for individual \eqn{i} for the given feature. 
#' We fit a simple linear regression of expression \eqn{y_{ij}(t_{ij})} on time \eqn{t_{ij}}. 
#' The intercept \eqn{\beta_0} and slope \eqn{\beta_1} are estimated via ordinary least squares:
#' \eqn{y_{ij}(t_{ij})= \beta_0 + \beta_1 t_{ij} + \epsilon_{ij}}, where \eqn{\epsilon_{ij} ~ N(0,\sigma^2_{\epsilon}).}
#' The second model (\code{modelsUsed}=1) is nonlinear where the straight line in regression replaced with a curve modelled using here for example a spline truncated line basis (\code{basis}="p-spline") as proposed Durban \emph{et al.} 2005:
#' 
#' \deqn{y_{ij}(t_{ij})= f(t_{ij}) +\epsilon_{ij},} 
#' 
#' where \eqn{\epsilon_{ij}~ N(0,\sigma_{\epsilon}^2).}
#' 
#' The penalized spline is represented by \eqn{f}, which depends on a set of knot positions \eqn{\kappa_1,...,\kappa_K} in the range of \eqn{{t_{ij}}}, some unknown coefficients \eqn{u_k}, an intercept \eqn{\beta_0} and a slope \eqn{\beta_1}. The first term in the above equation can therefore be expanded as:
#' \deqn{f(t_{ij})= \beta_0+ \beta_1t_{ij}+\sum\limits_{k=1}^{K}u_k(t_{ij}-\kappa_k)_+,}
#' with \eqn{(t_{ij}-\kappa_k)_+=t_{ij}-\kappa_k}, if \eqn{t_{ij}-\kappa_k  > 0, 0} otherwise.
#' 
#' The choice of the number of knots \eqn{K} and their positions influences the flexibility of the curve. 
#' If the argument \code{knots}=missing, we use a method proposed by Ruppert 2002 to estimate the number of knots given the measured number of time points \eqn{T}, so that the knots \eqn{\kappa_1 \ldots \kappa_K} are placed at quantiles of the time interval of interest: 
#'
#' \deqn{K= max(5,min(floor(\frac{T}{4}) , 40)).}
#' 
#' In order to account for individual variation, our third model (\code{modelsUsed}=2) adds a subject-specific random effect \eqn{U_i} to the mean response \eqn{f(t_{ij})}. 
#' Assuming \eqn{f(t_{ij})} to be a fixed (yet unknown) population curve, \eqn{U_i} is treated as a random realization of an underlying Gaussian process with zero-mean and variance \eqn{\sigma_U^2} and is independent from the random error \eqn{\epsilon_{ij}}:
#' 
#' \deqn{y_{ij}(t_{ij}) = f(t_{ij}) + U_i + \epsilon_{ij}}
#' 
#' with \eqn{U_{i} ~ N(0,\sigma_U^2)} and \eqn{\epsilon_{ij} ~ N(0,\sigma_{\epsilon}^2)}.

#' In the equation above, the individual curves are expected to be parallel to the mean curve as we assume the individual expression curves to be constant over time.
#' A simple extension to this model is to assume individual deviations are straight lines. The fourth model (\code{modelsUsed}=3) therefore fits individual-specific random intercepts \eqn{a_{i0}} and slopes \eqn{a_{i1}}:
#' 
#'  \deqn{y_{ij}(t_{ij}) = f(t_{ij}) + a_{i0} + a_{i1}t_{ij} + \epsilon_{ij}}
#'  
#' with \eqn{\epsilon_{ij} ~ N(0,\sigma_\epsilon^2)} and \eqn{(a_{i0},a_{i1})^T} ~ \eqn{ N(0,\Sigma).}
#' We assume independence between the random intercept and slope.
#'  @return lmmSpline returns an object of class \code{lmmspline} containing the following components:
#'  \itemize{
#' \item{predSpline}{\code{data.frame} containing predicted values based on linear model object or linear mixed effect model object.}
#' \item{modelsUsed}{\code{numeric} vector indicating the model used to fit the data. 0 = linear model, 1=linear mixed effect model spline (LMMS) with defined basis ('cubic' by default) 2 = LMMS taking subject-specific random intercept, 3 = LMMS with subject specific intercept and slope.}
#' \item{model}{\code{list} of models used to model time profiles.}
#' \item{derivative}{\code{logical} value indicating if the predicted values are the derivative information.}
#'  }
#' @references  Durban, M., Harezlak, J., Wand, M. P., & Carroll, R. J. (2005). \emph{Simple fitting of subject-specific curves for longitudinal data.} Stat. Med., 24(8), 1153-67.
#' @references  Ruppert, D. (2002). \emph{Selecting the number of knots for penalized splines.} J. Comp. Graph. Stat. 11, 735-757
#' @references  Verbyla, A. P., Cullis, B. R., & Kenward, M. G. (1999). \emph{The analysis of designed experiments and longitudinal data by using smoothing splines.} Appl.Statist, 18(3), 269-311.
#' @references  Straube J., Gorse A.-D., Huang B.E., Le Cao K.-A. (2015).  \emph{A linear mixed model spline framework for analyzing time course 'omics' data} PLOSONE, 10(8), e0134540.
#' @seealso \code{\link{summary.lmmspline}}, \code{\link{plot.lmmspline}}, \code{\link{predict.lmmspline}}, \code{\link{deriv.lmmspline}}
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' # running for samples in group 1
#' G1 <- which(kidneySimTimeGroup$group=="G1")
#' testLMMSpline<- lmmSpline(data=kidneySimTimeGroup$data[G1,],time=kidneySimTimeGroup$time[G1],
#'                  sampleID=kidneySimTimeGroup$sampleID[G1])
#' summary(testLMMSpline)
#' DerivTestLMMSplineTG<- lmmSpline(data=as.data.frame(kidneySimTimeGroup$data[G1,]),
#'                        time=kidneySimTimeGroup$time[G1],sampleID=kidneySimTimeGroup$sampleID[G1],
#'                        deri=TRUE,basis="p-spline")
#' summary(DerivTestLMMSplineTG)}
# setGeneric('lmmSpline',function(data,time,sampleID,timePredict,deri,basis,knots,keepModels,numCores){standardGeneric('lmmSpline')})
# setClassUnion("matrixOrFrame",c('matrix','data.frame'))
# setClassUnion("missingOrnumeric", c("missing", "numeric"))
# setClassUnion("missingOrcharacter", c("missing", "character"))
# setClassUnion("missingOrlogical", c("missing", "logical"))
# setClassUnion("factorOrcharacterOrnumeric", c("factor", "character","numeric"))
# # @rdname lmmSpline-methods
# # @aliases lmmSpline,matrixOrFrame,numeric,factorOrcharacterOrnumeric,
# # missingOrlogical,missingOrcharacter,missingOrnumeric,missingOrlogical,missingOrnumeric-method
# # @exportMethod lmmSpline
# 
# setMethod('lmmSpline',c(data="matrixOrFrame",time="numeric",sampleID="factorOrcharacterOrnumeric",timePredict="missingOrnumeric", deri="missingOrlogical", basis="missingOrcharacter",knots="missingOrnumeric",keepModels="missingOrlogical",numCores="missingOrnumeric"), function(data,time,sampleID,timePredict,deri,basis,knots,keepModels,numCores){
#   
#    lmmSplinePara(data=data,time=time,sampleID=sampleID,timePredict=timePredict,deri=deri,basis=basis, knots=knots,keepModels=keepModels,numCores=numCores)
# })
# @name lmmSpline

#' @docType methods
#' @rdname lmmSpline-methods
#' @export
lmmSpline <- function(data, time, sampleID, timePredict, deri, basis, knots,keepModels, numCores){

  if(missing(keepModels))
    keepModels <- F
  if(missing(timePredict))
    timePredict <- sort(unique(time))
  if(missing(basis))
    basis <- "cubic"
  
  if(missing(deri)){
    deri <- FALSE
  }else{
    deri <- deri
  }
 
  basis.collection <-  c("cubic","p-spline","cubic p-spline")
  if(!basis%in% basis.collection)
    stop(cat("Chosen basis is not available. Choose:", paste(basis.collection,collapse=', ')))
  if(diff(range(c(length(sampleID),length(time),nrow(data))))>0)
    stop("Size of the input vectors rep, time and nrow(data) are not equal")
  if(missing(knots)& (basis=="p-spline"|basis=='cubic p-spline'))
    warning("The number of knots is automatically estimated")
  if(deri & basis=='cubic')
    stop('To calculate the derivative choose either "p-spline" or "cubic p-spline" as basis')
  
  options(show.error.messages = TRUE) 
  
  i <- NULL
  fits <- NULL
  error <- NULL

  if(missing(numCores)){
    num.Cores <- detectCores()
  }else{
    num.Cores <- detectCores()
    if(num.Cores<numCores){
      warning(paste('The number of cores is bigger than the number of detected cores. Using the number of detected cores',num.Cores,'instead.'))
    }else{
      num.Cores <- numCores
    }
    
  }
  Molecule <- ''

  derivLme <- function(fit){ 
    #random slopes
    
    if(class(fit)=='lm'){
      beta.hat <- rep(fit$coefficients[2],length(unique(fit$model$time)))
      return(beta.hat)
      
    }else if(class(fit)=='lme'){
      u <- unlist(fit$coefficients$random$all)
      beta.hat <- fit$coefficients$fixed[2]
      Zt <-  fit$data$Zt[!duplicated(fit$data$Zt),]>0
      deriv.all <-    beta.hat + rowSums(Zt%*%t(u)) 
      return(deriv.all)
    }
  }
  
  #penalized cubic
  
  derivLmeCubic <- function(fit){ 
    #random slopes
    if(class(fit)=='lm'){
      beta.hat <- rep(fit$coefficients[2],length(unique(fit$model$time)))
      return(beta.hat)
      
    }else if(class(fit)=='lme'){
      u <- unlist(fit$coefficients$random$all)
      beta.hat <- fit$coefficients$fixed[2]
      PZ <-  fit$data$Zt[!duplicated(fit$data$Zt),]
      PZ <-PZ^(1/3)
      deriv.all <-    beta.hat + rowSums((PZ*PZ)%*%(t(u)*3)) 
      return(deriv.all)
    }
    
  }
  
  if(missing(knots))
    knots <-NULL
  nMolecules <- NULL
  nMolecules <- ncol(data)
  
  
  lme <- nlme::lme
  cl <- makeCluster(num.Cores,"SOCK")
  clusterExport(cl, list('data','lm','try','class','unique','anova','drop.levels','pdDiag','pdIdent','time','sampleID','melt','dcast','predict','derivLme','knots','derivLmeCubic','lme','keepModels','basis','data','other.reshape'),envir=environment())

  models <-list()

  
  new.data <- parLapply(cl,1:nMolecules,fun = function(i){
    expr <- data[,i]
    
    dataM <- as.data.frame(other.reshape(Rep=sampleID,Time=time,Data=unlist(expr)))
    dataM$all = rep(1, nrow(dataM))
    dataM$time = as.numeric(as.character(dataM$Time))
    dataM$Expr = as.numeric(as.character(dataM$Expr))
    
    
    #### CUBIC SPLINE BASIS ####
    if(basis=="cubic"){
      dataM$Zt <- lmeSplines::smspline(~ time, data=dataM)
      knots <- sort(unique(time))[2:(length(unique(time))-1)]
    }
    #### PENALIZED SPLINE BASIS#####
    if(basis%in%c("p-spline","cubic p-spline")){
      
      if(is.null(knots)){
        K <- max(6,min(floor(length(unique(dataM$time))/4),40))
      }else{
        K <- max(knots,6)
      }
      knots <- quantile(unique(dataM$time),seq(0,1,length=K+2))[-c(1,K+2)]
      if(min(knots)<=min(dataM$time) | max(knots)>=max(dataM$time))
        stop(cat('Make sure the knots are within the time range',range(dataM$time)[1],'to',range(dataM$time)[2]))
      PZ <- outer(dataM$time,knots,"-")
      if(basis=="cubic p-spline")
        PZ <- PZ^3
      PZ <- PZ *(PZ>0)
      dataM$Zt <- PZ 
      
    }
    

    
    if(deri){
      pred.spline = rep(NA,length(timePredict))
    }else{
      pred.spline =rep(NA,length(timePredict))
      pred.df <- data.frame(all=rep(1,length(timePredict)), time=timePredict)
      pred.df$Zt = lmeSplines::approx.Z(dataM$Zt, dataM$time, timePredict)
      
    }
    
    
    
    #library(nlme)
    fit0 <- NULL
    fit0  <- try(lm(Expr ~ time, data=dataM ))
    if(class(fit0) == 'try-error') {
      models <- list()
      error <- i
      pred.spline <- rep(NA,length(timePredict))
      fits <- NA
    }else{
    fit1 <- NULL
    fit1 <- try(lme(Expr ~ time, data=dataM, random=list(all=pdIdent(~Zt - 1)),
                    na.action=na.exclude, control=lmeControl(opt = "optim"))) 
    pvalue <-1
    if(class(fit1) != 'try-error') { 
      
      pvalue <- anova(fit1, fit0)$'p-value'[2]
      
    }
    
    if(pvalue <= 0.05){  
           
            fit2 <- NULL
            fit2 <- try(lme(Expr ~ time, data=dataM, 
                      random=list(all=pdIdent(~Zt - 1), Rep=pdIdent(~1)), 
                      na.action=na.exclude, control=lmeControl(opt = "optim")))
      
          if(class(fit2) != 'try-error') {  # to prevent errors stopping the loop
        
              pvalue = anova(fit1, fit2)$'p-value'[2]
          }else{ 
            pvalue=1
          }
      
          if(pvalue <= 0.05){  
                fit3 <-NULL
            fit3 <- try(lme(Expr ~ time, data=dataM, 
                        random=list(all=pdIdent(~Zt - 1), Rep=pdDiag(~time)), 
                        na.action=na.exclude, control=lmeControl(opt = "optim")))  
        
            if(class(fit3) != 'try-error') {  # to prevent errors stopping the loop
              pvalue = anova(fit2, fit3)$'p-value'[2]
        
            }else{
              pvalue=1
            }
              if(pvalue <= 0.05){
                    fits <- 3
                    models<- fit3
                    if(deri){
                      if(basis=='p-spline')
                        pred.spline = derivLme(fit3)
                      if(basis=='cubic p-spline')
                        pred.spline = derivLmeCubic(fit3)
                      
                    }else{
                      pred.spline = predict(fit3, newdata=pred.df, level=1, na.action=na.exclude)
                    }
                  }else{ # choose simpler model: fit2
                    fits <- 2
                    models <- fit2
                    if(deri){
                      if(basis=='p-spline')
                        pred.spline = derivLme(fit2)
                      if(basis=='cubic p-spline')
                        pred.spline = derivLmeCubic(fit2)
                    }else{
                      pred.spline = predict(fit2, newdata=pred.df, level=1, na.action=na.exclude)
                    }
                  } 
            
          }else{ 
            models <- fit1
            fits <- 1
            if(deri){
              if(basis=='p-spline')
                pred.spline = derivLme(fit1)
              if(basis=='cubic p-spline')
                pred.spline = derivLmeCubic(fit1)
            }else{
              pred.spline = predict(fit1, newdata=pred.df, level=1, na.action=na.exclude)
            }
          }
    }else{ 
     
      models <- fit0
      fits <-0
      if(deri){
        pred.spline = rep(fit0$coefficients[2],length(unique(dataM$time)))    
      }else{
        
        pred.spline = predict(fit0, newdata=pred.df, level=1, na.action=na.exclude)
      }
    }
    }
    if(!keepModels)
      keepModels <- list()
    
    return(list(pred.spl=pred.spline,fit=fits,models=models,error=error,knots=knots))
    
})
  stopCluster(cl)
  knots <- sort(unique(as.vector((sapply(new.data,'[[','knots')))))
  pred.spl <- matrix(sapply(new.data,'[[','pred.spl'),nrow=nMolecules,byrow=T)
  fits <-  unlist(sapply(new.data,'[[','fit'))
  error <-  unlist(sapply(new.data,'[[','error'))
  models <-list()
  if(keepModels){
    models <- sapply(new.data,'[[','models')

    if(is.matrix(models))
      models <- sapply(new.data,'[','models')
    }
  
  pred.spl = as.data.frame(pred.spl)
  MolNames <- as.character(unlist(colnames(data)))
  
  if(is.null(MolNames)| sum(is.na(MolNames))>0)
    MolNames <- 1:nrow(pred.spl)
  if(nrow(pred.spl)==length(MolNames))
    rownames(pred.spl)<-MolNames
  if(ncol(pred.spl)==length(timePredict))
    colnames(pred.spl) <- timePredict
  error2 <- "All features were modelled"
  if(length(error)>0){
    warning('The following features could not be fitted ',paste(MolNames[error],' ',sep='\n'))
    error2 <- ''
    error2 <- MolNames[error]
    
  }

  l <-new('lmmspline',predSpline=pred.spl,modelsUsed=fits,basis=basis,knots=knots,errorMolecules=error2,models=models, derivative=deri)
  return(l)
           
}

     
other.reshape <- function(Rep, Time, Data){
  lme.data<-NULL
  lme.data <- data.frame(Time=Time,Rep=Rep,as.matrix(Data))
  lme.data$Time = factor(drop.levels(lme.data$Time))
  lme.data$Rep = factor(drop.levels(lme.data$Rep))
  melt.lme.data <-NULL
  melt.lme.data <- melt(lme.data)
  cast.lme.data  <- NULL
  cast.lme.data <- dcast(melt.lme.data, variable+ Rep ~ Time)
  melt.lme.data2 <- NULL
  melt.lme.data2 <-  melt(data.frame(cast.lme.data))
  names(melt.lme.data2) <- c("Molecule",  "Rep", "Time", "Expr")
  melt.lme.data2$Time <- factor(gsub("^X", "", as.character(melt.lme.data2$Time)))
  return(as.data.frame(melt.lme.data2))
}
