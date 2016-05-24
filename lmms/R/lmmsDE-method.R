# Jasmin Straube, Queensland Facility of Advanced Bioinformatics
# Part of this script was borrowed from the lm function from the Stats package the lme function of the nlme package
# and functions from the lmeSpline, reshape, parallel, snow and gdata packages
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

#' Differential expression analysis using linear mixed effect model splines.
#' 
#' Function to fit a linear mixed effect model splines to perform differential expression analysis. The \code{\link{lmmsDE}} function fits LMM models with either a \code{cubic}, \code{p-spline} or \code{cubic p-spline} basis and compares the models to the null models. The type of basis to use is specified with the \code{basis} argument.
#' 
#' @import methods 
#' @importFrom nlme lme lmeControl pdIdent pdDiag
#' @importFrom parallel parLapply detectCores makeCluster clusterExport stopCluster
#' @importFrom gdata drop.levels
#' @importFrom lmeSplines smspline 
#' @importFrom reshape2 melt dcast
#' @importFrom  stats na.omit anova quantile fitted na.exclude predict predict.lm lm p.adjust
#' @usage lmmsDE(data, time, sampleID, group, type,
#' experiment, basis, knots,keepModels, numCores)
#' @param data \code{data.frame} or \code{matrix} containing the samples as rows and features as columns
#' @param time \code{numeric} vector containing the sample time point information.
#' @param sampleID \code{character}, \code{numeric} or \code{factor} vector containing information about the unique identity of each sample
#' @param group \code{character}, \code{numeric} or \code{factor} vector containing information about the group (or class) of each sample
#' @param type \code{character} indicating what type of analysis is to be performed. Options are \code{"time"} to identify differential expression over time, \code{"group"} to identify profiles with different baseline levels (intercepts), and \code{"time*group"} an interaction between these two . Use \code{"all"} to calculate all three types.
#' @param experiment \code{character} describing the experiment performed for correlation handling. Use \code{"all"} for data-driven selection of model; \code{"timecourse"} for replicated experiments with less variation in individual expression values (e.g. model organism, cell culture), \code{"longitudinal1"} for different intercepts and \code{"longitudinal2"} for different intercepts and slopes. 
#' @param basis \code{character} string. What type of basis to use, matching one of \code{"cubic"} smoothing spline as defined by Verbyla \emph{et al.} 1999, \code{"p-spline"} Durban \emph{et al.} 2005 or a \code{"cubic p-spline"}.
#' @param knots can take an integer value corresponding to the number of knots for the chosen basis or by default calculated as  in Ruppert 2002. Not in use for the 'cubic' smoothing spline basis.
#' @param numCores alternative \code{numeric} value indicating the number of CPU cores to be used for parallelization. Default value is automatically estimated.
#' @param keepModels alternative \code{logical} value if you want to keep the model output. Default value is FALSE
#' @details
#' lmmsDE extends the LMMS modelling framework to permit tests between groups, across time, and for interactions between the two implemented as described in Straube \emph{et al.} 2015. 
#' @return lmmsDE returns an object of class \code{lmmsde} containing the following components:
#' \item{DE}{\code{data.frame} returning p-values and adjusted p-values using Benjamini-Hochberg correction for multiple testing of the differential expression testing over time, group or their interaction.}
#' \item{modelsUsed}{\code{numeric} vector indicating the model used to fit the data. 1=linear mixed effect model spline (LMMS) with defined basis ('cubic' by default) 2 = LMMS taking subject-specific random intercept, 3 = LMMS with subject specific intercept and slope.}
#' \item{predTime}{\code{data.frame} containing predicted values based on linear model object or linear mixed effect model object.}
#' \item{predGroup}{\code{data.frame} containing predicted values based on linear model object or linear mixed effect model object.}
#' \item{predTime}{\code{data.frame} containing predicted values based on linear model object or linear mixed effect model object.}
#' \item{predTimeGroup}{\code{data.frame} containing predicted for the time*group model values based on linear model object or linear mixed effect model object.}
#' \item{modelTime}{a \code{list} of class \code{\link{lme}}, containing the models for every feature modelling the time effect.} 
#' \item{modelGroup}{a \code{list} of class \code{\link{lme}}, containing the models for every feature modelling group effect. }
#' \item{modelTimeGroup}{a \code{list} of class \code{\link{lme}}, containing the models for every feature modelling time and group interaction effect. }
#' \item{type}{an object of class \code{character}, describing the test performed either time, group, time*group or all. }
#' \item{experiment}{an object of class \code{character} describing the model used to perform differential expression analysis.}
#' @references  Durban, M., Harezlak, J., Wand, M. P., & Carroll, R. J. (2005). \emph{Simple fitting of subject-specific curves for longitudinal data.} Stat. Med., 24(8), 1153-67.
#' @references  Ruppert, D. (2002). \emph{Selecting the number of knots for penalized splines.} J. Comp. Graph. Stat. 11, 735-757
#' @references  Verbyla, A. P., Cullis, B. R., & Kenward, M. G. (1999). \emph{The analysis of designed experiments and longitudinal data by using smoothing splines.} Appl.Statist, 18(3), 269-311.
#' @references  Straube J., Gorse A.-D., Huang B.E., & Le Cao K.-A. (2015). \emph{A linear mixed model spline framework for analyzing time course 'omics' data} PLOSONE, 10(8), e0134540.

#' @seealso \code{\link{summary.lmmsde}}, \code{\link{plot.lmmsde}}
#' @examples 
#' \dontrun{
#' data(kidneySimTimeGroup)
#' lmmsDEtest <-lmmsDE(data=kidneySimTimeGroup$data,time=kidneySimTimeGroup$time,
#'               sampleID=kidneySimTimeGroup$sampleID,group=kidneySimTimeGroup$group)
#' summary(lmmsDEtest)}

# @export
#setGeneric('lmmsDE',function(data,time,sampleID,group,type,experiment,basis,knots,keepModels,numCores){standardGeneric('lmmsDE')})
#setClassUnion("missingOrnumeric", c("missing", "numeric"))
#setClassUnion("missingOrlogical", c("missing", "logical"))
#setClassUnion("missingOrcharacter", c("missing", "character"))
#setClassUnion("factorOrcharacterOrnumeric", c( "factor","character","numeric"))
#setClassUnion("matrixOrframe",c('matrix','data.frame'))
# @rdname lmmsDE-methods
# @aliases lmmsDE lmmsDE,matrixOrframe,numeric,factorOrcharacterOrnumeric,
# factorOrcharacterOrnumeric,missingOrcharacter,missingOrcharacter,missingOrcharacter,
# missingOrnumeric,missingOrlogical,missingOrnumeric-method
# @exportMethod lmmsDE

#setMethod('lmmsDE',c(data="matrixOrframe",time="numeric",sampleID="character",group="character",type="missing",experiment="missing",basis="missing",knots="missing",keepModels="missing",numCores="missing"), function(data,time,sampleID,group,type,experiment,basis,knots,keepModels,numCores){
#  lmmsDEPara(data=data,time=time,sampleID=sampleID,group=group,type=type,experiment=experiment,basis=basis,keepModels=keepModels,numCores=numCores)
#})


#setMethod('lmmsDE',c(data="matrixOrframe",time="numeric",sampleID="factorOrcharacterOrnumeric",group="factorOrcharacterOrnumeric",type="missingOrcharacter",experiment="missingOrcharacter",basis="missingOrcharacter",knots="missingOrnumeric",keepModels="missingOrlogical",numCores="missingOrnumeric"), function(data,time,sampleID,group,type,experiment,basis,knots,keepModels,numCores){
#  lmmsDEPara(data=data,time=time,sampleID=sampleID,group=group,type=type,experiment=experiment,basis=basis,keepModels=keepModels,numCores=numCores)
#})
# @method lmmsDE data.frame
#' @docType methods
#' @rdname lmmsDE-methods
#' @export
lmmsDE <- function(data, time, sampleID, group, type,experiment,basis,knots,keepModels,numCores){

 # require(lmeSplines)
 # require(parallel)
 # require(gdata)
 # require(reshape2)
  model.time <- list()
  model.time.group <- list()
  model.group <- list()
  if(missing(keepModels))
    keepModels <- F
  
  if(!is.logical(keepModels))
    stop('keepModels needs to be of type logical')
  
  if(missing(type))
    type <- 'all'
  
  if(missing(basis))
    basis <- 'cubic'
  
  if(missing(knots))
    knots <- NULL
  
  if(missing(experiment))
    experiment <- 'timecourse'
  
  if(missing(keepModels))
    keepModels <- F
  
  if(type=="time*group")
    type <- "grouptime"
  
  experiment.collection <-  c("timecourse","longitudinal1","longitudinal2","all")
  if(!experiment%in% experiment.collection)
    stop(paste("Chosen experiment is not available. Choose ", paste(experiment.collection,collapse=', ')))
  
  
  basis.collection <-  c("cubic","p-spline","cubic p-spline")
  if(!basis%in% basis.collection)
    stop(paste("Chosen basis is not available. Choose ",paste(basis.collection ,collapse=", ")))
  
  
  type.collection <-  c("all","time","group","grouptime")
  if(!type%in% type.collection)
    stop(paste("Chosen type is not available. Choose ",paste(type.collection,collapse=", ")))
  
  if(diff(range(c(length(sampleID),length(time),length(group),nrow(data))))>0)
    stop("Size of the input vectors sampleID, time, group and ncol(data) are not equal")
  if(is.null(knots)& (basis=="p-spline"|basis=='cubic p-spline'))
    warning("The number of knots is automatically estimated")
  
  
  if(length(unique(group))==1){
    warning("Only single group detected! Performing only time effect test.")
    type <- 'time'
  }
  
  
  nMolecules <- NULL
  nMolecules <- ncol(data)
  
  #levels.genotype <- unique(data$Group)
  pvals <- c()
  k<-0
  
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
  cl <- makeCluster(num.Cores,"SOCK")
  pred.time <- NA
  pred.time.group <- NA
  pred.group <- NA
  models <-list()
  ps <- rep(NA,nMolecules)
  pvalsdf <- matrix(NA,ncol=3,nrow=nMolecules)
  
  colnames(pvalsdf) <- c("Time","Group","Group_Time")
  fitsdf <- matrix(NA,ncol=3,nrow=nMolecules)
  
  colnames(fitsdf) <- c("Time","Group","Group_Time")
  clusterExport(cl, list('lme','sampleID','time','group','data','gsub','try','outer','pdDiag','melt','dcast','knots','basis','smspline','drop.levels','anova','type','pvals','keepModels','lmeControl','pdIdent','other.reshape.group'),envir=environment())
  
  #for(i in 1:nMolecules){
  new.data <- parLapply(cl,1:nMolecules,fun = function(i){
    
    if(type!="all")
      p2<-p1 <-p3 <- fitTime <-fitG <-fitTG <- NA
    expr <- data[,i]
    
    dataM <- other.reshape.group(Rep=sampleID,Time=time,Data=unlist(expr),Group=group)
    dataM$Group <- as.factor(dataM$Group)
    dataM$time = as.numeric(as.character(dataM$Time))
    dataM$Expr = as.numeric(as.character(dataM$Expr))
    dataM$all= rep(1,nrow(dataM))
    
    
    #### CUBIC SPLINE BASIS ####
    if(basis=="cubic"){
      dataM$Zt <- smspline(~ time, data=dataM)
      knots <- sort(unique(dataM$time))[-c(1,length(unique(dataM$time)))]
    }
    #### PENALIZED SPLINE BASIS#####
    if(basis%in%c("p-spline","cubic p-spline")){
      if(is.null(knots)){
        K <- max(5,min(floor(length(unique(dataM$time))/4),40))
        knots <- quantile(na.omit(unique(dataM$time)),seq(0,1,length=K+2))[-c(1,K+2)]
      }
      
      PZ <- outer(dataM$time,knots,"-")
      if(basis=="cubic p-spline"){
        PZ <- PZ^3
      }
      PZ <- PZ *(PZ>0)
      dataM$Zt <- PZ 
      
    }
    options(contrasts=c('contr.treatment','contr.poly'))
    tmp.data <- dataM
    p1 <- NA
    p2 <- NA
    p3 <- NA
    
    ###########group effect################# 
    
    if(type=="group"|type=="all"){
      fitG <- 1 
      fit0 <- NULL
      fit2 <- NULL
      nextF <- F
      p2TL1 <- 1
      
      if(experiment=="timecourse" | experiment=="all"){
        fit0T <-NULL
        fit0T <- try(lme(Expr ~1, data=tmp.data,
                         random=list(all=pdIdent(~Zt - 1)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        # print(fit0T)
        fit2T <- NULL
        fit2T<- try(lme(Expr ~as.factor(Group), data=tmp.data,
                        random=list(Group=pdIdent(~Zt - 1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        if(experiment=="timecourse" ){
          fit0 <- fit0T
          fit2 <- fit2T
          fitT <- 1
        }
        # print(fit2T)
      }
      
      
      if(experiment=='longitudinal1' | experiment=="all"){
        fit0L1 <- NULL
        fit0L1 <- try(lme(Expr ~1, data=tmp.data,
                          random=list(all=pdIdent(~Zt - 1),Rep=pdIdent(~1)),
                          na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit2L1 <-NULL
        fit2L1<- try(lme(Expr ~as.factor(Group), data=tmp.data,
                         random=list(Group=pdIdent(~Zt - 1), Rep=pdIdent(~1)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        if(experiment=="longitudinal1"){
          fit0 <- fit0L1
          fit2 <- fit2L1
          fitG <- 2
        }
        
      }
      
      if(experiment=="all"){
        if(class(fit2T) != 'try-error' & class(fit2L1) != 'try-error' & !is.null(fit2T) & !is.null(fit2L1)){
          
          p2TL1 <- anova(fit2T,fit2L1)$`p-value`[2][1]
          
          if(p2TL1<=0.05)
            nextF <- T
          
        }else{
          
          fit0 <- fit0T
          fit2 <- fit2T
          fitG <- 1
        }
        
        
        if(p2TL1>0.05){
          
          fit0 <- fit0T
          fit2 <- fit2T
          fitG <- 1
        }
        
      }
      if(experiment=="longitudinal2" | nextF ){
        fit0L2 <- try(lme(Expr ~1, data=tmp.data,
                          random=list(all=pdIdent(~Zt - 1),Rep=pdDiag(~time)),
                          na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        fit2L2<- try(lme(Expr ~as.factor(Group), data=tmp.data,
                         random=list(Group=pdIdent(~Zt - 1), Rep=pdDiag(~time)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        if(experiment=="longitudinal2"){
          fit0 <- fit0L2
          fit2 <- fit2L2
          fitG <- 3
        }
        
        if(nextF){
          
          if(class(fit2L2) != 'try-error' & class(fit2L1) != 'try-error' &!is.null(fit2L2) & !is.null(fit2L1)){
            
            p2L1L2 <- anova(fit2L1,fit2L2)$`p-value`[2][1]
            
            if(p2L1L2<0.05){
              fit0 <- fit0L2
              fit2 <- fit2L2
              fitG <- 3
            }else{
              fit0 <- fit0L1
              fit2 <- fit2L1
              fitG <- 2
            }
            
            
          }else{
            fit0 <- fit0L1
            fit2 <- fit2L1
            fitG <- 2
          }
        }
        
      }
      
      if(class(fit0) != 'try-error' & class(fit2) != 'try-error' & !is.null(fit2) & !is.null(fit0)){
        p1 <- anova(fit0,fit2)$`p-value`[2][1]        
      }
      if(keepModels){
        model.group <- fit2
      }
      
      tmp.time <- unique(tmp.data$time)
      tmp.group <- rep(unique(tmp.data$Group),length(tmp.time))
      tmp.time <- rep(tmp.time,each=2)
      pred.group <- rep(NA,length(tmp.time))
      if(class(fit2) != 'try-error'){
        tempdf <- data.frame(tmp.data$time,tmp.data$Group,fitted(fit2,level=1))
        tempdf <- tempdf[!duplicated(tempdf) & !is.na(tempdf[,3]),] 
  
        pred.group[paste(tmp.group,tmp.time)%in% paste(tempdf[,2],tempdf[,1])] <- tempdf[,3]
      }
    }
    
    ###########time effect#################
    if(type=="time"|type=="all"){
      fit0 <- NULL
      fit3 <- NULL
      fitTime <- 1
      nextF <- F
        
      if(experiment=="timecourse"|experiment=="all"){
        fit0T<- try(lme(Expr ~1, data=tmp.data,
                        random=list(all=pdIdent(~1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        fit3T<- try(lme(Expr ~time, data=tmp.data,
                        random=list(all=pdIdent(~Zt -1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        
        if(experiment=="timecourse"){
          fit0 <- fit0T
          fit3 <- fit3T
          fitTime <- 1
        }
      }
      
      
      if(experiment=="longitudinal1"|experiment=="all"){
        fit0L1 <-NULL
      fit0L1<- try(lme(Expr ~1, data=tmp.data,
                        random=list(all=pdIdent(~Zt - 1),Rep=pdIdent(~1)),
                        na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))

        
        fit3L1 <- NULL
        fit3L1<- try(lme(Expr ~time, data=tmp.data,
                         random=list(all=pdIdent(~Zt - 1),Rep=pdIdent(~1)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        if(experiment=="longitudinal1"){
          fit0 <- fit0L1
          fit3 <- fit3L1
          fitTime <- 2
        }
      }
      
      if(experiment=="all"){
        if(class(fit3L1) != 'try-error'& class(fit3T) != 'try-error' & !is.null(fit3L1) & !is.null(fit3T)){
          
          p1TL1 <- anova(fit3L1,fit3T)$`p-value`[2][1]
          
          if(p1TL1<=0.05){
            nextF <- T
          }else{
            fit0 <- fit0T
            fit3 <- fit3T
            fitTime <- 1
          }
          
        }else{
          fit0 <- fit0T
          fit3 <- fit3T
          fitTime <- 1
        }
        
      }
      if(experiment=="longitudinal2"|nextF){
        fit0L2 <- NULL 
        fit0L2<- try(lme(Expr ~1, data=tmp.data,
                         random=list(all=pdIdent(~1),Rep=pdDiag(~time)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit3L2 <- NULL
        fit3L2<- try(lme(Expr ~time, data=tmp.data,
                         random=list(all=pdIdent(~Zt - 1),Rep=pdDiag(~time)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        if(experiment=="longitudinal2"){
        fit0 <- fit0L2
        fit3 <- fit3L2
        fitTime <- 3
        }
     
        if(nextF){
        if(class(fit3L1) != 'try-error'& class(fit3L2) != 'try-error'& !is.null(fit3L1) & !is.null(fit3L2)){
          
          p3L1L2 <- anova(fit3L1,fit3L2)$`p-value`[2][1]
          
          if(p3L1L2<=0.05){
            fit0 <- fit0L2
            fit3 <- fit3L2
            fitTime <- 3
            
          }else{
            fit0 <- fit0L1
            fit3 <- fit3L1
            fitTime <- 2
          }
          
        }else{
          fit0 <- fit0L1
          fit3 <- fit3L1
          fitTime <- 2
        }
        
      } 
      }
      
      if(class(fit0) != 'try-error'& class(fit3) != 'try-error' & !is.null(fit0) & !is.null(fit3)){
        p2 <- anova(fit0,fit3)$`p-value`[2][1]
      }
      if(keepModels){
        model.time <- fit3
      }
      pred.time <- rep(NA,length(na.omit(unique(tmp.data$time))))
      if(class(fit3) != 'try-error')
        pred.time <- na.omit(unique(fitted(fit3,level=1)))
    }
    
    
    ###########group*time interaction#################
    if(type=="grouptime"|type=="all"){
      
      fitTG <- 0
      fit0 <- NULL
      fit5 <- NULL
      fit0T <- NULL
      fit5T <- NULL
      nextF <- F
      p3TL1 <- 1
      if(experiment=="timecourse"|  experiment=="all"){
        fit0T <- try(lme(Expr ~ as.factor(Group)+time,data=tmp.data,
                         random=list(all=pdIdent(~Zt - 1)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit5T <- try(lme(Expr ~ as.factor(Group)*time,data=tmp.data,
                         random=list(Group=pdIdent(~Zt - 1)),
                         na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        if(experiment=="timecourse"){
          fit0 <- fit0T
          fit5 <- fit5T
          fitTG <- 1
        }
        
      }  
      
      if(experiment=="longitudinal1"| experiment=="all") {
        fit0L1 <- try(lme(Expr ~ as.factor(Group)+time, data=tmp.data,
                          random=list(all=pdIdent(~Zt - 1), Rep=pdIdent(~1)),
                          na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit5L1 <- try(lme(Expr ~ as.factor(Group)*time,data=tmp.data,
                          random=list(Group=pdIdent(~Zt - 1), Rep=pdIdent(~1)),
                          na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
       if(experiment=="longitudinal1"){
        fit0 <- fit0L1
        fit5 <- fit5L1
        fitTG <- 2
       }
      }
      
      
      if(experiment=="all"){
        if(class(fit5L1) != 'try-error'& class(fit5T) != 'try-error' & !is.null(fit5L1) & !is.null(fit5T)){
          p3TL1 <- anova(fit5L1,fit5T)$`p-value`[2][1]
          
          if(p3TL1<=0.05){
            nextF <- T
          }else{
            fit0 <- fit0T
            fit5 <- fit5T
            fitTG <-1
          }
          
        }else{
          fit0 <- fit0T
          fit5 <- fit5T
          fitTG <- 1
        }
      }
      
      
      if(experiment=="longitudinal2"|nextF) {
        fit0L2 <- try(lme(Expr ~ as.factor(Group)+time, data=tmp.data,
                          random=list(all=pdIdent(~Zt - 1), Rep=pdDiag(~time)),
                          na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        fit5L2 <- try(lme(Expr ~ as.factor(Group)*time,data=tmp.data,
                          random=list(Group=pdIdent(~Zt - 1), Rep=pdDiag(~time)),
                          na.action=na.exclude, control=lmeControl(opt = "optim"),method="ML"))
        
        if(experiment=="longitudinal2"){
          fit0 <- fit0L2
          fit5 <- fit5L2
          fitTG <- 3
        }
        
        if(nextF){
        if(class(fit5L1) != 'try-error'& class(fit5L2) != 'try-error'& !is.null(fit5L1) & !is.null(fit5L2)){
          
          p3L1L2 <- anova(fit5L1,fit5L2)$`p-value`[2][1]
          
          if(p3L1L2<=0.05){
            fit0 <- fit0L2
            fit5 <- fit5L2
            fitTG <- 3
            
          }else{
            fit0 <- fit0L1
            fit5 <- fit5L1
            fitTG <- 2
          }
          
        }else{
          fit0 <- fit0L1
          fit5 <- fit5L1
          fitTG <- 2
        } 
      }
      }
      if(class(fit0) != 'try-error'& class(fit5) != 'try-error' & !is.null(fit0) & !is.null(fit5)){
        p3 <- anova(fit0,fit5)$`p-value`[2][1]
      }
      if(keepModels){
        model.time.group <- fit5
      }
      tmp.time <- sort(na.omit(unique(tmp.data$time)))
      tmp.group <- rep(unique(na.omit(tmp.data$Group)),length(tmp.time))
      tmp.time <- rep(tmp.time,each=2)
      pred.time.group <- rep(NA,length(tmp.time))
      
      if(class(fit5) != 'try-error'){
        
        tempdf <- data.frame(tmp.data$time,tmp.data$Group,fitted(fit5,level=1))
        tempdf <- tempdf[!duplicated(tempdf) & !is.na(tempdf[,3]),]   
        pred.time.group[paste(tmp.group,tmp.time)%in% paste(tempdf[,2],tempdf[,1])] <- tempdf[,3]
      }
    }
    
    # pvalsdf[i,] <- c(p2,p1,p3)
    # fitsdf[i,] <- c(fitTime,fitG,fitTG)
    
    pvals<-c(p2,p1,p3)
    
    
    
    if(keepModels){
      return(list(pvals=pvals,model.time=model.time,model.group=model.group,model.time.group=model.time.group,knots=knots,pred.time=pred.time,pred.time.group=pred.time.group,pred.group=pred.group,fits=c(fitTime,fitG,fitTG)))
    }else{
      return(list(pvals=pvals,knots=knots,pred.time=pred.time,pred.time.group=pred.time.group,pred.group=pred.group,fits=c(fitTime,fitG,fitTG)))
    }
    
    #     if(keepModels){
    #       return(list(pvals=pvals,model.time=model.time,model.group=model.group,model.time.group=model.time.group,knots=knots))
    #     }else{
    #       return(list(pvals=pvals,knots=knots,pred.time=pred.time,pred.time.group=pred.time.group,pred.group=pred.group))
    #     }
    
  })
  
  stopCluster(cl)
  new.df <- matrix(sapply(new.data,'[[','pvals'),nrow=nMolecules,ncol=3,byrow=T)
  knots <- unique(as.vector((sapply(new.data,'[[','knots'))))
  fits <- matrix(sapply(new.data,'[[','fits'),nrow=nMolecules,ncol=3,byrow=T)
  model.time <- list()
  model.time <- list()
  model.group <- list()
  model.time.group <- list()
  pred.time <- matrix()
  pred.group <- matrix()
  pred.group.time <- matrix()
  
  if(keepModels){
    model.time <- sapply(new.data,'[','model.time')
    model.group <- sapply(new.data,'[','model.group')
    model.time.group <- sapply(new.data,'[','model.time.group')
  }
    pred.time <- matrix(sapply(new.data,'[[','pred.time'),nrow=nMolecules,ncol=length(unique(time)),byrow=T)
    nameTime <-  na.omit(sort(unique(time)))
    if(ncol(pred.time)==length(nameTime))
      colnames(pred.time) <- nameTime
    pred.group <- matrix(sapply(new.data,'[[','pred.group'),nrow=nMolecules,ncol=(length(unique(time))*2),byrow=T)
    
    pred.group.time <- matrix(sapply(new.data,'[[','pred.time.group'),nrow=nMolecules,ncol=(length(unique(time))*2),byrow=T)
    nameGroup <- paste(na.omit(sort(unique(group))),rep(na.omit(sort(unique(time))),each=2))
    if(ncol(pred.group)==length(nameGroup) & ncol(pred.group.time)==length(nameGroup))
      colnames(pred.group) <- colnames(pred.group.time) <- nameGroup
    molnames <-  colnames(data)    

  if(is.null(molnames))
    molnames <- 1:ncol(data)
  
  if(type=="all"){
    
    df <- data.frame(Molecule=molnames,Time=new.df[,1],adj.Time=signif(p.adjust(new.df[,1],method="BH"),2), Group=new.df[,2],adj.Group=signif(p.adjust(new.df[,2],method="BH"),2),Group_Time=new.df[,3],adj.Group_Time=signif(p.adjust(new.df[,3],method="BH"),2))
  }
  if(type=="time"){
    df <- data.frame(Molecule=molnames,Time=new.df[,1],adj.Time=signif(p.adjust(new.df[,1],method="BH"),2))
  }
  if(type=="group"){
    df <- data.frame(Molecule=molnames, Group=new.df[,2],adj.Group=signif(p.adjust(new.df[,2],method="BH"),2))
  }
  if(type=="grouptime"){
    df <- data.frame(Molecule=molnames,Group_Time=new.df[,3],adj.Group_Time=signif(p.adjust(new.df[,3],method="BH"),2))
    
  }
    
  l <- new('lmmsde',DE=df, modelsUsed=fits,modelTime=model.time, modelGroup=model.group, modelTimeGroup=model.time.group,knots=knots,basis=basis,type=type,experiment=experiment,predTime=pred.time,predGroup=pred.group,predTimeGroup=pred.group.time)
  
      
  return(l)
  
}

other.reshape.group <- function(Rep, Time, Data,Group){
  lme.data<-NULL
  #require(reshape)
  if(sum(table(Rep,Time)>1)!=0)
    stop('Make sure you have one time point per sample ID.')
  lme.data <- data.frame(Time=Time,Rep=Rep,Group=Group,Data)
  
  lme.data$Time = factor(drop.levels(lme.data$Time))
  lme.data$Rep = factor(drop.levels(lme.data$Rep))
  lme.data$Group = factor(drop.levels(lme.data$Group))
  
  melt.lme.data <-NULL
 # melt <- reshape2::melt
  melt.lme.data <- melt(lme.data)
  #melt.lme.data <- melt(lme.data)
  cast.lme.data  <- NULL
  #cast.lme.data <- cast(melt.lme.data, variable+ Group+Rep ~ Time)
  cast.lme.data <- dcast(melt.lme.data, variable+ Group+Rep ~ Time)
  melt.lme.data2 <- NULL
 # melt.lme.data2 <-  melt(data.frame(cast.lme.data))
  melt.lme.data2 <-  melt(data.frame(cast.lme.data))
  
  names(melt.lme.data2) <- c("Molecule",  "Group","Rep", "Time", "Expr")
  melt.lme.data2$Time <- factor(gsub("^X", "", as.character(melt.lme.data2$Time)))
  return(as.data.frame(melt.lme.data2))
}
