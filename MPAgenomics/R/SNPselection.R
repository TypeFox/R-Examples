#'
#' This function selects, for each chromosome, the most relevant markers according to a response.
#'
#' @title markers selection
#' 
#' @param dataSetName The name of the data-set folder.
#' @param dataResponse A csv files or a data.frame with 2 columns : "files" and "response". The column "files" contains the filename to extract and the second column the response associated to the file.
#' @param chromosome A vector containing the number of the chromosomes for the SNPs selection.
#' @param signal either "CN" or "fracB". corresponding to which signal will be analyzed (default="CN").
#' @param normalTumorArray Only in the case of normal-tumor study. A csv file or a data.frame containing the mapping between normal and tumor files.
#' The first column contains the name of normal files and the second the names of associated tumor files.
#' @param onlySNP (only if signal="CN"). If TRUE, only the SNPs probes are used (default=FALSE).
#' @param nbFolds number of folds in the cross validation (default=10).
#' @param loss either "logistic" (binary response) or "linear" (quantitative response), default is "logistic"
#' @param plot If TRUE, cross-validation mean squared error is plotted (default=TRUE).
#' @param pkg Either "HDPenReg" or "spikeslab". Ued package in linear case.
#' @param ... Other parameters for HDlars, glmnet or spikeslab function.
#' 
#' @return a list containing length(chromosme) elements. Each element is a list containing
#' \describe{
#'   \item{chr}{chromosome corresponding to the signal.}
#'   \item{markers.index}{A vector containing the index of all selected markers.}
#'   \item{markers.position}{A vector containing the position of all selected markers.}
#'   \item{markers.names}{A vector containing the names of all selected markers.}
#'   \item{coefficient}{A vector containing the coefficients of all selected markers.}
#'   \item{intercept}{Intercept of the model.}
#   \item{fraction}{fraction of the l1 norm of coefficient selected by cross validation.}
#' }
#'
#' @details This function requires to use the aroma folder architecture. In your working directory, there must have the rawData folder and totalAndFracBData folder.
#' This function launches the lars algorithm on the CN or fracB data and uses a cross-validation to select the most appropriate solution.
#' 
#' 
#' @seealso HDPenReg, glmnet, spikeslab
#'
#' @author Quentin Grimonprez
#'
#' @export
#'
markerSelection=function(dataSetName,dataResponse,chromosome=1:22,signal=c("CN","fracB"),normalTumorArray,onlySNP=FALSE,nbFolds=10,loss=c("logistic","linear"),plot=TRUE,pkg=c("HDPenReg","spikeslab"),...)
{
  loss <- match.arg(loss)
  signal <- match.arg(signal)
  pkg <- match.arg(pkg)
  
  allpkg=TRUE
  if(!suppressPackageStartupMessages(require("aroma.affymetrix", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.affymetrix. For download it:\n")
    cat("source(\"http://www.braju.com/R/hbLite.R\")\n")
    cat(" hbLite(\"sfit\")\n")
    cat("source(\"http://bioconductor.org/biocLite.R\")\n")
    cat("biocLite(\"affxparser\")\n")
    cat("biocLite(\"DNAcopy\")\n")
    cat("biocLite(\"aroma.light\")\n")
    #     cat("source(\"http://aroma-project.org/hbLite.R\")\n")
    cat("install.packages(\"aroma.affymetrix\")\n")
    allpkg=FALSE
  }
  
  
  if(!suppressPackageStartupMessages(require("aroma.cn", quietly=TRUE) ) )
  {
    cat("Package not found: aroma.cn. For download it:\n")
    cat("install.packages(\"aroma.cn\")\n") 
    allpkg=FALSE
  }
  
  if(!allpkg)
    stop("You have to install some packages : Follow the printed informations.")
  
  ################## check parameters
  #check if the user is in the right place
  if(!("totalAndFracBData"%in%list.files()))
    stop("There is no \"totalAndFracBData\", check if you are in the good working directory or if you have run the signalPreProcess function before.")
  
  
  #dataSetName
  if(!is.character(dataSetName))
    stop("dataSetName must be the name of a folder in rawData.")
  if(length(grep(dataSetName,list.files("totalAndFracBData")))==0) #check if data are available for this dataSetName
    stop("The dataSetName you specify doesn't exist in the totalAndFracBData, run the ... function at first.")
  
  #dataResponse
  if(is.character(dataResponse))
    dataResponse=read.csv(dataResponse)
  else
  {
    if(!is.data.frame(dataResponse))
      stop("dataResponse must be either the path to the normalTumorArray csv file or a data.frame containing the response.") 
  }
  if(!("files"%in%names(dataResponse)))
    stop("dataResponse does not contain the column \"files\".")
  if(!("response"%in%names(dataResponse)))
    stop("dataResponse does not contain the column \"names\".")
  
  #loss
  #   if(!(loss%in%c("logistic","linear")))
  #     stop("loss must be either \"logistic\" or \"linear\".")
  
  #chromosome : vector of integer between 1 and 25
  if( !is.numeric(chromosome) || !is.vector(chromosome) )
    stop("chromosome must be a vector containing the different chromosomes to analyze.")
  chromsome=unique(chromosome)  #keep only unique values
  if( sum(sapply(chromosome,FUN=function(chr){!is.wholenumber(chr)}))>0 )
    stop("chromosome must be a vector containing integers corresponding to the different chromosomes to analyze.") 
  if(min(chromosome)<1)
    stop("chromosome must contain integer between 1 and 25.")
  if(max(chromosome)>25)
    stop("chromosome must contain integer between 1 and 25.")
  
  ########################
  
  #launch the right function depending of the signal
  res=switch(signal,
             "CN"=SNPselectionCNsignal(dataSetName,dataResponse,chromosome,normalTumorArray,onlySNP,nbFolds,loss,plot,pkg,...),
             "fracB"=SNPselectionFracBsignal(dataSetName,dataResponse,chromosome,normalTumorArray,nbFolds,loss,plot,pkg,...),
             "both"=stop("Not yet implemented."))
  
  return(res)
}


#
# This function selects, for each chromosome, the most relevant SNPs.
#
# @title SNPs selection from CN signal
# 
# @param dataSetName Name of the dataset folder in rawData 
# @param dataResponse response associated to the data
# @param chromosome chromosome used in the study
# @param normalTumorArray only if you have normal and tumor profile in your data folder. A csv file or a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param onlySNP (only if signal=\"CN\"). If TRUE, only the SNP markers will be used.
# @param nbFolds number of folds in the cross validation
# @param loss either \"logistic\" (binary response) or \"linear\" (quantitative response).
# @param plot if TRUE, plot the cross validation graphic
# @param pkg Either "HDPenReg" or "spikeslab". Ued package in linear case.
#
# @return a list containing length(chromosme) elements. Each element is a list containing
# \describe{
#   \item{chr}{chromosome corresponding to the signal.}
#   \item{variable}{A vector containing the index of all selected SNPs.}
#   \item{coefficient}{A vector containing the coefficeints of all selected SNPs.}
# }
#
# @details This function requires to use the aroma folder architecture. In the working directory, there must be the rawData folder and 
# totalAndFracBData folder.
#
# @author Quentin Grimonprez
#
SNPselectionCNsignal=function(dataSetName,dataResponse,chromosome,normalTumorArray,onlySNP,nbFolds=10,loss="logistic",plot=TRUE,pkg="HDPenReg",...)
{
  res=list()
  for(chr in chromosome)
  {
    #get the matrix with the copy number signal for the chromosome chr
    C=getCopyNumberSignal(dataSetName,chr,normalTumorArray,onlySNP,listOfFiles=as.character(dataResponse$files),verbose=FALSE)
    C=C[[paste0("chr",chr)]]
    gc()
    #extract the response in the right order
    ind=sapply(names(C)[3:(ncol(C)-1)],FUN=function(x){match(x,dataResponse$files)})    
    if(sum(is.na(ind))!=0)
      stop(paste0("A response is missing for the following files : ",paste(C$sampleNames[is.na(ind)],collapse=", ")))      
    
    response=dataResponse$response[ind]
    
    nbFolds=min(nbFolds,length(response))
    
    restemp=variableSelection(t(as.matrix(C[,3:(ncol(C)-1)])),response,nbFolds,loss,plot,pkg,...)
    
    res[[paste0("chr",chr)]]=list(chr=chr,markers.index=restemp$markers.index,markers.position=C$position[restemp$markers.index],
                                  markers.names=as.character(C$featureNames)[restemp$markers.index],coefficient=restemp$coefficient,
                                  intercept=restemp$intercept)
    
    #delete objects created during this loop
    rm(C,restemp)    
    gc()
  }
  
  return(res)
}

#
# This function selects, for each chromosome, the most relevant SNPs.
#
# @title SNPs selection from fracB signal
# 
# @param dataSetName Name of the dataset folder in rawData 
# @param dataResponse response associated to the data
# @param chromosome chromosome used in the study
# @param normalTumorArray only if you have normal and tumor profile in your data folder. A csv file or a data.frame with 2 columns: "normal" and "tumor".
# The first column contains the name of normal files and the second the names of associated tumor files.
# @param onlySNP (only if signal=\"CN\"). If TRUE, only the SNP markers will be used.
# @param nFolds number of folds in the cross validation
# @param loss either \"logistic\" (binary response) or \"linear\" (quantitative response).
# @param plot if TRUE, plot the cross validation graphic
# @param pkg Either "HDPenReg" or "spikeslab". Ued package in linear case.
# 
# @return a list containing length(chromosme) elements. Each element is a list containing
# \describe{
#   \item{chr}{chromosome corresponding to the signal.}
#   \item{variable}{A vector containing the index of all selected SNPs.}
#   \item{coefficient}{A vector containing the coefficeints of all selected SNPs.}
# }
#
# @details This function requires to use the aroma folder architecture. In the working directory, there must be the rawData folder and totalAndFracBData folder.
#
# @author Quentin Grimonprez
#
SNPselectionFracBsignal=function(dataSetName,dataResponse,chromosome,normalTumorArray,nbFolds=10,loss="logistic",plot=TRUE,pkg="HDPenReg",...)
{
  res=list()
  for(chr in chromosome)
  {
    #besoin que de la tumeur, rajouter onlyTumor comme parametre dans getFracB
    #ajouter listOfFile quand on veux ou on ne peux pas travailler sur ttes les donnees?
    
    #get the fracB signal for the chromosome chr
    
    fracB=getFracBSignal(dataSetName,chr,normalTumorArray,listOfFiles=as.character(dataResponse$files),verbose=FALSE)
    
    
    fracB=fracB[[paste0("chr",chr)]]$tumor
    gc()
    
    #extract the response in the right order
    ind=sapply(names(fracB)[3:(ncol(fracB)-1)],FUN=function(x){match(x,dataResponse$files)})    
    if(sum(is.na(ind))!=0)
      stop(paste0("A response is missing for the following files : ",paste(fracB$sampleNames[is.na(ind)],collapse=", ")))
    
    response=dataResponse$response[ind]
    
    nbFolds=min(nbFolds,length(response))
    
    restemp=variableSelection(t(as.matrix(fracB[,3:(ncol(fracB)-1)])),response,nbFolds,loss,plot,pkg,...)
    
    res[[paste0("chr",chr)]]=list(chr=chr,markers.index=restemp$markers.index,markers.position=fracB$position[restemp$markers.index],
                                  markers.names=as.character(fracB$featureNames)[restemp$markers.index],coefficient=restemp$coefficient,
                                  intercept=restemp$intercept)
    
    rm(fracB,restemp)
    gc()
  }
  
  return(res)
}


#'
#' This function selects the most relevant variables according to a response.
#'
#' @title SNPs selection
#' 
#' @param dataMatrix Matrix containing the data, each row is a different sample.
#' @param dataResponse response associated to the data.
#' @param nbFolds number of folds in the cross validation.
#' @param loss either "logistic" (binary response) or "linear" (quantitative response).
#' @param plot If TRUE plot cross-validation mean squared error (default=TRUE).
#' @param pkg Either "HDPenReg" or "spikeslab". Ued package in linear case.
#' @param ... spplementary arguments for cv.glmnet function in case of logistic loss or for HDlars or spikeslab function for linear loss.
#' 
#' @return a list containing 
#' \describe{
#'   \item{variable}{A vector containing the index of all selected variables.}
#'   \item{coefficient}{A vector containing the coefficients of all selected variables.}
#'   \item{intercept}{Intercept of the model.}
#'}
#'
#' @author Quentin Grimonprez
#' 
#' @export
variableSelection=function(dataMatrix,dataResponse,nbFolds=min(length(dataResponse),10),loss=c("logistic","linear"),plot=TRUE,pkg=c("HDPenReg","spikeslab"),...)
{
  loss <- match.arg(loss)
  pkg <- match.arg(pkg)
  
  #check plot (other parameters will be checcked in HDcvlars function)
  if(!is.logical(plot))
    stop("plot must be a logical.")
  
  if(loss=="linear")
  {
    if(pkg=="HDPenReg")
    {
      #lars algorithm for obtaining all the path
      reslars=HDlars(dataMatrix, dataResponse,...)
      
      #cross validation to choose the best lambda
      rescv=HDcvlars(dataMatrix, dataResponse, nbFolds,index = c(reslars@lambda,0), mode="lambda",...)
      
      if(plot)
      {
        plot(rescv)
      }
      
      
      
      #we compute the coefficients for the value given by the HDcvlars function
      #coeff=computeCoefficients(reslars,rescv$minIndex,mode="lambda")
      
      indKnee=Lmethod(rescv$cv)
      #var=reslars@variable[[which.min(rescv$cv)]]
      #coef=reslars@coefficient[[which.min(rescv$cv)]]
      var=reslars@variable[[indKnee]]
      coef=reslars@coefficient[[indKnee]]
      #var=coeff$variable
      #coef=coeff$coefficient
      intercept=reslars@mu
      if(length(var)!=0)
      {
        index=order(var)
        var=var[index]
        coef=coef[index]
        #index=order(coeff$variable)
        #var=coeff$variable[index]
        #coef=coeff$coefficient[index]
      }
      
      rm(reslars)
      #rm(reslars,coeff)
      gc()
    }
    else
    {
      if(pkg=="spikeslab")
      {
        ################ spikeslab
        #rescv=cv.spikeslab(x=dataMatrix, y=dataResponse, bigp.smalln = TRUE,K=nbFolds,... )
        rescv=spikeslab(x=dataMatrix, y=dataResponse, bigp.smalln = TRUE,... )
        
        intercept=rescv$y.center #rescv$spikeslab.obj$y.center #si cv.spikeslab
        var=which(rescv$gnet.scale!=0)
        coef=rescv$gnet.scale[var]
      }
    }
  }
  else
  {
    if(loss=="logistic")
    {
      rescv=cv.glmnet(dataMatrix, dataResponse, nfolds=nbFolds,family="binomial",...)
      if(plot)
      {
        plot(rescv)
      }
      coef=coef(rescv,s=rescv$lambda.min)
      
      ind=which(coef!=0)
      
      var=ind[-1]-1
      intercept=coef[1]
      coef=coef[ind[-1]]
      
      names(var)=NULL
      names(coef)=NULL
    }
  }
  
  
  res=list(markers.index=var,coefficient=coef,intercept=intercept)  
}


#Determining the Number of Clusters/Segments in Hierarchical Clustering/Segmentation Algorithms
#Stan Salvador and Philip Chan 
#
# This function searches the furthest point 
#
# @title Find a knee in a curve
#
# @param y ordiantes of the curve
# @param x abscissas of the curve
# @return index of the knee
#
knee=function(y,x=1:length(y))
{
  m=(y[length(y)]-y[1])/(x[length(x)]-x[1])
  p=y[1]-m*x[1]
  d=abs(m*x-y+p)/sqrt(2+m^2)
  knee=which.max(d)
  return(knee)
}


Lmethod=function(y,x=1:length(y))
{
  b=length(x)
  lrmse=rrmse=rmse=rep(NA,b)
  for(i in 2:(b-2))
  {
    l=lm(y[1:i]~x[1:i])
    r=lm(y[(i+1):b]~x[(i+1):b])
    lrmse[i]=sqrt(sum(l$residuals^2))
    rrmse[i]=sqrt(sum(r$residuals^2))
    
    rmse[i]=(i-1)/(b-1)*lrmse[i]+(b-i)/(b-1)*rrmse[i]
  }
  
  return(min(which.min(rmse)+1,b))
}