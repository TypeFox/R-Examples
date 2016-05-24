#' Simple descriptive statistics
#'
#' R Implementation of the SPSS Function \code{Descriptives}
#'
#' The xpssDescriptives function provides a set of descriptive statistic tools. 
#' 
#' \strong{\code{missing:}} 
#' \tabular{rll}{
#' 
#'\tab \code{variable} \tab removes user-, and system-missing data explicitly for every variable. 
#'\cr \tab \code{listwise} \tab performs a listwise-deletion.
#'\cr \tab \code{include} \tab includes all user-defined missing values.} 
#' 
#' \strong{\code{statistics:}}
#' \tabular{rll}{
#' 
#'\tab \code{kurtosis} \tab calculates the bulge of the variable. 
#'\cr \tab \code{max} \tab displays the maximum of the variable. 
#'\cr \tab \code{mean} \tab calculates the arithmetic mean, respectively the midpoint of the variable.
#'\cr \tab \code{min} \tab displays the minimum of the variable. 
#'\cr \tab \code{kurtosis} \tab calculates the bulge of the variable. 
#'\cr \tab \code{range} \tab displays the span between the minimum and the maximum value. 
#'\cr \tab \code{sekurtosis} \tab calculates the standrard error of the bulge of the variable. 
#'\cr \tab \code{semean} \tab displays the standard error of the arithmetic mean. 
#'\cr \tab \code{seskewness} \tab calculates the standrard error of the inclination of the variable. 
#'\cr \tab \code{skewness} \tab calculates the inclination of the variable. 
#'\cr \tab \code{stddev} \tab  displays the standard deviation of the variable. 
#'\cr \tab \code{sum} \tab calculates the sum of each observation within the variable. 
#'\cr \tab \code{variance} \tab displays the variance.}
#'
#'  \code{ztrans} input, is a list with elements varname and zname.  \code{varname} and  \code{zname} are either atomic characters or character vectors. \cr It is necessary that either both parameters are filled or blank.
#'  
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param variables atomic character or character vector with the name of the variables.
#' @param missing atomic character which specifiy the missing method. The method indicates what should happen when the data contains NAs. Default is \code{"variable"}.
#' @param statistics atomic chracter or character vector which determine the descriptiv statistics. Default are \code{"mean"}, \code{"max"}, \code{"min"}, \code{"stddev"}.
#' @param save logical indicator. TRUE adds the z-score of each variable to \code{x}. Default is FALSE.
#' @param ztrans list which specifies variables for z-transformation and name of z-transformed variables. Read Details for further information.
#' 
#' @return Output is a list object with descriptive statistic parameters. The specific outcomes of the selected variables are stored in a list object. Every variable is stored in a different list element.
#' 
#' If the parameter \code{save} is TRUE, a matrix with z-transformed values will be appended at the end of the list. If \code{ztrans} is blank, the name of the matrix will be Z*varname*. Otherwise whether \code{ztrans} is not empty the user specified description in \code{zname} will be the name of the z-transformed matrix of the variable \code{varname}.
#' 
#' @author Bastian Wiessner
#' @importFrom e1071 kurtosis skewness
#' @importFrom data.table as.data.table
#' @examples
#'
#' data(fromXPSS)
#' 
#' ## Analyzing Variable V5, Output contains default statistics
#' xpssDescriptives(x=fromXPSS,
#'                  variables="V5")
#'
#' ## Analyzing Variable V7_1, Output contains default statistics
#' ## and z-score of Variable V7_1
#' xpssDescriptives(x=fromXPSS,
#'                  variables="V7_1",
#'                  save = TRUE)
#'
#' ## Analyzing Variable V7_2, Output contains default statistics
#' ## and z-score of Variable V7_2 store in myZname
#' xpssDescriptives(x=fromXPSS,
#'                  variables="V7_2",
#'                  save = TRUE,
#'                  ztrans = list(varname = "V7_2", 
#'                                zname = "myZname"))
#'
#'  ## Analyzing Variable V7_2, Output contains kurtosis, skewness, semean and mean
#'  ## missing values are included
#'  ## z-score get calculated and store in myZname
#'  
#'   xpssDescriptives(x=fromXPSS,
#'                  variables="V7_2", 
#'                  statistics=c("kurtosis",
#'                                "skewness",
#'                                "semean",
#'                                "mean"),
#'                  missing="include",
#'                  save = TRUE,
#'                  ztrans = list(varname = "V7_2", 
#'                                zname = "myZname"))
#'                                

#' @export
xpssDescriptives <-  function(x, 
                              variables,
                              missing = "variable", 
                              statistics = c("mean",
                                             "max",
                                             "min",
                                             "stddev"),
                              save = FALSE,
                              ztrans = list(varname=NULL,
                                            zname=NULL))
{
  
  names(ztrans) <- c("varname", "zname")
  
  ####################################################################
  ####################################################################
  
  functiontype <- "AN"
  dataname <- eval(paste0(deparse(substitute(x))), envir = .GlobalEnv)
  x <- applyMetaCheck(x)
  options(warn=-1)
  ####################################################################
  ####################################################################
  ####################################################################
  for(i in 1:length(variables)){
    if(!(is.element(variables[i],names(x)))) {
      stop("The selected variables are not in the dataset")
    }
  }
  
  
  for(i in 1:length(variables))
  {
    if(class(x[,variables[i]]) != "numeric"){  
      stop("Variables are not numeric")
    }
  } 
  
  
  if(missing != "variable" && missing != "listwise" && missing != "include")  {
    stop("wrong 'missing' argument. Only the arguments 'variable', 'listwise', and 'include' are valid.")
  }
  if((statistics != "all") && (statistics != "default")&& (statistics != "n") && (statistics != "mean") && (statistics != "min") && (statistics != "max") && 
       (statistics != "stddev") && (statistics !="kurtosis") && (statistics != "range") && (statistics !="semean") && 
       (statistics != "skewness") && (statistics != "sum") && (statistics != "variance"))
  {
    stop("unknown statics command, only the following paramters are valid: 'all', 'default', 'kurtosis', 'mean', 'max', 'min', 'range', 'semean', 'skewness', 'stddev', 'sum', 'variance'")
  }
  
  if(!is.null(ztrans$zname) && is.null(ztrans$varname) || is.null(ztrans$zname) && !is.null(ztrans$varname))
  {
    stop("varname should contain the original variable name, zname the new name of the z-transformed variable")
  }
  
  if(!is.logical(save))
  {
    stop("the save argument has to be logical. f.e. TRUE or FALSE.")
  }
  
  #################################
  if(is.xpssFrame(x)){
    if(attributes(x)$SPLIT_FILE != FALSE){
      splitter <- unlist(str_split(attributes(x)$SPLIT_FILE,pattern=" "))
      splitter <- splitter[1:length(splitter)-1]
      variables <- c(variables,splitter)
    } 
  }
  descr <- list(variables) 
  #------------------ Missing FUnctions----------------------------------------------#
  if("variable" %in% missing)
  {
    for (i in 1:length(variables))
    {      
      pos <- which(x[,variables[i]] %in% attributes(x[,variables[i]])$defined.MIS)
      x[,variables[i]][pos] <- NA
      descr[[i]] <- na.omit(x[,variables[i]])
    }
  }
  if("listwise" %in% missing)
  {
    if(length(variables) >1){
      descr <- as.list(na.omit(x[,variables]))  
    } else {
      descr <- list(na.omit(x[,variables]))
    } 
  } 
  if("include" %in% missing)
  {
    variables <- c("V5","V6","V7_2")
    temp <- computeValue(x,variables)
    for(i in 1:length(variables)){
      descr[[i]] <- na.omit(temp[[i]])
    }
    
  }  
  names(descr) <- variables
  
  #----------------- Z - Transform -----------------------------------#
  if(save == T && is.null(ztrans$varname)) {
    for(i in 1:length(variables))
    {
      descr[[length(descr)+1]] <- scale(x[,variables[i]]) 
      names(descr)[length(descr)] <- paste("Z",variables[i], sep = "")  
    }    
  } 
  if(!is.null(ztrans$varname)) {
    if(is.element(ztrans$varname,variables)) {
      if((save==T) && (length(ztrans$varname) == (length(ztrans$zname))) && ((!is.null(ztrans$varname)) && (!is.null(ztrans$zname)))) {
        for(j in 1:length((intersect(ztrans$varname, variables)))){
          varnames <- names(descr)[which(names(descr)%in% ztrans$varname)]
          varnames <- paste("Z",varnames, sep = "")
          names(descr)[which(names(descr) %in% varnames[j])] <- ztrans$zname[j]           
        }
      } else if((save==F) && (length(ztrans$varname) == (length(ztrans$zname))) && ((!is.null(ztrans$varname)) && (!is.null(ztrans$zname)))){
        for(j in 1:length((intersect(ztrans$varname, variables)))){
          descr[[length(descr)+1]] <- scale(x[,ztrans$varname[j]])
          names(descr)[length(descr)] <- ztrans$zname[[j]]
        }
      }
    }
  }
  
  
  #------------------ Statistic FUnctions----------------------------------------------#
  tempmean <- NULL
  tempmin <- NULL
  tempmax <- NULL
  tempstddev  <- NULL
  tempkurtosis <- NULL
  tempsekurtosis <- NULL
  temprange <- NULL
  tempsemean <- NULL
  tempskewness <- NULL
  tempseskewness <- NULL
  tempsum <- NULL
  tempvariance <- NULL
  
  if(any(attributes(x)$SPLIT_FILE != FALSE)){
    descr <- as.data.table(descr)
    variables <- variables[length(variables)-(length(splitter))]
    splitter <- paste(splitter,collapse=",")
  }  
  if(any(attributes(x)$SPLIT_FILE != FALSE)){
    for(i in 1:length(variables))
    {
      express <- "list("
      if(is.element("n",statistics)){
        express <- paste0(express,"n=length(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("kurtosis",statistics)){
        express <- paste0(express,"kurt=kurtosis(get(variables[[i]]),type=2),")
      }
      if(is.element("max",statistics)){
        express <- paste0(express,"max=max(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("mean",statistics)){
        express <- paste0(express,"mean=mean(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("median",statistics)){
        express <- paste0(express,"median=median(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("min",statistics)){
        express <- paste0(express,"min=min(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("range",statistics)){
        express <- paste0(express,"range=diff(range(get(variables[[i]]),na.rm=T)),")
      }
      if(is.element("sekurt",statistics)){
        express <- paste0(express,"sekurt=2*(sqrt((6*length(na.omit(get(variables[[i]])))*(length(na.omit(get(variables[[i]])))-1)) / ((length(na.omit(get(variables[[i]])))-2)*(length(na.omit(get(variables[[i]])))+1)*(length(na.omit(get(variables[[i]])))+3))))* (sqrt(((length(na.omit(get(variables[[i]])))^2 -1)) / ((length(na.omit(get(variables[[i]])))-3)*(length(na.omit(get(variables[[i]])))+5)))),")
      }
      if(is.element("semean",statistics)){
        express <- paste0(express,"semean=sd(get(variables[[i]]),na.rm=T)/(sqrt(length(na.omit(get(variables[[i]]))))),")
      }
      if(is.element("skewness",statistics)){
        express <- paste0(express,"skew=skewness(get(variables[[i]]),type=2),")
      }
      if(is.element("seskewness",statistics)){
        express <- paste0(express,"seskew=sqrt((6*length(na.omit(get(variables[[i]])))*(length(na.omit(get(variables[[i]])))-1)) / ((length(na.omit(get(variables[[i]])))-2)*(length(na.omit(get(variables[[i]])))+1)*(length(na.omit(get(variables[[i]])))+3))),")
      }
      if(is.element("stddev",statistics)){
        express <- paste0(express,"stddev=sd(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("sum",statistics)){
        express <- paste0(express,"sum=sum(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("variance",statistics)){
        express <- paste0(express,"variance=var(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("default",statistics)){
        express <- paste0(express,"mean=mean(get(variables[[i]]),na.rm=T),sd=sd(get(variables[[i]]),na.rm=T),n=length(na.omit(get(variables[[i]]))),max=max(get(variables[[i]]),na.rm=T),min=min(get(variables[[i]]),na.rm=T),")
      }
      if(is.element("all",statistics)){
        express <- paste0(express,"n=length(na.omit(get(variables[[i]]))),kurt=kurtosis(get(variables[[i]]),type=2),max=max(get(variables[[i]]),na.rm=T),mean=mean(get(variables[[i]]),na.rm=T),median=median(get(variables[[i]]),na.rm=T),min=min(get(variables[[i]]),na.rm=T),range=diff(range(get(variables[[i]]),na.rm=T)),sekurt=2*(sqrt((6*length(na.omit(get(variables[[i]])))*(length(na.omit(get(variables[[i]])))-1)) / ((length(na.omit(get(variables[[i]])))-2)*(length(na.omit(get(variables[[i]])))+1)*(length(na.omit(get(variables[[i]])))+3))))* (sqrt(((length(na.omit(get(variables[[i]])))^2 -1)) / ((length(na.omit(get(variables[[i]])))-3)*(length(na.omit(get(variables[[i]])))+5)))),semean=sd(get(variables[[i]]),na.rm=T)/(sqrt(length(na.omit(get(variables[[i]]))))),skew=skewness(get(variables[[i]]),type=2),seskew=sqrt((6*length(na.omit(get(variables[[i]])))*(length(na.omit(get(variables[[i]])))-1)) / ((length(na.omit(get(variables[[i]])))-2)*(length(na.omit(get(variables[[i]])))+1)*(length(na.omit(get(variables[[i]])))+3))),stddev=sd(get(variables[[i]]),na.rm=T),sum=sum(get(variables[[i]]),na.rm=T),variance=var(get(variables[[i]]),na.rm=T),mean=mean(get(variables[[i]]),na.rm=T),sd=sd(get(variables[[i]]),na.rm=T),n=length(na.omit(get(variables[[i]]))),max=max(get(variables[[i]]),na.rm=T),min=min(get(variables[[i]]),na.rm=T),")
      }
      
      express <- str_sub(string=express,start=1,end=str_length(express)-1)  
      express <- paste0(express,")")
      descr <- descr[,eval(parse(text=express)),by=splitter]
      descr[order(splitter)]
    }
  } else {
    for(i in 1:length(variables))      {
      tempn <- N <- length(descr[[i]])
      if("mean" %in% statistics)
      {
        tempmean <- mean(descr[[i]])
      }
      if("min" %in% statistics)
      {
        tempmin <- min(descr[[i]])
      }
      if("max" %in% statistics)
      {
        tempmax<- max( descr[[i]])
      }
      if("stddev" %in% statistics)
      {
        tempstddev <- sd(descr[[i]])
      }
      if("kurtosis" %in% statistics)  {
        tempkurtosis <- kurtosis(descr[[i]], type=2)   
        tempsekurtosis <- 2*(sqrt((6*N*(N-1)) / ((N-2)*(N+1)*(N+3))))* (sqrt(((N^2 -1)) / ((N-3)*(N+5))))
      }  
      if("range" %in% statistics)  {
        temprange <- diff(range(descr[[i]]))
      } 
      if("skewness" %in% statistics)  {
        tempskewness <- skewness(descr[[i]], type=2)  
        tempseskewness <- sqrt((6*N*(N-1)) / ((N-2)*(N+1)*(N+3)))
      } 
      if("semean" %in% statistics)  {
        
        tempsemean <- sd(na.omit(descr[[i]])/(sqrt(N)))   
      }   
      if("sum" %in% statistics)  {
        tempsum <- sum( descr[[i]])    
      }   
      if("variance" %in% statistics)  {
        tempvariance <- var( descr[[i]])    
      } 
      if("all" %in% statistics)    {
        tempn <- length(descr[[i]])
        tempmean <- mean(descr[[i]])
        tempmin <- min(descr[[i]])
        tempmax <- max(descr[[i]])
        tempstddev <- sd( descr[[i]])
        tempkurtosis <- kurtosis( descr[[i]], type=2)
        tempsekurtosis <-  sqrt(((N^2 -1)) / ((N-3)*(N+5)))
        temprange <- diff(range(descr[[i]]))
        tempsemean <- sd(descr[[i]])/(sqrt(N))   
        tempskewness <-  skewness( descr[[i]], type=2)
        tempseskewness  <- sqrt((6*N*(N-1)) / ((N-2)*(N+1)*(N+3)))
        tempsum  <-  sum( descr[[i]])
        tempvariance  <- var( descr[[i]])
      } 
      if("default" %in% statistics) {
        tempn <- length(descr[[i]])
        tempmean <- mean(descr[[i]])
        tempmin <- min(descr[[i]])
        tempmax <- max(descr[[i]])
        tempstddev <- sd( descr[[i]])
      }
      descr[[i]] <- list("n" = tempn,
                         "mean" = tempmean,
                         "min" = tempmin, 
                         "max" = tempmax, 
                         "stddev" = tempstddev,
                         "kurtosis" = tempkurtosis,
                         "sekurtosis" = tempsekurtosis,
                         "range" = temprange,
                         "semean" = tempsemean,
                         "skewness" = tempskewness,
                         "seskewness" = tempseskewness,
                         "sum" = tempsum,
                         "variance" = tempvariance)
    }
  }
  
  
  
  
  for(i in 1:length(variables)){
    pos <- which(F==lapply(descr[[i]], is.null))
    descr[[i]] <- descr[[i]][pos]
  }
  
  options(warn=0)
  descr <- noquote(descr)
  return(descr)  
}
