#' Creates a xpssFrame Object
#'
#' xpssFrame creates a dataset as an xpssFrame object from a local file.
#'
#' @param x as character string with the name of the file.
#' @param \dots Arguments to pass on read.spss() from foreign.
#' @details  \code{x} the input data should be of the format \emph{.sav}.
#' 
#' The SPSS dataset contains the following attributes: \code{names}, \code{codepage}, \code{row.names}, \code{FILTER}, \code{TEMPORAY}, \code{SPLIT_FILE}, \code{DO_IF}, \code{SELECT_IF}, \code{WEIGHTS}, \code{class}
#' 
#' \tabular{rlll}{
#'  \tab Attribute \tab Type \tab Contain \cr
#' \tab \code{names} \tab atomic character or character vector \tab Name of the variable in the dataset \cr
#' \tab \code{codepage} \tab atomic numeric \tab ANSI code, which describe the encoding \cr
#' \tab \code{row.names}  \tab atomic numeric or numeric vector \tab row ID's  \cr
#' \tab \code{FILTER} \tab atomic logical or condition \tab Filtercondition or FALSE if \code{\link{xpssFilter}} is off\cr
#' \tab \code{TEMPORARY} \tab atomic logical \tab Logical statement if \code{\link{xpssTemporary}} is on or off \cr
#' \tab \code{SPLIT_FILE} \tab atomic logical \tab Split-File condition or FALSE if \code{\link{xpssSplitFile}} is off \cr
#' \tab \code{DO_IF} \tab atomic logical or condition \tab Do-If-condition or FALSE if \code{\link{xpssDoIf}} is off \cr
#' \tab \code{SELECT_IF} \tab atomic logical or condition \tab Select-If-condition or FALSE if \code{\link{xpssSelectIf}} is off \cr
#' \tab \code{WEIGHTS} \tab atomic character \tab  * will be implemented in a following update \cr
#' \tab \code{class} \tab atomic character or character \tab atomic or vector of names of classes the datasets inherits from.  \cr
#'}
#' 
#' Variable attributes are stored in the variable itself. A variable can have the following attributes: 
#' \code{defined.MIS}, \code{MIS}, \code{value.labels}, \code{variable.label}, \code{varname}
#'  \tabular{rlll}{
#'  \tab Attribute \tab Type \tab Contain \cr
#' \tab \code{defined.MIS} \tab Atomic numerics or atomic characters, respectively a numeric vector or character vector \tab Values which specify missing values \cr
#' \tab \code{MIS} \tab list with user-defined missings \tab POS contain the position of the user-defined missing, VAL the value of the user defined-missing \cr
#' \tab \code{value.labels}  \tab Named numeric or named character \tab  Value and label for a specific variable\cr
#' \tab \code{variable.label} \tab Atomic character \tab Label of the variable\cr
#' \tab \code{varname} \tab Atomic character \tab Name of the variable in the datasheet\cr
#'}
#' @author Andreas Wygrabek
#' @seealso \code{\link{read.spss}} \code{\link{as.xpssFrame}}
#' @importFrom foreign read.spss
#' @examples  \dontrun{
#' data <- xpssFrame(x="Testdata_1.sav")
#' }
#' 
#' @export
xpssFrame <- function(x, ...){
  
  
  options(warn=-1)
  if (!(is.character(x))){
    stop("Input has to be a string")}
  data <- suppressWarnings(read.spss(x,reencode="UTF-8",use.value.labels=F,to.data.frame=F))
  
  classBackUp <- sapply(data,class)
  for(xxx in colnames(data)){
    
    
    missings <- eval(parse(text=paste("attributes(data)$missings$",xxx,"$value", sep = "")))
    if(!(is.null(missings))){
      if(is.element(eval(parse(text=paste("attributes(data)$missings$",xxx ,"$type", sep = ""))),"range")) {
        eval(parse(text = paste("attr(data$",xxx,",'defined.MIS')$range <- c(from=missings[[1]],to=missings[[2]])",sep = ""))) 
        
      } else {
        eval(parse(text = paste("attr(data$",xxx,",'defined.MIS')$values <- missings",sep = ""))) 
      }
    }       
    
    if(length(missings > 0)){
      
      if(is.element(eval(parse(text=paste("attributes(data)$missings$",xxx ,"$type", sep = ""))),"range")) {
        misses <- vector()
        k <- 1
        for(j in 1:length(eval(parse(text=paste("data$",xxx, sep = ""))))){
          if((is.na(eval(parse(text=paste("data$",xxx, sep = "")))[j]) == FALSE) && (eval(parse(text=paste("attributes(data$",xxx,")$defined.MIS$range[[1]]", sep = ""))) <= eval(parse(text=paste("data$",xxx, sep = "")))[j]) && (eval(parse(text=paste("data$",xxx, sep = "")))[j] <= eval(parse(text=paste("attributes(data$",xxx,")$defined.MIS$range[[2]]", sep = "")))))
          {            
            misses[[k]] <- j
            k <- k+1
          }
        }
        POS <- misses
        eval(parse(text = paste("VAL <- data$",xxx,"[POS]", sep = "")))
        eval(parse(text = paste("missings <- cbind(POS, VAL)", sep = "")))
        
      } else {
        
        eval(parse(text = paste("POS <- which(is.element(data$",xxx,", missings))", sep = "")))
        eval(parse(text = paste("VAL <- data$",xxx,"[POS]", sep = "")))
        eval(parse(text = paste("missings <- cbind(POS, VAL)", sep = "")))
      }
      
      
      eval(parse(text = paste("attr(data$",xxx,",'MIS') <- missings",sep = ""))) 
    }
  }
  
  for(i in names(data)){
    
    if(eval(parse(text = paste("is.null(attributes(data)$label.table$",i,")", sep = "")))) {
      print(paste("No labels existing in", i, sep = " "))
    } else if (eval(parse(text = paste("length(attributes(data)$label.table$",i,")< ifelse(length(which(is.na(as.numeric(data$",i,"))))>0, length(unique(data$",i,"))-1,length(unique(data$",i,")))", sep = "")))){
      
      print(paste("Not every number is labelled in", i, sep = " "))
      
    } else {
      valLab <- character()
      
      #remove wihte spaces
      #eval(parse(text = paste("attributes(data)$label.table <- names(attributes(data)$label.table$",i,")", sep = "")))      
                 eval(parse(text = paste("valLab <- names(attributes(data)$label.table$",i,")", sep = "")))
                 eval(parse(text = paste("numbers <- attributes(data)$label.table$",i, sep = "")))
                 
                 
                 newVec <- numeric(length = nrow(data))
                 eval(parse(text = paste("newVec[which(is.na(data$",i,"))] <- NA", sep = "")))
                 
                 for(a in 1:length(valLab)){
                   
                   ########
                   
                   eval(parse(text = paste("logVec <- data$",i," == valLab[",a,"]", sep = "")))
                   eval(parse(text = paste("newVec[logVec] <- numbers[",a,"]", sep = "")))
                   
                   ########
                   
                 }
                 #eval(parse(text = paste("data$",i," <- newVec",sep = "")))
                 eval(parse(text = paste("attr(data$",i,", 'value.labels') <- numbers", sep = "")))
          }    
  }
  
  if(!is.null(attributes(data)$variable.labels))backup_varLabs <- attributes(data)$variable.labels
  if(!is.null(attributes(data)$names))backup_names <- attributes(data)$names
  if(!is.null(attributes(data)$class))backup_class <- attributes(data)$class
  if(!is.null(attributes(data)$rownames))backup_rownames <- attributes(data)$row.names
  
  if(exists("backup_class"))attr(data, "class") <- backup_class
  if(exists("backup_rownames")) attr(data, "rownames") <- backup_rownames
  
  if(exists("backup_names")){
    for(i in colnames(data)){
      eval(parse(text = paste("attr(data$",i,", 'varname') <- backup_names[which('",i,"' == colnames(data))]", sep = "")))
    }
  }
  
  if(exists("backup_varLabs")){
    for(i in colnames(data)){
      eval(parse(text = paste("attr(data$",i,", 'variable.label') <- backup_varLabs[which('",i,"' == colnames(data))][[i]]", sep = "")))
    }
  }
  
  
  
  
  for(i in colnames(data)){
    
    eval(parse(text = paste("data$",i,"[attributes(data$",i,")$MIS[,1]] <- NA", sep = "")))
    
  }
  
  ##############################################################
  ##############################################################
  
  
  attr(data, "FILTER") <- FALSE
  attr(data, "TEMPORARY") <- FALSE
  attr(data, "SPLIT_FILE") <- FALSE
  attr(data, "DO_IF") <- FALSE
  attr(data, "SELECT_IF") <- FALSE
  attr(data, "WEIGHTS") <- "none"
  
  ##############################################################
  ##############################################################
  ##############################################################
  
  attributes(data)$label.table <- NULL
  attributes(data)$missings <- NULL
  attributes(data)$variable.labels <- NULL
  
  
  attributeframe <- list()
  for(i in 1:length(data)){
    attributeframe[[i]] <- attributes(data[[i]])
  }
      
  for(i in 1:length(data)){
    if(!(any(grepl("[A-Za-z_]",data[[i]])))){
      data[[i]] <- as.numeric(data[[i]])
    }
    attributes(data[[i]]) <- attributeframe[[i]]
  }

  data <- as.data.frame(data)
  class(data) <- c("xpssFrame", "data.frame")
 
  options(warn=0)
  return(data)
}


