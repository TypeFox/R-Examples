#
# file: estimate.R
#
# This file is part of the R-package decisionSupport
#
# Authors:
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF)
#	http://www.worldagroforestry.org
#
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' @include rmvnorm90ci_exact.R
#' @include random.R
#' @include estimate1d.R
NULL
# Define global variables (ToDo: make local if possible):
if(getRversion() >= "2.15.1")  utils::globalVariables(c("variable",
                                                        "distribution"))
##############################################################################################
# estimate(distribution, lower, upper,..., correlation_matrix)
##############################################################################################
#' Create a multivariate estimate object.
#' 
#' \code{estimate} creates an object of \code{class estimate}. The concept of an estimate is 
#' extended from the 1-dimensional (cf. \code{\link{estimate1d}}) to the multivariate case. This 
#' includes the description of correlations between the different variables. An estimate of an 
#' n-dimensional variable is at minimum defined by each component being a 1-dimensional estimate. 
#' This means, that for each component, at minimum, the  type of its univariate parametric 
#' distribution, its 5\% - and 95\% quantiles must be provided. In probability theoretic terms, 
#' these are the marginal distributions of the components. Optionally, the individual median 
#' and the correlations between the components can be supplied.
#' @param distribution \code{character vector}: defining the types of the univariate parametric 
#' distributions. 
#' @param lower \code{numeric vector}: lower bounds of the 90\% confidence intervals, i.e the 5\%-quantiles 
#'   of this estimates components.
#' @param upper \code{numeric vector}: upper bounds of the 90\% confidence intervals, i.e the 95\%-quantiles
#'   of this estimates components.
#' @param ... in \code{estimate}: optional arguments that can be coerced to a data frame comprising
#'   further columns of the estimate (for details cf. below).\cr
#'           in \code{as.estimate}: arguments that can be coerced to a data frame comprising the 
#'           marginal distributions of the estimate components. Mandatory columns are \code{distribution}, \code{lower} and 
#'           \code{upper}. 
#' @param correlation_matrix \code{numeric matrix}: containing the correlations of the variables 
#'   (optional).
#' @details
#'   The input arguments inform the estimate about its marginal distributions and joint distribution, i.e.
#'    the correlation matrix.
#'   \subsection{The structure of the estimates marginal input information}{
#'     \describe{
#'       \item{in \code{estimate}}{
#'       The marginal distributions are defined by the arguments \code{distribution}, \code{lower}
#'       and \code{upper} and, optionally, by further columns supplied in \code{...} that can be 
#'       coerced to a \code{\link{data.frame}} with the same length as the mandatory arguments.
#'       }
#'       \item{in \code{as.estimate}}{
#'         The marginal distributions are completely defined in \code{...}. These arguments must be 
#'         coercible to a data.frame, all having the same length. Mandatory columns are 
#'        \code{distribution}, \code{lower} and  \code{upper}. 
#'       }
#'     }
#'     \subsection{Mandatory input columns}{
#'     \tabular{lll}{
#'       \bold{Column}       \tab  \bold{R-type}    \tab \bold{Explanation}\cr
#'       \code{distribution} \tab  \code{character vector} \tab  Marginal distribution types \cr
#'       \code{lower}        \tab  \code{numeric vector}   \tab  Marginal 5\%-quantiles \cr
#'       \code{upper}        \tab  \code{numeric vector}   \tab  Marginal 95\%-quantiles 
#'     }
#'     It must hold that \code{lower <= upper} for every component of the estimate.  
#'     }
#'     \subsection{Optional input columns}{
#'     The optional parameters in \code{...} provide additional characteristics of the marginal 
#'     distributions of the estimate. Frequent optional columns are:
#'     \tabular{lll}{
#'       \bold{Column}       \tab  \bold{R-type}                 \tab \bold{Explanation}\cr
#'       \code{variable}     \tab  \code{character vector}       \tab  Variable names\cr
#'       \code{median}       \tab  cf. below                     \tab  Marginal 50\%-quantiles \cr
#'       \code{method}       \tab  \code{character vector}       \tab  Methods for calculation of marginal distribution parameters
#'    }
#'    \subsection{The \code{median} column}{ 
#'      If supplied as input, any component of \code{median} can be either \code{NA}, \code{numeric}
#'      (and not \code{NA}) or the character string \code{"mean"}. If it equals \code{"mean"} it is
#'      set to \code{rowMeans(cbind(lower, upper))} of this component; if it is \code{numeric} it must
#'      hold that \code{lower <= median <= upper} for this component. In case that no element 
#'      \code{median} is provided, the default is \code{median=rep(NA, length(distribution))}.\cr
#'      The \code{median} is important for the different methods possible in generating the random 
#'      numbers (cf. \code{\link{random.estimate}}).
#'    } 
#'    }
#'  }
#'  \subsection{The structure of the estimates correlation input information}{
#'    The argument \code{correlation_matrix} is the sub matrix of the full correlation matrix of 
#'    the estimate containing all correlated elements. Thus, its row and column names must be a 
#'    subset of the variable names of the marginal distributions. This means, that the information  
#'    which variables are uncorrelated does not need to be provided explicitly.
#'    
#'    \code{correlation_matrix} must have all the properties of a correlation matrix, viz. symmetry, 
#'    all diagonal elements equal 1 and all of diagonal elements are between -1 and 1. 
#'  }
#'  
#' @return  An object of class \code{estimate} which is a list with components \code{$marginal} and 
#' \code{$correlation_matrix}:
#'   \describe{
#'     \item{\code{$marginal}}{
#'     is a \code{\link{data.frame}} with mandatory columns:  
#'     \tabular{lll}{
#'       \bold{Mandatory column}      \tab  \bold{R-type}                 \tab \bold{Explanation}\cr
#'       \code{distribution} \tab  \code{character vector} \tab  Distribution types \cr
#'       \code{lower}        \tab  \code{numeric vector}   \tab   5\%-quantiles\cr
#'       \code{median}       \tab  \code{numeric vector}   \tab  50\%-quantiles or \code{NA}\cr 
#'       \code{upper}        \tab  \code{numeric vector}   \tab  95\%-quantiles 
#'     }
#'     The \code{\link{row.names}} are the names of the variables. Each row has the properties of 
#'     an \code{\link{estimate1d}}. 
#'     
#'     Note that the \emph{\code{median}} is a mandatory element of an \code{estimate}, although it
#'     is not necessary as input. If a component of \code{median} is numeric and not \code{NA} it 
#'     holds that: \code{lower <= median <= upper}. In any case an \code{estimate} object has the
#'     property \code{any(lower <= upper)}.     
#'     }
#'     \item{\code{$correlation_matrix}}{
#'      is a symmetric matrix with row and column names being the subset of the variables supplied 
#'      in \code{$marginal} which are correlated. Its elements are the corresponding correlations.
#'      }
#'   }
#'   
#' @seealso  \code{\link{estimate1d}}, \code{\link{random.estimate}}, 
#' \code{\link{row.names.estimate}}, \code{\link{names.estimate}}, \code{\link{corMat}}, 
#' \code{\link{estimate_read_csv}} and \code{\link{estimate_write_csv}}.
#' @examples
#' # Create a minimum estimate (only mandatory marginal information supplied):
#' estimateMin<-estimate(c("posnorm", "lnorm"),
#'                       c(        4,       4),
#'                       c(       50,      10))
#' print(estimateMin) 
#' 
#' # Create an estimate with optional columns (only marginal information supplied):
#' estimateMarg<-estimate(           c("posnorm", "lnorm"),
#'                                   c(        4,       4),
#'                                   c(       50,      10),
#'                          variable=c("revenue", "costs"),
#'                          median = c(   "mean",      NA),
#'                          method = c(    "fit",      ""))
#' print(estimateMarg)
#' print(corMat(estimateMarg))
#' 
#' @export
estimate<-function(distribution, lower, upper, ..., correlation_matrix=NULL){
  # Check preconditions:
  ## Marginal mandatory arguments:
  ### Check argument types:
  if ( any(is.null(distribution) || !is.character(distribution)) )
    stop("\"distribution\" must be supplied as character string.")
  if ( any(is.null(lower) || is.na(lower<-as.numeric(lower))) )
    stop("\"lower\" must be supplied as numeric.")
  if ( any(is.null(upper) || is.na(upper<-as.numeric(upper))) )
    stop("\"upper\" must be supplied as numeric.")
  #### Check equality of dimension:
  if ( length(distribution) != length(lower)  || length(distribution) != length(upper) )
    stop("All input dimensions must be equal.")
  ### Check input semantics:
  if ( any(lower > upper) )
    stop("\"lower > upper\" for some variables")
  ## Marginal optional arguments:
  marginalOptional<-if(missing(...)) NULL else data.frame(..., stringsAsFactors=FALSE)
  median<-NULL
  if( !is.null(marginalOptional) )
    for ( i in names(marginalOptional) ){
      ### Check argument types:
      if (i == "variable" && any(!is.character(marginalOptional[[i]]))) 
        stop("Optional argument \"", i, "\" is not character for all entries.")
      if (i == "method" && any(!is.character(marginalOptional[[i]]))) 
        stop("Optional argument \"", i, "\" is not character for all entries.")
      #Process median:
      if (i == "median"){
        median<-marginalOptional[[i]]
        marginalOptional<-marginalOptional[!names(marginalOptional) %in% "median"]
        if ( !is.null(median) ){
          if ( length(distribution) != length(median))
            stop("\"median\" is not of the right length.")
          median[!is.na(median) & is.character(median) & median=="mean"]<-rowMeans(
            cbind(lower, upper))[!is.na(median) & is.character(median) & median=="mean"]
          #### Replace empty character with NA
          median[!is.na(median) & is.character(median) & median==""]<-
            rep(NA, length(distribution))[!is.na(median) & is.character(median) & median==""]
          #if( any(!is.na(median) & is.numeric(median) & (lower > median | median > upper)) )
          if ( any(is.na(median) & !is.na(median<-as.numeric(median))) ) 
            stop("\"median\": all values must be one of the following: \"numeric\", \"character\" 
                 with value \"mean\", \"NA\" or \"\".")
          else if( any(!is.na(median) & (lower > median | median > upper)) )
            stop("It must hold: \"lower <= median <= upper\" for all entries")
        }
        else 
          median<-rep(NA, length(distribution))
      }
    }
  if (is.null(median)) median<-rep(NA, length(distribution))
  ## Correlation matrix precondition check:
  if( !is.null(correlation_matrix)){
    if( !is.matrix(correlation_matrix) )
      correlation_matrix<-as.matrix(correlation_matrix)
    if( !identical( correlation_matrix, t(correlation_matrix) ) )
      stop("correlationMatrix must be a symmetric matrix.")
    if( !identical( as.vector(diag(correlation_matrix)), rep(1, nrow(correlation_matrix)) ) )
      stop("All diagonal elements of \"correlation_matrix\"  must be equal to 1.")
    # Check that all elements are between -1 and 1.
    if( any(abs(correlation_matrix) > 1 ) )
      stop("All values of \"correlation_matrix\" must be  >= -1 and <= 1.")
    #ToDo: Check that all rows are named
    #ToDo: check that rownames(correlation_matrix) is a subset of marginal names
  }
  
  # Create the marginal estimate:
  if( as.logical(length(marginalOptional)) )
    marginal<-data.frame(distribution=distribution, 
                         lower=lower, 
                         median=median, 
                         upper=upper, 
                         marginalOptional,
                         row.names=row.names(marginalOptional),
                         stringsAsFactors=FALSE) 
  else
    marginal<-data.frame(distribution=distribution, 
                         lower=lower, 
                         median=median, 
                         upper=upper, 
                         row.names=row.names(marginalOptional),
                         stringsAsFactors=FALSE) 
  if( !is.null(marginal$variable) ){
    rownames(marginal)<-marginal$variable
    marginal<-marginal[!colnames(marginal) %in% "variable"]
  }
  # Drop rows without variable name:
  marginal<-subset(marginal, row.names(marginal) != "")
  if( is.null(marginal$distribution) )
    stop("marginal must be supplied with a distribution column.")
  
  # Return object:
  returnObject=list(marginal=marginal , correlation_matrix=correlation_matrix)
  #ToDo: or  class(estimateObject)<-c("estimate", "data.frame") ???
  class(returnObject)<-"estimate"
  returnObject
}
##############################################################################################
# as.estimate( ..., correlation_matrix)
##############################################################################################
#' Coerce and transform to a multivariate estimate object.
#'
#' \code{as.estimate} tries to coerce a set of objects and transform them to \code{class estimate}.
#' @rdname estimate
#' @examples
#' # Create a minimum estimate from text (only mandatory marginal information supplied):
#' estimateTextMin<-"distribution, lower, upper
#'                   posnorm,      100,   1000
#'                   posnorm,      50,    2000
#'                   posnorm,      50,    2000
#'                   posnorm,      100,   1000"
#' estimateMin<-as.estimate(read.csv(header=TRUE, text=estimateTextMin, 
#'                           strip.white=TRUE, stringsAsFactors=FALSE))
#' print(estimateMin) 
#' 
#' # Create an estimate from text (only marginal information supplied):
#' estimateText<-"variable,  distribution, lower, upper, median, method
#'                revenue1,  posnorm,      100,   1000,  NA,        
#'                revenue2,  posnorm,      50,    2000,    ,     fit
#'                costs1,    posnorm,      50,    2000,  70,     calculate
#'                costs2,    posnorm,      100,   1000,  mean,             "
#' estimateMarg<-as.estimate(read.csv(header=TRUE, text=estimateText, 
#'                           strip.white=TRUE, stringsAsFactors=FALSE))
#' print(estimateMarg)
#' print(corMat(estimateMarg))
#' 
#' # Create an estimate from text (with correlated components): 
#' estimateTextMarg<-"variable,  distribution, lower, upper
#'                    revenue1,  posnorm,      100,   1000
#'                    revenue2,  posnorm,      50,    2000
#'                    costs1,    posnorm,      50,    2000
#'                    costs2,    posnorm,      100,   1000"
#' estimateTextCor<-",         revenue1, costs2
#'                   revenue1,        1,   -0.3
#'                   costs2,       -0.3,      1"
#' estimateCor<-as.estimate(read.csv(header=TRUE, text=estimateTextMarg, 
#'                           strip.white=TRUE, stringsAsFactors=FALSE),
#'                           correlation_matrix=data.matrix(read.csv(text=estimateTextCor, 
#'                                                                   row.names=1,
#'                                                                   strip.white=TRUE)))
#' print(estimateCor)
#' print(corMat(estimateCor))
#' @export
as.estimate<-function(..., correlation_matrix=NULL){
  # Coerce the marginal data.frame:
  # ToDo: if (...) contains an estimate: process correlation matrix!
  marginal<-data.frame(..., stringsAsFactors=FALSE)
  # Check preconditions:
  if (  is.null(marginal$distribution) )
    stop( "no \"distribution\" column!")
  if (  is.null(marginal$lower) )
    stop( "no \"lower\" column!")
  if (  is.null(marginal$upper)  )
    stop( "no \"upper\" column!")
  # Create and return the estimate:
  estimate(distribution=marginal[["distribution"]], 
           lower=marginal[["lower"]],
           upper=marginal[["upper"]],
           marginal[!names(marginal) %in% c("distribution", "lower", "upper")],
           correlation_matrix=correlation_matrix)
}
##############################################################################################
# row.names.estimate(x)
##############################################################################################
#' Get and set attributes of an \code{estimate} object.
#'
#' \code{row.names.estimate} returns the variable names of an \code{\link{estimate}} object which 
#' is identical to \code{row.names(x$marginal)}.
#' @param x an \code{\link{estimate}} object.
#' @seealso \code{\link{estimate}}, \code{\link{names.estimate}}, \code{\link{corMat.estimate}}, 
#'   \code{\link{corMat}}
#' @examples
#'  # Read the joint estimate information for the variables "sales", "productprice" and 
#'  # "costprice" from file:
#'  ## Get the path to the file with the marginal information:
#'  marginalFilePath=system.file("extdata","profit-4.csv",package="decisionSupport")
#'  ## Read marginal and correlation file into an estimate:
#'  parameterEstimate<-estimate_read_csv(fileName=marginalFilePath)
#'  print(parameterEstimate)
#'  ## Print the names of the variables of this estimate
#'  print(row.names(parameterEstimate))
#' @export
row.names.estimate<-function(x){
  row.names(x$marginal)
}
##############################################################################################
# names.estimate(x)
##############################################################################################
#' Return the column names of an \code{estimate} object.
#'
#' \code{names.estimate} returns the column names of an \code{\link{estimate}} object which is identical to
#' \code{names(x$marginal)}.
#' @rdname row.names.estimate
#' @examples
#'  ## Print the names of the columns of this estimate
#'   print(names(parameterEstimate))
#' @export
names.estimate<-function(x){
  names(x$marginal)
}
##############################################################################################
# generic: corMat(rho)
##############################################################################################
#' Return the Correlation Matrix.
#'
#' Return the correlation matrix of rho.
#' @param rho a distribution.
#' @export
corMat <- function(rho) UseMethod("corMat")
##############################################################################################
# corMat.estimate(rho)
##############################################################################################
#' Return the correlation matrix of an \code{estimate} object.
#'
#' \code{corMat.estimate} returns the full correlation matrix of an \code{\link{estimate}} object.
#' @param rho an \code{\link{estimate}} object.
#' @rdname row.names.estimate
#' @examples
#'  ## Print the full correlation matrix of this estimate
#'   print(corMat(parameterEstimate))
#' @export
corMat.estimate<-function(rho){
  # Create identity matrix:
  corMat<-diag(nrow=length(row.names(rho)))
  dimnames(corMat)<-list(row.names(rho), row.names(rho))
  # Replace the values for correlated elements:
  namesCorrelated<-row.names(rho$correlation_matrix)
  corMat[namesCorrelated, namesCorrelated]<-rho$correlation_matrix[namesCorrelated, namesCorrelated]
  # Return full correlation matrix:
  corMat
}
##############################################################################################
# generic: corMat<-(rho, correlationMatrix)
##############################################################################################
#' Replace correlation matrix.
#'
#' Replace the correlation matrix.
#' @param x a distribution.
#' @param value \code{numeric matrix}: new correlation matrix.
#' @export
`corMat<-`<-function(x, value) UseMethod("corMat<-")
##############################################################################################
# corMat.estimate(x)
##############################################################################################
#' Replace the correlation matrix of an \code{estimate} object.
#'
#' \code{'corMat<-.estimate'} replaces the correlation matrix of an \code{\link{estimate}} object.
#' @param value \code{numeric matrix}: new correlation matrix. For details cf. 
#'   \code{\link{estimate}}.
#' @rdname row.names.estimate
#' @seealso \code{\link{corMat<-}}
#' @export
`corMat<-.estimate` <- function(x, value){ 
  ## Correlation matrix precondition check:
  if( !is.null(value)){
    if( !is.matrix(value) )
      value<-as.matrix(value)
    if( !identical( value, t(value) ) )
      stop("value must be a symmetric matrix.")
    if( !identical( as.vector(diag(value)), rep(1, nrow(value)) ) )
      stop("All diagonal elements of \"value\"  must be equal to 1.")
    # Check that all elements are between -1 and 1.
    if( any(abs(value) > 1 ) )
      stop("All values of \"value\" must be  >= -1 and <= 1.")
    # Check that all rows are named
    if (  length(row.names(value)) != nrow(value) )
      stop("All rows of \"value\" must be named.")
    # Check that rownames(value) is a subset of marginal names
    if (  !any(row.names(value) %in% row.names(x)) )
      stop("Names of \"value\" must be a subset of \"row.name(x)\"") 
    # Check for uncorrelated variables and eliminate these rows and columns:
    correlatedVariables<-sapply(X=row.names(value), 
                                FUN=function(x) any((value!=0)[x, !row.names(value) %in% x]))
    x$correlation_matrix<-value[correlatedVariables, correlatedVariables]
    if( nrow(x$correlation_matrix) <= 1)
      x$correlation_matrix<-NULL
  } 
  else 
    x$correlation_matrix<-NULL
  # Return processed estimate:
  x
}
###############################################################################################
# estimate_read_csv(fileName, strip.white=TRUE, ...)
##############################################################################################
#' Read an Estimate from CSV - File.
#'
#' This function reads an \code{\link{estimate}} from the specified csv files. In this context, an 
#' estimate of several variables is defined by its marginal distribution types, its marginal
#' 90\%-confidence intervals \code{[lower,upper]} and, optionally, its correlations.
#' @param fileName Name of the file containing the marginal information of the estimate that 
#' should be read.
#' @inheritParams utils::read.csv
#' @param ... Further parameters to be passed to \code{\link[utils]{read.csv}}.
#' @return An object of type \code{\link{estimate}} which element \code{$marginal} is read from 
#'  file \code{fileName} and which element \code{$correlation_matrix} is read from file
#'  \code{gsub(".csv","_cor.csv",fileName)}.
#' @details An estimate might consists of uncorrelated and correlated variables. This is reflected
#'   in the input file structure, which is described in the following.
#'   \subsection{ CSV input file structures}{
#'   The estimate is read from one or two csv files: the marginal csv file which is mandatory and
#'   the correlation csv file which is optional. The marginal csv file contains the definition of
#'   the distribution of all variables ignoring potential correlations. The correlation csv file
#'   only defines correlations. \subsection{The structure of the marginal distributions input file
#'   (mandatory)}{
#'     File name structure: \code{<marginal-filename>.csv}
#'     
#'     Mandatory columns:
#'      \tabular{lll}{
#'       \bold{Column name}       \tab  \bold{R-type}    \tab \bold{Explanation}\cr
#'       \code{variable}          \tab  \code{character vector}       \tab  Variable names\cr
#'       \code{distribution}      \tab  \code{character vector} \tab  Marginal distribution types \cr
#'       \code{lower}             \tab  \code{numeric vector}   \tab  Marginal 5\%-quantiles \cr
#'       \code{upper}             \tab  \code{numeric vector}   \tab  Marginal 95\%-quantiles 
#'      }
#'      Frequent optional columns are:
#'      \tabular{lll}{
#'        \bold{Column name}       \tab  \bold{R-type}              \tab \bold{Explanation}\cr
#'        \code{description}       \tab  \code{character}           \tab  Short description of the variable.\cr
#'        \code{median}            \tab  cf. \code{\link{estimate}} \tab  Marginal 50\%-quantiles \cr
#'        \code{method}            \tab  \code{character vector}    \tab  Methods for calculation of marginal distribution parameters
#'      }
#'      Columns without names are ignored. Rows where the \code{variable} field is empty are also dropped.
#'   }
#'   \subsection{The structure of the correlation file (optional)}{
#'      File name structure: \code{<marginal-filename>_cor.csv}
#'      
#'      Columns and rows are named by the corresponding variables. Only those variables need to be 
#'      present which are correlated with others.
#'      
#'      The element \code{["rowname","columnname"]} contains the correlation between the variables 
#'      \code{rowname} and \code{columnname}. Uncorrelated elements have to be set to \code{0}. The
#'      diagonal element \code{["name","name"]} has to be set to \code{1}.
#'      
#'      The matrix must be given in symmetric form.
#'   }
#' }
#' @seealso \code{\link{estimate_write_csv}}, \code{\link[utils]{read.csv}}, \code{\link{estimate}}
#' @examples
#'  # Read the joint estimate information for the variables "sales", "productprice" and 
#'  # "costprice" from file:
#'  ## Get the path to the file with the marginal information:
#'  marginalFilePath=system.file("extdata","profit-4.csv",package="decisionSupport")
#'  ## Read the marginal information from file "profit-4.csv" and print it to the screen as
#'  ## illustration:
#'  read.csv(marginalFilePath, strip.white=TRUE)
#'  ## Read the correlation information from file "profit-4_cor.csv" and print it to the screen as
#'  ## illustration: 
#'  read.csv(gsub(".csv","_cor.csv",marginalFilePath), row.names=1)
#'  ## Now read marginal and correlation file straight into an estimate:
#'  parameterEstimate<-estimate_read_csv(fileName=marginalFilePath)
#'  print(parameterEstimate)
#' @export
estimate_read_csv <- function(fileName, strip.white=TRUE, ...){
  marginal<-NULL
  correlation_matrix<-NULL
  marginalFilename<-fileName
  # Read marginal data:
  #marginal<-read.csv(marginalFilename,row.names="variable", strip.white=strip.white, stringsAsFactors=FALSE, ...)
  marginal<-read.csv(marginalFilename, strip.white=strip.white, stringsAsFactors=FALSE, ...)
  # ToDo: replace subset() such that reference to global variable "variable" becomes obsolete:
  marginal<-subset(marginal,variable!="")
  marginal<-data.frame(marginal,row.names="variable")
  # Read correlation data:
  # Generate correlation filename:
  correlationFilename<-gsub(".csv","_cor.csv",marginalFilename)
  # Read correlation file if it exists:
  if(file.exists(correlationFilename))
    correlation_matrix<-data.matrix(read.csv(correlationFilename, row.names=1))
  
  # Return object
  as.estimate(marginal, correlation_matrix=correlation_matrix)
}
###############################################################################################
# estimate_write_csv(estimate, fileName, strip.white=TRUE, ...)
##############################################################################################
#' Write an Estimate to CSV - File.
#'
#' This function writes an \code{\link{estimate}} to the specified csv file(s).
#' @param fileName \code{character}: File name for the output of the marginal information of the 
#'   estimate. It must end with \code{.csv}.
#' @param estimate  \code{estimate}: Estimate object to write to file.
#' @param varNamesAsColumn \code{logical}: If \code{TRUE} the variable names will be written as a
#' separate column, otherwise as row names.
#' @param quote a \code{logical} value (TRUE or FALSE) or a numeric vector. If
#'   TRUE, any character or factor columns will be surrounded by double quotes.
#'   If a numeric vector, its elements are taken as the indices of columns to
#'   quote. In both cases, row and column names are quoted if they are written.
#'   If FALSE, nothing is quoted. Parameter is passed on to \code{\link{write.csv}}.
#' @param ... Further parameters to be passed to \code{\link[utils]{write.csv}}.
#' @details
#'   The marginal information of the \code{estimate} is written to file \code{fileName=<marginal-filename>.csv}. If 
#'   the estimate contains correlated variables, the correlation matrix is written to the separate
#'   file \code{<marginal-filename>_cor.csv}.
#' @seealso \code{\link{estimate_read_csv}}, \code{\link{estimate}}, \code{\link[utils]{write.csv}}
#' @export
estimate_write_csv <- function(estimate, fileName, varNamesAsColumn=TRUE, quote=FALSE, ...){
  marginalFilename=fileName
  # Write marginal data to file:
  if (varNamesAsColumn){
    marginal<-cbind(estimate$marginal,variable=row.names(estimate))
    row.names(marginal)<-NULL
  } else
    marginal<-estimate$marginal
  write.csv(x=marginal, file=marginalFilename, row.names=!varNamesAsColumn,  quote=FALSE, ...)
  # Write correlation data if exists:
  if( !is.null(estimate$correlation_matrix) ){
    # Generate correlation filename:
    correlationFilename<-gsub(".csv","_cor.csv",marginalFilename)
    # Wirte correlation file:
    write.csv(x=estimate$correlation_matrix, file=correlationFilename, quote=FALSE, ...)
  }
}
##############################################################################################
# random.estimate(rho,n,method, ...)
##############################################################################################
#' Generate random numbers for an estimate.
#'
#' This function generates random numbers for general multivariate
#' distributions that are defined as an \code{\link{estimate}}.
#' @param rho \code{estimate}: multivariate distribution to be randomly sampled.
#' @param n \code{integer}:Number of observations to be generated. 
#' @param method \code{character}: Particular method to be used for random number generation.
#' @param relativeTolerance \code{numeric}: the relative tolerance level of deviation of the
#'   generated confidence interval from the specified interval. If this deviation is greater than
#'   \code{relativeTolerance} a warning is given.
#' @param ... Optional arguments to be passed to the particular random number
#'  generating function.
#' @details
#' 	\subsection{Generation of uncorrelated components}{
#' 		Implementation: \code{\link{random.estimate1d}}
#'
#' 	}
#' 	\subsection{Generation of correlated components}{
#' 		Implementation: \code{\link{rmvnorm90ci_exact}}
#' 	}
#' @examples
#'	variable=c("revenue","costs")
#'	distribution=c("norm","norm")
#'  lower=c(10000,  5000)
#'  upper=c(100000, 50000)
#'  estimateObject<-as.estimate(variable, distribution, lower, upper)
#'  x<-random(rho=estimateObject, n=10000)
#'  apply(X=x, MARGIN=2, FUN=quantile, probs=c(0.05, 0.95))
#'  cor(x)
#'  colnames(x)
#'  summary(x)
#'  hist(x[,"revenue"])
#'  hist(x[,"costs"])
#'  
#'  # Create an estimate with median and method information:
#'  estimateObject<-estimate(         c("posnorm", "lnorm"),
#'                                    c(        4,       4),
#'                                    c(       50,      10),
#'                           variable=c("revenue", "costs"),
#'                           median = c(   "mean",      NA),
#'                           method = c(    "fit",      ""))
#'  # Sample random values for this estimate:
#'  x<-random(rho=estimateObject, n=10000)
#'  # Check the results 
#'  apply(X=x, MARGIN=2, FUN=quantile, probs=c(0.05, 0.95))
#'  summary(x)
#'  hist(x[,"revenue"], breaks=100)
#'  hist(x[,"costs"], breaks=100)
#'  
#' @seealso \code{\link{estimate}}, \code{\link{random.estimate1d}}, \code{\link{random}}
#' @export
random.estimate <- function(rho,n,method="calculate", relativeTolerance=0.05, ...){
  #ToDo: test
  x<-NULL
  if ( !is.null(rho$correlation_matrix) ){
    # Select correlated variables:
    rhoCorrelated<-list(marginal=NULL,correlation_matrix=NULL)
    class(rhoCorrelated)<-"estimateCorrelated"
    namesCorrelated<-row.names(rho$correlation_matrix)
    rhoCorrelated$marginal<-rho$marginal[namesCorrelated, ]
    rhoCorrelated$correlation_matrix<-rho$correlation_matrix
    # Generate correlated variables
    x<-random(rho=rhoCorrelated,
              n=n,
              method=method,
              relativeTolerance=relativeTolerance,
              ...)
    # Select uncorrelated variables if there are any:
    if( length(namesUnCorrelated
               <-row.names(rho$marginal[!(row.names(rho$marginal) %in% namesCorrelated ),]) ) ){
      rhoUnCorrelated<-rho$marginal[namesUnCorrelated, ]
      class(rhoUnCorrelated)<-c("estimateUnCorrelated", class(rhoUnCorrelated))
      x<-cbind(x, random(rho=rhoUnCorrelated, 
                         n=n, 
                         method=method,
                         relativeTolerance=relativeTolerance,
                         ...))
    }
  } else {
    class(rho$marginal)<-c("estimateUnCorrelated", class(rho$marginal))
    x<-random(rho=rho$marginal, 
              n=n, 
              method=method, 
              relativeTolerance=relativeTolerance,
              ...)
  }
  # Return the generated random variables:
  x
}

##############################################################################################
# random.estimateCorrelated(rho,n,method, relativeTolerance, ...)
##############################################################################################
# Generate the random numbers for the correlated subset of an estimate
#
random.estimateCorrelated <- function(rho, n, method, relativeTolerance=0.05, ...){
  x<-NULL
  if(method=="calculate"){
    if( identical( rho$marginal$distribution, rep("norm", nrow(rho$marginal)) ) ){
      x<-rmvnorm90ci_exact(n=n,
                           lower=data.matrix(rho$marginal["lower"]),
                           upper=data.matrix(rho$marginal["upper"]),
                           correlationMatrix=rho$correlation_matrix)
    }
    else
      stop("correlated variables must all be of type \"norm\".")
  }
  else
    stop ("method must be  \"calculate\".")
  # Return the generated random numbers:
  x
}
##############################################################################################
# random.estimateUnCorrelated(rho , n,method, relativeTolerance, ...)
##############################################################################################
# Generate the random numbers for the uncorrelated subset of an estimate
#
random.estimateUnCorrelated <- function(rho, n, method="calculate", relativeTolerance=0.05, ...){
  x<-NULL
  for(i in row.names(rho)){
    if(0){
      x<-cbind(x,matrix(withCallingHandlers(
        random(rho=as.estimate1d(rho[i,]), n=n, method=method, relativeTolerance=relativeTolerance, ...),
        warning=function(w) warning("Variable: ", i, "\t distribution: ", rho[i,"distribution"], "\n", 
                                    w$message, call. = FALSE, immediate.=TRUE),
        error=function(e) stop("Variable: ", i, "\n", e$message)
      ), nrow=n, ncol=1, dimnames=list(NULL,i)), deparse.level=1
      )
    }
    x<-withCallingHandlers(cbind(x,matrix(
      random(rho=as.estimate1d(rho[i,]), n=n, method=method, relativeTolerance=relativeTolerance, ...),
      nrow=n, ncol=1, dimnames=list(NULL,i)), 
      deparse.level=1),
      warning=function(w) {
        warning("Variable: ", i, "\t distribution: ", rho[i,"distribution"], "\n", 
                w$message, call. = FALSE, immediate.=TRUE)
      },
      #warning=function(w) warnings("Variable: ", i, "\t distribution: ", rho[i,"distribution"], "\n"),
      error=function(e) stop("Variable: ", i, "\n", e$message)
    )
  }
  #  Return the sampled multivariate values:
  x
}
