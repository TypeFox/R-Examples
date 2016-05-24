#############################################################################
#
#   This file is part of the R package "RSNNS".
#
#   Author: Christoph Bergmeir
#   Supervisor: José M. Benítez
#   Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Library General Public
#   License as published by the Free Software Foundation; either
#   version 2 of the License, or (at your option) any later version.
# 
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Library General Public License for more details.
# 
#   You should have received a copy of the GNU Library General Public License
#   along with this library; see the file COPYING.LIB.  If not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#   Boston, MA 02110-1301, USA.
#
#############################################################################


#' An S4 class that is the main class of RSNNS. Each instance of this class 
#' contains a pointer to a C++ object of type SnnsCLib, i.e. an instance 
#' of the SNNS kernel.
#' 
#' The only slot \code{variables} holds an environment with all member variables. 
#' Currently, there are two members (constructed by the object factory):
#' \describe{
#' \item{snnsCLibPointer}{A pointer to the corresponding C++ object}
#' \item{serialization}{a serialization of the C++ object, in SNNS .net format}
#' }
#' The member variables are not directly present as slots but wrapped in an environment
#' to allow for changing the serialization (by call by reference). 
#' 
#' An object of this class is used internally by all the models in the package. 
#' The object is always accessible by \code{model$snnsObject$...}
#' 
#' To make full use of the SNNS functionalities, you might want to use this class directly.
#' Always use the object factory \code{\link{SnnsRObjectFactory}} to construct an object, 
#' and the calling mechanism \code{\link{$}} to call functions. Through the calling mechanism,
#' many functions of SnnsCLib are present that are not documented here, but in the SNNS User 
#' Manual. So, if you choose to use the low-level interface, it is highly recommended to have
#' a look at the demos and at the SNNS User Manual. 
#' 
#' 
#' @title The main class of the package
# @slot variables holds an environment with all member variables of the class 
#' @references 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen. 
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#' 
#' @name SnnsR-class
#' @seealso \code{\link{$}}, \code{\link{SnnsRObjectFactory}}
#' @examples
#' \dontrun{demo(encoderSnnsCLib)} 
#' \dontrun{demo(art1_lettersSnnsR)}
#' \dontrun{demo(art2_tetraSnnsR)} 
#' \dontrun{demo(artmap_lettersSnnsR)} 
#' \dontrun{demo(eight_elmanSnnsR)}
#' \dontrun{demo(rbf_irisSnnsR)}
#' \dontrun{demo(rbf_sinSnnsR)}
#' \dontrun{demo(rbfDDA_spiralsSnnsR)}
#' \dontrun{demo(som_cubeSnnsR)}
#' 
#' 
#' #This is the demo eight_elmanSnnsR
#' #Here, we train an Elman network
#' #and save a trained and an untrained version
#' #to disk, as well as the used training data
#' 
#' basePath <- ("./")
#' 
#' data(snnsData)
#' 
#' inputs <- snnsData$eight_016.pat[,inputColumns(snnsData$eight_016.pat)]
#' outputs <- snnsData$eight_016.pat[,outputColumns(snnsData$eight_016.pat)]
#' 
#' snnsObject <- SnnsRObjectFactory()
#' 
#' snnsObject$setLearnFunc('JE_BP')
#' snnsObject$setUpdateFunc('JE_Order')
#' snnsObject$setUnitDefaults(1,0,1,0,1,'Act_Logistic','Out_Identity')
#' 
#' snnsObject$elman_createNet(c(2,8,2),c(1,1,1),FALSE)
#' 
#' 
#' patset <- snnsObject$createPatSet(inputs, outputs)
#' snnsObject$setCurrPatSet(patset$set_no)
#' 
#' snnsObject$initializeNet(c(1.0,  -1.0,  0.3,  1.0,  0.5) )
#' snnsObject$shufflePatterns(TRUE)
#' snnsObject$DefTrainSubPat()
#' 
#' snnsObject$saveNet(paste(basePath,"eight_elmanSnnsR_untrained.net",sep=""),
#'                                           "eight_elmanSnnsR_untrained")
#' 
#' parameters <- c(0.2, 0, 0, 0, 0)
#' maxit <- 1000
#' 
#' error <- vector()
#' for(i in 1:maxit) {
#'   res <- snnsObject$learnAllPatterns(parameters)
#'   if(res[[1]] != 0) print(paste("Error at iteration ", i, " : ", res, sep=""))
#'   error[i] <- res[[2]]
#' }
#' 
#' error[1:500]
#' plot(error, type="l")
#' 
#' snnsObject$saveNet(paste(basePath,"eight_elmanSnnsR.net",sep=""),
#'                                              "eight_elmanSnnsR")
#' snnsObject$saveNewPatterns(paste(basePath,"eight_elmanSnnsR.pat",sep=""), 
#'                                                          patset$set_no)
setClass( "SnnsR", representation( variables="environment" ))

#snnsCLibPointer = "externalptr",
#serialization="character")

#' Enable calling of C++ functions as methods of \code{SnnsR-class} objects.
#'
#' This function makes methods of SnnsR__ and SnnsCLib__ accessible via "$". If
#' no SnnsR__ method is present, then the according SnnsCLib__ method is
#' called. This enables a very flexible method handling. To mask a method from
#' SnnsCLib, e.g. to do some parameter checking or postprocessing, only a method
#' with the same name, but beginning with SnnsR__ has to be present in R.  See
#' e.g. \code{\link{SnnsRObject$initializeNet}} for such an implementation.
#' 
#' Error handling is also done within the method caller. If the result of a
#' function is a list with a member \code{err},  then \code{SnnsCLib__error} is
#' called to use the SNNS kernel function to get the corresponding error message
#' code and an R warning is thrown containing this message.
#' 
#' Furthermore, a serialization mechanism is implemented which all models
#' present in the package use to be able to be saved and loaded by R's normal
#' save/load mechanism (as RData files).
#' 
#' The completely trained object can be serialized with
#' 
#' \code{s <- snnsObject$serializeNet("RSNNS_untitled")}
#' 
#' \code{snnsObject@@variables$serialization <- s$serialization}
#' 
#' For the models implemented, this is done in \code{\link{SnnsRObject$train}}. If the S4 object is then saved and loaded, 
#' the calling mechanism will notice on the next use of a function that the pointer to the C++ SnnsCLib object is \code{nil}, 
#' and if a serialization is present, the object is restored from this serialization before the method is called.
#'
# @export
#' @title Method caller for SnnsR objects 
#' @rdname SnnsRObjectMethodCaller
#' @name SnnsRObjectMethodCaller
#' @param x object of class \link{SnnsR-class}
#' @param name function to call
#' @usage \S4method{$}{SnnsR}(x, name) 
#' @aliases $,SnnsR-method $
setMethod( "$", "SnnsR", function(x, name ){
      function(...) {
        #print(x)
        
        if(is.nil(x@variables$snnsCLibPointer)) {
          if( x@variables$serialization[1] != "") {
            
            eval.parent({x@variables$snnsCLibPointer <- .Call("SnnsCLib__new", package="RSNNS")}, 2 )
            SnnsR__deserialize(x, x@variables$serialization)
          } else {
            warning("The internal SnnsCLib object is not present, nor its serialization. Exiting..")
            return()
          }
        }
        
        myFunc <- mget(paste( "SnnsR", name, sep = "__" ), mode="any", 
            envir = as.environment(-1), 
            ifnotfound = list(FALSE), inherits=TRUE)
        
        #very useful for debugging. Everytime an SnnsR or SnnsCLib function is
        #called, its name is printed
        
        #print(name)
        
        if(is.function(myFunc[[1]])) {
          res <- myFunc[[1]](x, ... )
        }
        else {
          myFuncName <- paste( "SnnsCLib", name, sep = "__" )
          res <- .Call( myFuncName, x@variables$snnsCLibPointer , ... )
        }
        
        if(is.list(res))
          if(!is.null(res$err)) {
            err <- res$err
            if(err != 0) {
              msg <- .Call( "SnnsCLib__error", x@variables$snnsCLibPointer , err )
              warning(paste("SNNS error message in ", name, " : ", msg, sep=""))
            }
          }
        return(res)
      }     
    } )


#' Object factory to create a new object of type \code{\link{SnnsR-class}}.
#'
#' This function creates a new object of class \code{\link{SnnsR-class}}, initializes its only slot \code{variables}
#' with a new environment, generates a new C++ object of class SnnsCLib, and an empty object serialization. 
#' 
#' @title SnnsR object factory
#' @export
#' @seealso \code{\link{$}}, \code{\link{SnnsR-class}}
#' @examples 
#' mySnnsObject <- SnnsRObjectFactory()
#' mySnnsObject$setLearnFunc('Quickprop')
#' mySnnsObject$setUpdateFunc('Topological_Order') 
SnnsRObjectFactory <- function(){

  snnsObject <- new( "SnnsR")
  
  snnsObject@variables <- new.env()
  
  snnsObject@variables$snnsCLibPointer <- .Call("SnnsCLib__new", package="RSNNS")
  snnsObject@variables$serialization <- ""
  
  snnsObject
}

