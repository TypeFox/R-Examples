
## Function overwrite parameters of an existing MaxControl object using
## parameters supplied in a single list.
## We do not make it to a method: the signature would be indistinguishable
## from add(maxControl, ...) where ... is a single list
addControlList <- function(x, y, check=TRUE) {
   ## add list y to the control
   ## 
   ## check    only accept known control options.
   ##          useful if attaching known control list
   ##          if false, no checks performed and can add arbitrary list
   ##          
   setSlot <- function(openName,
                       slotName=openName[1],
                       convert=function(x) x
                       ) {
      ## Store potentially differently named value in slot
      ## 
      ## openName    vector of accepted name forms
      ## slotName    corresponding actual slot name
      ## convert     how to convert the value
      ## 
      if(!any(openName %in% names(y))) {
         return(NULL)
      }
      i <- tail(which(names(y) %in% openName), 1)
                           # pick the last occurrence: allow user to overwrite defaults
      slot(x, slotName) <- convert(y[[i]])
      assign("x", x, envir=parent.frame())
                           # save modified x into parent frame
   }
   if(!inherits(x, "MaxControl")) {
      stop("'x' must be of class 'MaxControl'")
   }
   if(is.null(y)) {
      return(x)
   }
   if(!inherits(y, "list")) {
      stop("Control arguments to 'maxControl' must be supplied in the form of a list")
   }
   if(check) {
      knownNames <- union(openParam(x), slotNames(x))
      if(any(uNames <- !(names(y) %in% knownNames))) {
         cat("Unknown control options:\n")
         print(names(y)[uNames])
         stop("Unknown options not accepted")
      }
   }
   ##
   setSlot("tol")
   setSlot("reltol")
   setSlot("gradtol")
   setSlot("lambdatol")
   setSlot("qrtol")
   ## QAC
   setSlot(c("qac", "QAC"), "qac")
   setSlot(c("marquardt_lambda0", "Marquardt_lambda0"))
   setSlot(c("marquardt_lambdaStep", "Marquardt_lambdaStep"))
   setSlot(c("marquardt_maxLambda", "Marquardt_maxLambda"))
   ## NM
   setSlot(c("nm_alpha", "NM_alpha", "alpha"))
   setSlot(c("nm_beta", "NM_beta", "beta"))
   setSlot(c("nm_gamma", "NM_gamma", "gamma"))
   ## SANN
   setSlot(c("sann_cand", "SANN_cand", "cand"))
   setSlot(c("sann_temp", "SANN_temp", "temp"))
   setSlot(c("sann_tmax", "SANN_tmax", "tmax"), convert=as.integer)
   setSlot(c("sann_randomSeed", "SANN_randomSeed", "random.seed"),
           convert=as.integer)
   ##
   setSlot("iterlim", convert=as.integer)
   setSlot(c("printLevel", "print.level"), convert=as.integer)
   ##
   validObject(x)
   return(x)
}

### Method for 'MaxControl' objects: add the second argument, list
setMethod("maxControl", signature("MaxControl"), addControlList)
