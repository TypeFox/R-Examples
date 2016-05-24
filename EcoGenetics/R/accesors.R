#' Generic accessors for EcoGenetics objects
#' 
#' @name EcoGenetics accessors
#' 
#' @description These functions provide basic access to the slots of any object
#' of the formal classes defined in EcoGenetics. In addition, 
#' for ecogen objects a replace method is defined for the accessor, which enables
#' basic processing of individual slots. 
#' 
#' @param X Any S4 object of the formal classes defined in EcoGenetics.
#'
#' @format <\strong{ecoslot.}> + <\strong{name of the slot}> + <\strong{(name of the object)}>
#' 
#' @details The accessor notation in EcoGenetics consists of the prefix 
#' "ecoslot." followed by the name of the slot of interest. E.g., 
#' to access the slot "IN" of the object "X", type ecoslot.IN(X).
#' 
#' For example, the class eco.correlog of the function \code{\link{eco.correlog}}
#' has the slots OUT, IN, BREAKS CARDINAL, etc. 
#' An object X of class eco.correlog, generated with this function, 
#' has access to these slots using:
#' ecoslot.OUT(X), ecoslot.IN(X), ecoslot.BREAKS(X) and ecoslot.CARDINAL(X), 
#' respectively. 
#' 
#' 
#'    -----------------------------------------------------------------------
#' 
#' \strong{FOR ECOGEN CLASS ONLY}:
#' 
#' - The accessors in objects of class \code{\link{ecogen}} 
#' have a double usage. First, the extraction of data included in a slot (\strong{get mode}).  
#' Second, the assignation of data (\strong{set mode}). The data assigned 
#' to ecogen objects using this method is properly pre-processed. 
#' 
#' - The get mode is defined for ecoslot.XY, ecoslot.P, ecoslot.G, 
#' ecoslot.A, ecoslot.E, ecoslot.S and ecoslot.OUT. For any ecogen object X, 
#' type ecoslot.SLOT(X), where SLOT is the slot of interest: ecoslot.XY(X), 
#' ecoslot.P(X), ecoslot.G(X), ecoslot.A(X), ecoslot.E(X), ecoslot.S(X) 
#' and ecoslot.OUT(X, ...). In the latter, the three dots (...) are
#' objects in the slot OUT. 
#' 
#' - The set mode is a replacement method, i.e.,  
#' the assignation ecoslot.SLOT(X) <- VALUE is defined for the ecogen class. 
#' 
#' For a generic ecogen object "eco", the defined replacement methods are:
#' 
#' 1. ecoslot.P(eco) <- P data frame
#' 
#' 2. ecoslot.G(eco, ...) <- G data frame
#' 
#' #' In this case the three dots (...) consist of the following 
#' variables that can be passed to the function:
#' 
#' - G.processed
#' 
#' - order.G
#' 
#' - type
#' 
#' - ploidy
#' 
#' - sep
#' 
#' - ncod
#' 
#' - missing
#' 
#' - NA.char
#' 
#' - poly.level
#' 
#' - rm.empty.ind
#' 
#' See the arguments of \code{\link{ecogen}} for details about these values, and the
#' Examples section for details about usage.
#' 
#' 
#' \strong{IMPORTANT}: The assignation of data in the slot G, creates automatically the slot A.
#' An accessor is defined for the slot A only in get mode (to get the data frame in this slot, 
#' but not for replacement purposes). The slot A cannot be replaced with accessors 
#' and is generated when a genetic data frame is assigned to the slot G. 
#' 
#' 3. ecoslot.E(eco) <- E data frame
#' 
#' 4. ecoslot.S(eco) <- S data frame
#' 
#' 5. ecoslot.C(eco) <- C data frame
#' 
#' 6. ecoslot.OUT(eco, value) <- results to store in the slot OUT.
#' Here, value means any object(s). 
#' Several objects can be passed as a list. 
#' See the section Examples.
#' 
#' @rdname EcoGenetics-accessors
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar} 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' #--------------
#' # GENERAL USE
#' #--------------
#' 
#' # Example 1
#' 
#' data(eco.test)
#' 
#' ## Test with phenotypic traits (example of ?eco.correlog)
#' moran <- eco.correlog(Z=eco[["P"]][,1], XY = eco[["XY"]], method = "I", smax=10, size=1000)
#  
#' # the slots are accesed with the generic format:
#' # ecoslot. + name of the slot + (name of the object)
#' 
#' ecoslot.OUT(moran)      # slot OUT
#' ecoslot.BREAKS(moran)   # slot BREAKS
#' 
#' 
#' #----------------------------------------------
#' # SPECIFIC USE OF ACCESSORS WITH ECOGEN OBJECTS
#' #----------------------------------------------
#' 
#' #1) GET MODE
#' 
#' # Example 2
#' 
#' # Example with G data of class "data.frame", corresponding to
#' # microsatellites of a diploid organism.
#' 
#' eco <- ecogen(XY = coordinates, P = phenotype, G = genotype,
#'               E = environment, S = structure, order.G = TRUE)
#' eco              
#' 
#' 
#' # Access to the slots   
#'
#' ecoslot.XY(eco) 
#' ecoslot.P(eco)
#' ecoslot.G(eco)
#' ecoslot.A(eco)
#' ecoslot.E(eco)
#' ecoslot.S(eco)
#' ecoslot.C(eco)
#' ecoslot.OUT(eco)         
#' 
#' # For ecogen objects, the double square brackets ("[[")
#' # are symbolic abbreviations of the accessors:
#' 
#' ecoslot.XY(eco) 
#' # is identical to:
#' eco[["XY"]]
#' 
#' #2) SET MODE (REPLACEMENT OF SLOTS)
#' 
#' # Example 3
#' 
#' eco <- ecogen(XY = coordinates, P = phenotype)
#' eco
#' 
#' ecoslot.G(eco, order.G = TRUE) <- genotype
#' 
#' # this is identical to
#' eco[["G", order.G=TRUE]] <- genotype
#' 
#' ecoslot.E(eco) <- environment
#' ecoslot.S(eco) <- structure
#' 
#' # Storing the data of Example 1 in the slot OUT
#' 
#' ecoslot.OUT(eco) <- moran
#'  
#' # Storing several data
#' singers <- c("carlos_gardel", "billie_holiday")
#' golden.number <- (sqrt(5) + 1) / 2
#' ecoslot.OUT(eco) <- list(singers, golden.number)    # several objects must be passed as a list
#'
#' # In get mode, ecoslot.OUT has a double behavior:
#' # when only the name of the ecogen object is passed to
#' # the accessor, it has an overview method, 
#' # returning a data frame with the objects
#' # stored and their classes
#' 
#' ecoslot.OUT(eco)
#' 
#' # ecoslot.OUT in get mode, has two arguments:
#' # ecoslot.OUT(object, ...)
#' # here, the three dots (...) means any object(s) stored in the slot OUT.
#' 
#' ecoslot.OUT(eco, "moran", "singers")
#' 
#' # In double square brackets notation, this is equivalent to
#' eco[["OUT"]][c("moran", "singers")]
#' # This is, it works as a list and have no overview method
#' eco[["OUT"]]
#' eco[["OUT"]]["moran"]
#' 
#' # ecoslot.OUT in set mode, has two arguments:
#' # ecoslot.OUT(object, value)
#' # here value means object(s) to store in the slot OUT. Several objects
#' # must be passed as a list.
#' # The names of the input data is recoded in case of
#' # duplicates
#' 
#' ecoslot.OUT(eco) <- singers
#' ecoslot.OUT(eco)
#' ecoslot.OUT(eco) <- list(singers, singers, singers)
#' ecoslot.OUT(eco)
#' 
#' # The set operation is equivalent with double square brackets
#' eco[["OUT"]] <- list(singers, c(1,2,3))
#' ecoslot.OUT(eco)
#' 
#' }
#' 
#' 
#' @export 


ecoslot.IN <- function(X) {
  X@IN
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.BREAK <- function(X) {
  X@BREAK
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.CARDINAL <- function(X) {
  X@CARDINAL
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.NAMES <- function(X) {
  X@NAMES
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.METHOD <- function(X) {
  X@METHOD
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.DISTMETHOD <- function(X) {
  X@DISTMETHOD
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.TEST <- function(X) {
  X@TEST
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.NSIM <- function(X) {
  X@NSIM
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.PADJUST <- function(X) {
  X@PADJUST
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.OBS<- function(X) {
  X@OBS
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.EXP <- function(X) {
  X@EXP
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.PVAL <- function(X) {
  X@PVAL
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.ALTER <- function(X) {
  X@ALTER
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.MULTI <- function(X) {
  X@MULTI
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.XY <- function(X) {
  X@XY
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.COND <- function(X) {
  X@COND
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.PAR <- function(X) {
  X@PAR
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.PAR.VAL <- function(X) {
  X@PAR.VAL
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.ROW.SD <- function(X) {
  X@ROW.SD
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.SELF <- function(X) {
  X@SELF
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.SP <- function(X) {
  X@SP
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.W <- function(X) {
  X@W
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.NONZERO <- function(X) {
  X@NONZERO
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.NONZEROIND <- function(X) {
  X@NONZEROIND
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.AVERAGE <- function(X) {
  X@AVERAGE
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.MEAN <- function(X) {
  X@MEAN
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.LOGMEAN <- function(X) {
  X@LOGMEAN
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.CARDINAL <- function(X) {
  X@CARDINAL
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.BREAKS <- function(X) {
  X@BREAKS
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.MLM <- function(X) {
  X@MLM
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.SUMMARY.MLM <- function(X) {
  X@SUMMARY.MLM
}


#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.ANOVA.MLM <- function(X) {
  X@ANOVA.MLM
}


#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.DF1 <- function(X) {
  X@DF1
}


#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.DF2 <- function(X) {
  X@DF2
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.TREES <- function(X) {
  X@TREES
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.PREDICTIONS <- function(X) {
  X@PREDICTIONS
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.FREQUENCIES <- function(X) {
  X@FREQUENCIES
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.POLY.DEG <- function(X) {
  X@POLY.DEG
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.RES <- function(X) {
  X@RES
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.MODEL <- function(X) {
  X@MODEL
}

#' @rdname EcoGenetics-accessors
#' @export 

ecoslot.ANALYSIS <- function(X) {
  X@ANALYSIS
}
