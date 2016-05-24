#' mod identify dType
#'
#' This function identify the data type of each curve  from a \code{\linkS4class{TLum.Analysis}} object.
#' It also add the new data type "testdose" and "preheat" based on the \code{comment} present in the \code{\linkS4class{TLum.Analysis}} object or a sequence vector.
#'
#' @param object
#'  \code{\linkS4class{TLum.Analysis}} (\bold{required}): object containing the initial TL curves.
#' @param method
#'  \link{character} (with default): Defines the methode use to identify the new data type ("comment", "sequence", "temperature+dose").
#' @param ref
#'  \link{list} (with default): Contains the reference values to identify the new data type.
#'
#' @return
#'  This function provides a new \code{\linkS4class{TLum.Analysis}} with the new dtype. \cr
#'
#' @author David Strebler, University of Cologne (Germany).
#'
#' @export mod_update.dType

mod_update.dType <- function(

  object,

  method ="comment",

  ref = list(sequence=c("Preheat", "Natural", "Background", "Testdose", "Background"),
             oneByOne = FALSE,
             protocol = "SAR",
             preheat = NA,
             testdose = NA)
){
  #Reference method
  C_COMMENT <- "comment"
  C_SEQUENCE <- "sequence"
  C_TEMPERATURE_DOSE <- "temperature+dose"

  #Reference dType
  C_NATURAL <- "Natural"
  C_NATURAL.DOSE <- "Natural+Dose"
  C_DOSE <- "Dose"
  C_BLEACH <- "Bleach"
  C_BLEACH.DOSE <- "Bleach+Dose"
  C_NATURAL.BLEACH <- "Nat.(Bleach)"
  C_NATURAL.DOSE.BLEACH <- "Nat.+Dose(Bleach)"
  C_BACKGROUND <- "Background"
  C_PREHEAT <- "Preheat" #new#
  C_TESTDOSE <- "Testdose" #new# (replace only "Dose")

  #Variation dType
  V_NATURAL <- c("natural",
                 "nat.")
  V_NATURAL.DOSE <- c("natural+dose",
                      "nat+dose",
                      "nat.+dose",
                      "natural + dose",
                      "nat + dose",
                      "nat. + dose")
  V_DOSE <- c("dose")
  V_BLEACH <- c("bleach")
  V_BLEACH.DOSE <- c("bleach+dose",
                     "bleach + dose")
  V_NATURAL.BLEACH <- c("nat.(bleach)",
                        "natural(bleach)",
                        "nat. (bleach)",
                        "natural (bleach)")
  V_NATURAL.DOSE.BLEACH <- c("nat.+dose(bleach)",
                             "natural+dose(bleach)",
                             "nat. + dose (bleach)",
                             "natural + dose (bleach)")
  V_BACKGROUND <- c("background")
  V_PREHEAT <- c("preheat",
                 "pre-heat",
                 "pre heat")
  V_TESTDOSE <- c("testdose",
                  "test-dose",
                  "test dose")

  V_DTYPE <- c(V_NATURAL,
               V_NATURAL.DOSE,
               V_DOSE,
               V_BLEACH,
               V_BLEACH.DOSE,
               V_NATURAL.BLEACH,
               V_NATURAL.DOSE.BLEACH,
               V_BACKGROUND,
               V_PREHEAT,
               V_TESTDOSE)

  V_NEW <- c(V_PREHEAT,
             V_TESTDOSE)

  # ------------------------------------------------------------------------------
  # Integrity Check
  # ------------------------------------------------------------------------------
  if (missing(object)){
    stop("[mod_update.dType] Error: Input 'object' is missing.")
  }else if (!is(object,"TLum.Analysis")){
    stop("[mod_update.dType] Error: Input 'object' is not of type 'TLum.Analysis'.")
  }

  if(!is.character(method)){
    stop("[mod_update.dType] Error: Input 'method' is not of type 'character'.")
  }

  # ------------------------------------------------------------------------------

  new.protocol <- object@protocol
  records <- object@records

  nRecords <- length(records)

  new.records <- list()
  comments <- vector()

  new.records <- list()

  if(method == C_COMMENT){

    for(i in 1: nRecords){
      temp.curve <- records[[i]]

      temp.metadata <- temp.curve@metadata

      temp.dtype <- tolower(temp.metadata$DTYPE)
      temp.comment <- tolower(temp.metadata$COMMENT)

      if(temp.comment %in% V_PREHEAT){
        new.dtype <- C_PREHEAT
      }else if( (temp.comment %in% V_TESTDOSE) && (temp.dtype %in% V_DOSE)){
        new.dtype <- C_TESTDOSE
      }else{
        new.dtype <- temp.metadata$DTYPE
      }

      new.metadata <- temp.metadata
      new.metadata$DTYPE <- new.dtype

      new.curve <- temp.curve
      new.curve@metadata <- new.metadata

      new.records <- c(new.records, new.curve)

      comments <- c(comments, temp.comment)
    }

    if(!(TRUE %in% (V_NEW %in% comments))){
      warning("[mod_identify.dType] warning: Comment from the input 'objet' do not includes any new 'data type'.")
    }

  }else if(method == C_SEQUENCE){
    sequence <- ref$sequence
    OneByOne <- ref$OneByOne
    protocol <- ref$protocol

    ## Value check
    if(is.null(sequence)){
      stop("[mod_update.dType] Error: Input 'ref$sequence' is missing.")
    }else if(!is.character(sequence)){
      stop("[mod_update.dType] Error: Input 'sequence' is not of type 'character'.")
    }

    if(is.null(OneByOne)){
      stop("[mod_update.dType] Error: Input 'ref$OneByOne' is missing.")
    }else if(!is.logical(OneByOne)){
      stop("[mod_update.dType] Error: Input 'sequence' is not of type 'logical'.")
    }

    if(is.null(protocol)){
      stop("[mod_update.dType] Error: Input 'ref$protocol' is missing.")
    }else if(!is.character(sequence)){
      stop("[mod_update.dType] Error: Input 'protocol' is not of type 'character'.")
    }

    if((nRecords%%length(sequence))>0){
      stop("[mod_update.dType] Error: The number of records should be a multiple of the sequence length.")
    }

    lSequence <- length(sequence)
    nSequence <- nRecords/lSequence

    if(!OneByOne){
      for (i in 1:(nSequence-1)){
        for(j in 1:lSequence){
          k <- (i-1)*lSequence + j

          temp.curve <- records[[k]]
          temp.metadata <- temp.curve@metadata

          temp.dtype <- temp.metadata$DTYPE
          temp.comment <- temp.metadata$COMMENT
          temp.new.dtype <- tolower(sequence[j])

          if(temp.new.dtype != tolower(temp.dtype)){
            if(temp.new.dtype %in% V_PREHEAT){
              new.dtype <- C_PREHEAT
              new.comment <- new.dtype

            }else if((temp.new.dtype %in% V_TESTDOSE) && (tolower(temp.dtype) %in% V_DOSE)){
              new.dtype <- C_TESTDOSE
              new.comment <- new.dtype

            }else{
              new.dtype <- temp.dtype
              new.comment <- temp.comment
            }

          }else{
            new.dtype <- temp.dtype
            new.comment <- temp.comment
          }

          new.metadata <- temp.metadata
          new.metadata$DTYPE <- new.dtype
          new.metadata$COMMENT <- new.comment

          new.curve <- temp.curve
          new.curve@metadata <- new.metadata

          new.records <- c(new.records, new.curve)
        }
      }
    }else{
      stop("[mod_update.dType] Error: this method is not yet implemented.")
    }

  }else if(method == C_TEMPERATURE_DOSE){
    preheat <- ref$preheat
    testdose <- ref$testdose

    ## Value check
    if(is.null(preheat)){
      stop("[mod_update.dType] Error: Input 'ref$preheat' is missing.")
    }else if(!is.numeric(preheat)){
      stop("[mod_update.dType] Error: Input 'preheat' is not of type 'numeric'.")
    }

    if(is.null(testdose)){
      stop("[mod_update.dType] Error: Input 'ref$testdose' is missing.")
    }else if(!is.numeric(testdose)){
      stop("[mod_update.dType] Error: Input 'testdose' is not of type 'numeric'.")
    }


    for(i in 1: nRecords){
      temp.curve <- records[[i]]

      temp.metadata <- temp.curve@metadata

      temp.dtype <- temp.metadata$DTYPE
      temp.Tmax <- temp.metadata$HIGH
      temp.dose <- temp.metadata$IRR_TIME
      temp.comment <- temp.metadata$COMMENT

      if(temp.Tmax <= preheat){
        new.dtype <- C_PREHEAT
        new.comment <- new.dtype

      }else if(temp.dose == testdose && (tolower(temp.dtype) %in% V_DOSE)){
        new.dtype <- C_TESTDOSE
        new.comment <- new.dtype

      }else{
        new.dtype <- temp.dtype
        new.comment <- temp.comment
      }

      new.metadata <- temp.metadata
      new.metadata$DTYPE <- new.dtype
      new.metadata$COMMENT <- new.comment

      new.curve <- temp.curve
      new.curve@metadata <- new.metadata

      new.records <- c(new.records, new.curve)
    }

    warning.message <- paste("[mod_update.dType] warning:",
                             testdose, "(s) can only be used as a testdose.",
                             "Every curves with Tmax <=", preheat, "(\u00b0C) is considered as a preheat curve.")

    warning(warning.message)

  }else{

    error.message <- paste("[mod_update.dType] Error: the method",
                           method,
                           "do not exist.",
                           "The method available are:",
                           c(C_COMMENT, C_SEQUENCE, C_TEMPERATURE_DOSE),
                           ".")
    stop(error.message)
  }

  # NEW TLUM.ANALYSIS
  new.TLum.Analysis <- set_TLum.Analysis(records= new.records,
                                         protocol=new.protocol)

  return(new.TLum.Analysis)
}
