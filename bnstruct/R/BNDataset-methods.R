#' initialize a \code{\link{BNDataset}} object.
#' 
#' @name BNDataset
#' @rdname BNDataset-class
#' @docType method
#' @aliases initialize,BNDataset-method
#' 
#' @param .Object an empty BNDataset.
#' 
#' @return a BNDataset object.
setMethod("initialize",
          "BNDataset", function(.Object)  
          {
            validObject(.Object)
            return(.Object)
          })


#' @export
BNDataset <- function(data, discreteness, variables = NULL, node.sizes = NULL, ...)
{
  dataset <- new("BNDataset")
  
  # # this is here for 2 purposes:
  # # 1. spare changes all over the package in order to remove name field
  # # 2. keep a suggestion on how to get variable name
  # # name(dataset) <- deparse(substitute(dataset))
  # above: seems heavy...
  name(dataset) <- "BNDataset"
  
  # The presence of ONLY data and discreteness, and them being 2 strings, mean that two files are passed:
  # - data file (data)
  # - header file (discreteness)
  if ( is.null(variables)        &&  is.null(node.sizes)        &&
      !is.null(data)             && !is.null(discreteness)      &&
       length(discreteness) == 1 &&  is.character(discreteness) &&
       length(data) == 1         &&  is.character(data)           ) {
    
    dataset <- read.dataset(dataset, data, discreteness, ...)
    validObject(dataset)
    return(dataset)
  }
  
  other.args <- list(...)
  if ("starts.from" %in% names(other.args))
    starts.from <- other.args$starts.from
  else
    starts.from <- 1
  
  if(length(variables) > 1)
  {
    variables(dataset) <- variables
    num.variables(dataset) <- length(variables)
  }
  
  if (length(node.sizes) > 1)
    node.sizes(dataset) <- node.sizes
  
  if (length(discreteness) > 1)
    discreteness(dataset) <- discreteness
  
  if (!is.null(data))
  {
    raw.data(dataset) <- as.matrix(data) + (1 - starts.from)
    if (is.null(variables)) {
      variables(dataset) <- rownames(data)
      warning("Variable names guessed from data. Please check for consistency with your actual data.")
    }
    if (is.null(node.sizes)) {
      node.sizes <- rep(0, length(variables))
      for (v in 1:length(variables))
      {
        node.sizes[v] <- max(data[,v][which(!is.na(data[,v]))]) - min(data[,v][which(!is.na(data[,v]))]) + 1
      }
      warning("Variable cardinalities guessed from data. Please check for consistency with your actual data. Otherwise, execution may terminate with errors later.")
    }
  }
  
  num.items(dataset) <- nrow(dataset@raw.data)
  
  validObject(dataset)
  
  if(length(dataset@variables) > 0 && has.raw.data(dataset))
    colnames(dataset@raw.data) <- dataset@variables
  
  return(dataset)
}



# validator
setValidity("BNDataset",
            function(object)
            {
              retval <- NULL
              if (object@num.variables > 0 && length(object@variables) > 0 && length(object@variables) != object@num.variables)
              {
                retval <- c(retval, "incoherent number of variable names")
              }
              if (object@has.raw.data && ncol(object@raw.data) != object@num.variables)
              {
                retval <- c(retval, "incoherent number of variables in raw dataset")
              }
              if (object@has.imputed.data && ncol(object@imputed.data) != object@num.variables)
              {
                retval <- c(retval, "incoherent number of variables in imputed dataset")
              }
              if(object@num.variables > 0 && length(object@discreteness) > 1 &&
                   length(object@discreteness) != object@num.variables)
              {
                retval <- c(retval, "incoherent number of variable statuses")
              }
                            
              if (object@num.variables > 0 && length(object@node.sizes) == object@num.variables && object@has.raw.data)
              {
                warn <- c()
                halt <- c()
                for (var in 1:object@num.variables)
                {
                  if ( object@discreteness[var] &&
                      (min(object@raw.data[,var][which(!is.na(object@raw.data[,var]))]) > 1 ||
                       max(object@raw.data[,var][which(!is.na(object@raw.data[,var]))]) < object@node.sizes[var]))
                  {
                    warn <- c(warn, var)
                  }
                  if ( object@discreteness[var] &&
                      (min(object@raw.data[,var][which(!is.na(object@raw.data[,var]))]) < 1 ||
                       max(object@raw.data[,var][which(!is.na(object@raw.data[,var]))]) > object@node.sizes[var]))
                  {
                    halt <- c(halt, var)
                  }
                }
                if (length(halt) > 0)
                {
                  wrongs <- strcat("Dataset contains values out of bounds for variables ", halt, sep=" ")
                  retval <- c(retval, wrongs)
                } else if (length(warn) > 0)
                {
                  wrongs <- strcat("Not all of the possible values have been observed for variables ", warn, sep = " ")
                  warning(wrongs)
                }
              }
              
              if (object@num.variables > 0 && length(object@node.sizes) == object@num.variables && object@has.imputed.data)
              {
                warn <- c()
                halt <- c()
                for (var in 1:object@num.variables)
                {
                  if ( object@discreteness[var] &&
                      (min(object@imputed.data[,var][which(!is.na(object@imputed.data[,var]))]) > 1 ||
                       max(object@imputed.data[,var][which(!is.na(object@imputed.data[,var]))]) < object@node.sizes[var]))
                  {
                    warn <- c(warn, var)
                  }
                  if ( object@discreteness[var] &&
                      (min(object@imputed.data[,var][which(!is.na(object@imputed.data[,var]))]) < 1 ||
                       max(object@imputed.data[,var][which(!is.na(object@imputed.data[,var]))]) > object@node.sizes[var]))
                  {
                    halt <- c(halt, var)
                  }
                }
                if (length(halt) > 0)
                {
                  wrongs <- strcat("Dataset contains values out of bounds for variables ", halt, sep=" ")
                  retval <- c(retval, wrongs)
                } else if (length(warn) > 0)
                {
                  wrongs <- strcat("Not all of the possible values have been observed for variables ", warn, sep= " ")
                  warning(wrongs)
                }
              }
              
              if (is.null(retval)) return (TRUE)
              return(retval)
            }
)

#' @rdname name
#' @aliases name,BNDataset
setMethod("name", "BNDataset", function(x) { return(slot(x, "name")) } )

#' @rdname num.variables
#' @aliases num.variables,BNDataset
setMethod("num.variables", "BNDataset", function(x) { return(slot(x, "num.variables")) } )

#' @rdname variables
#' @aliases variables,BNDataset
setMethod("variables", "BNDataset", function(x) { return(slot(x, "variables")) } )

#' @rdname discreteness
#' @aliases discreteness,BNDataset
setMethod("discreteness",
          "BNDataset",
          function(x)
          {
            return(slot(x, "discreteness"))
          })

#' @rdname node.sizes
#' @aliases node.sizes,BNDataset
setMethod("node.sizes", "BNDataset", function(x) { return(slot(x, "node.sizes")) } )


#' @rdname header.file
#' @aliases header.file,BNDataset
setMethod("header.file", "BNDataset", function(x) return(slot(x, "header.file")))


#' @rdname data.file
#' @aliases data.file,BNDataset
setMethod("data.file", "BNDataset", function(x) return(slot(x, "data.file")))

#' @rdname num.variables
#' @aliases num.variables,BNDataset
setMethod("num.variables","BNDataset", function(x) return(slot(x, "num.variables")))

#' @rdname num.items
#' @aliases num.items,BNDataset
setMethod("num.items", "BNDataset", function(x) return(slot(x, "num.items")))

#' @rdname has.boots
#' @aliases has.boots,BNDataset
setMethod("has.boots", "BNDataset", function(x) return(slot(x, "has.boots")))

#' @rdname has.imputed.boots
#' @aliases has.imputed.boots,BNDataset
setMethod("has.imputed.boots", "BNDataset", function(x) return(slot(x, "has.imputed.boots")))

#' @rdname boots
#' @aliases boots,BNDataset
setMethod("boots", "BNDataset", function(x) return(slot(x, "boots")))

#' @rdname imp.boots
#' @aliases imp.boots,BNDataset
setMethod("imp.boots", "BNDataset", function(x) return(slot(x, "imp.boots")))

#' @rdname num.boots
#' @aliases num.boots,BNDataset
setMethod("num.boots", "BNDataset", function(x) return(slot(x, "num.boots")))


#' @name name<-
#' @aliases name<-,BNDataset-method
#' @docType methods
#' @rdname name-set
setReplaceMethod("name",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "name") <- value
                   validObject(x)
                   return(x)
                 })


#' @name variables<-
#' @aliases variables<-,BNDataset-method
#' @docType methods
#' @rdname variables-set
setReplaceMethod("variables",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "variables")  <- value
                   num.variables(x)      <- length(value)
                   validObject(x)
                   return(x)
                 })


#' @name discreteness<-
#' @aliases discreteness<-,BNDataset-method
#' @docType methods
#' @rdname discreteness-set
setReplaceMethod("discreteness",
                 "BNDataset",
                 function(x, value)
                 {
                   if (is.logical(value))
                     slot(x, "discreteness") <- value
                   if (is.integer(value) || is.numeric(value))
                   {
                     d <- rep(F, num.variables(x))
                     d[value] <- T
                     slot(x, "discreteness") <- value
                   }
                   if (is.character(value))
                     slot(x, "discreteness") <- sapply(1:length(value), FUN=function(i){ !is.na(match(value[i],c('d',"D"))) })
                   validObject(x)
                   return(x)
                 })


#' @name node.sizes<-
#' @aliases node.sizes<-,BNDataset-method
#' @docType methods
#' @rdname node.sizes-set
setReplaceMethod("node.sizes",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "node.sizes") <- value
                   validObject(x)
                   return(x)
                 })


#' @rdname has.raw.data
#' @aliases has.raw.data,BNDataset
setMethod("has.raw.data",
          "BNDataset",
          function(x)
          {
            return(slot(x, "has.raw.data"))
          })


#' @rdname has.imputed.data
#' @aliases has.imputed.data,BNDataset
setMethod("has.imputed.data",
          "BNDataset",
          function(x)
          {
            return(slot(x, "has.imputed.data"))
          })


#' @rdname raw.data
#' @aliases raw.data,BNDataset
setMethod("raw.data",
          "BNDataset",
          function(x)
          {
            if (has.raw.data(x))
               return (x@raw.data)
            stop("The dataset contains no data.")
          })


#' @rdname imputed.data
#' @aliases imputed.data,BNDataset
setMethod("imputed.data",
          "BNDataset",
          function(x)
          {
            if (has.imputed.data(x))
              return (x@imputed.data)
            stop("The dataset contains no imputed data. ",
                 "Please impute data before learning.\nSee > ?impute for help.")
          })


#' @name header.file<-
#' @aliases header.file<-,BNDataset-method
#' @docType methods
#' @rdname header.file-set
setReplaceMethod("header.file",
                 "BNDataset",
                 function(x, value)
                  {
                    slot(x, "header.file") <- value
                    return(x)
                  })


#' @name data.file<-
#' @aliases data.file<-,BNDataset-method
#' @docType methods
#' @rdname data.file-set
setReplaceMethod("data.file",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "data.file") <- value
                   return(x)
                 })


#' @name num.variables<-
#' @aliases num.variables<-,BNDataset-method
#' @docType methods
#' @rdname num.variables-set
setReplaceMethod("num.variables",
                 "BNDataset",
                 function(x, value)
                  {
                    slot(x, "num.variables") <- value
                    validObject(x)
                    return(x)
                  })


#' @name num.items<-
#' @aliases num.items<-,BNDataset-method
#' @docType methods
#' @rdname num.items-set
setReplaceMethod("num.items",
                 "BNDataset",
                 function(x, value)
                 {
                  slot(x, "num.items") <- value
                  validObject(x)
                  return(x)
                 })


#' @name boots<-
#' @aliases boots<-,BNDataset-method
#' @docType methods
#' @rdname boots-set
setReplaceMethod("boots",
                 "BNDataset",
                 function(x, value)
                 {
                  slot(x, "boots")     <- value
                  slot(x, "num.boots") <- length(value)
                  slot(x, "has.boots") <- TRUE
                  validObject(x)
                  return(x)
                 })


#' @name num.boots<-
#' @aliases num.boots<-,BNDataset-method
#' @docType methods
#' @rdname num.boots-set
setReplaceMethod("num.boots",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "num.boots") <- value
                   validObject(x)
                   return(x)
                 })


#' @name imp.boots<-
#' @aliases imp.boots<-,BNDataset-method
#' @docType methods
#' @rdname imp.boots-set
setReplaceMethod("imp.boots",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "imp.boots")     <- value
                   slot(x, "num.boots")     <- length(value)
                   slot(x, "has.imputed.boots") <- TRUE
                   validObject(x)
                   return(x)
                 })


#' @name raw.data<-
#' @aliases raw.data<-,BNDataset-method
#' @docType methods
#' @rdname raw.data-set
setReplaceMethod("raw.data",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "raw.data")    <- value
                   slot(x, "has.raw.data") <- TRUE
                   num.items(x) <- nrow(value)
                   validObject(x)
                   return(x)
                 })


#' @name imputed.data<-
#' @aliases imputed.data<-,BNDataset-method
#' @docType methods
#' @rdname imputed.data-set
setReplaceMethod("imputed.data",
                 "BNDataset",
                 function(x, value)
                 {
                   slot(x, "imputed.data") <- value
                   slot(x, "has.imputed.data")  <- TRUE
                   num.items(x) <- nrow(value)
                   validObject(x)
                   return(x)
                 })

#' @rdname complete
#' @aliases complete,BNDataset
setMethod("complete",
          "BNDataset",
          function(x, complete.vars=seq_len(num.variables(x)))
          {
            y <- x
            rd <- raw.data(y)
            raw.data(y)  <- rd[complete.cases(rd[,complete.vars]),]
            num.items(y) <- nrow(raw.data(y))
            slot(y, "imputed.data")      <- matrix(c(0))
            slot(y, "has.imputed.data")  <- FALSE
            slot(y, "boots")             <- list(NULL)
            slot(y, "has.boots")         <- FALSE
            slot(y, "imp.boots")         <- list(NULL)
            slot(y, "has.imputed.boots") <- FALSE
            slot(y, "num.boots")         <- 0
            validObject(y)
            return(y)
          })


# redefinition of print() for BNDataset objects
#' @rdname print
#' @aliases print,BNDataset
setMethod("print",
          "BNDataset",
          function(x, show.raw.data = FALSE, show.imputed.data = FALSE, ...)
          {
            
            str <- "\nDataset: \n"
            #str <- paste(str, name(x), sep = '')
            #str <- paste(str, "\n", sep = '')
            cat(str)
            str <- "\nnum.variables "
            str <- paste(str, num.variables(x), sep = '')
            str <- paste(str, "\n", sep = '')
            cat(str)
            str <- "\nvariables\n"
            cat(str)
            cat(variables(x))
            str <- "\ndiscreteness\n"
            cat(str)
            cat(discreteness(x))
            str <- "\nnode.sizes\n"
            cat(str)
            cat(node.sizes(x))
            str <- "\nnum.items\n"
            cat(str)
            cat(num.items(x))
            str <- "\nimputation\n"
            cat(str)
            cat(has.imputed.data(x))
            str <- "\nhas.boots\n"
            cat(str)
            cat(has.boots(x))
            str <- "\nhas.imputed.boots\n"
            cat(str)
            cat(has.imputed.boots(x))
            str <- "\nnum.boots\n"
            cat(str)
            cat(num.boots(x))
            
            
            if (show.raw.data == TRUE && has.raw.data(x))
            {
              cat("\nRaw data:\n")
              print(raw.data(x))
            }
            if (show.imputed.data == TRUE && has.imputed.data(x))
            {
              cat("\nImputed data:\n")
              print(imputed.data(x))
            }
            
            cat("\n")
          })

#' @rdname impute
#' @aliases impute,BNDataset
setMethod("impute",
          "BNDataset",
          function(object, k.impute = 10)
          {
            # assumes raw data is ok
            bnstruct.start.log("performing imputation ...")
            object@imputed.data <- knn.impute(object@raw.data, k.impute,
                                              setdiff(1:length(object@node.sizes), c()))
            object@has.imputed.data  <- TRUE
            bnstruct.end.log("imputation finished.")
            return(object)
          })

#' @rdname bootstrap
#' @aliases bootstrap,BNDataset
setMethod("bootstrap",
          "BNDataset",
          function(object, num.boots = 100, seed = 0, imputation = FALSE, k.impute = 10)
          {
            if (imputation)
              bnstruct.start.log("Generating bootstrap samples with imputation ...")
            else
              bnstruct.start.log("Generating bootstrap samples ...")
            
            # assumes raw data is ok
            object@has.boots <- TRUE
            object@num.boots <- num.boots
            
            set.seed(seed)
            if (num.boots >= 1)
            {
              boot.sample <- matrix(sample.int(object@num.items,
                                               size = num.boots * object@num.items,
                                               replace=TRUE),
                                    object@num.items, num.boots)
              
              if (imputation)
                object@has.imputed.boots <- TRUE
              
              for (i in 1:num.boots)
              {
                object@boots[[i]] <- object@raw.data[boot.sample[,i],]
                
                if (imputation)
                  object@imp.boots[[i]] <- knn.impute(object@boots[[i]],
                                                      k.impute,
                                                      setdiff(1:length(object@node.sizes),c()) )
                
              }
            }
            
            bnstruct.end.log("Bootstrap samples generated.")
            return(object)
          })


#' @rdname boot
#' @aliases boot,BNDataset
setMethod("boot",
          c("BNDataset", "numeric"),
          function(dataset, index, use.imputed.data = FALSE)
          {
            if (!use.imputed.data && !dataset@has.boots)
              stop('No bootstrap samples available for dataset.')
            
            if (use.imputed.data && !dataset@has.imputed.boots)
              stop('No imputed bootstrap samples available for dataset. ',
                   "Please impute data before learning.\nSee > ?impute for help.")
            
            if (index <= 0 || index > dataset@num.boots)
              stop('Sample index out of range for dataset.\n')
            
            if (use.imputed.data)
              return(dataset@imp.boots[[index]])
            
            return(dataset@boots[[index]])
          })

