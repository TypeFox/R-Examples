######################################################################
## CODE AND DOCUMENTATION FOR PACKAGE sdcTarget
######################################################################

#' @name sdcTarget-package
#' @docType package
#' @import methods magic foreach parallel doParallel tuple
#' @description
#'   Provides functions to calculate target matrices for
#'   statistical disclosure control.
#' @details
#'   Parallel processing is now implemented by default
#'   to speed up large searches. The number of processors
#'   must be specified. If not specified, will default
#'   to linear processing.
#' @references
#'   \href{http://www.statslife.org.uk/events/annual-conference/conference-blog/1752-data-linkage-and-privacy}{Summary of talk at RSS 2014 Conference}
#'
NULL

#' @title
#'   Hashing Definition (S4 Class)
#' @description
#'   This class defines the data hash that is used to identify
#'   cells to be targeted for statistical disclosure control.
#' @details
#'   The hashing definition presently handles only categorical
#'   fields, which limits the applicability of this software
#'   accordingly.
#' @slot parts
#'   The number of hashing components.
#' @slot coverage
#'   The number of elements in each hashing component.
#' @slot lengths
#'   The number of levels of each field to be hashed.
#' @slot na.exists
#'   Whether there is an NA value in each field to be hashed.
#' @slot na.recode
#'   Whether NA values should be treated as levels.
#' @slot levels
#'   The number of levels of each field to be hashed.
#' @examples
#' new("sdcHashingDefinitionClass")
#' @export
setClass("sdcHashingDefinitionClass",
    slots = c(parts = "integer",
              coverage = "integer",
              lengths = "integer",
              na.exists = "logical",
              na.recode = "logical",
              levels = "list"),
    prototype = list(parts = integer(0),
                     coverage = integer(0),
                     lengths = integer(0),
                     na.exists = logical(0),
                     na.recode = TRUE,
                     levels = vector(mode = "list", length = 0) ),
    validity = function(object) {
      identical(length(coverage), parts) &
      identical(length(lengths), length(levels)) &
      identical(length(lengths), length(na.exists)) &
      identical(length(na.recode), length(1)) &
      identical(names(levels), names(lengths)) &
      identical(names(levels), names(na.exists))
    } )

#' @describeIn sdcHashingDefinitionClass
#' @param .Object
#'   A \code{\linkS4class{sdcHashingDefinitionClass}} object.
#' @param ...
#'   The optional parameters specifying the data
#'   and whether \code{NA} values are to be recoded and hashed.
#'   Recoding is \code{TRUE} by default.
#' @examples
#' new("sdcHashingDefinitionClass")
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' new("sdcHashingDefinitionClass", X = my.X)
#' new("sdcHashingDefinitionClass", X = my.X, na.recode = FALSE)
#' @export
setMethod("initialize", "sdcHashingDefinitionClass",
    function(.Object, ...) {

      # Get any specified constructor elemnts
      my.names <- names(match.call(expand.dots = TRUE))

      # Get the data, if any
      my.tmp <- pmatch(my.names, c("x", "data", "X"))
      if(sum(!is.na(my.tmp)) > 1)
        stop("Specify only one dataset on which the hashing is to be defined.")
      my.data_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.data <- NULL
      if(my.data_index)
        my.data <- eval(parse(text = match.call()[my.names[my.data_index]]))
      if(!is.null(my.data) && !is.data.frame(my.data))
        stop("The data should be in a data frame.")

      # Check if recoding is specified
      my.tmp <- pmatch(my.names, c("na.recode"))
      my.na_recode_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.na_recode <- TRUE
      if(my.na_recode_index)
        my.na_recode <- eval(parse(text = match.call()[my.names[my.na_recode_index]]))
      if(!is.logical(my.na_recode) || length(my.na_recode) != 1)
        stop("Recoding should be specified as a single logical value.")
      .Object@na.recode <- my.na_recode

      if(!is.null(my.data)) {
        .Object@levels <- lapply(lapply(my.data, unique), sort)
        .Object@na.exists <- sapply(my.data, function(XX) { any(is.na(XX)) })
        .Object@lengths <- sapply(.Object@levels, length) + .Object@na.recode
        my.tmp1 <- as.integer(rep(1, length(.Object@lengths)+1))
        my.number <- 1
        for(j in seq(.Object@lengths)) {
          my.tmp1[j+1] <- suppressWarnings(.Object@lengths[j]*my.tmp1[j])
          if(is.na(my.tmp1[j+1])) {
            my.number <- my.number + 1
            my.tmp1[j+1] <- .Object@lengths[j] }
        }
        .Object@parts <- as.integer(my.number)
        my.tmp2 <- my.tmp1[-1] > shift(my.tmp1[-1], -1)
        my.tmp2[length(my.tmp2)] <- TRUE
        my.tmp3 <- unique((seq(my.tmp1[-1]))*my.tmp2)
        .Object@coverage <- (my.tmp3 - shift(my.tmp3, 1))[-1]
      }

      return(invisible(.Object))
    })

#' @title
#'   S4 Hash Class
#' @description
#'   This class defines the data hash that is used to identify
#'   cells to be targeted for statistical disclosure control.
#' @details
#'   The hashing definition presently handles only categorical
#'   fields, which limits the applicability of this software
#'   accordingly.
#' @slot .Data
#'   A matrix.
#' @slot Hdef
#'   The hashing definition on which the
#'   \code{\linkS4class{sdcHashClass}} object is built.
#' @slot forward
#'   Logical indicating whether a forward hash or a backward
#'   hash, in terms of field order, has been used.
#' @slot hash
#'   A character string containing one hash representation
#'   for each data record.
#' @examples
#' new("sdcHashClass")
#' @export
setClass("sdcHashClass",
    contains = c("matrix"),
    slots = c(Hdef = "sdcHashingDefinitionClass",
              forward = "logical",
              hash = "character"),
    prototype = prototype(Hdef = new("sdcHashingDefinitionClass"),
                          forward = FALSE,
                          hash = character(0)),
    validity = function(object) {
      TRUE
    } )

#' @describeIn sdcHashClass
#' @details
#'   Does not currently permit use of an arbitrary hashing definition.
#' @param .Object
#'   A \code{\linkS4class{sdcHashingDefinitionClass}} object.
#' @param ...
#'   The optional parameters specifying the data
#'   and other options.
#' @examples
#' new("sdcHashClass")
#' my.X <- data.frame(matrix(ifelse(runif(5000)>.5, TRUE, FALSE), ncol = 50))
#' new("sdcHashClass", X = my.X)
#' new("sdcHashClass", X = my.X, na.recode = FALSE, which = 2:4,
#'     forwardHashing = TRUE)
#' @export
setMethod("initialize", "sdcHashClass",
    function(.Object, ...) {

      # Get any specified constructor elemnts
      my.names <- names(match.call(expand.dots = TRUE))

      # Get the data, if any
      my.tmp <- pmatch(my.names, c("x", "data", "X"))
      if(sum(!is.na(my.tmp)) > 1)
        stop("Specify only one dataset on which the hashing is to be defined.")
      my.data_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.data <- NULL
      if(my.data_index)
        my.data <- eval(parse(text = match.call()[my.names[my.data_index]]))
      if(!is.null(my.data) && !is.data.frame(my.data))
        stop("The data should be in a data frame.")

      # Check if forward hashing is specified
      my.tmp <- pmatch(my.names, c("forwardHash", "forward.hash"))
      my.forward_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.forward <- FALSE
      if(my.forward_index)
        my.forward <- eval(parse(text = match.call()[my.names[my.forward_index]]))
      if(!is.logical(my.forward) || length(my.forward) != 1)
        stop("Forward hashing should be specified as a single logical value.")
      .Object@forward <- my.forward

      if(!is.null(my.data)) {
      
        # Identify the fields to be hashed
        my.tmp <- pmatch(my.names, c("which", "fields", "variables"))
        my.fields_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
        # Work with fields by name
        my.fields <- names(my.data)
        if(my.fields_index)
          my.fields <- eval(parse(text = match.call()[my.names[my.fields_index]]))
        if(is.numeric(my.fields) && length(my.fields) < ncol(my.data))
          my.fields <- names(my.data)[my.fields]
        if(!is.character(my.fields))
          stop("Fields for hashing should be specified by name or index.")
        if(any(my.fields %!in% names(my.data)))
          stop(paste("Unidentified fields:",
                     paste(setdiff(my.fields, names(my.data)), collapse = " ")))
        # Convert to a logical indicator of fields to be hashed
        my.field_last <- 0
        my.fields <- ifelse(seq(ncol(my.data)) %in% match(my.fields, names(my.data)),
                           TRUE, FALSE)
        my.field_last <- max(which(my.fields))

        # Get the hashing definition object
        .Object@Hdef <- new("sdcHashingDefinitionClass", ...)
        
        # Prepare the vectors that will accept the calculated hashes
        #   * calculate greatest hash bin
        my.bin_max <- min(seq(.Object@Hdef@parts)[
            my.field_last <= cumsum(.Object@Hdef@coverage)])

        #   * create a list entry for each necessary hash bin
        my.list <- vector("list", my.bin_max)
        my.list_names <- vector("character", my.bin_max)

        i <- 0                                          # Hash base index
        #   * loop over the hash bins
        for(j in seq(my.bin_max)) {
          #   - get relevant lengths (number of category values per variable)
          my.lengths <- .Object@Hdef@lengths[(i+1):min(my.field_last,
                                             .Object@Hdef@coverage[j]+i)]
          #   - get the subset of fields over which to calculate, this being
          #     a vector of logicals that is used as a multiplier in the hash
          #     calculation, so that only the specified fields are included
          my.fields_subset <- my.fields[(i+1):min(my.field_last,
                                        .Object@Hdef@coverage[j]+i)]
          #   - calculate the multipliers per variable
          if(.Object@forward) my.multipliers <- c(1, my.lengths[-length(my.lengths)])
          else my.multipliers <- c(1, rev(my.lengths[-1]))
          #   - calculate hash in the event that it includes more than one variable
          if((my.field_last-i-1)>0) {
            if(.Object@forward)
              my.list[[j]] <- apply(my.data[,(i+1):min(my.field_last,
                                                       .Object@Hdef@coverage[j]+i)], 1,
                                    function(XX, YY = cumprod(my.multipliers),
                                             ZZ = my.fields_subset) {
                                      sum(ZZ*as.numeric(XX)*YY) } )
            else my.list[[j]] <- apply(my.data[,(i+1):min(my.field_last,
                                                          .Object@Hdef@coverage[j]+i)], 1,
                                       function(XX, YY = rev(cumprod(my.multipliers)),
                                                ZZ = my.fields_subset) {
                                         sum(ZZ*as.numeric(XX)*YY) } )
          } else my.list[[j]] <- my.data[, my.field_last]
          #   - name the hash bin
          my.list_names[j] <- paste("Hash", j, sep = "")
          #   - increment the variable index to start the next hash bin
          i <- i + .Object@Hdef@coverage[j]
        }

        #   * concatenate the hash bins to create the overall hash
        if(my.bin_max>1) {
          eee <- my.list[[1]]
          for(j in 2:my.bin_max) eee <- cbind(eee, my.list[[j]])
          .Object@hash <- apply(eee, 1, function(XX) {
            if(any(is.na(XX))) return(NA)
            else return(paste(XX, collapse="X")) })
        } else .Object@hash <- as.character(my.list[[1]])

        # Create and label the matrix
        .Object@.Data <- do.call(cbind, my.list)
        if(is.null(rownames(my.data)))
          rownames(.Object) <- 1:nrow(my.data)
        else rownames(.Object) <- rownames(my.data)
        colnames(.Object) <- my.list_names

      }
      
      return(invisible(.Object))
    })
    
#' @title
#'   Get a Hash List Subset
#' @description
#'   Returns a subset of a \code{\linkS4class{sdcHashClass}} object.
#' @details
#'   Does what it says on the package.
#' @param x
#'   An \code{\linkS4class{sdcHashClass}} object.
#' @param i
#'   An index of the hash.
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(5000)>.5, TRUE, FALSE), ncol = 50))
#' my.hc <- new("sdcHashClass", X = my.X)
#' my.hc[2:3]
#' @seealso \code{\linkS4class{sdcHashClass}},
#' @export
setMethod("[", signature = "sdcHashClass",
    function(x, i) {
      if(length(x@Hdef@parts) == 0)
        warning("Empty sdcHashClass object. No subset taken.")
      else {
        if(any(i %!in% seq(x@Hdef@parts)))
          stop(paste("Subset index not within 1:",
                     x@Hdef@parts, sep=""))
        x@.Data <- x@.Data[, i, drop = FALSE]
        x@hash <- apply(x@.Data, 1, paste, collapse = "X")
        my.in <- rep(ifelse(seq(x@Hdef@parts) %in% i, TRUE, FALSE),
                     x@Hdef@coverage)
        x@Hdef@levels <- x@Hdef@levels[which(my.in)]
        x@Hdef@na.exists <- x@Hdef@na.exists[which(my.in)]
        x@Hdef@lengths <- x@Hdef@lengths[which(my.in)]
        x@Hdef@coverage <- x@Hdef@coverage[i]
        x@Hdef@parts <- length(i)
      }
      return(invisible(x))
    })

#' @describeIn show
#' @title
#'   Show Method For \code{\linkS4class{sdcHashClass}} Objects
#' @description
#'   Returns an abbreviated summary of a \code{\linkS4class{sdcHashClass}} object.
#' @details
#'   As described in the documentation for \code{\link[methods]{show}},
#'   the default for printing out S4 objects is a call to
#'   a \code{\link[methods]{show}} method, whereas
#'   \code{\link[base]{print}} is the default for S3 objects.  This gives
#'   the possibility to have an alternative way to print the contents
#'   of a user-defined S4 class object.
#' @examples
#' show(new("sdcHashClass"))
#' @export
setMethod("show", signature = "sdcHashClass",
    function(object) {
      return(invisible(NULL))
    })

#' @describeIn show
#' @description
#'   The summary method runs and returns the output of
#'   some calculations on an \code{\linkS4class{sdcHashClass}}
#'   object, together with an abbreviated summary of the
#'   object itself.
#' @param object
#'   A \code{\linkS4class{sdcHashClass}} object.
#' @examples
#' summary(new("sdcHashClass"))
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' summary(my.hc)
#' @export
setMethod("summary", signature = "sdcHashClass",
    function(object) {

#      # Look at the inverse hash within one or more cells to check the data
#      #   X: A dataset (standardized or not)
#      #   Y: A hash list derived from the dataset X
#      get5TablesByHash <- function(X, Y) {
#        ccc <- Y@hash[duplicated(Y@hash)]   # Get the duplicate hash values
#        ddd <- sample.int(length(ccc), 5)         # Identify 5 cells at random
#        for(i in ddd) {
#          eee <- (1:length(Y@hash))[Y@hash==ccc[i]]
#          if(verbose) {
#            cat(paste("Data for Hash ", Y@hash[i], ":\n", sep=""))
#            print(X[eee,])
#          }
#        }
#      }

      Hdef <- object@Hdef
      cellsX <- table(object@hash)

      cat(paste("Number of Fields in Hash:", sum(Hdef@coverage)))
      cat(paste("\nNumber of Hash Parts:", Hdef@parts))
      cat(paste("\nNumber of Fields in Each Hash Part:", paste(Hdef@coverage, collapse=" ")))
      cat(paste("\nNumber of Records Hashed:", length(object@hash),
                "(check)", sum(cellsX),
                "(check uniques plus duplicates)",
                length(unique(object@hash)) + sum(duplicated(object@hash))))
      cat(paste("\nNumber of Records in Cells With >=2 Records:",
                sum(duplicated(object@hash) | duplicated(object@hash, fromLast=TRUE)),
                "(check 1)", 2*length(unique(object@hash[duplicated(object@hash)])) +
                                 sum(duplicated(object@hash[duplicated(object@hash)])),
                "(check 2)", sum(cellsX[cellsX>1])))
      cat(paste("\nNumber of Records in Cells With >=3 Records:",
                  3*length(unique(object@hash[triplicated(object@hash)])) +
                      sum(duplicated(object@hash[triplicated(object@hash)])),
                "(check)", sum(cellsX[cellsX>2])))
      cat(paste("\nNumber of Cells With >=2 Records:",
                length(unique(object@hash[duplicated(object@hash)])),
                "(check)", length(cellsX[cellsX>1])))
      cat(paste("\nNumber of Cells With >=3 Records:",
                length(unique(object@hash[triplicated(object@hash)])),
                "(check)", length(cellsX[cellsX>2])))
      cat(paste("\nNumber of Cells With One Record Alone:",
                sum(1-(duplicated(object@hash) |
                    duplicated(object@hash, fromLast=TRUE))),
                "(check)", length(cellsX[cellsX==1])))
      my.retval <- table(cellsX)

      cat("\nEmpirical Distribution of Cell Size:\n")
      invisible(return(my.retval))
    })

#' @name sdcProcessCells
#' @title
#'   Identify or Remove Small or Large Cells From a Hash List
#' @description
#'   Returns a hash list identifying the appropriate cells, or a
#'   subset of the data provided with the specified cells removed.
#' @param Y
#'   A hash list as returned by \code{\linkS4class{sdcHashClass}}.
#' @param X
#'   A data frame corresponding to \code{Y}.
#' @param cutoff
#'   The lower or upper cutoff of cell size.
#' @keywords big cells, small cells, identify, remove
#' @seealso \code{\linkS4class{sdcHashClass}},
NULL

#' @rdname sdcProcessCells
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' identifySmallCells(my.hc)
#' @export
identifySmallCells <- function(Y, cutoff=1) {
  cellsY <- table(Y@hash)   #  Note that by default NA's are excluded here
  setdiff(Y@hash, Y@hash[match(names(cellsY[cellsY>cutoff]), Y@hash)])
}

#' @rdname sdcProcessCells
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' identifyBigCells(my.hc)
#' @export
identifyBigCells <- function(Y, cutoff=1) {
  cellsY <- table(Y@hash)   #  Note that by default NA's are excluded here
  Y@hash[match(names(cellsY[cellsY>cutoff]), Y@hash)]
}

#' @rdname sdcProcessCells
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' identifySmallCellRecords(my.hc)
#' @export
identifySmallCellRecords <- function(Y, cutoff=1)
  matchAll(identifySmallCells(Y, cutoff=cutoff), Y@hash)
  
#' @rdname sdcProcessCells
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' identifyBigCellRecords(my.hc)
#' @export
identifyBigCellRecords <- function(Y, cutoff=1)
  matchAll(identifyBigCells(Y, cutoff=cutoff), Y@hash)
  
#' @rdname sdcProcessCells
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' removeSmallCellRecords(my.X, my.hc)
#' @export
removeSmallCellRecords <- function(X, Y, cutoff=1)
  X[identifyBigCellRecords(Y, cutoff=cutoff),]
  
#' @rdname sdcProcessCells
#' @examples
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.hc <- new("sdcHashClass", X = my.X)
#' removeBigCellRecords(my.X, my.hc)
#' @export
removeBigCellRecords <- function(X, Y, cutoff=1)
  X[identifySmallCellRecords(Y, cutoff=cutoff),]

#' @title
#'   S4 Target Definition
#' @description
#'   The SDC target definition is used to calculate a target matrix.
#' @details
#'   The hashing definition presently handles only categorical
#'   fields, which limits the applicability of this software
#'   accordingly.
#' @slot dim
#'   The dimensions of the data on which the target definition
#'   is calculated.
#' @slot dimnames
#'   A list of row and column names of the data on which the
#'   target definition is calculated, in the form returned
#'   by \code{\link[base]{dim}}.
#' @slot cutoff
#'   The minimum residual cell size.
#' @slot rownames
#'   The row names of the target elements.
#' @slot size
#'   The size of each target element.
#' @slot index
#'   The index of each target element.
#' @slot hash
#'   The target hash vector.
#' @examples
#' new("sdcTargetDefinitionClass")
#' @export
setClass("sdcTargetDefinitionClass",
    slots = c(dim = "integer",
              dimnames = "list",
              cutoff = "integer",
              rownames = "character",
              size = "integer",
              index = "character",
              hash = "character"),
    prototype = list(dim = as.integer(c(0, 0)),
                     dimnames = list(character(0),
                                     character(0)),
                     cutoff = as.integer(1),
                     rownames = character(0),
                     size = integer(0),
                     index = character(0),
                     hash = character(0) ),
    validity = function(object) {
      TRUE
    } )

#' @describeIn sdcTargetDefinitionClass
#' @param .Object
#'   An \code{\linkS4class{sdcTargetDefinitionClass}} object.
#' @param ...
#'   The optional parameters specifying the data
#'   and other options.
#'   Among these are
#' \tabular{ll}{
#' Name   \tab Description                                         \cr
#' X      \tab A data frame with the data, not necessarily         \cr
#'        \tab standardised.                                       \cr
#' Hdef   \tab A hashing definition, if generated externally,      \cr
#'        \tab perhaps because it comes from a larger set of data. \cr
#' cutoff \tab The minimum residual cell size.                     \cr
#' }
#' @examples
#' set.seed(256)
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' new("sdcTargetDefinitionClass", X = my.X)
#' @seealso \code{\linkS4class{sdcTargetMatrixClass}},
#' @export
setMethod("initialize", "sdcTargetDefinitionClass",
    function(.Object, ...) {
      # Variables starting with X derive from the data.
      # Variables starting with H relate to the hashes.
      #
      getStandardData <- function (X, H) {
        ddd <- X
        if(!all(class(X) %in% "factor"))
          ddd <- data.frame(lapply(X, function(Z) { factor(Z) }))
        my.special <- matchNone(ddd)
        for(i in seq(ddd)) {
          if(class(ddd[[i]]) %in% "factor") {
            for(j in seq(H@levels[[i]]))
              levels(ddd[[i]])[levels(ddd[[i]]) == H@levels[[i]][j]] <-
                  paste(my.special, j, sep = "")
            levels(ddd[[i]]) <- substr(levels(ddd[[i]]),
                                       start = nchar(my.special) + 1,
                                       stop = max(nchar(levels(ddd[[i]]))))
          }
          if(H@na.recode)
            levels(ddd[[i]])[is.na(X[[i]])] <- 0
        }
        ddd
      }

      # Get any specified constructor elements
      my.names <- names(match.call(expand.dots = TRUE))

      # Get the data
      my.tmp <- pmatch(my.names, c("x", "data", "X"))
      if(sum(!is.na(my.tmp)) > 1)
        stop("Specify only one dataset on which the hashing is to be defined.")
      my.data_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.data <- NULL
      if(my.data_index)
        my.data <- eval(parse(text = match.call()[my.names[my.data_index]]))
      if(!is.null(my.data) && !is.data.frame(my.data))
        stop("The data should be in a data frame.")

      # Check if recoding is specified
      my.tmp <- pmatch(my.names, c("na.recode"))
      my.na_recode_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.na_recode <- TRUE
      if(my.na_recode_index)
        my.na_recode <- eval(parse(text = match.call()[my.names[my.na_recode_index]]))
      if(!is.logical(my.na_recode) || length(my.na_recode) != 1)
        stop("Recoding should be specified as a single logical value.")

      # Get the hashing definition
      my.tmp <- pmatch(my.names, c("hdef", "Hdef"))
      if(sum(!is.na(my.tmp)) > 1)
        stop("Specify only one hashing definition.")
      my.Hdef_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.Hdef <- NULL
      if(my.Hdef_index)
        my.Hdef <- eval(parse(text = match.call()[my.names[my.Hdef_index]]))
      else my.Hdef <- new("sdcHashingDefinitionClass",
                          X = my.data, na.recode = my.na_recode)
      if(!is(my.Hdef, "sdcHashingDefinitionClass"))
        stop("Hashing definition must be of class sdcHashDefinitionClass.")

      # Check if cutoff is specified
      my.tmp <- pmatch(my.names, c("cutoff"))
      my.cutoff_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.cutoff <- 1
      if(my.cutoff_index)
        my.cutoff <- eval(parse(text = match.call()[my.names[my.cutoff_index]]))
      if(!is.numeric(my.cutoff) || length(my.cutoff) != 1)
        stop("Cutoff should be specified as a single integer value.")
      my.cutoff <- as.integer(my.cutoff)

      # Check if maximum number of cores to use is specified
      my.tmp <- pmatch(my.names, c("cores", "processors", "threads"))
      my.cores_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.cores <- 1
      if(my.cores_index)
        my.cores <- eval(parse(text = match.call()[my.names[my.cores_index]]))
      if(!is.numeric(my.cores) || length(my.cores) != 1)
        stop("Number of cores should be specified as a single integer value.")
      my.cores <- as.integer(my.cores)

      # Check if debug output is needed
      my.tmp <- pmatch(my.names, c("debugging"))
      my.debug_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.debug <- FALSE
      if(my.debug_index)
        my.debug <- eval(parse(text = match.call()[my.names[my.debug_index]]))
      if(!is.logical(my.debug) || length(my.debug) != 1)
        stop("Debugging output should be specified as a single logical value.")
        
      # Get the global temporary directory (used in debugging)
      my.tempdir <- tempdir()

      if(my.data_index) {
      
        # Step 1: Make a dataset with NA values set to 0 in each category.
        cat("Preparing the data...")
        flush.console()
        my.Xstrd <- getStandardData(my.data, my.Hdef)
        cat(paste(" done.", "\nCalculating statistics and",
                  "removing popular cell records...\n"))
        flush.console()

        # Step 2: Remove records in the popular cells (cutoff=1)
        #   * remove all records with the same actual _or_ missing values
        #   * does not currently permit use of an arbitrary hashing definition
        my.Xtest <- my.Xstrd
        my.H <- new("sdcHashClass", X = my.Xtest)
      
        #   * print a characterization of the root hash
        print(summary(my.H))
        flush.console()
      
        #   * identify big cell elements for removal first time through loop
        SUcurRecordSet <- identifyBigCellRecords(my.H, cutoff = my.cutoff)
        cat(paste(" done.  Records remaining to search:",
                    nrow(my.Xtest) - length(SUcurRecordSet), "out of",
                    nrow(my.Xtest), "\n"))
        flush.console()

        #   * initialize loop controling variables
        SUmaxSize <- ncol(my.data)
        SUcurSize <- SUindex <- 0

        #   * initialize return class
        my.target <- new("sdcTargetDefinitionClass")
        my.target@dim <- as.integer(dim(my.data))
        my.target@dimnames <- dimnames(my.data)
        my.target@cutoff <- as.integer(my.cutoff)

        # Step 3: Stop processing if left with a small enough final cell
        while(nrow(my.Xtest) - length(SUcurRecordSet) > my.cutoff) {
        
          #   * increment the size of record-discriminant sought
          SUcurSize <- SUcurSize + 1
          
          #   * whittle down the data for searching    
          if(length(SUcurRecordSet)>0) my.Xtest <- my.Xtest[-SUcurRecordSet,]

          #   * initialize counters
          SUcurIter <- 0
          SUmaxIter <- choose(length(my.Hdef@levels),
                              length(my.Hdef@levels) - SUcurSize) 

          #   * initialize collection variables for the outer loop
          SUcurVector <- SUcurLength <- vector(length=0, mode="integer")
          SUcurRecords <- vector(length=0, mode="integer")
          SUcurHashes <- SUcurIndex  <- vector(length=0, mode="character")

          # Step 3-1: Loop for a particular SU size
          cat(paste("\nSearching over ", SUmaxIter,
                    " target specifications of size ", SUcurSize,
                    "...", sep=""))
          flush.console()
        
          # Step 3-2: Allocate calculations to cores
          my.combinations <- combn(length(my.Hdef@levels),
                                   length(my.Hdef@levels) - SUcurSize,
                                   simplify=FALSE)
          my.last <- min(floor(length(my.combinations)/2), my.cores)
          ddd <- floor(quantile(seq(my.combinations),
                                seq(0, 1, length.out = my.last + 1)))
          my.index <- lapply(seq(ddd[-length(ddd)]),
                                          function(i) { ddd[i]:(ddd[i+1]-1) })
          my.index[[my.last]] <- c(my.index[[my.last]], length(my.combinations))

          # Step 3-3: Prepare for parallel processing
          my.cluster <- makeCluster(min(floor(length(my.combinations)/2), my.cores))
          registerDoParallel(my.cluster)
          
          SUobj <- foreach(i = my.index,
                           .packages = c("tuple", "sdcTarget")) %dopar% {

            # Get the subset of combinations for this core
            SUcurSubset <- my.combinations[i]

            # Initialize collection variables for the inner loop
            SUcurVector <- SUcurLength <- vector(length=0, mode="integer")
            SUcurRecords <- vector(length=0, mode="integer")
            SUcurHashes <- SUcurIndex  <- vector(length=0, mode="character")

            if(my.debug) {
              my.time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S_%Z")
              my.tmp <- paste(paste(rep("0", max(0, 3-nchar(SUcurSize))),
                                    collapse=""), SUcurSize, sep="")
              sink(file = tempfile(pattern = paste("sdcTarget_", my.tmp, "_", my.time, sep = ""),
                                   tmpdir = my.tempdir,
                                   fileext = ".txt"))
              print("==========")
              cat("Index:\n")
              print(i)
              cat("\nRemaining Data:\n")
              print(my.Xtest)
              cat("\nHash Definition:\n")
              print(my.Hdef)
              print("=====")
            }

            for(SUcurDefinition in SUcurSubset) {

              # Allow for tracking the progress of the search in real time
              SUcurIter <<- SUcurIter + 1
            
              # SUcurDefinition gives the set of variables included
              # SUcurInvDefinition gives the set of variables excluded
              SUcurInvDefinition <- setdiff(1:length(my.Hdef@levels), SUcurDefinition)
              # SUindex is the overall combination counter that gets incremented
              #   each time one passes through this loop (i.e., across set sizes)
              SUindex <- SUindex + 1
              # Htmp holds the variable subset specific hash for testing, with
              #   missing values propagated to hashes if na.hash is FALSE, given
              #   the recode of Hdef pursuant to na.hash (see above)

              Htmp <- new("sdcHashClass", X = my.Xtest, fields = SUcurDefinition)

              if(my.debug) {
                print("---")
                print(SUcurDefinition)
                print(Htmp)
              }

              # Identify the remaining popular cells
              ddd <- identifyBigCells(Htmp, cutoff = my.cutoff)
              
              if(my.debug) print(ddd)

              # If there are any big cells, identify their records for exclusion
              if(length(ddd)>0) {
                # Identify the big cell records (like identifyBigCellRecords)
                #   (call to the external function would be less efficient here)
                eee <- matchAll(unique(ddd), Htmp@hash)
                SUcurVector <- c(SUcurVector, ddd)
                SUcurLength <- c(SUcurLength, length(ddd))
                SUcurRecords <- c(SUcurRecords, eee)
                SUcurHashes <- c(SUcurHashes, Htmp@hash[eee])
                SUcurIndex <- c(SUcurIndex, rep(paste(SUcurInvDefinition, collapse=" "),
                                                length(eee)))
              }
            }

            if(my.debug) sink()

            return(list(SUcurVector = SUcurVector,
                        SUcurLength = SUcurLength,
                        SUcurRecords = SUcurRecords,
                        SUcurHashes = SUcurHashes,
                        SUcurIndex = SUcurIndex))
          }
          
          stopCluster(my.cluster)
          
          # Step 3-3: Put the lists together
          for(i in 1:length(SUobj)) {
            SUcurVector <- c(SUcurVector, SUobj[[i]]$SUcurVector)
            SUcurLength <- c(SUcurLength, SUobj[[i]]$SUcurLength)
            SUcurRecords <- c(SUcurRecords, SUobj[[i]]$SUcurRecords)
            SUcurHashes <- c(SUcurHashes, SUobj[[i]]$SUcurHashes)
            SUcurIndex <- c(SUcurIndex, SUobj[[i]]$SUcurIndex)
          }

          # Step 4: Save and report the findings
          SUcurRecordSet <- unique(SUcurRecords)
          my.target@rownames <- c(my.target@rownames,
                                  rownames(my.Xtest[SUcurRecords,]))
          my.target@size <- c(my.target@size,
                                 as.integer(rep(SUcurSize, length(SUcurRecords))))
          my.target@index <- c(my.target@index, SUcurIndex)
          my.target@hash <- c(my.target@hash, SUcurHashes)
          cat("\n ...")
          cat(paste(" found", length(SUcurVector), "targets of size", SUcurSize,
                    "in", length(SUcurRecordSet), "records."))
          flush.console()

          # * uniqueness is assured for residual unmatched data
          if(SUcurSize == SUmaxSize - 1) {
            if(length(SUcurRecordSet)>0)
              my.Xtest <- my.Xtest[-SUcurRecordSet,]
            if(nrow(my.Xtest) > 0) {
              my.target@rownames <- c(my.target@rownames, rownames(my.Xtest))
              my.target@size <- c(my.target@size,
                                  as.integer(rep(SUmaxSize, nrow(my.Xtest))))
              my.target@index <- c(my.target@index,
                                   rep(paste(seq(SUmaxSize), collapse=" "),
                                                nrow(my.Xtest)))
              Htmp <- new("sdcHashClass", X = my.Xtest)
              my.target@hash <- c(my.target@hash, Htmp@hash)
              cat(paste("\nAdded", nrow(my.Xtest),
                        "unique records to the target definition.\n"))
              flush.console()
            }
            break()
          } else {
            cat(paste("\nRecords remaining to search:",
                      nrow(my.Xtest)-length(SUcurRecordSet), "\n"))
            flush.console()
          }

        } # end of while()

        return(my.target)

      } else return(.Object)
      
    })

#' @title
#'   Estimate Time To Complete Calculation of Target Definition
#' @description
#'   This code provides a rough approximation for the number of
#'   seconds that a calculation of target definition may run.
#' @param X
#'   A data frame with the data (not necessarily standardised).
#' @param cutoff
#'   The minimum residual cell size.
#' @param cores
#'   The number of cores to be used in parallel processing.
#'   The default is linear processing: i.e., \code{core = 1}.
#' @param na.recode
#'   Whether NA values should be treated as levels.
#' @param Hdef
#'   A hashing definition (might be generated from a larger set of data).
#' @examples
#' set.seed(256)
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' estimateTimeToComputeTargetDefinition(X = my.X)
#' @export
estimateTimeToComputeTargetDefinition <- function(
  X,
  cutoff = 1,
  cores = 1,
  na.recode = TRUE,
  Hdef = new("sdcHashingDefinitionClass",
             X = X, na.recode = na.recode)
) {
  # Variables starting with X derive from the data.
  # Variables starting with H relate to the hashes.
  #
  getStandardData <- function (X, H) {
    ddd <- X
    if(!all(class(X) %in% "factor"))
      ddd <- data.frame(lapply(X, function(Z) { factor(Z) }))
    my.special <- matchNone(ddd)
    for(i in seq(ddd)) {
      if(class(ddd[[i]]) %in% "factor") {
        for(j in seq(H@levels[[i]]))
          levels(ddd[[i]])[levels(ddd[[i]]) == H@levels[[i]][j]] <-
              paste(my.special, j, sep = "")
        levels(ddd[[i]]) <- substr(levels(ddd[[i]]),
                                   start = nchar(my.special) + 1,
                                   stop = max(nchar(levels(ddd[[i]]))))
      }
      if(H@na.recode)
        levels(ddd[[i]])[is.na(X[[i]])] <- 0
    }
    ddd
  }

  # Get the data
  my.data <- X
  if(is.null(my.data) || !is.data.frame(my.data))
    stop("Data must be provided in a data frame.")

  # Check if recoding is specified
  my.na_recode <- na.recode
  if(!is.logical(my.na_recode) || length(my.na_recode) != 1)
    stop("Recoding should be specified as a single logical value.")

  # Get the hashing definition
  my.Hdef <- Hdef
  if(!is(my.Hdef, "sdcHashingDefinitionClass"))
    stop("Hashing definition must be of class sdcHashDefinitionClass.")

  # Check if cutoff is specified
  my.cutoff <- cutoff
  if(!is.numeric(my.cutoff) || length(my.cutoff) != 1)
    stop("Cutoff should be specified as a single integer value.")
  my.cutoff <- as.integer(my.cutoff)

  # Check how many cores to use in estimate
  my.cores <- cores
  if(!is.numeric(my.cores) || length(my.cores) != 1)
    stop("Number of cores should be specified as a single integer value.")
  my.cores <- as.integer(my.cores)

  SUallTime <- as.list(system.time({

    # Step 1: Make a dataset with NA values set to 0 in each category.
    cat("Preparing the data...")
    flush.console()
    my.Xstrd <- getStandardData(my.data, my.Hdef)
    cat(paste(" done.", "\nCalculating statistics and",
              "removing popular cell records..."))
    flush.console()

    # Step 2: Remove records in the popular cells (cutoff=1)
    #   * remove all records with the same actual _or_ missing values
    #   * does not currently permit use of an arbitrary hashing definition
    my.Xtest <- my.Xstrd
    my.H <- new("sdcHashClass", X = my.Xtest)

    #   * identify big cell elements for removal first time through loop
    SUcurRecordSet <- identifyBigCellRecords(my.H, cutoff = my.cutoff)
    cat(paste(" done.\nRecords remaining to search:",
              nrow(my.Xtest) - length(SUcurRecordSet), "out of",
              nrow(my.Xtest), "\n\n"))
    flush.console()

  }))$elapsed

  #   * initialize loop controling variables
  SUmaxSize <- ncol(my.data)
  SUcurSize <- SUindex <- 0
  SUcurTime <- vector(length=0, mode="numeric")
  SUcurLength <- vector(length=0, mode="integer")

  # Step 3: Stop processing if left with a small enough final cell
  while(nrow(my.Xtest) - length(SUcurRecordSet) > my.cutoff) {

    SUallTime <- SUallTime + as.list(system.time({

      #   * increment the size of record-discriminant sought
      SUcurSize <- SUcurSize + 1

      #   * whittle down the data for searching    
      if(length(SUcurRecordSet)>0) my.Xtest <- my.Xtest[-SUcurRecordSet, , drop = FALSE]

      #   * initialize counters
      my.last <- 3
      SUcurIter <- 0
      SUmaxIter <- min(choose(length(my.Hdef@levels),
                              length(my.Hdef@levels) - SUcurSize), my.last)

      # Step 3-1: Loop for a particular SU size
      cat(paste("Estimating over ", SUmaxIter,
                " target specifications of size ", SUcurSize,
                "...", sep=""))
      flush.console()

      # Step 3-2: Prepare to calculate times
      my.combinations <- sample(combn(length(my.Hdef@levels),
                                      length(my.Hdef@levels) - SUcurSize,
                                      simplify=FALSE), SUmaxIter)

      # Get the subset of combinations for this core
      SUcurSubset <- my.combinations

    }))$elapsed

    SUcurTime <- c(SUcurTime, as.list(system.time(

      for(SUcurDefinition in SUcurSubset) {

        # Allow for tracking the progress of the search in real time
        SUcurIter <<- SUcurIter + 1
            
        # SUcurDefinition gives the set of variables included
        # SUcurInvDefinition gives the set of variables excluded
        SUcurInvDefinition <- setdiff(1:length(my.Hdef@levels), SUcurDefinition)
        # SUindex is the overall combination counter that gets incremented
        #   each time one passes through this loop (i.e., across set sizes)
        SUindex <- SUindex + 1
        # Htmp holds the variable subset specific hash for testing, with
        #   missing values propagated to hashes if na.hash is FALSE, given
        #   the recode of Hdef pursuant to na.hash (see above)
        Htmp <- new("sdcHashClass", X = my.Xtest, fields = SUcurDefinition)

        # Identify the remaining popular cells
        ddd <- identifyBigCells(Htmp, cutoff = my.cutoff)
              
      }))$elapsed / SUmaxIter)

    SUcurLength <- c(SUcurLength, SUcurSize)

    print(SUcurTime)

    # Step 4: Save the findings
    cat(paste(" done.\n"))
    flush.console()

    my.decay <- 0.5
    SUremainder <- floor(nrow(my.Xtest) * (1 - my.decay))
    if(SUremainder < SUmaxIter) break()
    SUcurRecordSet <- sample(seq(nrow(my.Xtest)),
                             max(SUmaxIter, SUremainder))

  } # end of while()

  my.estimate <- choose(ncol(my.data), 1:(ncol(my.data)-1))[1:length(SUcurTime)] *
                     SUcurTime / my.cores
  SUallTime <- SUallTime + sum(my.estimate)
  
  cat(paste("\nEstimated Overall Seconds:", SUallTime, "\n"))
  return(invisible(list(overall = SUallTime, inner.loop = my.estimate, cores = cores)))
  
}

#' @title
#'   S4 Target Matrix
#' @description
#'   The taget matrix is a matrix with the same dimensions as the
#'   data from which it is derived, that indicates the number of
#'   combinations at a specific level of targetting for which the
#'   synthesis of a data element will make the record "sufficiently
#'   common".
#' @details
#'   Additional information is stored in the slots.
#' @slot .Data
#'   A matrix.
#' @slot Tdef
#'   A target definition class object.
#' @export
setClass("sdcTargetMatrixClass",
    contains = c("matrix"),
    slots = c(Tdef = "sdcTargetDefinitionClass"),
    prototype = prototype(Tdef = new("sdcTargetDefinitionClass")),
    validity = function(object) {
      TRUE
    } )

#' @describeIn sdcTargetMatrixClass
#' @param .Object
#'   An \code{\linkS4class{sdcTargetMatrixClass}} object.
#' @param ...
#'   The optional parameters specifying the basis of the
#'   target matrix (Tdef).
#' @examples
#' set.seed(256)
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.tdef <- new("sdcTargetDefinitionClass", X = my.X)
#' new("sdcTargetMatrixClass", Tdef = my.tdef)
#' @seealso \code{\linkS4class{sdcTargetDefinitionClass}},
#' @export
setMethod("initialize", "sdcTargetMatrixClass",
    function(.Object, ...) {

      # Get any specified constructor elements
      my.names <- names(match.call(expand.dots = TRUE))

      # Get the data
      my.tmp <- pmatch(my.names, c("Tdef", "tdef"))
      if(sum(!is.na(my.tmp)) > 1)
        stop("Specify only one dataset on which the hashing is to be defined.")
      my.tdef_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.tdef <- NULL
      if(my.tdef_index)
        my.tdef <- eval(parse(text = match.call()[my.names[my.tdef_index]]))
      if(!is.null(my.tdef) && !is(my.tdef, "sdcTargetDefinitionClass"))
        stop("Target definition not an sdcTargetDefinitionClass object.")

      if(!is.null(my.tdef) && (length(my.tdef@index) > 0)) {
        # Count of targets found by target specification, ordered by prevalence
        #   (descending), and within prevelance ordered by number of variables
        #   counted in definition:
        ddd <- table(my.tdef@index)[order(-table(my.tdef@index),
                   unlist(lapply(gregexpr(" ",names(table(my.tdef@index))),
                       length)))]
        # Count of targets found by number of associated variables (where it is
        #   assumed that the target definition corresponds to the values identified
        #   as risky, and not the opposite) times the number of such variables
        #   divided by total values in dataset (i.e., proportion of SDC data):
        eee <- sum(ddd*(unlist(lapply(lapply(gregexpr(" ",
                   dimnames(ddd)[[1]]), function(Z) {Z>0}), sum))+1)) /
                       prod(my.tdef@dim) * 100
        cat(paste("Overall proportion of values for SDC: ",
                  signif(eee,3), "%", sep=""))
      } else if(!is.null(my.tdef) && any(my.tdef@dim > 0))
        cat("Overall proportion of values for SDC: 0%\n")
      flush.console()

      if(!is.null(my.tdef) && (length(my.tdef@index) > 0)) {
        # Retrieve the row name associated with each target
        fff <- unlist(lapply(strsplit(my.tdef@rownames, "\\."), function(Z) { Z[1] }))
        fffu <- unique(fff)
        fffunum <- length(fffu)
        # Proportion of rows with data subject to SDC:
        ggg <- fffunum / my.tdef@dim[1] * 100
        cat(paste("\nProportion of data records subject to SDC: ",
                  signif(ggg,3), "%   (", fffunum,
                  " records out of ", my.tdef@dim[1], ")", sep=""))
        flush.console()

        # Proportion of records having multiple targets by target size:
        hhhn2 <- hhhn <- hhhd <- vector(length=my.tdef@dim[2], mode="integer")
        for(i in seq(my.tdef@dim[2])) {
          hhhi <- table(fff[my.tdef@size == i])
          hhhn[i] <- length(hhhi[hhhi > 1])
          hhhn2[i] <- length(unique(fff[my.tdef@size == i][which(unlist(lapply(
                          strsplit(my.tdef@rownames[my.tdef@size == i], "\\."),
                              length)) > 1)]))
          hhhd[i] <- length(unique(fff[my.tdef@size == i]))
        }
        cat(paste("\nProportion of records having multiple targets by target size:",
                  "\n", paste(signif(hhhn/hhhd*100,3), collapse="% "),
                  "\n   (NaN indicates that no targets of that size were found)",
                  "\n   Check two different calculations: ", identical(hhhn, hhhn2)))

        # Check that within record there is only one size of target
        cat("\nVerifying that each record associates with only one target size... ")
        flush.console()
        cat(length(table(apply(table(fff, my.tdef@size), 1, function(Z) {sum(Z>0)})))==1)

        # Create target substitution matrix for the data 
        cat("\nCreating target substitution matrix for the SDC data subset...")
        flush.console()
        yyy <- matrix(0, fffunum, my.tdef@dim[2]) 
        #   * open the target specifications
        zzz <- lapply(strsplit(my.tdef@index, " "), as.integer)
        #   * map the target records to their unique record values
        www <- match(fff, fffu)
        #   * add up the targets in each record
        for(i in 1:length(fff))
          yyy[www[i], zzz[[i]]] <- yyy[www[i], zzz[[i]]] + 1
        cat(" done.\nMapping to original data...")
        flush.console()

        yyyo <- matrix(0, my.tdef@dim[1], my.tdef@dim[2]) 
        rownames(yyyo) <- my.tdef@dimnames[[1]]
        colnames(yyyo) <- my.tdef@dimnames[[2]]

        www <- match(fffu, my.tdef@dimnames[[1]])
        for(i in 1:length(fffu)) yyyo[www[i],] <- yyy[i,]
        cat(" done.")
        flush.console()

        # Print out various checks on the target substitution matrix
        cat("\nMinimum number of records at each level of preservation:")
        print(table(apply(yyyo, 1, function(Z) {sum(Z==0)})))
        cat("\nPercent synthesis required by field for maximum partial synthesis:\n")
        print(signif(apply(yyyo, 2, function(Z) {sum(Z!=0)})/nrow(yyyo)*100, 3))
        flush.console()
      } else {
        yyyo <- matrix(nrow=0, ncol=0)
        my.tdef <- new("sdcTargetDefinitionClass")
      }
      
      # Prepare to return the sdcTargetMatrixClass object
      .Object@.Data <- yyyo
      .Object@Tdef <- my.tdef
      return(invisible(.Object))
    })

#' @title
#'   S4 Substitution Matrix
#' @description
#'   The substitution matrix is a matrix with the same dimensions as the
#'   data from which it is derived, that indicates which elements are to
#'   be subjected to a statistical disclosure control process.
#' @details
#'   The substitution matrix is calculated from a target matrix.
#'   Specification of forwards direction with a \code{cutoff} of 0
#'   results in a complete data synthesis substitution matrix.
#'   Specification of forwards direction with a \code{cutoff} of 1,
#'   or backwards direction with a sufficiently large \code{cutoff},
#'   results in the maximum partial synthesis substitution matrix.
#'   Specification of backwards direction with a \code{cutoff} of 1
#'   results in the minimum partial synthesis substitution matrix.
#'   Maximum partial synthesis is the default.
#' @slot .Data
#'   A matrix.
#' @slot forwards
#'   Indicates the direction in which to process the target matrix.
#'   Defaults to \code{TRUE}.
#' @slot cutoff
#'   Cutoff with respect to which an element in the target matrix
#'   should be indicated for substitution. Default value is 1.
#' @slot T
#'   A target matrix object.
#' @export
setClass("sdcSubstitutionMatrixClass",
    contains = c("matrix"),
    slots = c(forwards = "logical",
              cutoff = "numeric",
              T = "sdcTargetMatrixClass"),
    prototype = prototype(forwards = TRUE,
                          cutoff = 1,
                          T = new("sdcTargetMatrixClass",
                                  Tdef = new("sdcTargetDefinitionClass"))),
    validity = function(object) {
      TRUE
    } )
    
#' @describeIn sdcSubstitutionMatrixClass
#' @param .Object
#'   An \code{\linkS4class{sdcSubstitutionMatrixClass}} object.
#' @param ...
#'   The optional parameters specifying the basis of the
#'   substitution matrix. These include an indicator of
#'   forwards or backwards matrix calculation
#'   (\code{forwards}),
#'   a cutoff value relative to entries in the target matrix
#'   (\code{cutoff}), and the target matrix basis for the
#'   calculation (\code{T}).
#' @examples
#' set.seed(256)
#' my.X <- data.frame(matrix(ifelse(runif(500)>.5, TRUE, FALSE), ncol = 5))
#' my.smc <- new("sdcSubstitutionMatrixClass",
#'               T = new("sdcTargetMatrixClass",
#'                       Tdef = new("sdcTargetDefinitionClass", X = my.X)))
#' @seealso \code{\linkS4class{sdcTargetMatrixClass}},
#' @export
setMethod("initialize", "sdcSubstitutionMatrixClass",
    function(.Object, ...) {

      # Get any specified constructor elements
      my.names <- names(match.call(expand.dots = TRUE))

      # Get the target matrix
      my.tmp <- pmatch(my.names, c("Tmatrix", "tmatrix"))
      if(sum(!is.na(my.tmp)) > 1)
        stop("Specify only one dataset on which the hashing is to be defined.")
      my.T_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.T <- NULL
      if(my.T_index)
        my.T <- eval(parse(text = match.call()[my.names[my.T_index]]))
      if(is.null(my.T) || !is(my.T, "sdcTargetMatrixClass"))
        stop("No target matrix provided.")
    
      # Check if direction of calculation is specified
      my.tmp <- pmatch(my.names, c("forwards"))
      my.forwards_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.forwards <- TRUE
      if(my.forwards_index)
        my.forwards <- eval(parse(text = match.call()[my.names[my.forwards_index]]))
      if(!is.logical(my.forwards) || length(my.forwards) != 1)
        stop("Direction of calculation (forwards) should be specified as a single logical value.")

      # Check if cutoff is specified
      my.tmp <- pmatch(my.names, c("cutoff"))
      my.cutoff_index <- ifelse(all(is.na(my.tmp)), 0, min(which(!is.na(my.tmp))))
      my.cutoff <- 1
      if(my.cutoff_index)
        my.cutoff <- eval(parse(text = match.call()[my.names[my.cutoff_index]]))
      if(!is.numeric(my.cutoff) || length(my.cutoff) != 1)
        stop("Cutoff should be specified as a single integer value.")
      my.cutoff <- as.integer(my.cutoff)

      # Create a working set of target variable names
      my.target_names <- colnames(my.T)

      # Stop if the target columns are unnamed
      if((ncol(my.T) > 0) && is.null(my.target_names))
        stop("Columns of the target matrix have no names.")

      # Turn the target matrix into a data frame, so that it will keep
      #   column names and format even if it turns into a single column.
      my.Target <- data.frame(my.T)
      if(length(my.target_names)==1) names(my.Target) <- my.target_names

      # Apply the appropriate directional test
      if(my.forwards) {
        cat(paste("Applying forwards cut with cutoff ", my.cutoff, ".\n", sep=""))
        .Object@.Data <- apply(my.Target, 2, function(Z) { Z >= my.cutoff })
      } else {
        cat(paste("Applying backwards cut with cutoff ", my.cutoff, ".\n", sep=""))
        # Catch if target has been reduced to a vector from a matrix
        ddd <- apply(my.Target, 1, max)
        .Object@.Data <- apply(my.Target, 2, function(Z) { (Z > 0) & (Z > ddd - my.cutoff) })
      }
      rownames(.Object@.Data) <- rownames(my.T)

      # Fill and return the substitution matrix object
      .Object@forwards <- my.forwards
      .Object@cutoff <- my.cutoff
      .Object@T <- my.T
      
      return(invisible(.Object))
    })

#getSdcRows <- function(..., criterion="random",
#                       rows=NULL, patterns=NULL, like=NULL,
#                       size=NULL, quantity=NULL,
#                       X=NULL, TD=NULL,
#                       n=1) {
#  ddd <- charmatch(criterion,
#                   c("random", "rows", "patterns", "like",
#                     "size", "quantity"))
#  if(!is.null(rows)) ddd <- 2
#  if(!is.null(patterns)) ddd <- 3
#  if(!is.null(like)) ddd <- 4
#  if(ddd==1 && !is.null(TD)) {
#    fff <- unlist(lapply(strsplit(TD$RowID, "\\."), function(Y) { Y[1] }))
#    if(!is.null(size)) fff <- fff[TD$Size==size]
#    if(!is.null(quantity)) fff <- names(table(fff)[table(fff)==3])
#    if(is.numeric(n) && length(n)==1 && n>0)
#      fff <- sample(fff, min(length(fff),n))
#    if(!is.null(X)) return(X[match(fff, rownames(X)),])
#    else return(fff)
#  } else if(ddd==2 && !is.null(TD)) {
#    fff <- unlist(lapply(strsplit(TD$RowID, "\\."), function(Y) { Y[1] }))
#    if(!is.null(X)) return(X[match(fff, rownames(X)),])
#    else {
#      ggg <- matchAll(rows, fff)
#      ggg <- ggg[order(TD$RowID[ggg])]
#      return(data.frame(RowID=unlist(lapply(strsplit(TD$RowID[ggg], "\\."),
#                                            function(Y) { Y[1] })),
#                        Size=TD[[2]][ggg],
#                        Index=TD[[3]][ggg],
#                        Hash=TD[[4]][ggg]))
#    }
#  } else if(ddd==3 && !is.null(patterns) && !is.null(TD) && !is.null(X)) {
#    if(is.numeric(patterns))
#      eee <- matchAll(paste(sort(patterns), collapse=" "), TD$Index)
#    else eee <- matchAll(patterns, TD$Index)
#    fff <- unique(unlist(lapply(strsplit(TD$RowID[eee], "\\."),
#                         function(Y) { Y[1] })))
#    ggg <- X[match(fff[order(TD$Hash[fff])], rownames(X)),]
#    return(ggg)
#  } else if(ddd==4 && !is.null(like) && !is.null(TD) && !is.null(X)) {
#    fff <- unlist(lapply(strsplit(TD$RowID, "\\."), function(Y) { Y[1] }))
#    ggg <- matchAll(like, fff)
#    hhh <- getSdcRows(patterns=unique(TD$Index[ggg]), TD=TD, X=X)
#    return(hhh)
#  } else {         
#    stop("Please give a full specification.")
#  }
#}
