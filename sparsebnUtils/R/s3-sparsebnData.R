#
#  s3-sparsebnData.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/28/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#------------------------------------------------------------------------------#
# sparsebnData S3 Class for R
#------------------------------------------------------------------------------#

#
# sparsebnData S3 class skeleton
#
# Data
# * data.frame data     // data
# * character type      // either "continuous", "discrete", or "mixed"
# * list ivn            // list of nodes under intervention for each row (observation)
# * list levels         // names of levels for each variable
#

#
# A convenience class for storing metadata associated with discrete and continuous data. In addition
#  to the data.frame containing all of the data, this class keeps track of the data type, number
#  of levels in each variable (discrete data only), and possible interventions on the observations.
#
# ivn:
#  Each component of the list 'ivn' is an integer vector indicating which nodes (columns) were under
#  intervention for the corresponding row. Note that 'ivn' must be a list in order to accommodate
#  the possibility that multiple nodes are under intervention, and different nodes may have
#  different treatments (i.e. number and/or identity of manipulated nodes).
#
# Examples:
#  1) data$ivn[1] = NULL: No nodes were under intervention for the first observation
#  2) data$ivn[10] = c(1): The first node was under intervention for the tenth observation
#  3) data$ivn[120] = c(1,5,500): The 1st, 5th, and 500th nodes were under intervention for the 120th observation
#
# levels:
#  This is a list containing the different levels for each node / variable in the dataset.
#
# Examples:
#  1) data$levels = NULL: There are no levels - data is continuous.
#  2) data$levels[1] = c(1, 2, 3): The first node has three possible levels, 1, 2, or 3.
#  3) data$levels[50] = c("A", "B"): The 50th node has two possible levels, "A" or "B.
#

#' sparsebnData class
#'
#' This class stores data that may contain interventions on some or all of the observations. It also
#' allows for the degenerate case with no interventions, i.e. purely observational data.
#'
#' The structure of a \code{sparsebnData} object is very simple: It contains a \code{data.frame} object,
#' a type identifier (i.e. discrete or continuous), a list of factor levels, and a list of interventions.
#' The 'levels' list should be the same size as the number of nodes and consist of names of the different
#' levels for each node. The 'ivn' list should be the same size as the number of rows in the dataset,
#' and each component indicates which column(s) in the dataset is (are) under intervention. If an
#' observation has no interventions, then the corresponding component is \code{NULL}. Thus, if the data is
#' purely observational, this list should contain only \code{NULL} values.
#'
#' Also inherits from \code{\link{list}}.
#'
#' @param x a \code{\link{data.frame}} or \code{\link{matrix}} object.
#' @param type either '\code{discrete}' or '\code{continuous}'.
#' @param levels (optional) list of level for each node.
#' @param ivn (optional) list of interventions for each observation.
#' @param n (optional) number of rows from data matrix to print.
#' @param ... (optional) additional arguments.
#'
#' @section Slots:
#' \describe{
#' \item{\code{data}}{(data.frame) Dataset.}
#' \item{\code{type}}{(character) Type of data: Either "continuous", "discrete", or "mixed".}
#' \item{\code{levels}}{(list) List of levels for each column in \code{data}.}
#' \item{\code{ivn}}{(list) List of columns under intervention for each row in \code{data}.}
#' }
#'
#' @section Methods:
#' \code{\link{print}}
#' \code{\link{num.samples}}
#' \code{\link{is.obs}}
#' \code{\link{count.levels}}
#' \code{\link{count.interventions}}
#' \code{\link{as.data.frame}}
#'
#' @docType class
#' @name sparsebnData
NULL

#' @rdname sparsebnData
#' @export
is.sparsebnData <- function(x){
    inherits(x, "sparsebnData")
} # END IS.SPARSEBNDATA

# sparsebnData constructor
#' @export
sparsebnData.list <- function(x, ...){

    if( !is.list(x)){
        stop("Input must be a list!")
    } else if( length(x) != 4 || !setequal(names(x), c("data", "type", "levels", "ivn"))){
        stop("Input is not coercable to an object of type sparsebnFit, check list for the following elements: data (data.frame), type (character), levels (list), ivn (list)")
    } else if( !check_if_data_matrix(x$data)){
        stop(sprintf("Component 'data' must be a valid data.frame or numeric object! <Current type: %s>", class(x$data)))
    } else if(!(x$type %in% c("continuous", "discrete", "mixed"))){
        stop(sprintf("\'type\' must be one of the following: \'continuous\', \'discrete\', \'mixed\'."))
    } else if(!is.null(x$levels)){
        if(ncol(x$data) != length(x$levels)){
            stop("The length of the levels list must equal the number of columns in the data!")
        }
    } else if(nrow(x$data) != length(x$ivn)){
        stop("The length of the ivn list must equal the number of rows in the data!")
    }

    num_missing <- count_nas(x$data)
    if(num_missing > 0){
        warning(has_missing_values(num_missing))
    }

    ### Final output
    structure(x, class = "sparsebnData")
} # END SPARSEBNDATA.LIST

# sparsebnData constructor
#  Default constructor for data.frame input
#' @rdname sparsebnData
#' @export
sparsebnData.data.frame <- function(x, type, levels, ivn, ...){

    type_list <- c("continuous", "discrete")

    ### User must specify type
    if(missing(type)){
        stop("The data type (continuous or discrete?) was not specified: Must choose type = 'continuous' or type = 'discrete'.")
        ivn <- vector("list", length = nrow(x))
    } else{
        match_string <- pmatch(type, type_list) # use partial matching to select type
        if(is.na(match_string)){ # if there was no match, error
            stop("Invalid 'type' entered: Must match one of \'continuous\', \'discrete\', \'mixed\'.")
        } else{ # if match was found, use it
            type <- type_list[match_string]
        }
    }

    #
    # If the user fails to specify a list of interventions, ASSUME all rows are observational. If the data
    #  is experimental, the user needs to specify this by passing in 'ivn' (see also sparsebnData.list).
    #
    if(missing(ivn)){
        message("A list of interventions was not specified: Assuming data is purely observational.")
        ivn <- vector("list", length = nrow(x))
    }

    #
    # If the user fails to specify a list of levels, attempt to infer them automatically.
    #
    if(missing(levels)){
        # message("A list of levels was not specified: Assuming data is continuous.")
        if(type == "continuous"){
            levels <- NULL
        } else{
            levels <- auto_count_levels(x)
        }
    }

    ### Final output
    sparsebnData.list(list(data = x, type = type, levels = levels, ivn = ivn))
} # END SPARSEBNDATA.DATA.FRAME

# sparsebnData constructor
#  Default constructor for matrix input
#' @rdname sparsebnData
#' @export
sparsebnData.matrix <- function(x, type, levels, ivn, ...){
    sparsebnData.data.frame(as.data.frame(x), type, levels, ivn)
} # END SPARSEBNDATA.MATRIX

#' @describeIn num.samples Extracts the number of samples of \link{sparsebnData} object.
#' @export
num.samples.sparsebnData <- function(x){
    nrow(x$data)
} # END NUM.SAMPLES.SPARSEBNDATA

#' Check if data is observational
#'
#' Returns TRUE if the data contains no interventions, i.e. is purely observational
#'
#' @param data a \code{\link{sparsebnData}} object.
#'
#' @export
is.obs <- function(data){
    all(unlist(lapply(data$ivn, is.null)))
} # END IS.OBS

#' Count the number of rows under intervention
#'
#' Returns the number of rows with at least one intervention
#'
#' @param data a \code{\link{sparsebnData}} object.
#'
#' @export
count.interventions <- function(data){
    sum(unlist(lapply(data$ivn, function(x) !is.null(x))))
} # END COUNT.INTERVENTIONS

#' Count the number of levels per variable
#'
#' Returns the number of levels per variable as an ordered vector.
#'
#' @param data a \code{\link{sparsebnData}} object.
#'
#' @export
count.levels <- function(data){
    unlist(lapply(data$levels, length))
} # END COUNT.LEVELS

# Default print method
#' @rdname sparsebnData
#' @export
print.sparsebnData <- function(x, n = 5L, ...){
    # print(utils::head(data$data, n = n), row.names = FALSE)
    .print_data_frame(x$data, topn = n)

    cat(sprintf("\n%d total rows (%d rows omitted)\n", num.samples(x), num.samples(x) - 2*n))
    if(is.obs(x)){
        cat(sprintf("Observational data with %s observations", x$type))
    } else{

        cat(sprintf("%s data w/ interventions on %d/%d rows.", capitalize(x$type), count.interventions(x), num.samples(x)))
    }
    ### Add a message about the interventions as well / if purely obs, etc.
} # END PRINT.SPARSEBNDATA

#' Convert a sparsebnData object back to a data.frame
#'
#' @param x a \code{\link{sparsebnData}} object.
#' @param ... (optional) additional argument to \code{as.data.frame}.
#'
#' @method as.data.frame sparsebnData
#' @export
as.data.frame.sparsebnData <- function(x, ...){
    data.frame(x$data, ...)
} # END AS.DATA.FRAME.SPARSEBNDATA

### Internal method for picking the correct family for fitting parameters
pick_family.sparsebnData <- function(data){
    if(data$type == "continuous"){
        return("gaussian")
    } else if(data$type == "discrete"){
        return("binomial")
    } else{
        stop("'mixed' type not supported for inference yet!")
    }
}

### Borrow the print.data.table method from the 'data.table' package without needing to import the entire package
###  This is an experimental method!
.print_data_frame <- function(x,
                              topn=5,   # (5) print the top topn and bottom topn rows with '---' inbetween
                              nrows=5, # (100) under this the whole (small) table is printed, unless topn is provided
                              row.names = TRUE, ...){
    if (!is.numeric(nrows)) nrows = 100L
    if (!is.infinite(nrows)) nrows = as.integer(nrows)
    if (nrows <= 0L) return(invisible())   # ability to turn off printing
    if (!is.numeric(topn)) topn = 5L
    topnmiss = missing(topn)
    topn = max(as.integer(topn),1L)
    if (nrow(x) == 0L) {
        if (length(x)==0L)
           cat("Null data.table (0 rows and 0 cols)\n")  # See FAQ 2.5 and NEWS item in v1.8.9
        else
           cat("Empty data.table (0 rows) of ",length(x)," col",if(length(x)>1L)"s",": ",paste(utils::head(names(x),6),collapse=","),if(ncol(x)>6)"...","\n",sep="")
        return()
    }
    if (topn*2<nrow(x) && (nrow(x)>nrows || !topnmiss)) {
        toprint = rbind(utils::head(x, topn), utils::tail(x, topn))
        rn = c(seq_len(topn), seq.int(to=nrow(x), length.out=topn))
        printdots = TRUE
    } else {
        toprint = x
        rn = seq_len(nrow(x))
        printdots = FALSE
    }
    toprint=format(toprint, ...)
    # FR #5020 - add row.names = logical argument to print.data.table
    if (isTRUE(row.names)) rownames(toprint)=paste(format(rn,right=TRUE),":",sep="") else rownames(toprint)=rep.int("", nrow(x))
    if (is.null(names(x))) colnames(toprint)=rep("NA", ncol(toprint)) # fixes bug #4934
    if (printdots) {
        toprint = rbind(utils::head(toprint,topn),"---"="",utils::tail(toprint,topn))
        rownames(toprint) = format(rownames(toprint),justify="right")
        print(toprint,right=TRUE,quote=FALSE)
        return(invisible())
    }
    if (nrow(toprint)>20L)
        # repeat colnames at the bottom if over 20 rows so you don't have to scroll up to see them
        toprint=rbind(toprint,matrix(colnames(toprint),nrow=1)) # fixes bug #4934
    print(toprint,right=TRUE,quote=FALSE)
    invisible()
}
