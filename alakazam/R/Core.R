# Common input/output and data structure manipulation functions for Alakazam

#' @include Alakazam.R
NULL

#### File I/O functions ####

#' Read a Change-O tab-delimited database file
#' 
#' \code{readChangeoDb} reads a tab-delimited database file created by a Change-O tool 
#' into a data.frame.
#'
#' @param    file       tab-delimited database file output by a Change-O tool.
#' @param    select     columns to select from database file.
#' @param    drop       columns to drop from database file.
#' @param    seq_upper  if \code{TRUE} convert sequence columns to upper case;
#'                      if \code{FALSE} do not alter sequence columns. See Value 
#'                      for a list of which columns are effected.

#' 
#' @return   A data.frame of the database file. Columns will be imported as is, except for 
#'           the following columns which will be explicitly converted into character 
#'           values:
#'           \itemize{
#'             \item  SEQUENCE_ID
#'             \item  CLONE
#'             \item  SAMPLE
#'           }
#'           And the following sequence columns which will converted to upper case if
#'           \code{seq_upper=TRUE} (default).
#'           \itemize{
#'             \item  SEQUENCE_INPUT
#'             \item  SEQUENCE_VDJ
#'             \item  SEQUENCE_IMGT
#'             \item  JUNCTION
#'             \item  GERMLINE_IMGT
#'             \item  GERMLINE_IMGT_D_MASK
#'           }
#'                   
#' @seealso  Wraps \code{\link{read.table}}.
#' @family   file input and output functions
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
#' 
#' # Subset columns and convert sequence fields to upper case
#' df <- readChangeoDb(file, select=c("SEQUENCE_ID", "SEQUENCE_IMGT"), seq_upper=TRUE)
#' 
#' # Drop columns and do not alter sequence field case
#' df <- readChangeoDb(file, drop=c("D_CALL", "DUPCOUNT"), seq_upper=FALSE)
#' 
#' @export
readChangeoDb <- function(file, select=NULL, drop=NULL, seq_upper=TRUE) {
    # Define column data types
    seq_columns <- c("SEQUENCE_INPUT", "SEQUENCE_VDJ", "SEQUENCE_IMGT", "JUNCTION", 
                     "GERMLINE_IMGT", "GERMLINE_IMGT_D_MASK")
    text_columns <- c("SEQUENCE_ID", "CLONE", "SAMPLE")
    
    # Read file
    db_df <- read.delim(file, as.is=TRUE, na.strings=c("", "NA", "None"))
    
    # Select columns
    select_columns <- colnames(db_df)
    if(!is.null(select)) { select_columns <- intersect(select_columns, select) }
    if(!is.null(drop)) { select_columns <- setdiff(select_columns, drop) }
    db_df <- subset(db_df, select=select_columns)
    
    # Convert sequence fields to upper case
    if (seq_upper) {
        for (x in intersect(seq_columns, select_columns)) {
            db_df[, x] <- toupper(db_df[, x]) 
        }
    }
    
    # Convert text fields to character
    for (x in intersect(text_columns, select_columns)) {
        db_df[, x] <- as.character(db_df[, x])
    }
    
    return(db_df)
}


#' Write a Change-O tab-delimited database file
#' 
#' \code{writeChangeoDb} is a simple wrapper around \code{\link{write.table}} with defaults 
#' appropriate for writing a Change-O tab-delimited database file from a data.frame.
#'
#' @param    data  data.frame of Change-O data.
#' @param    file  output file name.
#' 
#' @return   NULL
#' 
#' @seealso  Wraps \code{\link{write.table}}.
#' @family   file input and output functions
#' 
#' @examples
#' \dontrun{
#'   # Write a database
#'   writeChangeoDb(data, "changeo_output.tab")
#' }
#' 
#' @export
writeChangeoDb <- function(data, file) {
    write.table(data, file=file, quote=FALSE, sep="\t", row.names=FALSE)
}


#' Create a temporary folder
#'
#' \code{makeTempDir} creates a randomly named temporary folder in the 
#' system temp location.
#' 
#' @param    prefix  prefix name for the folder.
#' 
#' @return   The path to the temporary folder.
#' 
#' @seealso  This is just a wrapper for \code{\link{tempfile}} and 
#'           \code{\link{dir.create}}.
#' @family   file input and output functions
#' 
#' @examples
#' makeTempDir("Clone50")
#' 
#' @export
makeTempDir <- function(prefix) {
    temp_path <- tempfile(paste0(prefix, "-temp-"))
    dir.create(temp_path)
    
    return(temp_path)
}


#### Data manipulation functions ####

#' Translate a vector of strings
#' 
#' \code{translateStrings} modifies a character vector by substituting one or more 
#' strings with a replacement string.
#'
#' @param    strings      vector of character strings to modify.
#' @param    translation  named character vector or a list of character vectors specifying 
#'                        the strings to replace (values) and their replacements (names).
#' 
#' @return   A modified \code{strings} vector.
#' 
#' @details
#' Does not perform partial replacements. Each translation value must match a complete 
#' \code{strings} value or it will not be replaced.  Values that do not have a replacement
#' named in the \code{translation} parameter will not be modified.
#' 
#' Replacement is accomplished using \code{\link{gsub}}.
#' 
#' @seealso  See \code{\link{gsub}} for single value replacement in the base package.
#' 
#' @examples
#' # Using a vector translation
#' strings <- LETTERS[1:5]
#' translation <- c("POSITION1"="A", "POSITION5"="E")
#' translateStrings(strings, translation)
#' 
#' # Using a list translation
#' strings <- LETTERS[1:5]
#' translation <- list("1-3"=c("A","B","C"), "4-5"=c("D","E"))
#' translateStrings(strings, translation)
#' 
#' @export
translateStrings <- function(strings, translation) {
    # TODO:  use match instead for complete matching?  Currently regex characters in values will mess up the matching.
    for (n in names(translation)) {
        rep_regex <- paste(translation[[n]], collapse='|')
        strings <- gsub(paste0("^(", rep_regex, ")$"), n, strings)
    }
    
    return(strings)
}


# Check data.frame for valid columns and issue message if invalid
#
# @param   data     data.frame to check
# @param   columns  vector of column names to check
# @param   logic    one of "all" or "any" controlling whether all or at least one of
#                   the columns must be valid
# @return  TRUE if columns are valid and a string message if not.
checkColumns <- function(data, columns, logic=c("all", "any")) {
    # Check arguments
    logic <- match.arg(logic)
    
    data_names <- names(data)
    if (logic == "all") {
        # Check that all columns exist
        for (f in columns) {
            if (!(f %in% data_names)) { 
                msg <- paste("The column", f, "was not found") 
                return(msg)
            }
        }        
        # Check that all values are not NA
        for (f in columns) {
            if (all(is.na(data[, f]))) { 
                msg <- paste("The column", f, "contains no data") 
                return(msg)
            }
        }
    } else if (logic == "any") {
        # Check that columns exist
        if (!any(columns %in% data_names)) {
            msg <- paste("Input must contain at least one of the columns:", paste(columns, collapse=", "))
            return(msg)
        }
        # Check that all values are not NA
        invalid <- sapply(columns, function(f) all(is.na(data_names[, f])))
        if (all(invalid)) { 
            msg <- paste("None of the columns", paste(columns, collapse=", "), "contain data") 
            return(msg)
        }
    }
    
    # Return TRUE if all checks pass
    return(TRUE)
}

#### Plotting functions ####

# Define universal plot settings
#
# @param    sizing  defines the style and sizing of the theme. One of 
#                   \code{c("figure", "window")} where \code{sizing="figure"} is appropriately
#                   sized for pdf export at 7 to 7.5 inch width, and \code{sizing="window"}
#                   is sized for an interactive session.
#
# @return   A ggplot2 theme object.
getBaseTheme <- function(sizing=c("figure", "window")) {
    # Check arguments
    sizing <- match.arg(sizing)
    
    # Define universal plot settings appropriate for PDF figures
    if (sizing == "figure") {
        base_theme <- theme_bw() + 
            theme(text=element_text(size=8)) +
            theme(plot.title=element_text(size=8)) +
            theme(plot.background=element_blank(),
                  panel.grid.major=element_blank(), 
                  panel.grid.minor=element_blank()) +
            theme(strip.background=element_blank(),
                  strip.text=element_text(size=7, face='bold')) +
            theme(axis.title=element_text(size=8, vjust=0.25),
                  axis.text.x=element_text(size=8, vjust=0.5, hjust=0.5),
                  axis.text.y=element_text(size=8)) +
            theme(legend.text=element_text(size=7),
                  legend.title=element_text(size=7),
                  legend.key.height=grid::unit(10, "points"), 
                  legend.key.width=grid::unit(10, "points"))
    } else if (sizing == "window") {
        # Define universal plot settings appropriate for an interactive session
        base_theme <- theme_bw() + 
            theme(text=element_text(size=14)) +
            theme(plot.title=element_text(size=16)) +
            theme(strip.background=element_rect(fill='white'), 
                  strip.text=element_text(size=14, face='bold')) +
            theme(axis.title=element_text(size=16, vjust=0.25),
                  axis.text.x=element_text(size=16, vjust=0.5, hjust=0.5),
                  axis.text.y=element_text(size=16)) +
            theme(legend.text=element_text(size=14),
                  legend.title=element_text(size=14),
                  legend.key.height=grid::unit(18, "points"), 
                  legend.key.width=grid::unit(18, "points"))
    }
    
    return(base_theme)
}


#' Plot multiple ggplot objects
#' 
#' Plots multiple ggplot objects in an equally sized grid.
#' 
#' @param   ...    ggplot objects to plot.
#' @param   ncol   number of columns in the plot.
#' @return  NULL
#' 
#' @seealso
#' \link{ggplot}.
#' 
#' @references
#' Modified from:
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)
#' 
#' @export
multiggplot <- function(..., ncol=1) {
    p <- list(...)
    n <- length(p)
    layout <- matrix(seq(1, ncol*ceiling(n/ncol)), ncol=ncol, nrow=ceiling(n/ncol))
    
    # Plot
    if (n == 1) {
        plot(p[[1]])
    } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:n) {
            idx <- as.data.frame(which(layout == i, arr.ind=T))
            plot(p[[i]], vp=grid::viewport(layout.pos.row = idx$row, layout.pos.col=idx$col))
        }
    }
}


#### Generic statistical functions ####

#' Weighted meta-analysis of p-values via Stouffer's method
#'
#' \code{stoufferMeta} combines multiple weighted p-values into a meta-analysis p-value
#' using Stouffer's Z-score method.
#' 
#' @param    p   numeric vector of p-values.
#' @param    w   numeric vector of weights.
#' 
#' @return   A named numeric vector with the combined Z-score and p-value in the form
#'           \code{c(Z, pvalue)}.
#' 
#' @examples
#' # Define p-value and weight vectors
#' p <- c(0.1, 0.05, 0.3)
#' w <- c(5, 10, 1)
#'
#' # Unweighted
#' stoufferMeta(p)
#' 
#' # Weighted
#' stoufferMeta(p, w)
#' 
#' @export
stoufferMeta <- function(p, w=NULL) {
    if (is.null(w)) {
        w <- rep(1, length(p))/length(p)
    } else {
        if (length(w) != length(p)) { stop("Length of p and w must equal.") }
        w <- w/sum(w)
    }
    x <- qnorm(1 - p)
    Z  <- sum(w*x) / sqrt(sum(w^2))
    pvalue <- 1 - pnorm(Z)
    
    return(c(Z=Z, pvalue=pvalue))
}

# Stirling's approximation of the binomial coefficient
# 
# Calculates Stirling's approximation of the binomial coefficient for large numbers.
#
# @param    n  a vector of n.
# @param    k  a vector of k.
#
# @return   The approximation of log(n choose k). For n < 100 \link{lchoose} is used.
lchooseStirling <- function(n, k) {
    if (any(n < k)) {
        stop("n must be >= k")
    }
    
    n_len <- length(n)
    k_len <- length(k)
    nCk <- rep(NA, max(n_len, k_len))
    nCk[n == k] <- 0
    
    # a = index n_small
    # i = index k_small
    # x = index nCk_small
    #
    # b = index n_large
    # j = index k_large
    # y = index nCk_large
    #
    # Check for vector inputs and assign indexing
    if (n_len >= 1 & k_len >= 1 & n_len == k_len) {
        a <- i <- x <- (n < 100 & n != k)
        b <- j <- y <- (n >= 100 & n != k)
    } else if (n_len > 1 & k_len == 1) {
        a <- x <- (n < 100 & n != k)
        b <- y <- (n >= 100 & n != k)
        i <- j <- TRUE
    } else if (n_len == 1 & k_len > 1) {
        a <- (n < 100)
        b <- !a
        i <- j <- (n != k)
        x <- if (n < 100) { i } else { NULL }
        y <- if (n >= 100) { i } else { NULL }
    } else {
        stop("Inputs are wrong. n and k must have the same length or be length one.")
    } 
    
    
    # Small n
    nCk[x] <-  lchoose(n[a], k[i])
    
    # Large n indices
    nCk[y] <- n[b]*log(n[b]) - k[j]*log(k[j]) - (n[b] - k[j])*log(n[b] - k[j]) + 
        0.5*(log(n[b]) - log(k[j]) - log(n[b] - k[j]) - log(2*pi))
    
    #     .nCk <- function(n, k) {
    #         n*log(n) - k*log(k) - (n - k)*log(n - k) + 
    #         0.5*(log(n) - log(k) - log(n - k) - log(2*pi))
    #     }
    
    return(nCk)
}

