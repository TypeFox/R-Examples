#' @import graphics
#' @import grDevices
#' @import methods
#' @importFrom stats ave cor dist lm median model.frame model.response na.fail pchisq pnorm prcomp predict predict.lm quantile sd terms 
#' @import utils
NULL

#' Wrapper for several methods to test if a variable is empty
#'
#' This is mainly an internal function but as other dependent packages also
#' use it sometimes and it generally is quite handy to have it is exported for
#' public use.
#' 
#' @param x A variable.
#' @param false_triggers Whether \code{FALSE} should be considered as empty.
#' @return Logical telling if variable is blank.
#' @examples
#' is_blank(NULL)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
is_blank <- function(x, false_triggers=FALSE){
    if(is.function(x)) return(FALSE) # Some of the tests below trigger warnings when used on functions
    is.null(x) ||
    length(x) == 0 ||
    all(is.na(x)) ||
    all(x=="") ||
    (false_triggers && all(!x))
}

#' Detect if modeling results contains multiple procedures
#' 
#' @param result Modeling results, as returned by \code{\link{evaluate}}.
#' @return Logical scalar.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
is_multi_procedure <- function(result){
    belongs_to_method <- function(r){
        !is.null(names(r)) &&
        identical(sort(names(r)),
                  c("error", "importance", "model", "prediction"))
    }
    if(inherits(result, "modeling_result")){
        if(all(sapply(result, belongs_to_method))){
            return(FALSE)
        } else if(all(sapply(result, sapply, belongs_to_method)) &&
                  all(sapply(result, length) == length(result[[1]]))){
            return(TRUE)
        }
    }
    stop("Invalid modeling result.")
}

#' Get names for modeling procedures
#' 
#' @param procedure List of modeling procedures.
#' @return A character vector of suitable non-duplicate names.
#' @author Christofer \enc{Bäcklin}{Backlin}
name_procedure <- function(procedure){
    if(inherits(procedure, "modeling_procedure"))
        stop("`name_procedure` only takes lists of modeling procedures.")
    if(is.null(names(procedure))){
        name <- subtree(procedure, TRUE, "method", error_value=NA, warn=-1)
        name <- ifelse(is.na(name), "model", name)
    } else {
        name <- names(procedure)
    }
    ave(name, name, FUN=function(x){
        if(length(x) == 1) x else sprintf("%s (%i)", x, seq_along(x))
    })
}

#' Replace values with something else
#'
#' @param x Variable containing NAs.
#' @param pattern The values in \code{x} to be replaced. Can also be a
#'   function.
#' @param replacement The value which is to replace the values matching
#'   \code{pattern}.
#' @param invert Whether to fill all values except the ones matching
#'   \code{pattern}.
#' @return An imputed version of \code{x}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @examples
#' fill(1:10, function(x) x %% 2 == 1, 0)
#' na_fill(c(1,2,NA,4,5), 3)
#' @export
fill <- function(x, pattern, replacement, invert=FALSE){
    if(is.function(pattern)){
        x[xor(invert, pattern(x))] <- replacement
    } else {
        x[xor(invert, x %in% pattern)] <- replacement
    }
    x
}
#' @rdname fill
#' @export
na_fill <- function(x, replacement){
    fill(x, is.na, replacement)
}

#' Load a package and offer to install if missing
#' 
#' If running R in interactive mode, the user is prompted for installing
#' missing packages. If running in batch mode an error is thrown.
#'
#' @param pkg Package name.
#' @param reason A status message that informs the user why the package is
#'   needed.
#' @return Nothing
#' @examples
#' nice_require("base", "is required to do anything at all")
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
nice_require <- function(pkg, reason){
    pkg.loaded <- sapply(pkg, requireNamespace)
    if(!all(pkg.loaded) && interactive()){
        pkg <- pkg[!pkg.loaded]
        if(missing(reason))
            reason <- paste(if(length(pkg) > 1) "are" else "is", "required but not installed")
        cat(sprintf("Package%s %s %s. Install now? (Y/n) ",
                    if(length(pkg) > 1) "s" else "",
                    paste("`", pkg, "`", sep="", collapse=", "),
                    reason))
        if(grepl("^(y(es)?)?$", tolower(readline()))){
            install.packages(pkg)
            pkg.loaded <- sapply(pkg, requireNamespace)
            if(!all(pkg.loaded))
                stop("Cannot load package `%s`.", pkg[!pkg.loaded])
        } else {
            stop(sprintf("%s was not installed.",
                         paste("`", pkg, "`", sep="", collapse=", ")))
        }
    }
}

#' Trapezoid rule numerical integration
#' 
#' Only intended for internal use.
#'
#' @param x Integrand.
#' @param y Function values.
#' @return Area under function.
#' @examples
#' x <- seq(0, pi, length.out=100)
#' trapz(x, sin(x))
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
trapz <- function(x,y){
    idx <- order(x)
    x <- x[idx]
    y <- y[idx]
    idx <- !is.na(x) & !is.na(y)
    x <- x[idx]
    y <- y[idx]

    n <- length(x)
    if(n != length(y)) stop("x and y must have same length")
    sum((y[-1]+y[-n])/2 * (x[-1] - x[-n]))
}

#' List all available methods
#' 
#' This function searches all attached packages for methods compatible with the
#' \pkg{emil} framework.
#' 
#' @param pos Location to search in, see \code{\link{ls}}.
#' @return A data frame.
#' @examples
#' list_method()
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
list_method <- function(pos=search()){
    method <- lapply(pos, function(p){
        grep("^(fit|predict|importance)_", ls(p), value=TRUE)
    })
    data.frame(location = rep(pos, sapply(method, length)),
               method = unlist(method)) %>%
        extract_("method", c("plugin", "method"), c("^(fit|predict|importance)_(.*)$")) %>%
        mutate_(plugin = "factor(plugin, c('fit', 'predict', 'importance'))", dummy = TRUE) %>%
        spread_("plugin", "dummy") %>%
        mutate_(fit = "!is.na(fit)", predict = "!is.na(predict)", importance = "!is.na(importance)")
}

#' Convert subsetting vectors.
#'
#' @param y Response vector.
#' @param subset A subsetting vector on arbitrary form.
#' @return A subsetting vector.
#' @examples
#' y <- runif(20)
#' ind <- sort(sample(20, 10))
#' identical(y[ind], y[logical_subset(ind)])
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
positive_integer_subset <- function(y, subset){
    # Not quite sure why this line is needed,
    # but if not present the function will use subset=TRUE if subset is missing!
    # Must be some funky evaluation problem...
    stopifnot(!missing(y) && !missing(subset))
    UseMethod("positive_integer_subset", subset)
}
#' @method positive_integer_subset default
#' @export
positive_integer_subset.default <- function(y, subset){
    seq_along(y)[subset]
}
#' @method positive_integer_subset fold
#' @export
positive_integer_subset.fold <- function(y, subset){
    index_fit(subset)
}

#' @noRd
logical_subset <- function(y, subset){
    stopifnot(!missing(y) && !missing(subset))
    UseMethod("logical_subset", subset)
}
#' @method logical_subset default
#' @export
logical_subset.default <- function(y, subset){
    seq_along(y) %in% seq_along(y)[subset]
}
#' @method logical_subset fold
#' @export
logical_subset.fold <- function(y, subset){
    subset > 0
}

#' Convert factors to logicals
#' 
#' Factors are converted to logical vectors or matrices depending on the number
#' of levels. Ordered factors are converted to matrices where each column
#' represent a level, coded \code{TRUE} for observations that match the level
#' and \code{FALSE} otherwise.
#' Unordered factors are converted in a similar way but coded \code{TRUE} for
#' observations that match the level \emph{or a higher level}.
#' Interpred in words, the star rating example below returns a matrix containing 
#' a column named \dQuote{3 stars} that contains \code{TRUE} for observations
#' with at least three stars and \code{FALSE} for observations with fewer than
#' three stars.
#'
#' @param x Factor.
#' @param base Level to consider as the basis for comparison. Can be either
#'   integer or character.
#'   Note that \code{base = 4} is interpreted as a level named "4",
#'   but \code{base = 4L} is interpreted as the fourth level.
#' @param drop Whether to keep the base level. The base level column never holds
#'   any information that cannot be deduced from the remaining columns.
#' @examples
#' # Binary factor
#' email <- factor(sample(2, 20, TRUE), labels=c("unverified", "verified"))
#' factor_to_logical(email)
#' 
#' # Unordered multi-level factors
#' wine_preferences <- factor(sample(3, 20, TRUE), 
#'                            labels=c("red", "white", "none"))
#' factor_to_logical(wine_preferences, base="none")
#' 
#' fruit <- factor(sample(4, 20, TRUE),
#'                 labels = c("apple", "banana", "cantaloup", "durian"))
#' fruit[sample(length(fruit), 3)] <- NA
#' factor_to_logical(fruit, drop=FALSE)
#' 
#' # Ordered factor
#' rating <- factor(1:5, labels = paste(1:5, "stars"), ordered=TRUE)
#' factor_to_logical(rating)
#' 
#' # Ordered factor with custom base
#' tie_break <- factor(1:5, 
#'                     labels=c("SetAlice", "AdvAlice", "Deuce", "AdvBob", "SetBob"),
#'                     ordered = TRUE)
#' tie_status <- as.data.frame(
#'     factor_to_logical(tie_break, base="Deuce", drop=FALSE)
#' )
#' print(tie_status)
#' tie_break[tie_status$AdvAlice]
#' tie_break[tie_status$SetBob]
#' tie_break[tie_status$Deuce]
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
factor_to_logical <- function(x, base=1L, drop=TRUE){
    stopifnot(length(base) <= 1)
    stopifnot(is.factor(x))
    my_drop <- drop
    if(length(base) == 0){
        if(is.ordered(x)) 
            stop("Ordered factors cannot be converted without a base level.")
        baseL <- 1L
        my_drop <- FALSE
    } else if(is.integer(base)){
        if(base < 1L)
            stop("Invalid base level (non-positive).")
        if(base > nlevels(x))
            stop("Invalid base level (larger than number of levels).")
        baseL <- base
    } else {
        if(is.numeric(base))
            warning("`base` is numeric but not an integer, interpreting it as a character.")
        if(!base %in% levels(x))
            stop("Invalid base level.")
        baseL <- match(base, levels(x))
    }
    levelsL <- 1:nlevels(x) - baseL
    if(is.ordered(x)){
        topleft <- diag(baseL-1)
        topleft[upper.tri(topleft)] <- 1
        bottomleft <- matrix(0, nrow = nlevels(x)-baseL+1, ncol = baseL-1)

        bottomright <- diag(nlevels(x)-baseL)
        bottomright[lower.tri(bottomright)] <- 1L
        topright <- matrix(0, nrow = baseL, ncol = nlevels(x)-baseL)

        newx <- cbind(rbind(topleft, bottomleft), 1L, rbind(topright, bottomright))
    } else {
        newx <- diag(nlevels(x))
    }
    mode(newx) <- "logical"
    colnames(newx) <- levels(x)
    if(my_drop) newx <- newx[, -baseL, drop=FALSE]
    newx[as.integer(x), , drop=my_drop]
}

#' Get the most common value
#' 
#' @param x Vector.
#' @param na.rm Whether to ignore missing values when calculating the mode.
#'   Note that modes may be identified even if \code{x} contains missing values
#'   as long as they are too few to affect the result.
#' @param allow_multiple Controls what is returned if \code{x} contains more
#'   than one mode. If \code{TRUE} all modes are returned, if \code{FALSE}
#'   \code{NA} is returned.
#' @return The most common values or values in \code{x} or \code{NA} if could
#'   not be determined.
#' @examples
#' mode(mtcars$cyl)
#' mode(chickwts$feed)
#' mode(unlist(strsplit("Hello Dolly!", "")))
#' 
#' # Multiple modes
#' mode(iris$Species)
#' 
#' # Missing values
#' x <- rep(1:4, 4)
#' x[2:4] <- NA
#' mode(x)
#' mode(x, na.rm=TRUE)
#' 
#' x <- c(rep(1:3, c(4,2,1)), NA)
#' mode(x, na.rm=FALSE)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
mode <- function(x, na.rm=FALSE, allow_multiple=TRUE){
    if(na.rm){
        count <- table(x, useNA="no")
        count_max <- count == max(count)
    } else {
        count <- table(x, useNA="always")
        na_count <- tail(count, 1)
        count <- count[-length(count)]
        count_max <- count == max(count)
        if(any(count[!count_max] + na_count >= count[count_max]))
            return(NA)
    }
    if(is.factor(x)){
        x_mode <- factor(unname(which(count_max)), seq_len(nlevels(x)), levels(x))
    } else {
        x_mode <- names(count)[count_max]
        if(!is.character(x)) x_mode <- as(x_mode, class(x))
    }
    if(!allow_multiple && length(x_mode) > 1) NA
    else x_mode
}

#' Summarize a potentially long character vector to a compact form
#' 
#' @param type The type of entities listed in \code{names}.
#' @param name Names of entities that should be summarized.
#' @return Character scalar.
#' @examples
#' example_string("letter", LETTERS)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @noRd
example_string <- function(type, names){
    stopifnot(length(names) > 0)
    sprintf("%s `%s`%s", type, names[1],
        if(length(names) > 1)
            sprintf(" and %i other %s%s", length(names)-1, type,
                if(length(names) > 2) "s" else ""))
}

