
# TODO:
# - methods for common "data analysis" classes
#   e.g., "lm" ?
# - complete support for lists
# - NAMESPACE
# - help files

# NOTE that all.equal does something similar, but in a less
# organised way and with less flexible reporting of the
# transformations attempted.

# The "comparison" class
comparison <- function(model, comparison, result, transform,
                       partialM=model, partialC=comparison,
                       partialT=transform) {
    comp <- list(result=result,
                 transform=transform,
                 tM=model,
                 tC=comparison,
                 tMpartial=partialM,
                 tCpartial=partialC,
                 partialTransform=partialT)
    class(comp) <- "comparison"
    comp
}

print.comparison <- function(x, ...) {
    tr <- x$transform
    if (length(tr) == 0)
        cat(x$result, "\n", sep="")
    else
        cat(x$result, "\n",
            paste(paste("  ", tr, sep=""), collapse="\n"),
            "\n", sep="")
}

# Accessors
# Make isTRUE() generic then add a method for "comparison" objects
isTRUE <- function(x) {
    UseMethod("isTRUE")
}

isTRUE.default <- base::isTRUE

isTRUE.comparison <- function(x) {
    x$result
}

transforms <- function(comp) {
    comp$transform
}

# The "multipleComparison" class
multipleComparison <- function(model, comparison, result,
                               detailedResult, transform,
                               partialM=model, partialC=comparison,
                               partialT=transform) {
    # FIXME: this is a data-frame-specific multiple-comparison
    # because it does the check for identical rownames and colnames here
    comp <- list(result=result && all(detailedResult),
                 detailedResult=detailedResult,
                 transform=transform,
                 tM=model,
                 tC=comparison,
                 tMpartial=partialM,
                 tCpartial=partialC,
                 partialTransform=partialT)
    class(comp) <- c("multipleComparison", "comparison")
    comp
}

print.multipleComparison <- function(x, ...) {
    tr <- x$transform
    if (x$result)
        detail <- ""
    else
        detail <- paste("[", paste(x$detailedResult, collapse=", "),
                        "]", sep="")
    if (length(tr) == 0)
        cat(x$result, " ", detail, "\n", sep="")
    else
        cat(x$result, " ", detail, "\n",
            paste(paste("  ", tr, sep=""), collapse="\n"),
            "\n", sep="")
}

###############
# Utilities
###############

coercedMsg <- function(from, to) {
    paste("coerced from <", class(from)[1], "> to <", to, ">", sep="")
}

nocomparisonMsg <- function(from, to) {
    paste("No comparison available between <",
          class(from), "> and <", class(to),
          ">", sep="")
}

# DO NOT restore class OR levels (for factors)
# OR dim OR dimnames (for matrices/arrays)
# If restoring names after a coercion, it is possible
# for length(names) to be wrong (e.g., coerce a
# vector with names attribute to a data frame
# so expecting only 1 name).  In such cases,
# (i.e., if length(names) is wrong) just rep()
# names to the right length to avoid error messages
# (unlikely to produce anything that compares nicely)
restoreAttrs <- function(x, attrs) {
    for (i in names(attrs)) {
        if (i %in% "names" && length(attrs[[i]]) != length(x))
            attrs[[i]] <- rep(attrs[[i]], length.out=length(x))
        if (!(i %in% c("class", "levels", "dim", "dimnames")))
            attr(x, i) <- attrs[[i]]
    }
    x
}

###############
# Check identity
###############
                                                                                compareIdentical <- function(model, comparison, transform=character(), ...) {
    UseMethod("compareIdentical")
}

compareIdentical.default <- function(model, comparison,
                                     transform=character(), ...) {
    result <- identical(model, comparison)
    comparison(model, comparison, result, transform)
}

###############
# Allow "minor" differences
###############

# ALWAYS ASSUME THAT compareEqual() IS CALLED AFTER compareIdentical()
# SO CAN ASSUME THAT model AND comparison ARE NOT IDENTICAL!

# all.equal() reports degree of difference to some extent
# for character, list, POSIXct, ...
# It also looks at attributes
# THIS function is just trying to look for equality
# and allowing for tiny differences
# In many cases, compareEqual() is the same as compareIdentical()
# but an example where they differ is in floating-point values
compareEqual <- function(model, comparison, transform=character(), ...) {
    UseMethod("compareEqual")
}

# Following all.equal()'s lead ...
compareEqual.default <- function(model, comparison, 
                                 transform=character(),
                                 ...) {
    if (is.language(model) || is.function(model) || is.environment(model)) {
        compareEqual.language(model, comparison, transform, ...)
    } else if (is.recursive(model)) {
        transform <- c(transform, "model treated as list")
        # unlist() may not remove class
        compareEqual.list(as.list(unclass(model)),
                          as.list(unclass(comparison)),
                          transform, ...)
    } else {
        transform <- c(transform, nocomparisonMsg(model, comparison))
        comparison(model, comparison, FALSE, transform)        
    }
}

# Following all.equal()'s lead ...
compareEqual.language <- function(model, comparison,
                                 transform=character(),
                                 ...) {
    if (mode(model) == "expression" &&
        mode(comparison) == "expression") {
        transform <- c(transform, "model treated as list")
        compareEqual.list(as.list(unclass(model)),
                          as.list(unclass(comparison)),
                          transform, ...)
    } else if (mode(model) == mode(comparison)) {
        transform <- c(transform, "model treated as character")
        compareEqual(deparse(model), deparse(comparison),
                     transform, ...)
    } else {
        transform <- c(transform, nocomparisonMsg(model, comparison))
        comparison(model, comparison, FALSE, transform)        
    }
}
                                  
compareEqual.logical <- function(model, comparison,
                                 transform=character(),
                                 ...) {
    compareIdentical(model, comparison, transform, ...)
}

compareEqual.numeric <- function(model, comparison,
                                 transform=character(),
                                 round=FALSE,
                                 ...) {
    # NOTE if all.equal() fails, the result may be character
    # HENCE the use of isTRUE()
    result <- isTRUE(all.equal(model, comparison))
    comp <- comparison(model, comparison, result, transform)
    if (is.numeric(comparison)) {
        # Allow rounding
        if (is.function(round)) {
            # Allow function (of single argument)
            if (!result) {
                roundedM <- round(model)
                roundedC <- round(comparison)
                roundedT <- c(transform, "rounded")
                result <- isTRUE(all.equal(roundedM, roundedC))
                comp <- comparison(roundedM, roundedC, result, roundedT,
                                   model, comparison, transform)
            }            
        } else if ((is.numeric(round) || round)) {
            # Allow round=TRUE or round=<number of digits to round>
            if (is.logical(round)) {
                round <- 0
            }
            if (!result) {
                roundedM <- round(model, round)
                roundedC <- round(comparison, round)
                roundedT <- c(transform, "rounded")
                result <- isTRUE(all.equal(roundedM, roundedC))
                comp <- comparison(roundedM, roundedC, result, roundedT,
                                   model, comparison, transform)
            }
        }
    }
    comp
}

# Allow check to drop whitespace at start and end
compareEqual.character <- function(model, comparison,
                                   transform=character(),
                                   ignoreCase=FALSE,
                                   trim=FALSE,
                                   ...) {
    result <- identical(model, comparison)
    comp <- comparison(model, comparison, result, transform)
    if (is.character(comparison)) {
        transM <- model
        transC <- comparison
        trans <- transform
        # Allow whitespace tolerance
        if (trim &&
            length(grep("^ | $", transC)) > 0) {
            transM <- gsub("^ +", "",
                           gsub(" +$", "", transM))
            transC <- gsub("^ +", "",
                           gsub(" +$", "", transC))
            trans <- c(trans, "trimmed whitespace")
        }
        # Allow case-insensitivity
        if (ignoreCase &&
            !identical(transM, transC)) {
            transM <- toupper(transM)
            transC <- toupper(transC)
            trans <- c(trans, "ignored case")
        }
        result <- identical(transM, transC)
        comp <- comparison(transM, transC, result, trans,
                           model, comparison, transform)
    }
    comp
}

compareEqual.factor <- function(model, comparison,
                                transform=character(),
                                dropLevels=FALSE,
                                ignoreLevelOrder=FALSE,
                                ...) {    
    result <- identical(model, comparison)
    comp <- comparison(model, comparison, result, transform)
    if (is.factor(comparison)) {
        transM <- model
        transC <- comparison
        trans <- transform
        # Allow dropping of unused levels
        if (dropLevels &&
            (length(levels(transM)) > length(unique(transM)) ||
             length(levels(transC)) > length(unique(transC)))) {
            transM <- factor(transM)
            transC <- factor(transC)
            trans <- c(trans, "dropped [unused] levels")            
        }
        # Allow order of levels to differ
        if (ignoreLevelOrder &&
            !identical(levels(transM), levels(transC)) &&
            all(transC %in% levels(transM))) {
            transC <- factor(transC, levels=levels(transM))
            trans <- c(trans, "reordered levels")
        }
        result <- identical(transM, transC)
        comp <- comparison(transM, transC, result, trans,
                           model, comparison, transform)
    }
    comp
}

compareEqual.Date <- function(model, comparison,
                              transform=character(),
                              ...) {
    compareIdentical(model, comparison, transform, ...)
}

# 'x' and 'y' are lists of names
# The question is:  for each component of 'x', where
# is an identical component in 'y' ?
dimnamesINdimnames <- function(x, y) {
    lapply(x,
           function(a, b) {
               which(sapply(b, identical, a))
           },
           y)
}

compareEqual.array <- function(model, comparison,
                               transform=character(),
                               ignoreDimOrder=FALSE,
                               ...) {
    # If both model and comparison have the same class
    # and the same number of dimensions and share the
    # same names for their dimensions ...
    if (ignoreDimOrder &&
        identical(class(model), class(comparison)) &&
        length(dim(model)) == length(dim(comparison)) &&
        all(names(dimnames(model)) %in% names(dimnames(comparison)))) {
        dimnamesMatches <- dimnamesINdimnames(dimnames(model),
                                              dimnames(comparison))
        # ... and the dimensions are just in a different order
        if (all(sapply(dimnamesMatches, length) == 1)) {
            # Retain and restore class
            # (because this handles matrices, tables, and arrays)
            cl <- class(comparison)
            # put the dimensions in the same order.
            comparison <- aperm(comparison, unlist(dimnamesMatches))
            class(comparison) <- cl
            transform <- c(transform, "reordered dimensions")
        }
    }
    # Allow for either character or numeric or logical matrix
    if (is.logical(model) && is.logical(comparison)) {
        comp <- compareEqual(as.logical(model),
                             as.logical(comparison),
                             transform, ...)
    } else if (is.numeric(model) && is.numeric(comparison)) {
        comp <- compareEqual(as.numeric(model),
                             as.numeric(comparison),
                             transform, ...)
    } else if (is.character(model) && is.character(comparison)) {
        comp <- compareEqual(as.character(model),
                             as.character(comparison),
                             transform, ...)
    } else {
        comp <- compareIdentical(model, comparison, transform, ...)
        warning("Only logical or numeric or character arrays currently supported")
    }
    # If the underlying data are the same, still need to check
    # that matrix has the correct dimensions and dimnames
    if (isTRUE(comp)) {
        comp$result <- identical(dim(model), dim(comparison)) &&
                       identical(dimnames(model), dimnames(comparison))
    }
    comp    
}

compareEqual.matrix <- compareEqual.array

compareEqual.table <- compareEqual.array

multipleTransform <- function(x, name, partial=FALSE) {
    if (partial)
        trans <- x$partialTransform
    else
        trans <- x$transform
    if (length(trans) > 0)
        paste("[", name, "] ", trans, sep="")
    else
        trans
}

# Names are in common, but names in common are not in same order
similarNames <- function(model, comparison, ignoreNameCase) {
    namesM <- names(model)
    namesC <- names(comparison)
    if (ignoreNameCase) {
        namesM <- toupper(namesM)
        namesC <- toupper(namesC)
    }
    all(namesM %in% namesC) &&
    !identical(order(namesM), order(namesC[1:length(namesM)]))
}

# Put columns/components of comparison that are in common with model
# first and in the model order
# This function assumes that previous checks have been done
# to ensure that ALL model names are in comparison
# (modulo case of names)
reorderComparison <- function(model, comparison, transform, ignoreNameCase,
                              rebuildFun) {
    # Only uppercase if necessary
    if (ignoreNameCase &&
        !all(names(model) %in% names(comparison))) {
        names(model) <- toupper(names(model))
        names(comparison) <- toupper(names(comparison))
        transform <- c(transform, "renamed")
    }
    namesM <- names(model)
    namesC <- names(comparison)
    namesInCommon <- namesC %in% namesM
    newComp <- comparison[match(namesM, namesC)]
    theRest <- comparison[!namesInCommon]
    list(model=model,
         comparison=rebuildFun(newComp, theRest),
         transform=transform)
}

compareEqualDF <- function(model, comparison, compareFun,
                           transform,
                           # There are certain transforms that should not
                           # be considered for individual columns
                           # so intercept them here ...
                           shorten=FALSE,
                           ignoreOrder=FALSE,
                           ignoreColOrder=FALSE,
                           ignoreNameCase=ignoreNameCase,
                           ...) {
    # If model and comparison have different number of
    # columns, we can give up straight away
    if (length(model) != length(comparison)) {
        return(comparison(model, comparison, FALSE, transform))
    }
    # The partial transform is just the original (incoming)
    # version (i.e., ANY transformation of the columns is NOT
    # recorded in the partial result, which allows these to be
    # tried again following subsequent transformations of the
    # overall data frame).
    partialModel <- model
    partialComparison <- comparison
    partialTransform <- transform
    # Reorder columns by name (if allowed and necessary)
    if (ignoreColOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  cbind)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered columns")
    }
    # Only bother with columns in common
    multipleResult <- mapply(compareFun,
                             model,
                             comparison[seq_along(model)],
                             # Begin each col check with pristine transform
                             # There are certain transforms that should not
                             # be considered for individual columns
                             # so explicitly disallow them here
                             MoreArgs=list(shorten=FALSE, ignoreOrder=FALSE,
                               ignoreNameCase=ignoreNameCase,
                               ...),
                             SIMPLIFY=FALSE)
    detailedResult <- sapply(multipleResult, "[[", "result")
    # Build a fully-transformed and a partially-transformed
    # version of the result
    # Check for ANY transforms
    if (any(sapply(multipleResult,
                   function(x) { length(x$transform) > 0 }))) {
        # model may have had its names stripped by this point,
        # in which case, just use column numbers
        colNames <- names(multipleResult)
        if (is.null(colNames))
            colNames <- seq_along(model)
        transform <- c(transform,
                       unlist(mapply(multipleTransform,
                                     multipleResult,
                                     colNames)))
    }
    # Replace each column with its tM version
    for (i in seq_along(model)) {
        model[[i]] <- multipleResult[[i]]$tM
    }
    # Replace each column (that is in common with model) with its tC version
    for (i in seq_along(model)) {
        comparison[[i]] <- multipleResult[[i]]$tC
    }
    multipleComparison(model, comparison,
                       identical(attributes(model), attributes(comparison)),
                       detailedResult, transform,
                       partialModel, partialComparison,
                       partialTransform)
}

compareEqual.data.frame <- function(model, comparison,
                                    transform=character(),
                                    ignoreColOrder=FALSE,
                                    ignoreNameCase=FALSE,
                                    ...,
                                    recurseFun=compareEqual) {
    # comparison needs to be a data frame first, then ...
    if (is.data.frame(comparison)) {
        # Two data frames are equal if all columns are equal ...
        comp <- compareEqualDF(model, comparison, recurseFun, transform,
                               ignoreColOrder=ignoreColOrder,
                               ignoreNameCase=ignoreNameCase,
                               ...)
    } else {
        comp <- comparison(model, comparison, FALSE, transform)
    }
    comp
}

compareListComponents <- function(model, comparison,
                                  transform=character(),
                                  partialM, partialC, partialT,
                                  ignoreComponentOrder=FALSE,
                                  ignoreNameCase=FALSE,
                                  ...,
                                  recurseFun=compareEqual) {
        # Compare all components
        multipleResult <-
            mapply(recurseFun,
                   model,
                   # Only makes sense to compare components in common
                   comparison[seq_along(model)],
                   # transform is not specified;
                   # each component comparison starts with blank transform
                   MoreArgs=list(ignoreComponentOrder=ignoreComponentOrder,
                     ignoreNameCase=ignoreNameCase,
                     ...),
                   SIMPLIFY=FALSE)
        
        detailedResult <- sapply(multipleResult, "[[", "result")
        # Build a fully-transformed and a partially-transformed
        # version of the result
        # Check for ANY transforms
        if (any(sapply(multipleResult,
                       function(x) { length(x$transform) > 0 }))) {
            # model may not have names,
            # in which case, just use column numbers
            cNames <- names(multipleResult)
            if (is.null(cNames))
                cNames <- seq_along(model)
            transform <- c(transform,
                           unlist(mapply(multipleTransform,
                                         multipleResult,
                                         cNames)))
        }
        # Replace each column with its tM version
        for (i in seq_along(model)) {
            model[[i]] <- multipleResult[[i]]$tM
        }
        # Replace each column (that is in common with mode) with its tC version
        for (i in seq_along(model)) {
            comparison[[i]] <- multipleResult[[i]]$tC
        }
        # Overall result is TRUE if all children are equal
        # AND all attributes of the model and comparison are the same
        multipleComparison(model, comparison,
                           identical(attributes(model),
                                     attributes(comparison)),
                           detailedResult, transform,
                           partialM, partialC,
                           partialT)
}

compareEqual.list <- function(model, comparison,
                              transform=character(),
                              ignoreComponentOrder=FALSE,
                              ignoreNameCase=FALSE,
                              ...,
                              recurseFun=compareEqual) {
    if (length(model) == 0) { 
        # Special case of model zero-length
        if (length(comparison) == 0) {
            comp <- comparison(model, comparison,
                               identical(attributes(model),
                                         attributes(comparison)),
                               transform)
        } else {
            comp <- comparison(model, comparison, FALSE, transform)
        }
    } else if (length(comparison) > 0 &&
               is.list(comparison)) {
        # No matter what happens, the partial result is
        # what was originally passed in
        partialM <- model
        partialC <- comparison
        partialT <- transform
        # Allow for reordering by name (possibly ignoring name case)
        if (ignoreComponentOrder &&
            !is.null(names(model)) &&
            similarNames(model, comparison, ignoreNameCase)) {
            temp <- reorderComparison(model, comparison,
                                      transform, ignoreNameCase,
                                      c)
            model <- temp$model
            comparison <- temp$comparison
            transform <- c(temp$transform, "reordered components")
        }
        comp <- compareListComponents(model, comparison, transform,
                                      partialM, partialC, partialT,
                                      ignoreComponentOrder,
                                      ignoreNameCase,
                                      ...,
                                      recurseFun=recurseFun)
    } else {
        comp <- comparison(model, comparison, FALSE, transform)
    }
    comp
}

# Do a compareIdentical()
# If that fails, do a compareEqual (if allowed)
# Add untransformed model, comparison, transform to result
# so that anyone calling can ignore any transforms done by compareEqual()
same <- function(model, comparison, transform, equal, ...) {
    comp <- compareIdentical(model, comparison, transform, ...)
    if (!comp$result && equal) {
        comp <- compareEqual(model, comparison, transform, ...)
    }
    if (inherits(comp, "multipleComparison"))
        multipleComparison(comp$tM, comp$tC, comp$result,
                           comp$detailedResult,
                           comp$transform,
                           model, comparison, transform)
    else
        comparison(comp$tM, comp$tC, comp$result,
                   comp$transform,
                   model, comparison, transform)
}

###############
# Allow type coercion
###############

compareCoerce <- function(model, comparison,
                          transform=character(),
                          equal=TRUE, ...) {
    UseMethod("compareCoerce")
}

# Following all.equal()'s lead ...
compareCoerce.default <- function(model, comparison, 
                                  transform=character(),
                                  ...) {
    if (is.language(model) || is.function(model) || is.environment(model)) {
        compareCoerce.language(model, comparison, transform, ...)
    } else if (is.recursive(model)) {
        transform <- c(transform, "model treated as list")
        compareCoerce.list(as.list(unclass(model)),
                           as.list(unclass(comparison)),
                           transform, ...)
    } else {
        transform <- c(transform, nocomparisonMsg(model, comparison))
        comparison(model, comparison, FALSE, transform)        
    }
}

# Following all.equal()'s lead ...
compareCoerce.language <- function(model, comparison,
                                   transform=character(),
                                   ...) {
    if (mode(model) == "expression" &&
        mode(comparison) == "expression") {
        transform <- c(transform, "model treated as list")
        compareCoerce.list(as.list(unclass(model)),
                           as.list(unclass(comparison)),
                           transform, ...)
    } else if (mode(model) == mode(comparison)) {
        transform <- c(transform, "model treated as character")
        compareCoerce(deparse(model), deparse(comparison),
                      transform, ...)
    } else {
        transform <- c(transform, nocomparisonMsg(model, comparison))
        comparison(model, comparison, FALSE, transform)        
    }
}
                                  
compareCoerce.logical <- function(model, comparison,
                                  transform=character(),
                                  equal=TRUE, ...) {
    if (is.logical(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        # Retain all attributes so we can check for those separately
        attrs <- attributes(comparison)
        coerced <- as.logical(comparison)
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "logical"))
        comp <- same(model, coerced, transform, equal, ...)
    }
    comp    
}

compareCoerce.integer <- function(model, comparison,
                                  transform=character(),
                                  equal=TRUE, ...) {
    # If we already have an integer, do not record any transformation
    # (possibly does a redundant check for identical or equal
    #  if those check have been done previously)
    if (is.integer(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        # Retain all attributes so we can check for those separately
        attrs <- attributes(comparison)
        # Special case for coming from a factor
        if (is.factor(comparison)) {
            coerced <- as.integer(as.character(comparison))
        } else {
            coerced <- try(suppressWarnings(as.integer(comparison)))
        }
        if (inherits(coerced, "try-error")) {
            transform <- c(transform, nocomparisonMsg(model, comparison))
            comp <- comparison(model, comparison, FALSE, transform)
        } else {
            coerced <- restoreAttrs(coerced, attrs)
            transform <- c(transform, coercedMsg(comparison, "integer"))
            comp <- same(model, coerced, transform, equal, ...)
        }
    }
    comp    
}

compareCoerce.numeric <- function(model, comparison,
                                  transform=character(),
                                  equal=TRUE, ...) {
    # The comparison has to be an atomic or a factor, otherwise just fail
    if (!(is.atomic(comparison) || is.factor(comparison))) {
        comp <- comparison(model, comparison, FALSE, transform)
    # If we already have a numeric, do not record any transformation
    # (possibly does a redundant check for identical or equal
    #  if those check have been done previously)
    # NOTE that integer vector is a special case
    } else if (is.numeric(comparison) && !is.integer(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        # Retain all attributes so we can check for those separately
        attrs <- attributes(comparison)
        # Special case for coming from a factor
        if (is.factor(comparison)) {
            coerced <- as.numeric(as.character(comparison))
        } else {
            coerced <- suppressWarnings(as.numeric(comparison))
        }
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "numeric"))
        comp <- same(model, coerced, transform, equal, ...)
    }
    comp    
}

compareCoerce.character <- function(model, comparison,
                                    transform=character(),
                                    equal=TRUE, ...) {
    # If we already have character, do not record any transformation
    # (possibly does a redundant check for identical or equal
    #  if those check have been done previously)
    if (is.character(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        attrs <- attributes(comparison)
        # Retain all attributes so we can check for those separately
        coerced <- as.character(comparison)
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "character"))
        comp <- same(model, coerced, transform, equal, ...)
    }
    comp    
}

compareCoerce.factor <- function(model, comparison,
                                 transform=character(),
                                 equal=TRUE, ...) {
    # If we already have factor, do not record any transformation
    # (possibly does a redundant check for identical or equal
    #  if those check have been done previously)
    if (is.factor(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        if (!is.atomic(comparison)) {
            transform <- c(transform, nocomparisonMsg(model, comparison))
            comp <- comparison(model, comparison, FALSE, transform)            
        } else {
            attrs <- attributes(comparison)
            # Retain all attributes so we can check for those separately
            coerced <- as.factor(comparison)
            coerced <- restoreAttrs(coerced, attrs)
            transform <- c(transform, coercedMsg(comparison, "factor"))
            comp <- same(model, coerced, transform, equal, ...)
        }
    }
    comp    
}

compareCoerce.Date <- function(model, comparison,
                                  transform=character(),
                                  equal=TRUE, ...) {
    if (inherits(comparison, "Date")) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        # Retain all attributes so we can check for those separately
        attrs <- attributes(comparison)
        coerced <- as.Date(comparison)
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "Date"))
        comp <- same(model, coerced, transform, equal, ...)
    }
    comp    
}

compareCoerce.array <- function(model, comparison,
                                transform=character(),
                                equal=TRUE, ...) {
    if (is.array(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        attrs <- attributes(comparison)
        # Retain all attributes so we can check for those separately
        coerced <- as.array(comparison)
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "array"))
        comp <- same(model, coerced, transform, equal, ...)        
    }
    comp
}

compareCoerce.matrix <- function(model, comparison,
                                 transform=character(),
                                 equal=TRUE, ...) {
    if (is.matrix(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        attrs <- attributes(comparison)
        # Retain all attributes so we can check for those separately
        coerced <- as.matrix(comparison)
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "matrix"))
        comp <- same(model, coerced, transform, equal, ...)        
    }
    comp
}

compareCoerce.table <- function(model, comparison,
                                transform=character(),
                                equal=TRUE, ...) {
    if (is.table(comparison)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        attrs <- attributes(comparison)
        # Retain all attributes so we can check for those separately
        coerced <- as.table(comparison)
        coerced <- restoreAttrs(coerced, attrs)
        transform <- c(transform, coercedMsg(comparison, "table"))
        comp <- same(model, coerced, transform, equal, ...)        
    }
    comp
}

# This differs from compareEqualDF in two ways:
# - it does NOT have a 'compareFun' argument,
#   we ONLY call compareCoerce() recursively
# - any transformations on the columns are persistent
compareCoerceDF <- function(model, comparison, transform,
                            ...) {
    # If model and comparison have different number of
    # columns, we can give up straight away
    if (ncol(model) != ncol(comparison)) {
        return(comparison(model, comparison, FALSE, transform))
    }
    # Only bother with columns in common
    multipleResult <- mapply(compareCoerce,
                             model,
                             comparison[seq_along(model)],
                             # Begin each col check with pristine transform
                             MoreArgs=list(...),
                             SIMPLIFY=FALSE)
    detailedResult <- sapply(multipleResult, "[[", "result")
    # Build a fully-transformed and a partially-transformed
    # version of the result
    # The partial transform retains any transformations
    # by compareCoerce(), but not any from the compareEqual()
    # within the compareCoerce()
    # (i.e., build partial result out of partial results
    #  from compareCoerce() on columns)
    partialModel <- model
    partialComparison <- comparison
    partialTransform <- transform
    # Check for ANY transforms
    if (any(sapply(multipleResult,
                   function(x) { length(x$transform) > 0 }))) {
        # model may have had its names stripped by this point,
        # in which case, just use column numbers
        colNames <- names(multipleResult)
        if (is.null(colNames))
            colNames <- seq_along(model)
        transform <- c(transform,
                       unlist(mapply(multipleTransform,
                                     multipleResult,
                                     colNames)))
        partialTransform <- c(partialTransform,
                              unlist(mapply(multipleTransform,
                                            multipleResult,
                                            colNames, partial=TRUE)))
    }
    # Replace each column with its tM version
    for (i in seq_along(model)) {
        model[[i]] <- multipleResult[[i]]$tM
    }
    for (i in seq_along(model)) {
        partialModel[[i]] <- multipleResult[[i]]$tMpartial
    }
    # Replace each column (that is in common with model) with its tC version
    for (i in seq_along(model)) {
        comparison[[i]] <- multipleResult[[i]]$tC
    }
    for (i in seq_along(model)) {
        partialComparison[[i]] <- multipleResult[[i]]$tCpartial
    }
    multipleComparison(model, comparison,
                       identical(attributes(model), attributes(comparison)),
                       detailedResult, transform,
                       partialModel, partialComparison,
                       partialTransform)
}

compareCoerce.data.frame <- function(model, comparison,
                                     transform=character(),
                                     equal=TRUE,
                                     ignoreColOrder=FALSE,
                                     ignoreNameCase=FALSE,
                                     ...) {
    # If we don't have a data frame, try to coerce
    if (!is.data.frame(comparison)) {
        attrs <- attributes(comparison)
        # Retain all attributes so we can check for those separately
        comparison <- as.data.frame(comparison)
        comparison <- restoreAttrs(comparison, attrs)
        transform <- c(transform, coercedMsg(comparison, "data frame"))
    }
    # Order columns before trying coercion
    if (ignoreColOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        # Permanently modify model and comparison
        # (to avoid other code redoing the transform)
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  cbind)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered columns")
    } 
    # Now see if any columns can be coerced
    comp <- compareCoerceDF(model, comparison, transform, ...)
    comp
}

coerceListComponents <- function(model, comparison,
                                 transform=character(),
                                 ignoreComponentOrder=FALSE,
                                 ignoreNameCase=FALSE,
                                 ...) {
    # If comparison has fewer components than model, we can give
    # up straight away
    if (length(model) > length(comparison)) {
        return(comparison(model, comparison, FALSE, transform))
    }
    # Compare all components
    multipleResult <-
        mapply(compareCoerce,
               model,
               # Only makes sense to compare components in common
               comparison[seq_along(model)],
               # transform is not specified;
               # each component comparison starts with blank transform
               MoreArgs=list(ignoreComponentOrder=ignoreComponentOrder,
                 ignoreNameCase=ignoreNameCase,
                 ...),
               SIMPLIFY=FALSE)
    detailedResult <- sapply(multipleResult, "[[", "result")
    # Build a fully-transformed and a partially-transformed
    # version of the result
    # The partial transform retains any transformations
    # by compareCoerce(), but not any from the compareEqual()
    # within the compareCoerce()
    # (i.e., build partial result out of partial results
    #  from compareCoerce() on columns)
    partialModel <- model
    partialComparison <- comparison
    partialTransform <- transform
    # Check for ANY transforms
    if (any(sapply(multipleResult,
                   function(x) { length(x$transform) > 0 }))) {
        # model may have had its names stripped by this point,
        # in which case, just use column numbers
        cNames <- names(multipleResult)
        if (is.null(cNames))
            cNames <- seq_along(model)
        transform <- c(transform,
                       unlist(mapply(multipleTransform,
                                     multipleResult,
                                     cNames)))
        partialTransform <- c(partialTransform,
                              unlist(mapply(multipleTransform,
                                            multipleResult,
                                            cNames, partial=TRUE)))
    }
    # Replace each column with its tM version
    for (i in seq_along(model)) {
        model[[i]] <- multipleResult[[i]]$tM
    }
    for (i in seq_along(model)) {
        partialModel[[i]] <- multipleResult[[i]]$tMpartial
    }
    # Replace each column (that is in common with model) with its tC version
    for (i in seq_along(model)) {
        comparison[[i]] <- multipleResult[[i]]$tC
    }
    for (i in seq_along(model)) {
        partialComparison[[i]] <- multipleResult[[i]]$tCpartial
    }
    # Overall result is TRUE if all children are equal
    # AND all attributes of the model and comparison are the same
    multipleComparison(model, comparison,
                       identical(attributes(model),
                                 attributes(comparison)),
                       detailedResult, transform,
                       partialModel, partialComparison,
                       partialTransform)
}

compareCoerce.list <- function(model, comparison,
                               transform=character(),
                               equal=TRUE,
                               ignoreComponentOrder=FALSE,
                               ignoreNameCase=FALSE,
                               ...) {
    if (length(model) == 0) { 
        # Special case of model zero-length
        if (length(comparison) == 0) {
            comp <- comparison(model, comparison,
                               identical(attributes(model),
                                         attributes(comparison)),
                               transform)
        } else {
            comp <- comparison(model, comparison, FALSE, transform)
        }
    } else if (length(comparison) > 0) {
        if (!inherits(comparison, "list")) {
            transform <- c(transform, coercedMsg(comparison, "list"))
            # Retain all attributes so we can check for those separately
            attrs <- attributes(comparison)
            comparison <- as.list(comparison)
            class(comparison) <- NULL
            comparison <- restoreAttrs(comparison, attrs)
        }
        # Allow for reordering by name (possibly ignoring name case)
        if (ignoreComponentOrder &&
            !is.null(names(model)) &&
            similarNames(model, comparison, ignoreNameCase)) {
            temp <- reorderComparison(model, comparison,
                                      transform, ignoreNameCase,
                                      c)
            model <- temp$model
            comparison <- temp$comparison
            transform <- c(temp$transform, "reordered components")
        }
        comp <- coerceListComponents(model, comparison, transform,
                                     ignoreComponentOrder,
                                     ignoreNameCase,
                                     ...)
    } else {
        comp <- comparison(model, comparison, FALSE, transform)
    }
    comp
}

###############
# Allow shortening
###############

compareShorten <- function(model, comparison,
                        transform=character(),
                        equal=TRUE, ...) {
    UseMethod("compareShorten")
}

compareShorten.default <- function(model, comparison,
                                   transform=character(),
                                   equal=TRUE, ...) {
    nM <- length(model)
    nC <- length(comparison)
    if (nM != nC) {
        if (nM > nC) {
            model <- model[1:nC, drop=FALSE]
            transform <- c(transform, "shortened model")
        } else {
            comparison <- comparison[1:nM, drop=FALSE]
            transform <- c(transform, "shortened comparison")
        }
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        comp <- same(model, comparison, transform, equal, ...)
    }
    comp
}

dropDimensions <- function(model, comparison, transform) {
    UseMethod("dropDimensions")
}

# For a matrix, comparison must be forced to 2D
dropDimensions.matrix <- function(model, comparison, transform) {
    dimM <- dim(model)
    dimC <- dim(comparison)
    ndimM <- length(dimM)
    ndimC <- length(dimC)
    if (ndimM > ndimC) {
        comparison <- matrix(comparison, ncol=1)
        transform <- c(transform, "added comparison dimension")
    } else {
        # Keep the first two dimensions
        dimnamesC <- dimnames(comparison)
        comparison <- array(comparison, dimC[1:2])
        dimnames(comparison) <- dimnamesC[1:2]
        transform <- c(transform, "dropped comparison dimension(s)")        
    }
    list(model=model, comparison=comparison,
         transform=transform)
}

# For an array, either model or comparison will drop dimensions to match
dropDimensions.array <- function(model, comparison, transform) {
    dimM <- dim(model)
    dimC <- dim(comparison)
    ndimM <- length(dimM)
    ndimC <- length(dimC)
    dimnamesM <- dimnames(model)
    dimnamesC <- dimnames(comparison)
    if (ndimM > ndimC) {
        model <- array(model, dimM[1:ndimC])
        dimnames(model) <- dimnamesM[1:ndimC]
        transform <- c(transform, "dropped model dimension(s)")
    } else {
        comparison <- array(comparison, dimC[1:ndimM])
        dimnames(comparison) <- dimnamesC[1:ndimM]
        transform <- c(transform, "dropped comparison dimension(s)")
    }
    list(model=model, comparison=comparison,
         transform=transform)
}

# When dropping table dimensions, "collapse" instead
# (i.e., sum across the dimension being dropped)
dropDimensions.table <- function(model, comparison, transform) {
    dimM <- dim(model)
    dimC <- dim(comparison)
    ndimM <- length(dimM)
    ndimC <- length(dimC)
    dimnamesM <- dimnames(model)
    dimnamesC <- dimnames(comparison)
    if (ndimM > ndimC) {
        model <- apply(model, 1:ndimC, sum)
        transform <- c(transform, "collapsed model dimension(s)")
    } else {
        comparison <- apply(comparison, 1:ndimM, sum)
        transform <- c(transform, "collapsed comparison dimension(s)")
    }
    list(model=model, comparison=comparison,
         transform=transform)    
}

# Shortening an array only drops entire dimensions
compareShorten.array <- function(model, comparison,
                                 transform=character(),
                                 equal=TRUE, ...) {
    dimM <- dim(model)
    dimC <- dim(comparison)
    ndimM <- length(dimM)
    ndimC <- length(dimC)
    if (!identical(dimM, dimC)) {
        # Special case of comparison has a zero dim
        # (i.e., it is zero-extent)
        # which means that both model and comparison need
        # to be reduced to zero.
        # This will terminate any further comparisons because
        # model and comparison become identical.
        if (any(dimC == 0)) {
            model <- comparison <- vector()
            transform <- c(transform, "model reduced to zero extent")
        } else {
            # Drop any extra dimensions
            if (ndimM != ndimC) {
                dropResult <- dropDimensions(model, comparison, transform)
                model <- dropResult$model
                comparison <- dropResult$comparison
                transform <- dropResult$transform
                dimM <- dim(model)
                dimC <- dim(comparison)
                ndimM <- length(dimM)
                ndimC <- length(dimC)
            }
        }
    }
    same(model, comparison, transform, equal, ...)
}

compareShorten.matrix <- compareShorten.array

compareShorten.table <- compareShorten.array

compareShorten.data.frame <- function(model, comparison,
                                      transform=character(),
                                      equal=TRUE,
                                      colsOnly=TRUE,
                                      ignoreColOrder=FALSE,
                                      ignoreNameCase=FALSE,
                                      ...) {
    # Order columns before shortening
    if (ignoreColOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  cbind)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered columns")
    }
    # Drop cols
    nM <- length(model)
    nC <- length(comparison)
    if (nM != nC) {
        if (nM > nC) {
            model <- model[, 1:nC, drop=FALSE]
            transform <- c(transform, "shortened model")
        } else {
            comparison <- comparison[, 1:nM, drop=FALSE]
            transform <- c(transform, "shortened comparison")
        }
    }
    # Also allow for checking number of rows
    if (!colsOnly) {
        nrM <- nrow(model)
        nrC <- nrow(comparison)
        # If the number of rows are already the same, don't transform
        if (nrM != nrC) {
            if (nrM > nrC) {
                model <- model[1:nrC, , drop=FALSE]
                transform <- c(transform, "shortened model rows")
            } else {
                comparison <- comparison[1:nrM, , drop=FALSE]
                transform <- c(transform, "shortened comparison rows")
            }
        }
    }
    same(model, comparison, transform, equal, ...)
}

compareShorten.list <- function(model, comparison,
                                transform=character(),
                                equal=TRUE,
                                ignoreComponentOrder=FALSE,
                                ignoreNameCase=FALSE,
                                ...) {
    # Order components before shortening
    if (ignoreComponentOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  c)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered components")
    }
    # Drop components
    nM <- length(model)
    nC <- length(comparison)
    if (nM > nC) {
        model <- model[1:nC]
        transform <- c(transform, "shortened model")        
    }
    if (nM < nC) {
        comparison <- comparison[1:nM]
        transform <- c(transform, "shortened comparison")
    }
    same(model, comparison, transform, equal, ...)
}

###############
# Ignore order
###############

compareIgnoreOrder <- function(model, comparison,
                               transform=character(),
                               equal=TRUE, ...) {
    UseMethod("compareIgnoreOrder")
}

# More complex objects will need to write specific methods
compareIgnoreOrder.default <- function(model, comparison,
                                       transform=character(),
                                       equal=TRUE, ...) {
    if (is.language(model) || is.function(model) || is.environment(model)) {
        comp <- compareIgnoreOrder.language(model, comparison, transform, ...)
    } else if (is.recursive(model)) {
        transform <- c(transform, "model treated as list")
        comp <- compareIgnoreOrder.list(as.list(unclass(model)),
                                        as.list(unclass(comparison)),
                                        transform, ...)
    } else if (!is.atomic(comparison)) {
        # Model is vector but comparison is not
        transform <- c(transform, nocomparisonMsg(model, comparison))
        comp <- comparison(model, comparison, FALSE, transform)
    } else { 
        # Handle any sort of vector ?
        # Only try reordering if it is necessary
        modelOrder <- order(model)
        comparisonOrder <- order(comparison)
        if (!identical(modelOrder, comparisonOrder)) {        
            # Retain all attributes so we can check for those separately
            attrsM <- attributes(model)
            attrsC <- attributes(comparison)
            sortedM <- sort(model)
            sortedC <- sort(comparison)
            sortedM <- restoreAttrs(sortedM, attrsM)
            sortedC <- restoreAttrs(sortedC, attrsC)
            transform <- c(transform, "sorted")
            comp <- same(sortedM, sortedC, transform, equal, ...)
        } else {
            comp <- same(model, comparison, transform, equal, ...)
        }
    }
    comp    
}

# Following all.equal()'s lead ...
compareIgnoreOrder.language <- function(model, comparison,
                                        transform=character(),
                                        ...) {
    if (mode(model) == "expression" &&
        mode(comparison) == "expression") {
        transform <- c(transform, "model treated as list")
        compareIgnoreOrder.list(as.list(unclass(model)),
                                as.list(unclass(comparison)),
                                transform, ...)
    } else if (mode(model) == mode(comparison)) {
        transform <- c(transform, "model treated as character")
        compareIgnoreOrder(deparse(model), deparse(comparison),
                           transform, ...)
    } else {
        transform <- c(transform, nocomparisonMsg(model, comparison))
        comparison(model, comparison, FALSE, transform)        
    }
}
                                  
# Put the each component of dimnames in alphabetical order
compareIgnoreOrder.array <- function(model, comparison,
                                     transform=character(),
                                     equal=TRUE, ...) {
    dimnamesM <- dimnames(model)
    dimnamesC <- dimnames(comparison)
    if (length(dimnamesM) == length(dimnamesC)) {
        if (!all(unlist(mapply(identical, dimnamesM, dimnamesC,
                               SIMPLIFY=FALSE)))) {
            dimOrderM <- lapply(dimnamesM, order)
            dimOrderC <- lapply(dimnamesC, order)
            model <- do.call("[", c(list(model), dimOrderM))
            comparison <- do.call("[", c(list(comparison), dimOrderC))
        }
    }
    same(model, comparison, transform, equal, ...)
}

compareIgnoreOrder.data.frame <- function(model, comparison,
                                          transform=character(),
                                          equal=TRUE, 
                                          ignoreColOrder=FALSE,
                                          ignoreNameCase=FALSE,
                                          ...) {
    # Order columns before ordering rows
    if (ignoreColOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  cbind)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered columns")
    } 
    # Only try reordering if it is necessary
    modelOrder <- do.call("order", model)
    comparisonOrder <- do.call("order", comparison)
    if (!identical(modelOrder, comparisonOrder)) {
        # Try ordering the rows
        model <- model[modelOrder, , drop=FALSE]    
        comparison <- comparison[comparisonOrder, , drop=FALSE]
        transform <- c(transform, "sorted")
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        comp <- same(model, comparison, transform, equal, ...)
    }
    comp
}

compareIgnoreOrder.list <- function(model, comparison,
                                    transform=character(),
                                    equal=TRUE,
                                    ignoreNameCase=FALSE,
                                    ...) {
    # Order components
    if (!is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  c)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "sorted by name")
    }
    same(model, comparison, transform, equal, ...)
}

###############
# Ignore name case
###############

compareIgnoreNameCase <- function(model, comparison,
                                  transform=character(),
                                  equal=TRUE, ...) {
    UseMethod("compareIgnoreNameCase")
}

# Handle any sort of vector ?
# More complex objects will need to write specific methods
compareIgnoreNameCase.default <- function(model, comparison,
                                          transform=character(),
                                          equal=TRUE, ...) {
    namesM <- names(model)
    namesC <- names(comparison)
    # If the names are already identical, don't transform
    if (identical(namesM, namesC)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        names(model) <- toupper(namesM)
        names(comparison) <- toupper(namesC)
        transform <- c(transform, "renamed")
        comp <- same(model, comparison, transform, equal, ...)
    }
    comp    
}

compareIgnoreNameCase.data.frame <- function(model, comparison,
                                             transform=character(),
                                             equal=TRUE, 
                                             colsOnly=TRUE,
                                             ignoreColOrder=FALSE,
                                             ignoreNameCase=FALSE,
                                             ...) {
    # Order columns before comparing name case
    if (ignoreColOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  cbind)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered columns")
    } 
    # The default method handles the column names
    comp <- NextMethod()
    # Also allow for checking row names
    if (!colsOnly) {
        model <- comp$tM
        comparison <- comp$tC
        transform <- comp$transform
        rownamesM <- rownames(model)
        rownamesC <- rownames(comparison)
        # If the rownames are already identical, don't transform
        if (!identical(rownamesM, rownamesC)) {
            rownames(model) <- toupper(rownamesM)
            rownames(comparison) <- toupper(rownamesC)
            transform <- c(transform, "renamed rows")
            comp <- same(model, comparison, transform, equal, ...)
        }
    }
    comp    
}

compareIgnoreNameCase.list <- function(model, comparison,
                                       transform=character(),
                                       equal=TRUE,
                                       ignoreComponentOrder=FALSE,
                                       ...) {
    # Order components 
    if (ignoreComponentOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase=TRUE)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase=TRUE,
                                  c)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered components")
    }
    same(model, comparison, transform, equal, ...)
}

###############
# Ignore names
###############

compareIgnoreNames <- function(model, comparison,
                               transform=character(),
                               equal=TRUE, ...) {
    UseMethod("compareIgnoreNames")
}

# Handle any sort of vector ?
# More complex objects will need to write specific methods
compareIgnoreNames.default <- function(model, comparison,
                                       transform=character(),
                                       equal=TRUE, ...) {
    namesM <- names(model)
    namesC <- names(comparison)
    # If neither object has names or the names are identical, don't transform
    if (identical(namesM, namesC)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        model <- unname(model)
        comparison <- unname(comparison)
        transform <- c(transform, "dropped names")
        comp <- same(model, comparison, transform, equal, ...)
    }
    comp    
}

compareIgnoreNames.data.frame <- function(model, comparison,
                                          transform=character(),
                                          equal=TRUE, 
                                          colsOnly=TRUE,
                                          ignoreColOrder=FALSE,
                                          ignoreNameCase=FALSE,
                                          ...) {
    # Order columns before dropping names
    if (ignoreColOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  cbind)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered columns")
    } 
    # The default method handles the column names
    comp <- NextMethod()
    # Also allow for dropping row names
    # Only bother to check if the number of rows are the same
    if (!colsOnly &&
        dim(model)[1] == dim(comparison)[1]) {
        model <- comp$tM
        comparison <- comp$tC
        transform <- comp$transform
        rownamesM <- rownames(model)
        rownamesC <- rownames(comparison)
        # If the rownames are already identical, don't transform
        if (!identical(rownamesM, rownamesC)) {
            rownames(model) <- 1:(dim(model)[1])
            rownames(comparison) <- 1:(dim(model)[1])
            transform <- c(transform, "dropped row names")
            comp <- same(model, comparison, transform, equal, 
                         colsOnly=colsOnly, ...)
        }
    }
    comp    
}

compareIgnoreNames.list <- function(model, comparison,
                                    transform=character(),
                                    equal=TRUE,
                                    ignoreComponentOrder=FALSE,
                                    ignoreNameCase=FALSE,
                                    ...) {
    # Order components (that are in common)
    if (ignoreComponentOrder &&
        !is.null(names(model)) &&
        similarNames(model, comparison, ignoreNameCase)) {
        temp <- reorderComparison(model, comparison,
                                  transform, ignoreNameCase,
                                  c)
        model <- temp$model
        comparison <- temp$comparison
        transform <- c(temp$transform, "reordered components")
    }
    # Drop names
    if (identical(names(model), names(comparison))) {        
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        model <- unname(model)
        comparison <- unname(comparison)
        transform <- c(transform, "dropped names")
        comp <- same(model, comparison, transform, equal, ...)
    }
    comp
}

###############
# Ingore attributes
###############

compareIgnoreAttrs <- function(model, comparison,
                               transform=character(),
                               equal=TRUE, ...) {
    UseMethod("compareIgnoreAttrs")
}

# Handle any sort of vector ?
# More complex objects will need to write specific methods
compareIgnoreAttrs.default <- function(model, comparison,
                                       transform=character(),
                                       equal=TRUE, ...) {
    attrsM <- attributes(model)
    attrsC <- attributes(comparison)
    # If neither object has attributes
    # or the attributes are identical, don't transform
    if (identical(attrsM, attrsC)) {
        comp <- same(model, comparison, transform, equal, ...)
    } else {
        attributes(model) <- NULL
        attributes(comparison) <- NULL
        transform <- c(transform, "dropped attributes")
        comp <- same(model, comparison, transform, equal, ...)
    }
    comp    
}

###############
# Some pre-packaged higher-level functions
###############

# The ellipsis allows for, e.g.,
# ignoreCase=FALSE for calls to compareEqual() on strings
# These need to be passed to all comparisons so that
# they can pass them on to all subcomparisons

compare <- function(model, comparison,
                    # Which comparisons to perform
                    equal=TRUE,
                    coerce=allowAll,
                    shorten=allowAll,
                    ignoreOrder=allowAll,
                    ignoreNameCase=allowAll,
                    ignoreNames=allowAll,
                    ignoreAttrs=allowAll,
                    # Toggles for various comparisons
                    round=FALSE, # Only allow rounding if specified
                                 # A value of TRUE means round to 0 d.p.
                                 # so using round=allowAll is dangerous
                    ignoreCase=allowAll,
                    trim=allowAll,
                    dropLevels=allowAll,
                    ignoreLevelOrder=allowAll,
                    ignoreDimOrder=allowAll,
                    ignoreColOrder=allowAll,
                    ignoreComponentOrder=allowAll,
                    colsOnly=!allowAll,
                    # Turn all options on
                    # (still can override per option)
                    allowAll=FALSE) {
    comp <- compareIdentical(model, comparison)
    if (!comp$result && equal) {
        comp <- compareEqual(comp$tM,
                             comp$tC,
                             comp$transform,
                             # For all compareEqual() calls
                             equal=equal,
                             # These for recursive calls to compare()
                             coerce=coerce,
                             shorten=shorten,
                             ignoreOrder=ignoreOrder,
                             ignoreNameCase=ignoreNameCase,
                             ignoreNames=ignoreNames,
                             ignoreAttrs=ignoreAttrs,
                             # These for various compare*() calls
                             # (the standalone transforms)
                             round=round,
                             ignoreCase=ignoreCase,
                             trim=trim,
                             dropLevels=dropLevels,
                             ignoreLevelOrder=ignoreLevelOrder,
                             ignoreDimOrder=ignoreDimOrder,
                             ignoreColOrder=ignoreColOrder,
                             ignoreComponentOrder=ignoreComponentOrder,
                             colsOnly=colsOnly,
                             recurseFun=compare)
    }
    if (!comp$result && coerce) {
        # NOTE: use tMpartial, tCpartial, and partialTransform
        # (i.e., WITHOUT any compareEqual() transformations
        #  because these will be REAPPLIED after the transform)
        comp <- compareCoerce(comp$tMpartial, 
                              comp$tCpartial,
                              comp$partialTransform,
                              # For all compareEqual() calls
                              equal=equal,
                              # These for recursive calls to compare()
                              coerce=coerce,
                              shorten=shorten,
                              ignoreOrder=ignoreOrder,
                              ignoreNameCase=ignoreNameCase,
                              ignoreNames=ignoreNames,
                              ignoreAttrs=ignoreAttrs,
                              # These for various compare*() calls
                              # (the standalone transforms)
                              round=round,
                              ignoreCase=ignoreCase,
                              trim=trim,
                              dropLevels=dropLevels,
                              ignoreLevelOrder=ignoreLevelOrder,
                              ignoreDimOrder=ignoreDimOrder,
                              ignoreColOrder=ignoreColOrder,
                              ignoreComponentOrder=ignoreComponentOrder,
                              colsOnly=colsOnly,
                              recurseFun=compare)
    }
    if (!comp$result && shorten) {
        # NOTE: use tMpartial, tCpartial, and partialTransform
        # (i.e., WITHOUT any compareEqual() transformations
        #  because these will be REAPPLIED after the transform)
        comp <- compareShorten(comp$tMpartial, 
                               comp$tCpartial,
                               comp$partialTransform,
                               # For all compareEqual() calls
                               equal=equal,
                               # These for recursive calls to compare()
                               coerce=coerce,
                               shorten=shorten,
                               ignoreOrder=ignoreOrder,
                               ignoreNameCase=ignoreNameCase,
                               ignoreNames=ignoreNames,
                               ignoreAttrs=ignoreAttrs,
                               # These for various compare*() calls
                               # (the standalone transforms)
                               round=round,
                               ignoreCase=ignoreCase,
                               trim=trim,
                               dropLevels=dropLevels,
                               ignoreLevelOrder=ignoreLevelOrder,
                               ignoreDimOrder=ignoreDimOrder,
                               ignoreColOrder=ignoreColOrder,
                               ignoreComponentOrder=ignoreComponentOrder,
                               colsOnly=colsOnly,
                               recurseFun=compare)
    }
    if (!comp$result && ignoreOrder) {
        # NOTE: use tMpartial, tCpartial, and partialTransform
        # (i.e., WITHOUT any compareEqual() transformations
        #  because these will be REAPPLIED after the transform)
        comp <- compareIgnoreOrder(comp$tMpartial, 
                                   comp$tCpartial,
                                   comp$partialTransform,
                                   # For all compareEqual() calls
                                   equal=equal,
                                   # These for recursive calls to compare()
                                   coerce=coerce,
                                   shorten=shorten,
                                   ignoreOrder=ignoreOrder,
                                   ignoreNameCase=ignoreNameCase,
                                   ignoreNames=ignoreNames,
                                   ignoreAttrs=ignoreAttrs,
                                   # These for various compare*() calls
                                   # (the standalone transforms)
                                   round=round,
                                   ignoreCase=ignoreCase,
                                   trim=trim,
                                   dropLevels=dropLevels,
                                   ignoreLevelOrder=ignoreLevelOrder,
                                   ignoreDimOrder=ignoreDimOrder,
                                   ignoreColOrder=ignoreColOrder,
                                   ignoreComponentOrder=ignoreComponentOrder,
                                   colsOnly=colsOnly,
                                   recurseFun=compare)
    }
    if (!comp$result && ignoreNameCase) {
        comp <- compareIgnoreNameCase(comp$tMpartial, 
                                      comp$tCpartial,
                                      comp$partialTransform,
                                      # For all compareEqual() calls
                                      equal=equal,
                                      # These for recursive calls to compare()
                                      coerce=coerce,
                                      shorten=shorten,
                                      ignoreOrder=ignoreOrder,
                                      ignoreNameCase=ignoreNameCase,
                                      ignoreNames=ignoreNames,
                                      ignoreAttrs=ignoreAttrs,
                                      # These for various compare*() calls
                                      # (the standalone transforms)
                                      round=round,
                                      ignoreCase=ignoreCase,
                                      trim=trim,
                                      dropLevels=dropLevels,
                                      ignoreLevelOrder=ignoreLevelOrder,
                                      ignoreDimOrder=ignoreDimOrder,
                                      ignoreColOrder=ignoreColOrder,
                                      ignoreComponentOrder=ignoreComponentOrder,
                                      colsOnly=colsOnly,
                                      recurseFun=compare)
    }
    if (!comp$result && ignoreNames) {
        comp <- compareIgnoreNames(comp$tMpartial, 
                                   comp$tCpartial,
                                   comp$partialTransform,
                                   # For all compareEqual() calls
                                   equal=equal,
                                   # These for recursive calls to compare()
                                   coerce=coerce,
                                   shorten=shorten,
                                   ignoreOrder=ignoreOrder,
                                   ignoreNameCase=ignoreNameCase,
                                   ignoreNames=ignoreNames,
                                   ignoreAttrs=ignoreAttrs,
                                   # These for various compare*() calls
                                   # (the standalone transforms)
                                   round=round,
                                   ignoreCase=ignoreCase,
                                   trim=trim,
                                   dropLevels=dropLevels,
                                   ignoreLevelOrder=ignoreLevelOrder,
                                   ignoreDimOrder=ignoreDimOrder,
                                   ignoreColOrder=ignoreColOrder,
                                   ignoreComponentOrder=ignoreComponentOrder,
                                   colsOnly=colsOnly,
                                   recurseFun=compare)
    }
    if (!comp$result && ignoreAttrs) {
        comp <- compareIgnoreAttrs(comp$tMpartial, 
                                   comp$tCpartial,
                                   comp$partialTransform,
                                   # For all compareEqual() calls
                                   equal=equal,
                                   # These for recursive calls to compare()
                                   coerce=coerce,
                                   shorten=shorten,
                                   ignoreOrder=ignoreOrder,
                                   ignoreNameCase=ignoreNameCase,
                                   ignoreNames=ignoreNames,
                                   ignoreAttrs=ignoreAttrs,
                                   # These for various compare*() calls
                                   # (the standalone transforms)
                                   round=round,
                                   ignoreCase=ignoreCase,
                                   trim=trim,
                                   dropLevels=dropLevels,
                                   ignoreLevelOrder=ignoreLevelOrder,
                                   ignoreDimOrder=ignoreDimOrder,
                                   ignoreColOrder=ignoreColOrder,
                                   ignoreComponentOrder=ignoreComponentOrder,
                                   colsOnly=colsOnly,
                                   recurseFun=compare)
    }
    comp
}

# Given just the name of the model object,
# look for a comparison object, allowing for
# differences in case
compareName <- function(model, compName, ...,
                        ignore.case=TRUE,
                        # The environment containing comparison objects
                        compEnv=.GlobalEnv) {
    # Just search the given environment
    if (!exists(compName, envir=compEnv, inherits=FALSE)) {
        compName <- grep(paste("^", compName, "$", sep=""),
                         ls(envir=compEnv),
                         ignore.case=ignore.case,
                         value=TRUE)
        if (length(compName) > 0) {
            # Only ignore case if necessary
            if (ignore.case) {
                transform <- "renamed object"
            } else {
                transform <- character()
            }
            comparison <- get(compName,
                              envir=compEnv, inherits=FALSE)
            comp <- compare(model, comparison, ...)
            # Prepend any transform done here
            comp$transform <- c(transform, comp$transform)
        } else {
            comp <- comparison(model, NULL, FALSE, "object not found")
        }
    } else {
        comparison <- get(compName,
                          envir=compEnv, inherits=FALSE)
        comp <- compare(model, comparison, ...)
    }
    comp    
}

