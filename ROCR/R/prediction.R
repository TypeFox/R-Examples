prediction <- function(predictions, labels, label.ordering=NULL) {

    ## bring 'predictions' and 'labels' into list format,
    ## each list entry representing one x-validation run

    ## convert predictions into canonical list format
    if (is.data.frame(predictions)) {
        names(predictions) <- c()
        predictions <- as.list(predictions)
    } else if (is.matrix(predictions)) {
        predictions <- as.list(data.frame(predictions))
        names(predictions) <- c()
    } else if (is.vector(predictions) && !is.list(predictions)) {
        predictions <- list(predictions)
    } else if (!is.list(predictions)) {
        stop("Format of predictions is invalid.")
    } 
    ## if predictions is a list -> keep unaltered
  
    ## convert labels into canonical list format
    if (is.data.frame(labels)) {
        names(labels) <- c()
        labels <- as.list( labels)
    } else if (is.matrix(labels)) {
        labels <- as.list( data.frame( labels))
        names(labels) <- c()
    } else if ((is.vector(labels) ||
                is.ordered(labels) ||
                is.factor(labels)) &&
               !is.list(labels)) {
        labels <- list( labels)
    } else if (!is.list(labels)) {
        stop("Format of labels is invalid.")
    }
    ## if labels is a list -> keep unaltered

    ## Length consistency checks
    if (length(predictions) != length(labels))
      stop(paste("Number of cross-validation runs must be equal",
                 "for predictions and labels."))
    if (! all(sapply(predictions, length) == sapply(labels, length)))
      stop(paste("Number of predictions in each run must be equal",
                 "to the number of labels for each run."))
    
    ## only keep prediction/label pairs that are finite numbers
    for (i in 1:length(predictions)) {
        finite.bool <- is.finite( predictions[[i]] )
        predictions[[i]] <- predictions[[i]][ finite.bool ]
        labels[[i]] <- labels[[i]][ finite.bool ]
    }

    ## abort if 'labels' format is inconsistent across
    ## different cross-validation runs
    label.format=""  ## one of 'normal','factor','ordered'
    if (all(sapply( labels, is.factor)) &&
        !any(sapply(labels, is.ordered))) {
        label.format <- "factor"
    } else if (all(sapply( labels, is.ordered))) {
        label.format <- "ordered"
    } else if (all(sapply( labels, is.character)) || 
               all(sapply( labels, is.numeric)) ||
               all(sapply( labels, is.logical))) {
        label.format <- "normal"
    } else {
        stop(paste("Inconsistent label data type across different",
                   "cross-validation runs."))
    }
    
    ## abort if levels are not consistent across different
    ## cross-validation runs
    if (! all(sapply(labels, levels)==levels(labels[[1]])) ) {
        stop(paste("Inconsistent factor levels across different",
                   "cross-validation runs."))
    }
        
    ## convert 'labels' into ordered factors, aborting if the number
    ## of classes is not equal to 2.
    levels <- c()
    if ( label.format == "ordered" ) {
        if (!is.null(label.ordering)) {
            stop(paste("'labels' is already ordered. No additional",
                       "'label.ordering' must be supplied."))
        } else {
            levels <- levels(labels[[1]])
        }
    } else {
        if ( is.null( label.ordering )) {
            if ( label.format == "factor" ) levels <- sort(levels(labels[[1]]))
            else levels <- sort( unique( unlist( labels)))
        } else {
          ## if (!setequal( levels, label.ordering)) {
          if (!setequal( unique(unlist(labels)), label.ordering )) {
            stop("Label ordering does not match class labels.")
          }
          levels <- label.ordering
        }
        for (i in 1:length(labels)) {
            if (is.factor(labels))
              labels[[i]] <- ordered(as.character(labels[[i]]),
                                     levels=levels)
            else labels[[i]] <- ordered( labels[[i]], levels=levels)
        }

    }

    if (length(levels) != 2) {
        message <- paste("Number of classes is not equal to 2.\n",
                         "ROCR currently supports only evaluation of ",
                         "binary classification tasks.",sep="")
        stop(message)
    }

    ## determine whether predictions are continuous or categorical
    ## (in the latter case stop; scheduled for the next ROCR version)
    if (!is.numeric( unlist( predictions ))) {
        stop("Currently, only continuous predictions are supported by ROCR.")
    }

    ## compute cutoff/fp/tp data

    cutoffs <- list()
    fp <- list()
    tp <- list()
    fn <- list()
    tn <- list()
    n.pos <- list()
    n.neg <- list()
    n.pos.pred <- list()
    n.neg.pred <- list()
    for (i in 1:length(predictions)) {
        n.pos <- c( n.pos, sum( labels[[i]] == levels[2] ))
        n.neg <- c( n.neg, sum( labels[[i]] == levels[1] ))
        ans <- .compute.unnormalized.roc.curve( predictions[[i]], labels[[i]] )
        cutoffs <- c( cutoffs, list( ans$cutoffs ))
        fp <- c( fp, list( ans$fp ))
        tp <- c( tp, list( ans$tp ))
        fn <- c( fn, list( n.pos[[i]] - tp[[i]] ))
        tn <- c( tn, list( n.neg[[i]] - fp[[i]] ))
        n.pos.pred <- c(n.pos.pred, list(tp[[i]] + fp[[i]]) )
        n.neg.pred <- c(n.neg.pred, list(tn[[i]] + fn[[i]]) )
    }


    return( new("prediction", predictions=predictions,
                labels=labels,
                cutoffs=cutoffs,
                fp=fp,
                tp=tp,
                fn=fn,
                tn=tn,
                n.pos=n.pos,
                n.neg=n.neg,
                n.pos.pred=n.pos.pred,
                n.neg.pred=n.neg.pred))
}


## fast fp/tp computation based on cumulative summing
.compute.unnormalized.roc.curve <- function( predictions, labels ) {
    ## determine the labels that are used for the pos. resp. neg. class :
    pos.label <- levels(labels)[2]
    neg.label <- levels(labels)[1]

    pred.order <- order(predictions, decreasing=TRUE)
    predictions.sorted <- predictions[pred.order]
    tp <- cumsum(labels[pred.order]==pos.label)
    fp <- cumsum(labels[pred.order]==neg.label)

    ## remove fp & tp for duplicated predictions
    ## as duplicated keeps the first occurrence, but we want the last, two
    ## rev are used.
    ## Highest cutoff (Infinity) corresponds to tp=0, fp=0
    dups <- rev(duplicated(rev(predictions.sorted)))
    tp <- c(0, tp[!dups])
    fp <- c(0, fp[!dups])
    cutoffs <- c(Inf, predictions.sorted[!dups])
    
    return(list( cutoffs=cutoffs, fp=fp, tp=tp ))
}
