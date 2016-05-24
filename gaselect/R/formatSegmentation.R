#' Format the raw segmentation list returned from the C++ code into a usable list
#'
#' @param object The Evaluator object.
#' @param segments The raw segmentation list.
#' @return A list of the form replication -> outerSegment -> (calibration, validation, inner -> (test, train))
#' @include Evaluator.R
#' @rdname formatSegmentation-methods
setGeneric("formatSegmentation", function(object, segments) { standardGeneric("formatSegmentation"); });

#' @rdname formatSegmentation-methods
setMethod("formatSegmentation", signature(object = "GenAlgPLSEvaluator", segments = "list"), function(object, segments) {
    replInd <- rep(seq_len(object@numReplications), each = object@outerSegments * 2 * (object@innerSegments + 1));
    segByRepl <- split(segments, replInd);

    names(segByRepl) <- NULL;

    lapply(segByRepl, function(s) {
        segByOuter <- split(s, rep(seq_len(object@outerSegments), each = 2 * (object@innerSegments + 1)));

        names(segByOuter) <- NULL;

        lapply(segByOuter, function(s) {
            names(s) <- rep.int(c("train", "test"), object@innerSegments + 1);

            ls <- length(s);

            return(list(
                inner = unname(split(s[seq_len(ls - 2)], rep(seq_len(object@innerSegments), each = 2))),
                calibration = unname(s[[ls - 1]]),
                validation = unname(s[[ls]])
            ));
        });
    });
});

#' @rdname formatSegmentation-methods
setMethod("formatSegmentation", signature(object = "GenAlgUserEvaluator", segments = "list"), function(object, segments) {
    return(vector("list"));
});

#' @rdname formatSegmentation-methods
setMethod("formatSegmentation", signature(object = "GenAlgLMEvaluator", segments = "list"), function(object, segments) {
    return(vector("list"));
});

#' @rdname formatSegmentation-methods
setMethod("formatSegmentation", signature(object = "GenAlgFitEvaluator", segments = "list"), function(object, segments) {
    names(segments) <- rep.int(c("train", "test"), object@numSegments);
    seg <- split(segments, rep(seq_len(object@numSegments), each = 2));
    names(seg) <- NULL;

    return(list(list(list(inner = seg))));
});
