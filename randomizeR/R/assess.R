#' @include issue.R
#' @include randSeq.R
#' @include util.R
#' @include endpoint.R
NULL

###############################################
# --------------------------------------------#
# Class assessment                            #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the assessment class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateAssessment <- function(object) {
  errors <- character()
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for assessment
# --------------------------------------------

# Randomization paramters generic
setClass("assessment",
         slots = c(D = "data.frame", design = "character", N = "numeric", K = "numeric",
                   groups = "character"),
         validity = validateAssessment)
         

# --------------------------------------------
# Accesssor functions for asssessment
# --------------------------------------------

#' Method defining the $ operator for the assessemnt class
#' 
#' @inheritParams overview
setMethod("$", "assessment",
          function(x, name) slot(x, name))
        

# --------------------------------------------
# Show function for assessment
# --------------------------------------------

setMethod("show", "assessment", function(object) {
    # headline
    #cat("\nObject of class \"", class(object)[1],"\"\n\n", sep="")
    #cat("\nAssessment of ",slot(object, "design"),"\n\n", sep="")
    cat("\nAssessment of a randomization procedure \n\n", sep="")
    # iterate through all slots of the object
    names <- slotNames(object)
    names <- names[!(names == "D")] # without D
    for(name in names) {
      cat(name, "=", slot(object, name), "\n")
    }
    cat("\n") 
    # The data.frame D is printed seperately dependent on its size.
    if (dim(object@D)[1] <= 3) {
      if (nchar(as.character(object@D[1, 1])) >= 10)
        object@D[ ,1] <- paste(substr(object@D[, 1], 1, 9), "...")
      object@D[ ,-1] <- round(object@D[ ,-1], digits = 3)
      print(object@D)
    } else {
      cat("\nThe first 3 rows of", dim(object@D)[1], "rows of D: \n\n")
      obj <- object@D[1:3, ]
      if (nchar(as.character(obj[1, 1])) >= 10)
        obj[ ,1] <- paste(substr(obj[, 1], 1, 9), "...")
      obj[ ,-1] <- round(obj[ ,-1],digits = 3)
      print(obj)
      cat("...")
    }
    cat("\n") 
  }  
)


# --------------------------------------------
# Generic functions for using objects of type issue and randSeq
# --------------------------------------------

#' Assessing randomization sequences
#'
#' Assesses randomization sequences based on specified issues 
#' in clinical trials.
#'
#' @param randSeq object of class \code{randSeq}.
#' @param endp object of class \code{endpoint}, or \code{missing}.
#' @param ... at least one object of class \code{issue} or just a list of objects 
#' of the class \code{issue}.
#'
#' @details
#' Randomization sequences behave differently with respect to issues
#' like selection bias, chronological bias, or loss in power estimation.
#' The \code{assess} function evaluates the behaviour of randomization 
#' sequences with respect to these issues. 
#' The first argument should be a result of one of the functions 
#' \code{\link{genSeq}} or \code{\link{getAllSeq}}.
#' The second argument should be any number of \code{\link{issues}} arising 
#' in a clinical trial. The last argument \code{endp} may be provided if 
#' the assessment should take the distribution of the treamtent groups
#' into account, e.g. for power evaluation.
#'
#' @examples 
#' # assess the full set of Random Allocation Rule for N=4 patients
#' sequences <- getAllSeq(rarPar(4))
#' issue1 <- corGuess("CS")
#' issue2 <- corGuess("DS")
#' issue3 <- imbal("imb")
#' issue4 <- imbal("maxImb")
#' assess(sequences, issue1, issue2, issue3, issue4)
#'
#' # assess one sequence of the Big Stick Design with respect to correct guesses
#' sequence <- genSeq(bsdPar(10, 2), seed = 1909)
#' assess(sequence, issue1)
#'
#' # assess the same sequence with respect to selection bias
#' endp <- normEndp(c(2, 2), c(1, 1))
#' issue5 <- selBias("CS", 4, "exact")
#' issue6 <- setPower(2, "exact")
#' assess(sequence, issue1, issue5, issue6, endp = endp)
#'
#' # recommended plot for the assessment of rejection probabilities
#' RP <- getAllSeq(crPar(6))
#' cB <- chronBias(type = "linT", theta = 1/6, method = "exact")
#' sB <- selBias(type=  "CS", eta = 1/4, method = "exact")
#' normEndp <- normEndp(c(0, 0), c(1, 1))
#' A <- assess(RP, cB, sB, endp = normEndp)
#' D <- A$D
#' desiredSeq <- round(sum(D[,2][D[,3] <= 0.05 & D[,4] <= 0.05]), digits = 4)
#' colnames(D) <- c("Seq", "Prob", "SB", "linT")
#' g <- ggplot(D, aes(x = SB, y = linT))
#' g <- g + annotate("rect", xmin = 0, xmax = 0.05, ymin = 0, ymax = 0.05,
#' alpha=0.2, fill="green") 
#' g <- g + geom_point(alpha = 1/10, size = 3, col = "orange")
#' g <- g <- g + geom_vline(xintercept = 0.05, col = "red")
#' g <- g + geom_hline(yintercept = 0.05, col = "red")
#' g  <- g + geom_text(data = NULL, x = 0, y = 0,
#' label = paste("Proportion:", desiredSeq), hjust=0, vjust=0, size = 7)
#' g
#'
#' @return
#' \code{S4} object of class \code{assessment} summarizing the assessment of the 
#' randomization procedure.
#'
#' @seealso Representation of randomization procedures: \code{\link{randPar}}
#' @seealso Generation of randomization sequences: \code{\link{genSeq}}
#' @seealso \code{\link{issues}} for the assessment of randomization sequences
#' 
#' @name assess
NULL

#' @rdname assess
#'
#' @export
setGeneric("assess", function(randSeq, ..., endp) standardGeneric("assess"))

#' Summary of assessments of a randomization procedure
#' 
#' @param object assessment object.
#' @param ... additional arguments affecting the summary that will be produced.
#'
#' @details
#' For each issue the assessment of the sequences is summarized to permit a design-based 
#' assessment of the randomization procedure.
#' This approach uses the sequence-wise values of the assessment and the probabilities
#' in order to give an overall summary.
#'
#' @return 
#' Data frame with a summary of the assessment object. 
#' 
#' @examples 
#' # assess the full set of PBR(4)
#' seq <- getAllSeq(pbrPar(4))
#' issue <- corGuess("CS")
#' A <- assess(seq, issue)
#' summary(A)
#' 
#' @name summary
NULL

#' @rdname summary
#'
#' @export
setGeneric("summary")



# --------------------------------------------
# Methods for assessment
# --------------------------------------------

#' @rdname assess
setMethod("assess", signature(randSeq = "randSeq", endp = "missing"),
          function(randSeq, ...) {
            L <- list(...)
            if (length(L) == 1 && is.list(L[[1]])) {
              L <- c(...)
            }
            stopifnot(randSeq@K == 2, all(sapply(L, function(x)  is(x, "issue"))))
            stopifnot(randSeq@ratio == c(1, 1))
            D <- data.frame("Sequence" = apply(getRandList(randSeq), 1, function(x) paste(x, sep = "", collapse = "")))
            if (.hasSlot(randSeq, "seed")) { 
              D$Relative_Frequency <- 1/dim(randSeq@M)[1]
            } else {
              D$Probability <- getProb(randSeq)
            }
            
            D <- cbind(D, do.call(cbind, lapply(L, function(x)  getStat(randSeq, x))))
            
            new("assessment",
                D = D, design = getDesign(randSeq),
                N = randSeq@N, K = randSeq@K, groups = randSeq@groups)   
          }
)

#' @rdname assess
setMethod("assess", signature(randSeq = "randSeq", endp = "endpoint"),
          function(randSeq, ..., endp) {
            L <- list(...)
            if (length(L) == 1 && is.list(L[[1]])) {
              L <- c(...)
            }
            stopifnot(randSeq@K == 2, all(sapply(L, function(x) is(x, "issue"))))
            stopifnot(randSeq@ratio == c(1, 1))
            D <- data.frame("Sequence" = apply(getRandList(randSeq), 1, function(x) paste(x, sep = "", collapse = "")))
            if (.hasSlot(randSeq, "seed")) { 
              D$Relative_Frequency <- 1/dim(randSeq@M)[1]
            } else {
              D$Probability <- getProb(randSeq)
            }
            
            D <- cbind(D, do.call(cbind, lapply(L, function(x)  getStat(randSeq, x, endp = endp))))
            
            new("assessment",
                D = D, design = getDesign(randSeq),
                N = randSeq@N, K = randSeq@K, groups = randSeq@groups)   
          }
)

#' @rdname summary
setMethod("summary", signature(object = "assessment"), function(object) {
  D <- object@D
  colnames(D)[2] <- "Probability"
  probs <- D$Probability
  # if (dim(D)[1] == 1) stop("Selected randomization procedure(s) should have more than one generated sequence.")
  D$Probability <- D$Sequence <- NULL
  stat <- apply(D, 2, function(x) {
    ## weighted mean value
    x1 <- sum(x*probs)
    ## weighted standard deviation
    if (dim(D)[1] == 1) x2 <- NA else x2 <- sqrt(sum(probs*(x-x1)^2)/(1 - sum(probs^2))) 
    ## weighted quantiles
    sA <- data.frame(cbind(x, probs))
    d <- x # Save unordered vector x for the computation of max and min (later)
    sA <- sA[order(x),]
    wv <- cumsum(sA[ ,2])
    x <- sA[,1]
    x05 <- x[wv >= 0.05][1]
    x25 <- x[wv >= 0.25][1]
    x50 <- x[wv >= 0.5][1]
    x75 <- x[wv >= 0.75][1]
    x95 <- x[wv >= 0.95][1]
    c(x1, x2, max(d[probs>0]), min(d[probs>0]), x05, x25, x50, x75, x95)
  }) 
  rownames(stat) <- c("mean", "sd", "max", "min", "x05", "x25", "x50", "x75", "x95")
  round(stat, digits=3)
}
)




