# samediff.R
#
# last mod: Oct/16/2014, FW

# Calculate discrimination probabilities, P("different"), from same-different
# judgments
psi <- function(data, oa1 = "s1", oa2 = "s2", resp = "resp"){
    
    ## Check data
    if (!(oa1 %in% names(data) && oa2 %in% names(data)))
        stop("Stimulus names need to be given to oa1 and/or oa2.")
    if (!resp %in% names(data))
        stop("response variable not defined.")
    if (!all(sort(unique(data$resp)) == c("d", "s")))
        stop("response variable does not consist of 'd' and 's' answers only")
    
    ## Frequency different
    formula <- as.formula(paste("~", paste(oa1, oa2, sep=" + ")))
    freq <- as.matrix(unclass(xtabs(formula, data[data$resp == "d",])))
    attr(freq, "call") <- NULL

    ## Probability different
    n <- as.matrix(unclass(xtabs(formula, data)))
    attr(n, "call") <- NULL
    prob <- freq/n
    prob[which(is.na(prob))] <- 1

    x <- if(anyNA(suppressWarnings(
            nnam <- as.numeric(snam <- rownames(freq))))) snam else nnam
    y <- if(anyNA(suppressWarnings(
            nnam <- as.numeric(snam <- colnames(freq))))) snam else nnam

    retval <- list(prob=prob, ntrials=n, freq=freq, x=x, y=y)
    class(retval) <- "psi"
    retval
}

