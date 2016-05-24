weightassess <-
function(inputter, dataframe, weightvec, 
    prevec = NULL) {
    require(Hmisc)
    prx <- "Unweighted"
    if(!is.null(prevec))
      prx <- "Old Weights"
    if(sum(prevec==1)==length(prevec))
      prx <- "Unweighted"
    out <- list()
    if (is.null(prevec)) {
        prevec <- rep(1, length(weightvec))
    }
    for (i in names(inputter)) {
        target <- eval(parse(text = paste("inputter$", i, sep = "")))
        orign <- wtd.table(eval(parse(text = paste("dataframe$", 
            i, sep = ""))), weight = prevec)$sum.of.weights
        if (is.factor(eval(parse(text = paste("dataframe$", i, 
            sep = ""))))) {
          target <- target[match(names(orign), names(target))]
        }
        if (is.logical(eval(parse(text = paste("dataframe$", 
            i, sep = ""))))) {
          names(orign) <- c("TRUE", "FALSE")
          if(sum(names(target) %in% names(orign))==2)
            target <- target[match(names(orign), names(target))]
          else{
            names(target) <- c("FALSE", "TRUE")
            target <- target[match(names(orign), names(target))]
          }
        }
        origpct <- wtd.table(eval(parse(text = paste("dataframe$", 
            i, sep = ""))), weight = prevec)$sum.of.weights/sum(wtd.table(eval(parse(text = paste("dataframe$", 
            i, sep = ""))), weight = prevec)$sum.of.weights)
        newn <- wtd.table(eval(parse(text = paste("dataframe$", 
            i, sep = ""))), weights = weightvec)$sum.of.weights
        newpct <- wtd.table(eval(parse(text = paste("dataframe$", 
            i, sep = ""))), weights = weightvec)$sum.of.weights/sum(wtd.table(eval(parse(text = paste("dataframe$", 
            i, sep = ""))), weights = weightvec)$sum.of.weights)
        chpct <- newpct - origpct
        rdisc <- target - newpct
        odisc <- target - origpct
        nout <- cbind(target, orign, origpct, newn, newpct, chpct, 
            rdisc, odisc)
        colnames(nout) <- c("Target", paste(prx, "N"), paste(prx, "%"), 
            "Wtd N", "Wtd %", "Change in %", "Resid. Disc.", 
            "Orig. Disc.")
        eval(parse(text = paste("out$", i, "<- nout", sep = "")))
    }
    out
}

