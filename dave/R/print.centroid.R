print.centroid <-
function(x,...) {
    cnames<- strtrim(rownames(x$prob.table),25)
    rnames<- strtrim(colnames(x$prob.table),25)
    cnames<- paste(cnames,cnames[26:27]<- "  ")
    rnames<- paste(rnames,rnames[26:27]<- "  ")
    mm<- matrix(round(as.numeric(t(x$prob.table)),digits=2),nrow=x$nspec)
    mm<- matrix(as.character(mm),ncol=ncol(mm))
    lf<- mm == "0"
    mm[lf == TRUE] <- "."
    rownames(mm) <- rnames
    colnames(mm) <- cnames
    print(mm,quote = FALSE)
  }
