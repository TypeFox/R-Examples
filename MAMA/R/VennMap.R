fold.change<-function (x, varname) 
{
  for (i in 1:length(GEDM(x)))
  {
  group <- grep(varname,names(clinical(x)[[i]]))
    mean1j <- apply(GEDM(x)[[i]][, as.numeric(as.factor(clinical(x)[[i]][, group])) == 
        1], 1, mean)
    mean2j <- apply(GEDM(x)[[i]][, as.numeric(as.factor(clinical(x)[[i]][, group])) == 
        2], 1, mean)
    fc <- mean1j - mean2j
    if (i == 1) FCH <- fc else FCH <- cbind(FCH,fc)
     
    }
    colnames(FCH) <- datanames(x)    
    return(FCH)
}

gene.select.FC <- function (fch, cutoff) 
{
    lists <- list()
    for (i in 1:ncol(fch)) lists[[i]] <- rownames(fch)[abs(fch[, 
        i]) >= cutoff]
    names(lists) <- colnames(fch)
    return(lists)
}
