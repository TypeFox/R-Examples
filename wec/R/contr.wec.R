contr.wec <-
function(x, ref)
{
    frequencies <- table(x)
    n.cat       <- length(table(x))
    ref         <- which(levels(x) == ref)
    
    new.contrasts <- contr.treatment(n.cat, base=ref)
    new.contrasts[ref,] <- -1 * frequencies[-ref] / frequencies[ref]
    
    colnames(new.contrasts) <- names(frequencies[-ref])
    
    contrasts(x) <- new.contrasts
    return(x)
}
