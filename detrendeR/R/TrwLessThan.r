TrwLessThan = function (rwl, TRW = 0) 
{
    value <- TRW
    PrintMissingRings = function(x, ...) {
        x <- cbind(x, NA)
        n.c <- ncol(x)
        n.r <- nrow(x)
        for (r in 1:n.r) {
            cat(rownames(x)[r])
            for (c in 1:n.c) {
                if (is.na(x[r, c])) 
                  break
                cat(x[r, c])
                if (is.na(x[r, c + 1])) 
                  break
                if (c%%10 != 0) 
                  cat(" ")
                if (c%%10 == 0) {
                  cat("\n",format("", width=nChar-1))
        
                  }
            }
            cat("\n\n")
        }
    }
	names <- colnames(rwl)
	nChar <- max(sapply(names,nchar))+1
    #colnames(rwl) <- substr(paste(colnames(rwl), "       "), 
    #    1, 8)
    	if (nChar < 5) nChar <- 5
    colnames(rwl) <- format(names, width=nChar, justify="left")
    if (value < 0) 
        return(cat("Tree-ring width should be >=0! "))
    if (value == 0) 
        cat("\nMISSING RINGS\n", sep = "")
    first.year <- min(as.numeric(rownames(rwl)))
    Missing = function(rw, value, first.year) {
        years <- which(rw <= value)
        years <- years + first.year - 1
        Years = c(years, rep(NA, length(rw)))
        Years = Years[1:length(rw)]
        Years
    }
    a <- t(apply(rwl, 2, Missing, value = value, first.year = first.year))
    {
        if (all(is.na(a))) {
            if (value == 0) 
                cat("There are no missing rings.\n")
            else {
                cat("There is no ring narrower than ", value, 
                  ".", sep = "")
            }
        }
        else {
            a <- a[apply(a, 1, sum, na.rm = T) > 0, apply(a, 
                2, sum, na.rm = T) > 0]
            number.of.rings <- length(a[!is.na(a)])
            if (number.of.rings == 1 & value == 0) 
                cat("There is one missing ring.\n", sep = "")
            if (number.of.rings == 1 & value > 0) 
                cat("There is one ring narrower than ", value, 
                  ".\n", sep = "")
            if (value == 0 & number.of.rings > 1) 
                cat("There are ", number.of.rings, " missing rings.\n", 
                  sep = "")
            if (number.of.rings > 1 & value > 0) 
                cat("There are ", number.of.rings, " ring narrower than ", 
                  value, ".\n", sep = "")
            PrintMissingRings(a, nChar)
        }
    }
    cat("\n")
}




