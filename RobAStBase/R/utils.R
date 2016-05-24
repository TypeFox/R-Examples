.eq <-function(x,y = 0*x, tol = 1e-7) abs(x-y)<tol

.getDistr <- function(L2Fam){
        slots <- slotNames(L2Fam@distribution@param)
        slots <- slots[slots != "name"]
        nrvalues <- length(slots)
        if (nrvalues > 0) {
            values <- numeric(nrvalues)
            for (i in 1:nrvalues) 
                values[i] <- attributes(attributes(L2Fam@distribution)$param)[[slots[i]]]

            paramstring <- paste("(", paste(values, collapse = ", "), ")", sep = "")
        }else{
            paramstring <- NULL
        }
        distr <- paste(class(L2Fam@distribution)[1], paramstring, sep = "")
}

