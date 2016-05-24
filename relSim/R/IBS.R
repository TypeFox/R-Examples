IBS = function(prof1, prof2, nLoci = length(prof1) / 2, bPrint = FALSE){
    results =  .IBS(prof1, prof2, nLoci)

    if(bPrint){
        for(i in 1:nLoci){
            x = matrix(c(prof1[c(2 *i - 1, 2 * i)], prof2[c(2 * i - 1, 2 * i)]), nrow = 1)
            m = locusIBS(x)

            if(bPrint){
                cat(paste(prof1[c(2 * i - 1, 2 * i)], collapse = "/"),
                    "\t",
                    paste(prof2[c(2 * i - 1, 2 * i)], collapse = "/"),
                    "\t",
                    ifelse(m > 0, "TRUE", "FALSE"),
                    "\t",
                    m,"\n", sep = "")
            }
        }
    }

    return(results)
}
