"epi.offset" <- function(id.names)
    {   total <- length(id.names)
        counts <- as.vector(table(id.names))
        offset <- c(1)
        for (i in 2:length(counts)-1)
                {var <- counts[i] + offset[i]
                 offset <- c(offset, var)   
                    }
                offset <- c(offset, total)
                return(offset)
                }   
