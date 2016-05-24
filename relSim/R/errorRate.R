errorRate = function(simResults, bIBS = TRUE, bKI = FALSE,
                     rel = "UN", IBSthresh = 14:17,
                     KIthresh = c(1e3,1e4,1e5,1e6), nLoci = 13){

    rel = toupper(rel)
    if(!grepl("(UN|FS|PC)", rel)){
        stop("rel must be one of 'UN', 'FS', or 'PC'")
    }

    if(bIBS & !bKI){
        o = tabulate(simResults$ibs+1, nbins = 2*nLoci + 1)
        p = cumsum(o)/sum(o)

        if(rel == "UN"){
            return(1-p[IBSthresh])
        }else{
            return(p[IBSthresh])
        }
    }else if(!bIBS & bKI){
        if(rel == "UN"){
            results = matrix(rep(0, 2*length(KIthresh)), ncol = 2)
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                results[i,1] = sum(simResults$pc >= ki)/N
                results[i,2] = sum(simResults$sib >= ki)/N
            }

            return(results)
        }else if(grepl("(FS|PC)",rel)){
            results = rep(0, length(KIthresh))
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                if(rel == "FS"){
                    results[i] = sum(simResults$sib < ki)/N
                }else{
                    results[i] = sum(simResults$pc < ki)/N
                }
            }

            return(results)
        }
    }else if(bIBS & bKI){ ## both
        ## in this situation IBSthresh and KIBSthresh will have
        ## repeated elements
        cat("Both\n")

        if(length(IBSthresh) != length(KIthresh)){
            stop("IBSthresh and KIthresh must be of equal length when using both criteria")
        }

        if(rel == "UN"){
            results = matrix(rep(0, 2*length(KIthresh)), ncol = 2)
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                ibs = IBSthresh[i]
                results[i,1] = sum(simResults$ibs >= ibs
                                   & simResults$pc >= ki)/N
                results[i,2] = sum(simResults$ibs >= ibs
                                   & simResults$sib >= ki)/N
            }

            return(results)
        }else if(grepl("(FS|PC)",rel)){
            results = rep(0, length(KIthresh))
            N = nrow(simResults)

            for(i in 1:length(KIthresh)){
                ki = KIthresh[i]
                ibs = IBSthresh[i]
                if(rel == "FS"){
                    results[i] = sum(simResults$ibs < ibs
                                     | simResults$sib < ki)/N
                }else{
                    results[i] = sum(simResults$ibs < ibs
                                     | simResults$pc < ki)/N
                }
            }

            return(results)
        }
    }
}


