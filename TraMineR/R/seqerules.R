seqerules <- function(fsubseq, sortv=NULL, decreasing=FALSE) {
	
        res <- list("rule", "support", "confidence", "lift", "slift", "jmeasure")
        nrule <- 1
        subseqs <- fsubseq$subseq
	matcount <- seqeapplysub(fsubseq, method="presence")
	ntotal <- length(fsubseq$seqe)
        nsubseq <- 0
        for(i in 1:length(subseqs)) {
          #message("sseq = ", subseqs[i])
          sseq <- subseqs[i]
          splitted <- strsplit(as.character(sseq), "-")[[1]]
          if(length(splitted)>1) {
            nsubseq <- nsubseq + 1
            for(z in 1:(length(splitted)-1)) {
              ab <- as.character(subseqs[i])
              a <- paste(splitted[1:z], collapse="-")
              b <- paste(splitted[-(1:z)], collapse="-")
            
#              message("A=", a)
              counta <- sum(matcount[,a])
#              message("B=", b)
              countb <- sum(matcount[,b])
#              message("AB=", ab)
              countab <- sum(matcount[,ab])
              
              pa <- counta/ntotal
              pb <- countb/ntotal
              pab <- countab/ntotal
              confidence <- pab/pa
                                        
              lift <- pab/(pa*pb)
              ## standard lift
              lmb <- max((pa+pb-1), 1/ntotal)
              ups <- 1/max(pa,pb)
              slift <- (lift - lmb)/(ups - lmb)

              ## j-measure
              invconf <- 1-confidence
              invpb <- 1-pb
              jmeasure <- (confidence * log((confidence/pb), base=2)) + (invconf * log((invconf/invpb), base=2)) 

              ## implicative stat
              icstatmat <- implicativestat(matcount[,a], matcount[,b], type="indice")
              if(dim(icstatmat)[1]>1) {
                icstat <- icstatmat[2,2]
              }
              else {
                icstat <- NA
              }
              icstatp <- 1-pnorm(-icstat)
            
              res[["rule"]][nrule] <- paste(a, "=>", b) 
              res[["support"]][nrule] <- countab
              res[["confidence"]][nrule] <- confidence
              res[["lift"]][nrule] <- lift
              res[["slift"]][nrule]  <- slift
              res[["implicative"]][nrule] <- icstat
              res[["implicativep"]][nrule] <- icstatp
              res[["jmeasure"]][nrule] <- jmeasure
             
	      nrule <- nrule+1
              
            }
                                 
          }

       }

        
        
        res[["implicativepb1"]] <- sapply(res[["implicativep"]], function(x) (1-(1-x)^nrule))
        res[["implicativepb2"]] <- sapply(res[["implicativep"]], function(x) (1-(1-x)^nsubseq))
        
        resdata <- data.frame("Rules" = res$rule, "Support"= res$support, "Conf" = res$confidence, "Lift" = res$lift, "Standardlift" = res$slift, "JMeasure"=res$jmeasure, "ImplicStat" = res$implicative, "p-value" = res$implicativep, "p.valueB1" = res$implicativepb1, "p.valueB2"=res$implicativepb2)
        
	if(!is.null(sortv)) {
		or <- order(res[[sortv]],decreasing=decreasing)
		return(resdata[or,])
	}
	
	return(resdata)
          
}


