summaryChain = function(chain, probs = c(0.005, 0.025, 0.05, 0.5)) {
   thenames = names(chain)
  thenames = thenames[! thenames %in% "deviance"]
   
   themat = unlist(lapply(chain[thenames], is.matrix))
   
   
   probs = unique(signif(sort(c(probs, 1-probs)),5))
   
   getRes = function(qq) {
      thep = mean(qq>0)
      thep = min(thep, 1-thep)
   
      res = c(mean = mean(qq), 
        pval = thep,
        sd = sd(c(qq)),
        quantile(qq, probs = probs))

		names(res) = gsub("\\%", "pct", names(res))
		res
   }

    result = list(scalars = matrix(NA, sum(themat), length(getRes(1)),
      dimnames = list(thenames[themat], names(getRes(1)) ) )  )
    
    
   for(D in thenames[themat]) {
     result$scalars[D,] = getRes(chain[[D]])
   }
                
   vectorParams = thenames[!themat]
   vectorParams = vectorParams[!vectorParams %in% grep("Grid$", vectorParams, value=T)]
   
   for(D in vectorParams) {
    result[[D]] = t(apply(chain[[D]], 3, getRes))
   
   }
   
   # if this is a spatial model, get summaries of exponentials of fitted values
   
   if(length(grep("Spatial$",names(result)))) {
     for(D in grep("^Fitted", names(result), value=T)) {
        result[[paste("FittedRate", gsub("^FittedR", "", D), sep="")]]=
                t(apply(exp(chain[[D]]), 3, getRes))

     } 
   }
   

   return(result)

}
                                            