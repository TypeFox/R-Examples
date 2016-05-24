setMethod("simplifyD", "UnivarMixingDistribution",
           function(object){
                if(length(object@mixDistr)==1)
                   return(object@mixDistr[[1]])
                ep <- getdistrOption("TruncQuantile")
                mixC.old <- mixCoeff(object)
                if(any(mixC.old > (1 - ep)))
                   {
                    object@mixCoeff <- 1
                    object@mixDistr <- mixDistr[mixC.old > (1 - ep)]
                    return(flat.mix(object))
                   }
                if(all(mixC.old< ep))
                   return(flat.mix(object))
                mixC.new <- mixC.old[mixC.old>ep]
                mixC.new <- mixC.new/sum(mixC.new)
                mixD.new <- mixDistr(object)[mixC.old>ep]
                object@mixCoeff <- mixC.new
                object@mixDistr <- new("UnivarDistrList", mixD.new)
                return(flat.mix(object))}
           )

setMethod("simplifyD", "UnivarLebDecDistribution",
           function(object){
                ep <- getdistrOption("TruncQuantile")
                if(acWeight(object)> (1 - ep))
                   return(acPart(object))
                if(acWeight(object)< ep)
                   return(discretePart(object))
                return(object)}
           )

setMethod("simplifyD", "AbscontDistribution", 
           function(object)object)
setMethod("simplifyD", "DiscreteDistribution", 
           function(object)object)
