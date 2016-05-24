
# poisson.random.sampling
# returns sample of data using poisson sampling

setGeneric("poisson.random.sampling",
           function(x, rate, n)
           standardGeneric("poisson.random.sampling")
           )
setMethod("poisson.random.sampling", "yuima.data",
          function(x, rate, n){
			Data <- get.zoo.data(x)
            if(missing(rate)){
              rate <- rep(0.01,length(Data))
            }
            
			for(i in 1:(length(Data))) {
              T <- end(Data[[i]])-start(Data[[i]])
              Num <- length(Data[[i]])
              deltat <- T /(Num-1)
			  
              # expectation value of length
              Elength <- T*n*rate[i]

              ## make random numbers following exponential distribution
              Time <- diffinv(rexp(Elength,rate=deltat*n*rate[i]))
              # make up a deficit
              while(Time[length(Time)]< Num){
                # adding random numbers (by 10% of Elength)
                Time2 <- diffinv(rexp(trunc(Elength/10+1),rate=deltat*n*rate[i]))+Time[length(Time)]
                # restrain duplication
                Time <- append(Time,Time2[-1])
              }
              
              ## making time index
              # get rid of first element of X and elements that over n
              # round Time to unit
              xTime<-trunc(Time[0<Time & Num>Time])
              # get rid of elements of value = 0
              xTime<-xTime[0<xTime]

              # time index : (xTime-1)*deltat(X)
              idx <- unique(xTime)
              Data[[i]] <- Data[[i]][idx]
            }
             
            return(setData(original.data=Data))
          }
          )

setMethod("poisson.random.sampling","yuima",function(x,rate,n) return(poisson.random.sampling(x@data,rate,n)))
