setMethod("summary","Dataclass",
          function(object, dims0=1:obsDim(object), runs0 = 1:runs(object), ..., 
                   NOT.A.SIMULATION = TRUE){
            if(NOT.A.SIMULATION) 
               {cat(gettextf("name of Dataclass: %s\n",name(object)))
                cat(gettextf("filename of Dataclass: %s\n", filename(object)))}
            cat(gettextf("dimension of the observations: %d\n", obsDim(object)))
            cat(gettextf("number of runs: %d\n", runs(object)))
            cat(gettextf("size of sample: %d\n", samplesize(object)))
            ISOLD <- isOldVersion(object)
            if(is.null(Data(object))) simulate(object)
 
            lrun0=min(getdistrSimOption("MaxNumberofSummarizedRuns"), 
                      length(runs0))           
            ldim0=min(getdistrSimOption("MaxNumberofSummarizedObsDims"), 
                      length(dims0))           
            
            z0<-runs0[1:lrun0]
            y0<-dims0[1:ldim0]
            x0<-Data(object)[, y0, z0, drop = FALSE]
            apply(x0, c(2,3), summary)
          })

setMethod("summary","Simulation",
          function(object,...){
            if(is.null(Data(object)))
              stop("No Data found -> simulate first")            
            cat(gettextf("filename of simulation: %s\n",filename(object)))
            summary(as(object,"Dataclass"), dims0=1:obsDim(object), 
                    runs0 = 1:runs(object), ..., NOT.A.SIMULATION = FALSE)            
          })


setMethod("summary","Contsimulation",
          function(object,...){
            if(is.null(Data(object)))
              stop("No Data found -> simulate first")
            
            cat(gettextf("name of simulation: %s\n", filename(object)))
            cat(gettextf("rate of contamination: %f\n", rate(object)))
            cat(gettextf("real Data:\n"))
            summary(as(object,"Dataclass"), dims0=1:obsDim(object), 
                    runs0=1:runs(object), ..., NOT.A.SIMULATION = FALSE)            
          })

