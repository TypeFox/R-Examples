setMethod("initModelClass", signature(model="kin"),                       
          function (model) 
          {
            model@mirf <-  length(model@measured_irf) != 0
            model@irf <-  length(model@irfpar) != 0 || model@mirf
            model@dispmu <-  length(model@parmu) != 0
            model@disptau <-  length(model@partau) != 0
            model@fullk <-  length(dim(model@kmat)) != 1
            if(model@fullk) 
              model@seqmod <- FALSE
            model@wavedep <- (model@dispmu || model@disptau || model@weight || 
                                model@lclp0 || model@lclpequ || length(model@parmu) > 
                                0)
            model@getX <- !(model@dispmu || model@disptau || model@weight || 
                              length(model@parmu) > 0)
            model@clpdep <- model@wavedep
            if(model@fullk) 
              model@ncomp <- nrow(model@kmat) + length(model@kinpar2)
            else if (model@numericalintegration)
              model@ncomp <- length(model@initialvals)
            else 
              model@ncomp <- length(model@kinpar) + length(model@kinpar2)
            if(length(model@nl)==0) {
              model@ncolc <- array(model@ncomp, 1)
            } else {
              model@ncolc <- array(model@ncomp, model@nl)
            }
            
            if (length(model@cohspec$type) == 0) # TODO: find a better test
              model@cohspec$type <- ""
            if (length(model@oscspec$type) == 0) # TODO: find a better test
              model@oscspec$type <- ""
            if(length(model@speckin2$seqmod) == 0)
              model@speckin2$seqmod <- FALSE
            if(length(model@speckin2$jvec) == 0) {
              model@speckin2$fullk <- FALSE
            } else { 
              model@speckin2$fullk <- TRUE
            }
            model@usekin2 <- if( length(model@kinpar2) == 0) FALSE
            else TRUE
            if (length(model@cohspec) != 0) 
              model <- getCoh(model)
            if (length(model@oscspec) != 0) 
              model <- getOsc(model)
            model <- getAnisotropy(model)  
            
            model
          }) 

