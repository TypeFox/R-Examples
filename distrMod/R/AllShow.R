### from Matthias' thesis / ROptEst
setMethod("show", "ParamFamParameter", 
    function(object){
        cat(gettextf("An object of class \"%s\"\n", class(object)))
        cat(gettextf("name:\t%s\n", object@name))
        if(length(object@main) > 0){
            if(length(object@main) > 1){
                if(is.null(names(object@main)))
                    cat(paste(gettextf("element %s of main:\t%s\n", 
                              1:length(object@main), object@main), 
                              collapse = ""))
                else{
                    cat(paste(gettextf("%s:\t%s\n", names(object@main), 
                              object@main), collapse = ""))
                }
            }else{
                if(is.null(names(object@main)))
                    cat(gettextf("main:\t%s\n", object@main))
                else{
                    cat(gettextf("%s:\t%s\n", names(object@main), object@main))
                }
            }
        }
        if(!is.null(object@nuisance) && length(object@nuisance)){
            if(length(object@main) > 1){
                if(is.null(names(object@nuisance)))
                    cat(paste(gettextf("element %s of nuisance:\t%s\n", 
                              1:length(object@nuisance), object@nuisance), 
                              collapse = ""))
                else{
                    cat(gettext("nuisance:\n"))
                    cat(paste(gettextf("\t%s:\t%s\n", names(object@nuisance), 
                        object@nuisance), collapse = ""))
                }
            }else{
                if(is.null(names(object@nuisance)))
                    cat(gettextf("nuisance:\t%s\n", object@nuisance))
                else{
                    cat(gettext("nuisance:\n"))
                    cat(gettextf("%s:\t%s\n", names(object@nuisance), 
                        object@nuisance))
                }
            }
        }
        if(!is.null(object@fixed) && length(object@fixed)){
            if(length(object@main) > 1){
                if(is.null(names(object@fixed)))
                    cat(paste(gettextf("element %s of fixed part of param.:\t%s\n", 
                              1:length(object@fixed), object@fixed), 
                              collapse = ""))
                else{
                    cat(gettext("fixed part of param.:\n"))
                    cat(paste(gettextf("%s:\t%s\n", names(object@fixed), 
                        object@fixed), collapse = ""))
                }
            }else{
                if(is.null(names(object@fixed)))
                    cat(gettextf("fixed part of param.:\t%s\n", object@fixed))
                else{
                    cat(gettext("fixed part of param.:\n"))
                    cat(gettextf("\t%s:\t%s\n", names(object@fixed), 
                        object@fixed))
                }
            }
        }
        if(!identical(all.equal(trafo(object), diag(length(object)), 
                            tolerance = .Machine$double.eps^0.5), TRUE)){
            if(getdistrModOption("show.details")!="minimal"){
               if(is.function(object@trafo)){
                  if(getdistrModOption("show.details")=="maximal"){
                     cat(gettext("trafo:\n"))
                     print(object@trafo, quote = FALSE)                   
                  }else     
                     cat(gettext("slot trafo is a non-trivial function\n"))
               }else{
                  cat(gettext("trafo:\n"))
                  print(object@trafo, quote = FALSE)
               }
            }
        } 
    })

setMethod("show", "ParamWithShapeFamParameter",
    function(object){
       show(as(object,"ParamFamParameter"))
       if(object@withPosRestr)
          cat(gettext("Shape parameter must not be negative.\n"))
})
setMethod("show", "ParamWithScaleAndShapeFamParameter",
    getMethod("show", "ParamWithShapeFamParameter"))

setMethod("show", "ParamFamily", 
    function(object){
        cat(gettextf("An object of class \"%s\"\n", class(object)))
        cat(gettextf("### name:\t%s\n", object@name))
        cat(gettext("\n### distribution:\t"))
        print(object@distribution, quote = FALSE)
        cat(gettext("\n### param:\t"))
        show(object@param)
        if(length(object@props) != 0){
            cat(gettext("\n### props:\n"))
            show(object@props)
        }
    })
setMethod("show", "RiskType", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("risk type:\t", object@type, "\n")
    })
setMethod("show", "asUnOvShoot", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("risk type:\t", object@type, "\n")
        cat("width:\t", object@width, "\n")
    })
setMethod("show", "asHampel", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("risk type:\t", object@type, "\n")
        cat("bound:\t", object@bound, "\n")
    })
setMethod("show", "fiUnOvShoot", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("risk type:\t", object@type, "\n")
        cat("width:\t", object@width, "\n")
    })
setMethod("show", "fiHampel", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("risk type:\t", object@type, "\n")
        cat("bound:\t", object@bound, "\n")
    })

setMethod("show", "Estimate", 
    function(object){
        title <- gettextf("Evaluations of %s:\n", name(object))
        linet <- paste(paste(rep("-", nchar(title)-1),sep="",collapse=""),
                       "\n",sep="",collapse="") 
        cat(title)
        cat(linet)

        if(getdistrModOption("show.details")!="minimal")
           cat(gettextf("An object of class %s \n", dQuote(class(object))))

        if(getdistrModOption("show.details")!="minimal"){

            cat(gettextf("generated by call\n  "))
            print(estimate.call(object), quote = FALSE)
            if(length(object@samplesize) > 0)
               cat(gettextf("samplesize:   %d\n",object@samplesize))
           }
        
        trafo.mat <- object@trafo$mat
        trafo.fct <- object@trafo$fct
        
        asvar0 <- if(!is.null(object@asvar)) asvar(object) else NULL

        if(!is.null(asvar0)){
           sd0 <- sqrt(diag(asvar0)/object@samplesize)
           if(!is.null(untransformed.asvar(object)) &&
                    all(!is.na(untransformed.asvar(object))))
                untransformed.sd0 <- sqrt(diag(untransformed.asvar(object))/object@samplesize)
           else untransformed.sd0 <- NULL
           
           if(getdistrModOption("show.details")!="minimal")
              cat(gettextf("estimate:\n"))

           dim.est <- dim(object@estimate)
           if(is.null(dim.est))
              .show.with.sd(object@estimate,sd0)
           else{
              if(length(dim.est) >2) stop("not yet implemented")
              c.nms <- colnames(object@estimate)
              r.nms <- rownames(object@estimate)
              rn <- dim.est[1]; cn <- dim.est[2]
              if(rn == 1){
                 dim(object@estimate) <- NULL
                 names(object@estimate) <- c.nms
                 .show.with.sd(object@estimate,sd0)
              }else{
                 cni <- (1:cn)-1
                 for(k in 1:rn){
                     cat("Row [", r.nms[k], ",]:\n", sep="")
                     oe <- object@estimate[k,,drop=TRUE]
                     names(oe) <- paste("[",r.nms[k],",",c.nms,"]",sep="")
                     sd1 <- sd0[cni*rn+k]
                     .show.with.sd(oe,sd1)
                 }
              }
           }

           if(!is.null(object@nuis.idx)){
              cat(gettextf("nuisance parameter:\n"))
              print(nuisance(object), quote = FALSE)        
           }
           
           if(!is.null(object@fixed) && 
               getdistrModOption("show.details")!="minimal"){
              cat(gettextf("fixed part of the parameter:\n"))
              print(fixed(object), quote = FALSE)        
           }

           if(getdistrModOption("show.details")!="minimal"){
               cat(gettextf("asymptotic (co)variance (multiplied with samplesize):\n"))
               print(asvar(object)[,])
              }

           if(getdistrModOption("show.details")=="maximal"){
              if(!.isUnitMatrix(trafo.mat)){
                   if(!is.null(untransformed.sd0) && all(!is.na(untransformed.sd0))){
                      cat(gettextf("untransformed estimate:\n"))
                      .show.with.sd(object@untransformed.estimate,untransformed.sd0)
                   }else{
                      cat(gettextf("untransformed estimate:\n"))
                      print(object@untransformed.estimate, quote = FALSE)
                   }
                   if(!is.null(untransformed.asvar(object))){
                      cat(gettextf("asymptotic (co)variance of untransformed estimate (multiplied with samplesize):\n"))
                      print(untransformed.asvar(object)[,])
                     }
                   }
            }
        }else{

           cat("estimate:\n")
           print(object@estimate, quote = FALSE)

           if(!is.null(object@nuis.idx)){
              cat(gettextf("nuisance parameter:\n"))
              print(nuisance(object), quote = FALSE)        
           }

           if(!is.null(object@fixed) && 
               getdistrModOption("show.details")!="minimal"){
              cat(gettextf("fixed part of the parameter:\n"))
              print(fixed(object), quote = FALSE)        
           }
        } 


        if(getdistrModOption("show.details")=="maximal"){
           if(!.isUnitMatrix(trafo.mat)){
              cat("Transformation of main parameter:\n")
              print(trafo.fct)   
              cat("Trafo / derivative matrix:\n")
              print(trafo.mat[,], quote = FALSE)                
              }
        } 
        
        if(getdistrModOption("show.details")!="minimal"){
           if(nrow(object@Infos) > 0){
             cat("Infos:\n")
             print(object@Infos)
           }
       }
})
   
setMethod("show", "MCEstimate", 
    function(object){
       digits <- getOption("digits")
       show(as(object,"Estimate"))
       if(getdistrModOption("show.details")!="minimal"){
        cat("Criterion:\n")
        print(criterion(object), quote = FALSE)}
    })


setMethod("show", "Confint", 
    function(object){
        if (length(type(object))<2)
            cat(gettextf("A[n] %s confidence interval:\n",type(object)))
        else{
            cat(gettextf("A[n] %s confidence interval:\n",type(object)[1]))
            cat(gettextf("%s\n",type(object)[-1]), sep = "  ")
        }
        print(confint(object), quote = FALSE)
        if(getdistrModOption("show.details")!="minimal"){
            cat(gettextf("Type of estimator: %s\n", name.estimate(object)))

            if(length(object@samplesize.estimate) > 0)
               cat(gettextf("samplesize:   %d\n", samplesize.estimate(object)))

            cat(gettextf("Call by which estimate was produced:\n"))
            print(call.estimate(object), quote = FALSE)

            if(!is.null(nuisance.estimate(object))){
               cat(gettext("Nuisance parameter at which estimate was produced:\n"))
               print(nuisance.estimate(object), quote = FALSE)
            }
        }
        if(getdistrModOption("show.details")=="maximal"){        
            trafo.mat <- object@trafo.estimate$mat
            trafo.fct <- object@trafo.estimate$fct
        
            if(!is.null(fixed.estimate(object))){
               cat(gettext("Fixed part of the parameter at which estimate was produced:\n"))
               print(fixed.estimate(object), quote = FALSE)
            }

            if(!.isUnitMatrix(trafo.mat)){
               cat("Transformation of main parameter by which estimate was produced:\n")
               print(trafo.fct)   
               cat("Trafo / derivative matrix at which estimate was produced:\n")
               print(trafo.mat[,])                
            }
        }
})


setMethod("print", "ShowDetails", 
    function(x, digits = getOption("digits"), 
                show.details = c("maximal", "minimal", "medium")){

        ### match arg show.details
        if(missing(show.details))
           show.details <- getdistrModOption("show.details")
        else 
           show.details <- match.arg(show.details)

        oldDigits <- getOption("digits")
        options("digits" = digits)
        on.exit(options("digits" = oldDigits))
        # unfortunately, within this methods, the straightforward
        #          distrModOptions("show.details"=show.details)
        # /does not work/ --- our work-around
         old.distrModOptions <- distrModOptions()
         new.distrModOptions <- old.distrModOptions
         new.distrModOptions$"show.details" <- show.details
         env <- asNamespace("distrMod")
         assign.error <- FALSE
         if(is(try(assign(".distrModOptions", new.distrModOptions, envir = env),
                   silent = TRUE), "try-error")){
            assign.error <- TRUE
            old.show.details <- getdistrModOption("show.details")
            distrModOptions("show.details" = show.details)
         }
                   
        # end workaround
        show(object = x)
         
#        instead of distrModOptions("show.details"=old.show.details) ::
         if(assign.error) distrModOptions("show.details" = old.show.details)
         else  assign(".distrModOptions", old.distrModOptions, envir = env)
        })

#setMethod("print", "Confint", 
#    function(x, digits = getOption("digits")){
#        oldDigits <- getOption("digits")
#        options("digits" = digits)
#        show(object = x)
#        options("digits" = oldDigits)
#        })        
