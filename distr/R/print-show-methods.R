setMethod("show", "Symmetry", 
    function(object){ 
        cat(gettextf("type of symmetry:\t%s\n", object@type))
        if(!is.null(object@SymmCenter))
            cat(gettext("center of symmetry:\n"))
            print(object@SymmCenter, quote = FALSE)
    })

#------  UnivariateDistribution ---------- #


setMethod("show", "UnivariateDistribution",
          function(object){
            cls <- class(object)[1]
            cat(showobj(object, className = cls))
            ws <- .IssueWarn(object@.withArith, object@.withSim)
            if(!is.null(ws$msgA)) warning(ws$msgA)
            if(!is.null(ws$msgS)) warning(ws$msgS)
            }
          )

setMethod("show", "LatticeDistribution",
          function(object){
            cls <- class(object)[1]
            cat(showobj(object, className = cls))
            ws <- .IssueWarn(object@.withArith, object@.withSim)
            if(!is.null(ws$msgA)) warning(ws$msgA)
            if(!is.null(ws$msgS)) warning(ws$msgS)
          })


setMethod("showobj", "UnivariateDistribution",
          function(object, className = class(object)[1]){
#            x <- object
            txt <- gettextf("Distribution Object of Class: %s\n", className)
            parameter = param(object)
            Names = slotNames(parameter)
            if(length(Names) > 1){
              for(i in Names[Names != "name"])
                txt <- c(txt,
                gettextf("%s: %s\n", i, slot(parameter, i)))                
            }
            #cat(txt)
            return(txt)
          })

setMethod("print", "UnivariateDistribution",
          function(x, ...)show(x))

#------ DistributionList  ---------- #

setMethod("show", "DistrList", 
    function(object){
       cls <- class(object)[1]
       cat(showobj(object, className = cls))
       ws <- .IssueWarn(any(unlist(lapply(object,function(x)x@.withArith))), 
                        any(unlist(lapply(object,function(x)x@.withSim))))
       if(!is.null(ws$msgA)) warning(ws$msgA)
       if(!is.null(ws$msgS)) warning(ws$msgS)
    
      }
    )

setMethod("showobj", "DistrList", 
    function(object,className = class(object)[1]){
        txt <- gettextf("An object of class \"%s\"\n", className)
        for(i in 1:length(object)){
            s <- showobj(object[[i]])
            les <- length(s)
            if(les >1)
                 st <- c(gettextf("[[%i]] ",i), rep("      :",les-1))
            else st <- gettextf("[[%i]] ",i)
            sts <- paste(st, gettextf("%s", s), sep = "")
            txt <- c(txt, sts)
        }
        return(txt)
    })

#------ UnivarMixingDistribution ---------- #


setMethod("show", "UnivarMixingDistribution",
          function(object){
            cls <- class(object)[1]
            cat(showobj(object, className = cls))
            ws <- .IssueWarn(object@.withArith, object@.withSim)
            if(!is.null(ws$msgA)) warning(ws$msgA)
            if(!is.null(ws$msgS)) warning(ws$msgS)
          }
         )


setMethod("showobj", "UnivarMixingDistribution",
          function(object, className = class(object)[1]){
              txt <- gettextf("An object of class \"%s\"\n", className)
              l <- length(mixCoeff(object))
              txt <- c(txt,
                     "---------------------------------------------\n")
              txt <- c(txt,
                  gettextf("It consists of  %i components \n", l))
              txt <- c(txt,
                     gettextf("Components: \n"))
 
              for(i in 1:l){
                 s <- showobj(mixDistr(object)[[i]])
                 les <- length(s)
                 st <- if (les>1) 
                     c(gettextf("[[%i]]", i), rep("      :",les-1))
                 else  gettextf("[[%i]]", i)
                 txt <- c(txt,
                          paste(st, gettextf("%s",s),sep=""))}
              txt <- c(txt,
                     "---------------------------------------------\n")
              txt <- c(txt, 
                          gettextf("Weights: \n"))
              txt <- c(txt, gettextf("%f",
                          round(mixCoeff(object),3)))
            txt<-c(txt,"\n ---------------------------------------------\n")
            return(txt)
            }
          )

setMethod("showobj", "CompoundDistribution",
          function(object, className = class(object)[1]){
              txt <- gettextf("An object of class \"%s\"\n\n", className)
              txt <- c(txt,
                     gettextf("The frequency distribution is:\n%s",
                              paste(showobj(NumbOfSummandsDistr(object)), collapse="")))
              txt <- c(txt,
                     gettextf("The summands distribution is/are:\n%s",
                              paste(showobj(SummandsDistr(object)), collapse="")))
              txt <- c(txt,
                     gettextf("\nThis Distribution is:\n%s",
                              paste(showobj(simplifyD(object)), collapse="")))                                            
            return(txt)
            }
          )


#------ UnivarLebDecDistribution ---------- #

setMethod("show", "UnivarLebDecDistribution",
          function(object){
           sc <- sys.call(sys.parent(1))
           if(identical(sc[[1]], as.name("show"))){
              objs <- as.character(deparse(match.call(call = sc)$object))
           }else{
              objs <- "<obj>"
           }
            cls <- class(object)[1]
            cat(showobj(object, className = cls, objs = objs))
           ws <- .IssueWarn(object@.withArith, object@.withSim)          
           if(!is.null(ws$msgA)) warning(ws$msgA)
           if(!is.null(ws$msgS)) warning(ws$msgS)
          })

setMethod("showobj", "UnivarLebDecDistribution",
          function(object, className = class(object)[1], objs){
              if(missing(objs)) objs <- "<obj>" 
              else if(length(grep("<S4 object ", objs))) objs <- "<obj>"
              txt <- gettextf("An object of class \"%s\"\n", className)
              txt <- c(txt,
                       gettextf("--- a Lebesgue decomposed distribution:\n\n"))
              l <- length(mixCoeff(object))
              txt <- c(txt,
                gettextf("   Its discrete part (with weight %f) is a\n", 
                  round(discreteWeight(object),3)))
              txt <- c(txt, showobj(discretePart(object)))
              txt <- c(txt,
              gettextf("This part is accessible with 'discretePart(%s)'.\n\n", 
                  objs))
              txt <- c(txt,
              gettextf("   Its absolutely continuous part (with weight %f) is a\n", 
                  round(acWeight(object),3)))
              txt <- c(txt, showobj(acPart(object)))
              txt <- c(txt,
              gettextf("This part is accessible with 'acPart(%s)'.\n", 
                  objs))
             return(txt)
          }          
          )

