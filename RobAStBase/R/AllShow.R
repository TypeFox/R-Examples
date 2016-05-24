setMethod("show", "Neighborhood", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("type:\t", object@type, "\n")
        cat("radius:\t", object@radius, "\n")
    })
setMethod("show", "FixRobModel", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("###### center:\t")
        show(object@center)
        cat("\n###### neighborhood:\t")
        show(object@neighbor)
    })
setMethod("show", "InfRobModel", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("###### center:\t")
        show(object@center)
        cat("\n###### neighborhood:\t")
        show(object@neighbor)
    })
setMethod("show", "InfluenceCurve",
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("### name:\t", object@name, "\n")
        cat("\n### 'Curve':\t")
        show(object@Curve)
#        cat("\n### Risks:\n")
#        print(object@Risks)
        cat("\n### Infos:\n")
        print(object@Infos)
    })
setMethod("show", "IC",
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("### name:\t", object@name, "\n")
        L2Fam <- eval(object@CallL2Fam)
        cat("### L2-differentiable parametric family:\t", L2Fam@name, "\n")
        cat("\n### 'Curve':\t")
        show(object@Curve)
#        cat("\n### Risks:\n")
#        print(object@Risks)
        cat("\n### Infos:\n")
        print(object@Infos)
    })
setMethod("show", "ContIC", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("### name:\t", object@name, "\n")
        L2Fam <- eval(object@CallL2Fam)
        cat("\n### L2-differentiable parametric family:\t", L2Fam@name, "\n")
        cat("### param:\t")
        show(L2Fam@param)
        cat("\n### neighborhood radius:\t", object@neighborRadius, "\n")
#        cat("\n### 'Curve':\t")
#        show(object@Curve)
        cat("\n### clip:\t")
        show(object@clip)                
        cat("### cent:\t")
        show(object@cent)                
        cat("### stand:\n")
        show(object@stand)   
        if(!is.null(object@lowerCase)){
            cat("### lowerCase:\t")
            show(object@lowerCase)   
        }
#        cat("\n### Risks:\n")
#        show(object@Risks)
        cat("\n### Infos:\n")
        show(object@Infos)
    })
setMethod("show", "TotalVarIC", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("### name:\t", object@name, "\n")
        L2Fam <- eval(object@CallL2Fam)
        cat("\n### L2-differentiable parametric family:\t", L2Fam@name, "\n")
        cat("### param:\t")
        show(L2Fam@param)
        cat("\n### neighborhood radius:\t", object@neighborRadius, "\n")
#        cat("\n### 'Curve':\t")
#        show(object@Curve)
        cat("\n### clipLo:\t")
        show(object@clipLo)                
        cat("### clipUp:\t")
        show(object@clipUp)                
        cat("### stand:\n")
        show(object@stand)   
#        cat("\n### Risks:\n")
#        show(object@Risks)
        cat("\n### Infos:\n")
        show(object@Infos)
    })
setMethod("show", "ALEstimate", 
    function(object){
        digits <- getOption("digits")
        show(as(object,"Estimate"))
        if(getdistrModOption("show.details") != "minimal"){
            cat("asymptotic bias:\n")
            print(asbias(object), quote = FALSE)
        }
        if(getdistrModOption("show.details") == "maximal" && !is.null(pIC(object))){
            cat("(partial) influence curve:\n")
            show(pIC(object))
        }
    })
setMethod("show", "kStepEstimate", 
    function(object){
        digits <- getOption("digits")
        show(as(object,"ALEstimate"))
        if(getdistrModOption("show.details") != "minimal"){
            cat("steps:\n")
            print(steps(object), quote = FALSE)
        }
    })
setMethod("show", "MEstimate", 
    function(object){
        digits <- getOption("digits")
        show(as(object,"ALEstimate"))
        if(getdistrModOption("show.details") != "minimal"){
            cat("value of M equation:\n")
            print(Mroot(object), quote = FALSE)
        }
    })
setMethod("show", "OptionalpICList", function(object){
  if(is.null(object)) return(invisible(NULL))
  getMethod("show","pICList")(as(object,"pICList"))
})
setMethod("show", "pICList", function(object){
  if(!length(object)) return(invisible(NULL))
  cat("List of intermediate [p]IC's\n")
  for(i in 1:length(object)){
    oI <- object[[i]]
    cat("[p]IC number", i,":\n")

    if(is(oI,"IC"))
       show(oI)
    else{oIC <- object[[i]]@Curve
         for(j in 1:length(oIC))
             show(oIC[[j]]@Map)
    }
  }
})