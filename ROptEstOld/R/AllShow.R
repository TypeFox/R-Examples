setMethod("show", "ParamFamParameter", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("name:\t", object@name, "\n")
        cat("main:\t", object@main, "\n")
        if(!is.null(object@nuisance))
            cat("nuisance:\t", object@nuisance, "\n")
        if(!identical(all.equal(object@trafo, diag(length(object)), 
                            tolerance = .Machine$double.eps^0.5), TRUE)){
            cat("trafo:\n")
            print(object@trafo)
        }
    })
setMethod("show", "ParamFamily", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("### name:\t", object@name, "\n")
        cat("\n### distribution:\t")
        print(object@distribution)
        cat("\n### param:\t")
        show(object@param)
        if(length(object@props) != 0){
            cat("\n### props:\n")
            show(object@props)
        }
    })
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
