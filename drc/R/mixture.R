"mixture" <- function(object, model = c("CA", "Hewlett", "Voelund"), start, startm, control = drmc())
{ 
    model <- match.arg(model)

    if (!is.null(object$"call"$"pmodels"[[2]]))
    {
        bListEl <- as.formula(object$"call"$"pmodels"[[2]])
    } else {
        stop("'collapse' argument should be a formula")    
    }

    ## Constructing new collapse argument ... in common for Hewlett and Voelund
    curveid <- as.character(object$"call"$"curve")
    eListEl1 <- as.formula(paste("~ I(1/(", curveid, "/100))-1", sep = ""))
    eListEl2 <- as.formula(paste("~ I(1/(1-", curveid, "/100))-1", sep = ""))

    blstr <- paste(bListEl, collapse = "")
    el1str <- paste(eListEl1, collapse = "")
    el2str <- paste(eListEl2, collapse = "")
#
#    print(bListEl)
#    print(eListEl1) 
#    print(eListEl2)

    fct <- object$"fct"
    if (model == "CA")
    {
        if (length(fct$names) == 2)  # LL.2 
        {
            mixtfct <- hewlett(fixed = c( NA, 0, 1, NA, NA, 1))
            collapseNew2 <- list(bListEl, eListEl1, eListEl2)
            pmText <- paste(blstr, el1str, el2str, sep = ", ")
        
            eName <- fct$names[2]
            noLim <- 0  # number of c and d parameters
        } else if (length(fct$names) == 3)  # l3 function 
        {
            mixtfct <- hewlett(fixed = c( NA, 0, NA, NA, NA, 1))
            collapseNew2 <- list(bListEl, ~1, eListEl1, eListEl2)
            pmText <- paste(blstr, "~1", el1str, el2str, sep = ", ")
        
            eName <- fct$names[3]
            noLim <- 1  # number of c and d parameters
        } else if (length(fct$names) == 4)  # l4 function
        {
            mixtfct <- hewlett(fixed = c( NA, NA, NA, NA, NA, 1))
            collapseNew2 <- list(bListEl, ~1, ~1, eListEl1, eListEl2)
            pmText <- paste(blstr, "~1", "~1", el1str, el2str, sep = ", ")
            
            eName <- fct$names[4]
            noLim <- 2  # number of c and d parameters  
        } else {
            stop("Does not work for LL.5()")
        }
        
        if (missing(startm)) {startm <- NULL}
        
        class(mixtfct) <- "CA"  # otherwise it would be "Hewlett"
        mixtfct$"name" <- "ca"  # overruling the default name "hewlett"
    }
    
    if (model == "Hewlett")
    {
        if (length(fct$names) == 2)  # LL.2 
        {
            mixtfct <- hewlett(fixed = c(NA, 0, 1, NA, NA, NA))
            collapseNew2 <- list(bListEl, eListEl1, eListEl2, ~1)
            pmText <- paste(blstr, el1str, el2str, "~1", sep = ", ")
        
            eName <- fct$names[2]
            noLim <- 0  # number of c and d parameters
        } else if (length(fct$names) == 3)  # LL.3
        {
            mixtfct <- hewlett(fixed = c(NA, 0, NA, NA, NA, NA))
            collapseNew2 <- list(bListEl, ~1, eListEl1, eListEl2, ~1)
            pmText <- paste(blstr, "~1", el1str, el2str, "~1", sep = ", ")
        
            eName <- fct$names[3]
            noLim <- 1  # number of c and d parameters
        } else if (length(fct$names) == 4)  # LL.4
        {
            mixtfct <- hewlett()
            collapseNew2 <- list(bListEl, ~1, ~1, eListEl1, eListEl2, ~1)
            pmText <- paste(blstr, "~1", "~1", el1str, el2str, "~1", sep = ", ")
            
            eName <- fct$names[4]
            noLim <- 2  # number of c and d parameters  
        } else {
            stop("Does not work for LL.5()")
        }
        if (missing(startm)) {startm <- 1}
    }

    if (model == "Voelund")
    {
        if (length(fct$names) == 2)  # LL.2 
        {
            mixtfct <- voelund(fixed = c( NA, 0, 1, NA, NA, NA, NA))
            collapseNew2 <- list(bListEl, eListEl1, eListEl2, ~1, ~1)
            pmText <- paste(blstr, el1str, el2str, "~1", "~1", sep = ", ")
        
            eName <- fct$names[2]
            noLim <- 0  # number of c and d parameters
        } else if (length(fct$names) == 3)  # LL.3
        {
            mixtfct <- voelund(fixed = c( NA, 0, NA, NA, NA, NA, NA))
            collapseNew2 <- list(bListEl, ~1, eListEl1, eListEl2, ~1, ~1)
            pmText <- paste(blstr, "~1", el1str, el2str, "~1", "~1", sep = ", ")
        
            eName <- fct$names[3]
            noLim <- 1  # number of c and d parameters
        } else if (length(fct$names) == 4)  # LL.4
        {
            mixtfct <- voelund()
            collapseNew2 <- list(bListEl, ~1, ~1, eListEl1, eListEl2, ~1, ~1)
            pmText <- paste(blstr, "~1", "~1", el1str, el2str, "~1", "~1", sep = ", ")
            
            eName <- fct$names[4]
            noLim <- 2  # number of c and d parameters  
        } else {
            stop("Does not work for LL.5()")
        }
        if (missing(startm)) {startm <- c(3, 0.3)}
    }    
        
#    assign("collapseNew2", collapseNew2, envir = .GlobalEnv)    


    ## Checking if levels 0 and 100 are present 
    assayNo <- object$"dataList"$"curveid"   
    if (all(regexpr("0", as.character(unique(assayNo))) < 0 ))
    {
        stop("Level 0 is missing")   
    }
    if (all(regexpr("100", as.character(unique(assayNo))) < 0 ))
    {
        stop("Level 100 is missing")   
    }  
    
    ## Constructing starting values 
    sv <- coef(object)
    
    parNames1 <- object$"parNames"[[1]]
    parNames2 <- object$"parNames"[[3]]
#    eNames <- as.character(parNames2[regexpr(paste(eName, ":", sep=""), parNames1, fixed = TRUE) > 0])
#    eInd <- grep(paste(eName, ":", sep = ""), parNames1)
#    eNames <- as.character(parNames2[eInd])

#    pos0 <- match(paste("factor(", curveid, ")0", sep = ""), eNames)
#    pos1 <- match(paste("factor(", curveid, ")100", sep = ""), eNames)
#    if (is.na(pos0)) {pos0 <- 1}  # it is the intercept
#    if (is.na(pos1)) {pos1 <- 1}
#
#    noED50 <- length(eNames)
#    noB <- length(coef(object)) - noED50 - noLim
#    sv2 <- sv[c(1:(noB + noLim), noB + noLim + pos0, noB + noLim + pos1)]
#    sv2[noB+noLim+2] <- sv2[noB+noLim+1] + sv2[noB+noLim+2]
    
    if (missing(start))
    {
        sv2 <- c(sv[-grep(paste(eName, ":", sep = ""), parNames1)], 
        sv[grep(paste(eName, ":", "100", sep = ""), parNames1)],
        sv[grep(paste(eName, ":", "0", sep = ""), parNames1)])
        sv3 <- as.vector(c(sv2, startm))  # removing the element names
    } else {
        sv3 <- start
    } 
   
#    print(mixtfct)
#    print(collapseNew2)
#    print(sv3)
    mModel <- update(object, fct = mixtfct, pmodels = collapseNew2, start = sv3, control = control)
    
#    rm(collapseNew2, envir = .GlobalEnv)    

    mModel$deviance <- object$"fit"$"value"
    
    mModel$"anova"$"test" <- "F"
    mModel$"anova"$"anovaFit" <- object$"anova"$"anovaFit"  # model1
    mModel$"text" <- paste(model, "model")
    mModel$"pmodelsText" <- pmText

    return(mModel)
}
