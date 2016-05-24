`ensemble.formulae` <- function(
    x, factors=NULL, dummy.vars=NULL
)
{
# in older version of raster used layerNames instead of names
    vars <- names(x)
# exclude column for pb for data.frames
    vars <- vars[which(vars != "pb")]
    gamscope <- as.list(vars)
    names(gamscope) <- vars
    nv <- length(vars)
    nf <- length(factors)
    nd <- length(dummy.vars)
    if (is.null(factors) == F) {
        factors <- as.character(factors)
        for (i in 1:nf) {
            if (any(vars==factors[i])==FALSE) {
                cat(paste("\n", "WARNING: categorical variable '", factors[i], "' not among grid layers", "\n", "\n", sep = ""))
            }
        }
    }
    if (is.null(dummy.vars) == F) {
        dummy.vars <- as.character(dummy.vars)
        for (i in 1:nd) {
            if (any(vars==dummy.vars[i])==FALSE) {
               cat(paste("\n", "WARNING: dummy variable '", dummy.vars[i], "' not among grid layers", "\n", "\n", sep = ""))
            }
        }
    }
    results <- list(GBM.formula=NULL, RF.formula=NULL, GLM.formula=NULL, STEP.formula=NULL, GLMSTEP.scope=NULL,  
        GAM.formula=NULL, GAMSTEP.scope=NULL, MGCV.formula=NULL, MGCVFIX.formula=NULL,
        EARTH.formula=NULL, RPART.formula=NULL, NNET.formula=NULL,
        FDA.formula=NULL, SVM.formula=NULL, SVME.formula=NULL)
    numpb <- paste("pb ~ ")
    catpb <- paste("as.factor(pb) ~ ")
    stepvars <- paste(vars[1])
    allvars <- paste0("allvars", c(1:nv))
    for (i in 1:nv) {allvars[i] <- paste(vars[i])}
    glmvars <- gamvars <- mgcvvars <- mgcvfixvars <- explicitcatvars <- allvars
    for (i in 1:nv) {
        if (any(vars[i]==factors) == T) {
            explicitcatvars[i] <- paste("as.factor(", vars[i], ")", sep="")
            gamvars[i] <- paste("as.factor(", vars[i], ")", sep="")
            gamscope[[as.name(vars[i])]] <- as.formula(paste("~1 + as.factor(", vars[i], ")", sep=""))
        }else{
            if (any(vars[i]==dummy.vars) == T) {            
                gamscope[[as.name(vars[i])]] <- as.formula(paste("~1 +", vars[i], sep=""))
            }else{
                glmvars[i] <- paste(vars[i], "+ I(", vars[i], "^2) + I(", vars[i], "^3)", sep="")
                gamvars[i] <- paste("gam::s(", vars[i], ", 4)", sep="")
                mgcvvars[i] <- paste("s(", vars[i], ", k=4)", sep="")
                mgcvfixvars[i] <- paste("s(", vars[i], ", k=4, fx=T)", sep="")
                gamscope[[as.name(vars[i])]] <- as.formula(paste("~1 + ", vars[i], "+ gam::s(", vars[i], ", 4)", sep=""))
            }
        }
    }
    ne <- nv-nf
    earthvars <- NULL
    if (ne > 0) {
        earthvars <- paste0("earthvars", c(1:ne))        
        j <- 0
        for (i in 1:nv) {
            if (any(vars[i]==factors) == F) { 
                j <- j+1
                earthvars[j] <- paste(vars[i])
            }
        }
    }
    results$GBM.formula <- as.formula(paste(numpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    results$RF.formula <- as.formula(paste(numpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    results$GLM.formula <- as.formula(paste(numpb, paste(glmvars, sep="", collapse="+"), sep="", collapse="+"))
    results$STEP.formula <- as.formula(paste(numpb, stepvars, sep="", collapse="+"))
    results$GLMSTEP.scope <- list(upper=as.formula(paste("~", paste(glmvars, sep="", collapse="+"), sep="", collapse="+")), lower=as.formula(paste("~1")))
    results$GAM.formula <- as.formula(paste(numpb, paste(gamvars, sep="", collapse="+"), sep="", collapse="+"))
    results$GAMSTEP.scope <- gamscope
    results$MGCV.formula <- as.formula(paste(numpb, paste(mgcvvars, sep="", collapse="+"), sep="", collapse="+"))
    results$MGCVFIX.formula <- as.formula(paste(numpb, paste(mgcvfixvars, sep="", collapse="+"), sep="", collapse="+"))
# no categorical variables for earth
    results$EARTH.formula <- as.formula(paste(catpb, paste(earthvars, sep="", collapse="+"), sep="", collapse="+"))
    results$RPART.formula <- as.formula(paste(catpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    results$NNET.formula <- as.formula(paste(catpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    results$FDA.formula <- as.formula(paste(numpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    results$SVM.formula <- as.formula(paste(numpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    results$SVME.formula <- as.formula(paste(catpb, paste(allvars, sep="", collapse="+"), sep="", collapse="+"))
    if (is.null(factors) == F) {
        cat(paste("\n", "Note that categorical variables were not included by ensemble.formulae for EARTH", "\n", sep = ""))
    }
    return(results)
}


