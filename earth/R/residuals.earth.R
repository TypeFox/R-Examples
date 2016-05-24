# residuals.earth.R:

residuals.earth <- function(object = stop("no 'object' argument"), type = NULL, warn=TRUE, ...)
{
    glm.resids <- function(object, type)
    {
        g <- object$glm.list
        if(is.null(g))
            stop0("residuals.earth: type \"", type, "\" can be used ",
                  "only on earth-glm models")
        colnames <- ""
        for(imodel in seq_along(g)) {
            rval1 <- residuals(g[[imodel]], type) # invokes residuals.glm
            if(imodel == 1)
                rval <- rval1
            if(NROW(rval1) != NROW(rval)) # should never happen
                stop0("residuals.earth: glm.list[[", imodel, "]] does ",
                      "not conform to glm.list[[", 1, "]] ",
                      "(residuals have a different length)")
            if(imodel > 1) {
                colnames <- c(colnames)
                rval <- cbind(rval, rval1)
            }
        }
        rval
    }
    #--- residuals.earth starts here ---
    warn.if.dots(...)
    if(warn && is.null(type) && !is.null(object$glm.list))
        warning0("residuals.earth: returning earth (not glm) residuals")
    if(is.null(type))
        type <- "earth"
    types <- c("earth", "standardize", "delever", "deviance",
               "glm.pearson", "glm.working", "glm.response", "glm.partial",
               "response") # "response" same as "earth", for plotmo::plotmo_resids
    if(is.null(object$residuals)) # I think this can only happen for cv models
        stop0("earth object has no residuals field.\n",
              "       Use keepxy=TRUE in the call to earth.")
    resids <- switch(match.choices(type, types, "type"),
        earth        = object$residuals,
        standardize  = plotmo::plotmo_standardizescale(object) * object$residuals,
        delever      = object$residuals / sqrt(1 - hatvalues(object)),
        deviance     = if(is.null(object$glm.list))
                           object$residuals
                       else
                           glm.resids(object, "deviance"),
        glm.pearson  = glm.resids(object, "pearson"),
        glm.working  = glm.resids(object, "working"),
        glm.response = glm.resids(object, "response"),
        glm.partial  = glm.resids(object, "partial"),
        response     = object$residuals)

    if(!is.matrix(resids))
        resids <- matrix(resids, ncol = 1)
    if(type != "glm.partial")
        colnames(resids) <- colnames(object$residuals)
    rownames(resids) <- case.names(object)
    resids
}
resid.earth <- function(object = stop("no 'object' argument"), type = NULL, warn=TRUE, ...)
{
    residuals.earth(object, type, warn, ...)
}
