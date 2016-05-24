################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Formulation of an endemic-only twinstim as a Poisson-GLM with response the
### number of events per space-time cell of stgrid and offset log(dt*ds)
###
### Copyright (C) 2013-2014 Sebastian Meyer
### $Revision: 967 $
### $Date: 2014-08-25 10:42:07 +0200 (Mon, 25. Aug 2014) $
################################################################################


utils::globalVariables("area")  # in glm(), the 'offset' is evaluated in 'data'

glm_epidataCS <- function (formula, data, ...)
{
    if (missing(formula)) {
        covariates <- c("start", setdiff(names(data$stgrid), c(
            reservedColsNames_stgrid, obligColsNames_stgrid)))
        formula <- as.formula(paste0("~", paste0(covariates, collapse=" + ")))
    }

    ## for a type-specific model, we really have to set up the full
    ## "stkappagrid", i.e. with nBlocks*nTiles*nTypes rows
    typeSpecificModel <- "type" %in% all.vars(formula)
    typeNames <- levels(data$events@data$type)
    nTypes <- length(typeNames)
    
    ## aggregated number of events in each cell of the stgrid
    if (typeSpecificModel) {
        .stgrid <- do.call("rbind", lapply(typeNames, function (type) {
            cbind(data$stgrid, type=type, deparse.level=0)
        }))
        eventsByCell <- c(table(with(data$events@data, {
            interaction(tile, BLOCK, type, drop=FALSE, sep=".", lex.order=FALSE)
        })))
        .stgrid$nEvents <- eventsByCell[paste(
            .stgrid$tile, .stgrid$BLOCK, .stgrid$type, sep=".")]
    } else {
        .stgrid <- data$stgrid
        eventsByCell <- c(table(with(data$events@data, {
            interaction(tile, BLOCK, drop=FALSE, sep=".", lex.order=FALSE)
        })))
        .stgrid$nEvents <- eventsByCell[paste(
            .stgrid$tile, .stgrid$BLOCK, sep=".")]
    }
    .stgrid$nEvents[is.na(.stgrid$nEvents)] <- 0L
    stopifnot(sum(.stgrid$nEvents) == nobs(data))

    ## Fit corresponding Poisson-GLM
    environment(formula) <- environment() # to see typeSpecificModel and nTypes
    glm(update.formula(formula, nEvents ~ .),
        family = poisson(link="log"),
        data = .stgrid,
        offset = log((if(typeSpecificModel) 1 else nTypes)*(stop-start)*area),
        ...)
}
