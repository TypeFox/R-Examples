## this is an exact copy of John Fox's anovaTable from Rcmdr/R/model-menu.R

anovaTableII.HH <- function(){
    .activeModel <- ActiveModel()
    if (is.null(.activeModel)) return()
    if (!checkMethod("Anova", .activeModel)) {
        errorCondition(message=gettextRcmdr("There is no appropriate Anova method for a model of this class."))
        return()
        }
    doItAndPrint(paste("Anova(", .activeModel, ")", sep=""))
    }
