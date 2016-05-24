"anovaTableI.HH" <-
function(){
    .activeModel <- ActiveModel()
    if (!checkMethod("anova", .activeModel)) {
        errorCondition(message=gettextRcmdr("There is no appropriate anova method for a model of this class."))
        return()
        }
    doItAndPrint(paste("anova(", .activeModel, ")", sep=""))
    }

