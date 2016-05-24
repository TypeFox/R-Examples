effectsStdErr <-
function (reg) 
{
    if ((class(reg) != "lm") && (class(reg) != "nls")) {
        stop("Object of class \"lm\" expected\n")
    }
    summary.effects <- summary(reg)
    aliased <- summary.effects$aliased
    stderr <- rep("NA", length(aliased))
    stderr[!aliased] <- summary.effects$coefficients[, "Std. Error"]
    stderr[stderr == "NA"] <- NA
    return(as.numeric(stderr))
}
