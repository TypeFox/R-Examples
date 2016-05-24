printANODEVTable = function(xtbl,  sanitize.text.function = function(x){x},
                            test = NULL, ...){

    if(!is.null(test)){
        if(test=="Chisq"){
            colnames(xtbl) = c("Df", "Deviance", "Resid. Df","Resid. Dev.",
                    "$\\Pr(>\\chi^2)$")
        }else if(test=="F"){
            colnames(xtbl) = c("Df","Deviance","Resid. Df","Resid. Dev.",
                               "$F$ value","$\\Pr(>F)$")
        }else{
            cat(paste("Unknown option", test, "\n"))
        }
    }else{
        colnames(xtbl) = c("Df","Deviance","Resid. Df","Resid. Dev.")
    }

    print(xtbl, sanitize.text.function = sanitize.text.function, ...)
}
