printANOVATable = function(xtbl,  sanitize.text.function = function(x){x}, ...){
    colnames(xtbl) = c("Df","Sum Sq","Mean Sq","$F$ value","$\\Pr(>F)$")
    print(xtbl, sanitize.text.function = sanitize.text.function, ...)
}
