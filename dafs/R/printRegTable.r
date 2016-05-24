printRegTable = function(xtbl, sanitize.text.function = function(x){x},
                         test = 't',...){
    colnames(xtbl) = c("Estimate","Std. Error","$t$ value", "$\\Pr(>|t|)$")
    if(test=='z'){
        colnames(xtbl) = c("Estimate","Std. Error","$z$ value", "$\\Pr(>|Z|)$")
    }
    print(xtbl, sanitize.text.function = sanitize.text.function, ...)
}
