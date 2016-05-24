finco <-
structure(function (data, level) 
{
    n <- dim(data)[1]
    f <- dim(data)[2] - 1
    indic <- rep(0, f)
    output <- indic
    varia <- f + 1
    for (k in 1:f) {
        rate <- rep(1, f)
        if (k > 1) {
            varia <- c(where, varia)
        }
        for (m in 1:f) {
            if (indic[m] == 0) {
                which <- c(m, varia)
                rate[m] = inconsist(data[, which])
            }
        }
        where <- order(rate)[1]
        output[k] <- rate[where]
        indic[where] <- 1
        if (k > 1) {
            if ((output[k] > output[k - 1]) || (output[k] <= 
                level)) {
                indic <- rep(1, f)
            }
        }
    }
    which <- rev(which)
    which <- which[-1]
    cat("features selected and their inconsistency  rates")
    cat("\n")
    which1 <- which[1:(length(which) - 1)]
    list(varselec = which1, inconsis = output[1:length(which1)])
}, source = c("function(data, level)", "{", "# ***************************************************************", 
"# This function selects features using the FINCO algorithm", 
"# data: name of the dataset", "# level: minumum inconsistency level", 
"# Edgar Acuna, 2003", "# **************************************************************", 
"n <- dim(data)[1]", "f <- dim(data)[2] - 1", "indic <- rep(0, f)", 
"output <- indic", "varia <- f + 1", "for(k in 1:f) {", "#initializing the classification rate", 
"rate <- rep(1, f)", "if(k > 1) {", "varia <- c(where, varia)", 
"}", "for(m in 1:f) {", "if(indic[m] == 0) {", "which <- c(m, varia)", 
"rate[m]=inconsist(data[, which])", "}", "}", "#print(rate)", 
"#print(prov)", "where <- order(rate)[1]", "#cat(\"The entering feature\\n\")", 
"#print(where)", "#correct classfication rate of the entering feature", 
"output[k] <- rate[where]", "indic[where] <- 1", "#Avoiding ties", 
"if(k > 1) {", "if((output[k] > output[k - 1]) || (", "output[k] <= level)) {", 
"indic <- rep(1, f)", "}", "}", "}", "#print(which)", "which <- rev(which)", 
"which <- which[-1]", "cat(\"features selected and their inconsistency  rates\"", 
")", "cat(\"\\n\")", "which1 <- which[1:(length(which) - 1)]", 
"list(varselec = which1, inconsis = output[1:", "length(which1)])", 
"}"))
