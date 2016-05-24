perfect.trees<-function (TP, FN, TN, FP, study, data) 
{
    na.act <- getOption("na.action")
    if (missing(data)) 
        data <- NULL
    no.data <- is.null(data)
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }
    mf <- match.call()
    mf.TP <- mf[[match("TP", names(mf))]]
    mf.FN <- mf[[match("FN", names(mf))]]
    mf.TN <- mf[[match("TN", names(mf))]]
    mf.FP <- mf[[match("FP", names(mf))]]
    mf.study <- mf[[match("study", names(mf))]]
    mf.mods <- mf[[match("mods", names(mf))]]
    TP <- eval(mf.TP, data, enclos = sys.frame(sys.parent()))
    FN <- eval(mf.FN, data, enclos = sys.frame(sys.parent()))
    TN <- eval(mf.TN, data, enclos = sys.frame(sys.parent()))
    FP <- eval(mf.FP, data, enclos = sys.frame(sys.parent()))
    study <- eval(mf.study, data, enclos = sys.frame(sys.parent()))
    mods <- eval(mf.mods, data, enclos = sys.frame(sys.parent()))
    freq <- NULL
    treeid <- rep(seq(1:length(study)), each = 4)
    study_id <- rep(study, each = 4)
    for (i in 1:length(study)) {
        freq <- c(freq, TP[i], FN[i], FP[i], TN[i])
    }
    salida <- data.frame(study_id, treeid, freq)
    formultexA <- rep("a", length(study))
    formultexO <- rep("a", length(study))
    formultexF <- rep("a", length(study))
    formultexR <- rep("a", length(study))
    for (i in 1:length(study)) {
        formultexA[i] <- paste0("P", i, "*sensitivity,")
        formultexO[i] <- paste0("P", i, "*(1-sensitivity),")
        formultexF[i] <- paste0("(1-P", i, ")*(1-specificity),")
        formultexR[i] <- paste0("(1-P", i, ")*specificity,")
    }
    formultexR[length(study)] <- paste0("(1-P", i, ")*specificity )")
    formultex <- "specmodel<-mptspec("
    for (i in 1:length(study)) {
        formultex <- paste0(formultex, formultexA[i], formultexO[i], 
            formultexF[i], formultexR[i])
    }
    mpt_pf <- mpt(spec = eval(parse(text = formultex)), data = salida, 
        start = c(0.7, 0.5, rep(0.5, length(study))))
    resultado <- summary(mpt_pf)
    rownames(resultado$coefficients)[1] <- paste0("Prevalence_", 
        study[1])
    for (i in 4:(length(study) + 2)) {
        rownames(resultado$coefficients)[i] <- paste0("Prevalence_", 
            study[i - 2])
    }
    output <- resultado
    class(output) <- "perfect.trees"
    output
}