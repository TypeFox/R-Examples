imperfect.trees<-function(TP, FN, TN, FP, study, data, ...){
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
        formultexA[i] <- paste0("p", i, "*sy*sx+(1-p", i, ")*(1-ey)*(1-ex),")
        formultexO[i] <- paste0("p", i, "*sy*(1-sx)+(1-p", i, 
            ")*(1-ey)*ex,")
        formultexF[i] <- paste0("p", i, "*(1-sy)*sx+(1-p", i, 
            ")*ey*(1-ex),")
        formultexR[i] <- paste0("p", i, "*(1-sy)*(1-sx)+(1-p", 
            i, ")*ey*ex,")
    }
    formultexR[length(study)] <- paste0("p", i, "*(1-sy)*(1-sx)+(1-p", 
        i, ")*ey*ex)")
    formultex <- "specmodel<-mptspec("
    for (i in 1:length(study)) {
        formultex <- paste0(formultex, formultexA[i], formultexO[i], 
            formultexF[i], formultexR[i])
    }
    sy <- 0
    sx <- 0
    ex <- 0
    ey <- 0




    mpt222 <- mpt(spec=eval(parse(text = formultex)), data = salida, 
        start = c(0.7, 0.5, 0.5, 0.5, rep(0.5, length(study))))
    results <- summary(mpt222)
    coeff <- results$coefficients
    for (i in 1:dim(coeff)[1]) {
        eval(parse(text = paste0(rownames(coeff)[i], "<-", coeff[i, 
            1])))
    }
    if ((sy < 0.5) && (sx < 0.5) && (ex < 0.5) && (ey < 0.5)) {
        for (i in 1:(dim(coeff)[1] - 4)) {
            eval(parse(text = paste0("Pvl", i, "<-1-p", i)))
        }
        Seny <- 1 - ey
        Senx <- 1 - ex
        Espy <- 1 - sy
        Espx <- 1 - sx
        output <- summary(mpt222)
        output$Seny <- Seny
        output$Senx <- Senx
        output$Espy <- Espy
        output$Espx <- Espx
        output$Prevl <- rep(0, (dim(coeff)[1] - 4))
        for (i in 1:(dim(coeff)[1] - 4)) {
            eval(parse(text = paste0("output$Prevl[", i, "]<-Pvl", 
                i)))
        }
        output$study <- study
    }
    else {
        output <- summary(mpt222)
        output$Seny <- sy
        output$Senx <- sx
        output$Espy <- ey
        output$Espx <- ex
        for (i in 1:(dim(coeff)[1] - 4)) {
            eval(parse(text = paste0("output$Prevl[", i, "]<-p", 
                i)))
        }
        output$study <- study
    }
    class(output) <- "imperfect.trees"
    return(output)
}


