DiffAnalysis <-
function (fileIN = "resNorm.txt", n = 3, ind.array = 1:2, name.A = "A", 
    name.M = "M.norm", fileOUT = "ListOfGenes.txt", fileDelete = "GenesOutOfAnalysis.txt", 
    procs = c("bonferroni", "BH"), alpha = c(0.05, 0.05), dyeswap = TRUE, 
    indDS = c(2), fileID = NULL, function.trt = NULL, by.var = "ID", 
    varmixt.meth = TRUE, header = TRUE, sep = "\t", sep.write = "\t", 
    dec.write = ".", ...) 
{
    if (length(alpha) != length(procs)) 
        break
    if (dyeswap) {
        if (length(indDS) != length(ind.array)/2) 
            stop("Verify indDS argument. \n Its length is not equal to length(ind.array)/2")
    }
    cat("####################################### \n")
    cat("DIFFERENTIAL ANALYSIS \n")
    cat("Arrays used in the analysis : ", ind.array, "\n")
    cat("Outfile for lists of genes : ", fileOUT, "\n")
    cat("Outfile for genes not analysed : ", fileDelete, "\n")
    cat("Multiple testing procedures : ", procs, "\n")
    cat("Nominal Type I error rate for each procedure : ", alpha, 
        "\n")
    cat("####################################### \n")
    if (is.character(fileIN)) {
        GenFic <- read.table(fileIN, header = header, sep = sep, 
            ...)
    }
    else {
        GenFic <- fileIN
    }
    if (is.null(function.trt) == FALSE) 
        GenFic <- function.trt(GenFic, n, name.A, name.M, by.var)
    indA <- match(paste(name.A, ind.array, sep = ""), names(GenFic))
    indM <- match(paste(name.M, ind.array, sep = ""), names(GenFic))
    index <- c(1:n, indA, indM)
    GenFic <- GenFic[, index]
    if (!is.null(fileID)) {
        fichier2 <- read.table(fileID, header = header, sep = sep, 
            ...)
        n <- length(union(names(fichier2), names(GenFic)[1:n]))
        GenFic <- merge(fichier2, GenFic, by = names(fichier2)[names(fichier2) %in% 
            names(GenFic)])
    }
    indA <- match(paste(name.A, ind.array, sep = ""), names(GenFic))
    indM <- match(paste(name.M, ind.array, sep = ""), names(GenFic))
    if (dyeswap) {
        col.DS <- match(paste(name.M, indDS, sep = ""), names(GenFic))
        GenFic[, col.DS] <- (-GenFic[, col.DS])
    }
    names(GenFic) <- sub("M.norm", "Delta", names(GenFic))
    Int.Mean <- rowMeans(GenFic[, indA], na.rm = TRUE)
    LogRatio.Mean <- rowMeans(GenFic[, indM], na.rm = TRUE)
    All.Mean <- cbind(Int.Mean, LogRatio.Mean)
    tmp.LR <- GenFic[, indM]
    NbObs <- length(indM) - rowSums(1 * is.na(tmp.LR))
    good.NbObs <- which(NbObs > 1)
    if (length(good.NbObs) < length(NbObs)) {
        tmp.LR <- GenFic[-good.NbObs, ]
        write.table(tmp.LR, file = fileDelete, row.names = FALSE, 
            sep = sep.write, dec = dec.write, append = FALSE)
    }
    else {
        nom.colonne = names(GenFic)
        write.table(t(nom.colonne), file = fileDelete, row.names = FALSE, 
            sep = sep.write, dec = dec.write, col.names = FALSE)
    }
    don <- GenFic[good.NbObs, ]
    All.Mean <- All.Mean[good.NbObs, ]
    NbObs <- NbObs[good.NbObs]
    VarianceByGene <- rowVars(don[, indM], na.rm = TRUE)
    IndVarNulle <- which(VarianceByGene == 0)
    cat("Number of genes with a null variance : ", length(IndVarNulle), 
        "\n")
    if (length(IndVarNulle) != 0) {
        tmp.LR <- don[IndVarNulle, ]
        write.table(tmp.LR, file = fileDelete, row.names = FALSE, 
            sep = sep.write, dec = dec.write, append = FALSE)
        don <- don[-IndVarNulle, ]
        good.NbObs <- good.NbObs[-IndVarNulle]
        All.Mean <- All.Mean[-IndVarNulle, ]
        NbObs <- NbObs[-IndVarNulle]
        VarianceByGene <- VarianceByGene[-IndVarNulle]
    }
    Delta <- rowMeans(don[, indM], na.rm = TRUE)
    StatOfTest <- sqrt(NbObs) * (Delta/sqrt(VarianceByGene))
    PvalueRaw <- 2 * (1 - pt(abs(StatOfTest), NbObs - 1))
    PvalueAdj <- sapply(procs, function(meth) p.adjust(PvalueRaw, 
        meth))
    PvalueAdj <- data.frame(don, All.Mean, VarianceByGene, StatOfTest, 
        PvalueRaw, PvalueAdj)
    cat("##########################################################\n")
    cat("                  One variance by gene                    \n")
    cat("##########################################################\n")
    different <- apply(as.matrix(1:length(procs)), 1, FUN = function(x) PvalueAdj[which(PvalueAdj[, 
        grep(procs[x], names(PvalueAdj))] <= alpha[x]), ])
    for (compteur in 1:length(different)) {
        fichierRes <- different[[compteur]]
        cat("Multiple testing procedure : ", procs[compteur], 
            paste("(", alpha[compteur] * 100, "%)", sep = ""), 
            " \n")
        cat("Number of differentially expressed genes : ", nrow(fichierRes), 
            "\n")
        if (nrow(fichierRes) != 0) {
            write.table(fichierRes, file = paste("VarianceByGene-", 
                procs[compteur], "-", alpha[compteur], "-", fileOUT, 
                sep = ""), row.names = FALSE, sep = sep.write, 
                dec = dec.write)
        }
    }
    cat("##########################################################\n")
    cat("       A common variance for all the genes                \n")
    cat("##########################################################\n")
    Variance <- mean(VarianceByGene)
    echantillon <- ((NbObs - 1) * VarianceByGene)/Variance
    G <- nrow(don)
    liste1 <- which(echantillon < qchisq(1e-04, NbObs - 1))
    liste2 <- which(echantillon > qchisq(0.9999, NbObs - 1))
    GeneEnleve <- c(liste1, liste2)
    if (length(GeneEnleve) != 0) {
        tmp.LR <- don[GeneEnleve, ]
        write.table(tmp.LR, file = fileDelete, row.names = FALSE, 
            append = TRUE, col.names = FALSE, sep = sep.write, 
            dec = dec.write)
        cat("Number of dropped genes : ", length(GeneEnleve), 
            "\n")
        cat("Number of dropped genes with a too small variance : ", 
            length(liste1), "\n")
        cat("Number of dropped genes with a too high variance : ", 
            length(liste2), "\n")
        Delta2 <- Delta[-GeneEnleve]
        VarianceByGene2 <- VarianceByGene[-GeneEnleve]
        NbObs2 <- NbObs[-GeneEnleve]
        don2 <- don[-GeneEnleve, ]
        All.Mean2 <- All.Mean[-GeneEnleve, ]
    }
    else {
        Delta2 <- Delta
        VarianceByGene2 <- VarianceByGene
        don2 <- don
        All.Mean2 <- All.Mean
        NbObs2 <- NbObs
    }
    Variance2 <- mean(VarianceByGene2)
    echantillon2 <- ((NbObs2 - 1) * VarianceByGene2)/Variance2
    StatOfTest2 <- (sqrt(NbObs2) * Delta2)/sqrt(Variance2)
    PvalueRaw <- 2 * (1 - pnorm(abs(StatOfTest2)))
    PvalueAdj <- sapply(procs, function(meth) p.adjust(PvalueRaw, 
        meth))
    PvalueAdj <- data.frame(don2, All.Mean2, VarianceByGene = VarianceByGene2, 
        StatOfTest = StatOfTest2, PvalueRaw, PvalueAdj)
    write.table(PvalueAdj, file = paste("Complete", fileOUT, 
        sep = ""), row.names = FALSE, sep = sep.write, dec = dec.write)
    Sd2 <- sqrt(var(VarianceByGene2))
    cat(" Estimated variance (standard error): ", Variance2, 
        "(", Sd2, ")", "\n")
    different <- apply(as.matrix(1:length(procs)), 1, FUN = function(x) PvalueAdj[which(PvalueAdj[, 
        grep(procs[x], names(PvalueAdj))] <= alpha[x]), ])
    for (compteur in 1:length(different)) {
        fichierRes <- different[[compteur]]
        cat("Multiple testing procedure : ", procs[compteur], 
            paste("(", alpha[compteur] * 100, "%)", sep = ""), 
            " \n")
        cat("Number of differentially expressed genes : ", nrow(fichierRes), 
            "\n")
        if (nrow(fichierRes) != 0) 
            write.table(fichierRes, file = paste(procs[compteur], 
                "-", alpha[compteur], "-", fileOUT, sep = ""), 
                row.names = FALSE, sep = sep.write, dec = dec.write)
    }
    if (varmixt.meth == TRUE) {
        cat("############################################################\n")
        cat(" Clusters of genes with equal variance (Delmar et al. 2005)\n")
        cat("############################################################\n")
        resvm <- est.varmixt(VAR = VarianceByGene, Kmax = 50, 
            dfreedom = NbObs - 1)
        groups <- as.vector(table(max.col(resvm$tau)))
        cat("----------------------------------------- \n")
        if (length(groups) == 1) 
            groups[1] = length(VarianceByGene)
        for (i in 1:length(groups)) {
            cat("Group ", i, " \n")
            cat("   Number of genes ", groups[i], "\n")
            cat("   Variance        ", resvm$vars[i], "\n")
        }
        cat("----------------------------------------- \n")
        teststatVM = sqrt(NbObs) * Delta/sqrt(resvm$VM)
        if (length(groups) == 1) {
            pval.vm = pnorm(abs(teststatVM), mean = 0, sd = sqrt(resvm$vars/resvm$VM))
        }
        else {
            pval.tmp = apply(as.matrix(1:resvm$nmixt), 1, FUN = function(x) pnorm(abs(teststatVM), 
                mean = 0, sd = sqrt(resvm$vars[x]/resvm$VM)) * 
                resvm$tau[, x])
            pval.vm = rowSums(pval.tmp)
        }
        pval.vm <- 2 * (1 - pval.vm)
        resultatVM <- sapply(procs, function(meth) p.adjust(pval.vm, 
            meth))
        PValueAdj <- data.frame(GenFic[good.NbObs, ], Delta, 
            VarianceByGene, NbObs, resvm$tau, VarianceVM = resvm$VM, 
            StatOfTestVM = teststatVM, PValueVM = pval.vm, resultatVM)
        names(PValueAdj) <- sub("X", "tau", names(PValueAdj))
        PvalueAdj <- PValueAdj
        write.table(PValueAdj, file = paste("VM-Complete", fileOUT, 
            sep = ""), row.names = FALSE, sep = sep.write, dec = dec.write)
        diffrtVM <- apply(as.matrix(1:length(procs)), 1, FUN = function(x) PValueAdj[which(PValueAdj[, 
            grep(procs[x], names(PValueAdj))] <= alpha[x]), ])
        for (compteur in 1:length(diffrtVM)) {
            fichierResVM <- diffrtVM[[compteur]]
            cat("Multiple testing procedure : ", procs[compteur], 
                paste("(", alpha[compteur] * 100, "%)", sep = ""), 
                " \n")
            cat("Number of differentially expressed genes - VM : ", 
                nrow(fichierResVM), "\n")
            cat("-------------------------------------------------------------\n")
            if (nrow(fichierResVM) != 0) {
                names(fichierResVM) <- sub(name.M, "Delta", names(fichierResVM))
                write.table(fichierResVM, file = paste(procs[compteur], 
                  "-", alpha[compteur], "-VM-", fileOUT, sep = ""), 
                  row.names = FALSE, sep = sep.write, dec = dec.write)
            }
        }
    }
    invisible(PvalueAdj)
}

