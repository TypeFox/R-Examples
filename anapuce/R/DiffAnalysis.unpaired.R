DiffAnalysis.unpaired <-
function (fileIN = "resNorm.txt", n = 3, cond1 = "cond1.", cond2 = "cond2.", 
    fileOUT = "ListOfGenes.txt", fileDelete = "GenesOutOfAnalysis.txt", 
    procs = c("bonferroni", "BH"), alpha = c(0.05, 0.05), fileID = NULL, 
    function.trt = NULL, by.var = "ID", varmixt.meth = TRUE, 
    header = TRUE, sep = "\t", sep.write = "\t", dec.write = ".", 
    ...) 
{
    if (length(alpha) != length(procs)) 
        break
    if (is.character(fileIN)) {
        GenFic <- read.table(fileIN, header = header, sep = sep, 
            ...)
    }
    else {
        GenFic <- fileIN
    }
    indC1 <- grep(cond1, names(GenFic))
    indC2 <- grep(cond2, names(GenFic))
    index <- c(1:n, indC1, indC2)
    cat("####################################### \n")
    cat("DIFFERENTIAL ANALYSIS \n")
    cat("Arrays - first condition : ", names(GenFic)[indC1], 
        "\n")
    cat("Arrays - second condition : ", names(GenFic)[indC2], 
        "\n")
    cat("Outfile for lists of genes : ", fileOUT, "\n")
    cat("Outfile for genes not analysed : ", fileDelete, "\n")
    cat("Multiple testing procedures : ", procs, "\n")
    cat("Nominal Type I error rate for each procedure : ", alpha, 
        "\n")
    cat("####################################### \n")
    GenFic <- GenFic[, index]
    if (is.null(function.trt) == FALSE) 
        GenFic <- function.trt(GenFic, n, cond1, cond2, by.var)
    if (!is.null(fileID)) {
        fichier2 <- read.table(fileID, header = header, sep = sep, 
            ...)
        n <- length(union(names(fichier2), names(GenFic)[1:n]))
        GenFic <- merge(fichier2, GenFic, by = names(fichier2)[names(fichier2) %in% 
            names(GenFic)])
    }
    Fic1 <- GenFic[, grep(cond1, names(GenFic))]
    Fic2 <- GenFic[, grep(cond2, names(GenFic))]
     NbObs1 <- length(indC1) - rowSums(1 * is.na(Fic1))
     NbObs2 <- length(indC2) - rowSums(1 * is.na(Fic2))
     don12 <- GenFic[unique(c(which(NbObs1 < 2), which(NbObs2 < 
        2))), ]
    cat("--------------------------------------------------------------------- \n")
    cat("Number of genes with less than two observations \n")
    cat("Condition 1 :    ", length(which(NbObs1 < 2)), "\n")
    cat("Condition 2 :    ", length(which(NbObs2 < 2)), "\n")
    cat("At least in one of the two conditions :    ", nrow(don12), 
        "\n")
    if ((nrow(GenFic) - length(don12)) > 0) {
        write.table(don12, file = "GenesOutOfAnalysis_NotEnoughObs.txt", 
            row.names = FALSE, sep = "\t", append = FALSE)
    }
    
   
     VarByGene1 <- rowVars(as.data.frame(Fic1), na.rm = TRUE)
    VarByGene2 <- rowVars(as.data.frame(Fic2), na.rm = TRUE)    
    IndVarNulle <- GenFic[unique(c(which(VarByGene1 == 0), which(VarByGene2 == 
        0))), ]
    cat("--------------------------------------------------------------------- \n")
    cat("Number of genes with a null variance :  \n")
    cat("Condition 1 : ", length(which(VarByGene1 == 0)), "\n")
    cat("Condition 2 : ", length(which(VarByGene2 == 0)), "\n")
    cat("At least in one of the two conditions : ", nrow(IndVarNulle), 
        "\n")
    if (length(IndVarNulle) != 0) {
        write.table(IndVarNulle, file = "GenesOutOfAnalysis_NullVar.txt", 
            row.names = FALSE, , append = FALSE, sep = "\t")
    }
    tmp <- unique(c(which(NbObs1 < 2), which(NbObs2 < 2), c(which(VarByGene1 == 
        0), which(VarByGene2 == 0))))
    cat("--------------------------------------------------------------------- \n")
    cat("Total number of genes out of the analysis : ", length(tmp), 
        "\n")
    cat("--------------------------------------------------------------------- \n")
    if (length(tmp) > 0) {
        GenFic <- GenFic[-tmp, ]
        Fic1 <- Fic1[-tmp, ]
        Fic2 <- Fic2[-tmp, ]
        NbObs1 <- NbObs1[-tmp]
        NbObs2 <- NbObs2[-tmp]
        VarByGene1 <- VarByGene1[-tmp]
        VarByGene2 <- VarByGene2[-tmp]
    }
    Moy1 <- rowMeans(GenFic[, grep(cond1, names(GenFic))], na.rm = TRUE)
    Moy2 <- rowMeans(GenFic[, grep(cond2, names(GenFic))], na.rm = TRUE)
    Delta <- Moy1 - Moy2
    pvalVartest = apply(as.matrix(1:nrow(GenFic)), 1, FUN = function(x) bartlett.test(list(as.numeric(GenFic[x, 
        grep(cond1, names(GenFic))]), as.numeric(GenFic[x, grep(cond2, 
        names(GenFic))])))$p.value)
    ou.bad = which(pvalVartest * nrow(GenFic) <= 0.05)
    if (length(ou.bad) > 0) {
        BadTest = data.frame(GenFic[ou.bad, ], Mean1 = Moy1[ou.bad], 
            Mean2 = Moy2[ou.bad], VarByGene1 = VarByGene1[ou.bad], 
            VarByGene2 = VarByGene2[ou.bad])
        cat(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")
        cat(" Warning : See VarianceByGeneNotEqual-Bonf0.05.txt              \n")
        cat(" Equality of variances rejected (Bonferroni 5%) :", 
            length(ou.bad), "\n")
        cat(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")
        write.table(BadTest, file = "VarianceByGeneNotEqual-Bonf0.05.txt", 
            row.names = FALSE)
    }
    PooledVar <- (VarByGene1 * (NbObs1 - 1) + VarByGene2 * (NbObs2 - 
        1))/(NbObs1 + NbObs2 - 2)
    denom <- sqrt(PooledVar * (1/NbObs1 + 1/NbObs2))
    StatOfTest <- Delta/denom
    PvalueRaw <- 2 * (1 - pt(abs(StatOfTest), (NbObs1 + NbObs2 - 
        2)))
    res <- sapply(procs, function(meth) p.adjust(PvalueRaw, meth))
    PvalueAdj <- data.frame(GenFic, Variance = PooledVar, VarianceByGene1 = VarByGene1, 
        VarianceByGene2 = VarByGene2, StatOfTest, PvalueRaw, 
        res)
    cat(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")
    cat("The variances of each group are treated as being equal and the pooled variance \n")
    cat("is used to estimate the variance of the difference between the two groups. \n")
    cat(" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")
    cat("##########################################################\n")
    cat("         One variance by gene                             \n")
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
    cat("      A common variance for all the genes                \n")
    cat("      Variances treated as equal using a pooled variance \n")
    cat("##########################################################\n")
    CommonVar1 <- mean(VarByGene1)
    CommonVar2 <- mean(VarByGene2)
    Sample1 <- ((NbObs1 - 1) * VarByGene1)/CommonVar1
    Sample2 <- ((NbObs2 - 1) * VarByGene2)/CommonVar2
    seuil <- 1/(length(NbObs1))
    ou1p <- which(Sample1 < qchisq(seuil, 2))
    ou1g <- which(Sample1 > qchisq(1 - seuil, 2))
    ou2p <- which(Sample2 < qchisq(seuil, 2))
    ou2g <- which(Sample2 > qchisq(1 - seuil, 2))
    cat("Number of dropped genes -- condition 1: ", length(c(ou1p, 
        ou1g)), "\n")
    cat("Number of dropped genes with a too small variance : ", 
        length(ou1p), "\n")
    cat("Number of dropped genes with a too high variance : ", 
        length(ou1g), "\n")
    cat("Number of dropped genes -- condition 2: ", length(c(ou2p, 
        ou2g)), "\n")
    cat("Number of dropped genes with a too small variance : ", 
        length(ou2p), "\n")
    cat("Number of dropped genes with a too high variance : ", 
        length(ou2g), "\n")
    GeneEnleve <- GenFic[unique(c(ou1p, ou1g, ou2p, ou2g)), ]
    write.table(GeneEnleve, fileDelete, row.names = FALSE, sep = "\t")
    cat("Total number of dropped genes :    ", nrow(GeneEnleve), 
        "\n")
    cat("--------------------------------------------------------------------- \n")
    GenFicm <- GenFic[-unique(c(ou1p, ou1g, ou2p, ou2g)), ]
    CommonVar1b <- mean(VarByGene1[-c(ou1p, ou1g)])
    CommonVar2b <- mean(VarByGene2[-c(ou2p, ou2g)])
    N1 <- NbObs1[-unique(c(ou1p, ou1g, ou2p, ou2g))]
    N2 <- NbObs2[-unique(c(ou1p, ou1g, ou2p, ou2g))]
    sp2 <- ((N1 - 1) * CommonVar1b + (N2 - 1) * CommonVar2b)/(N1 + 
        N2 - 2)
    Denom <- sqrt(sp2 * (N1 + N2)/(N1 * N2))
    Moy1 <- rowMeans(GenFicm[, grep(cond1, names(GenFic))], na.rm = TRUE)
    Moy2 <- rowMeans(GenFicm[, grep(cond2, names(GenFic))], na.rm = TRUE)   
    DeltaMoy <- Moy1 - Moy2
    StatDeTest <- DeltaMoy/Denom
    PvalueRaw <- 2 * (1 - pnorm(abs(StatDeTest)))
    res <- sapply(procs, function(meth) p.adjust(PvalueRaw, meth))
    PvalueAdj <- data.frame(GenFicm, VarianceByGene1 = VarByGene1[-unique(c(ou1p, 
        ou1g, ou2p, ou2g))], VarianceByGene2 = VarByGene2[-unique(c(ou1p, 
        ou1g, ou2p, ou2g))], StatOfTest = StatDeTest, PvalueRaw, 
        res)
    write.table(PvalueAdj, file = paste("Complete", fileOUT, 
        sep = ""), row.names = FALSE, sep = sep.write, dec = dec.write)
    cat(" Estimated variance (standard error) - condition 1: ", 
        CommonVar1b[1], "(", sqrt(var(VarByGene1[-c(ou1p, ou1g)])), 
        ")", "\n")
    cat(" Estimated variance (standard error) - condition 2: ", 
        CommonVar2b[1], "(", sqrt(var(VarByGene2[-c(ou2p, ou2g)])), 
        ")", "\n")
    cat(" Estimated pooled variance : ", sp2[1], "\n")
    cat("--------------------------------------------------------------------- \n")
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
            write.table(fichierRes, file = paste(procs[compteur], 
                "-", alpha[compteur], "-", fileOUT, sep = ""), 
                row.names = FALSE, sep = sep.write, dec = dec.write)
        }
    }
    if (varmixt.meth == TRUE) {
        cat("############################################################\n")
        cat(" Clusters of genes with equal variance (Delmar et al. 2005)\n")
        cat("############################################################\n")
        resvm <- est.varmixt(VAR = PooledVar, Kmax = 50, dfreedom = NbObs1 + 
            NbObs2 - 2)
        groups <- as.vector(table(max.col(resvm$tau)))
        cat("----------------------------------------- \n")
        if (length(groups) == 1) 
            groups[1] = length(PooledVar)
        for (i in 1:length(groups)) {
            cat("Group ", i, " \n")
            cat("   Number of genes ", groups[i], "\n")
            cat("   Variance        ", resvm$vars[i], "\n")
        }
        cat("----------------------------------------- \n")
        teststatVM = sqrt(NbObs1 * NbObs2) * Delta/sqrt(resvm$VM * 
            (NbObs1 + NbObs2))
        if (length(groups) == 1) {
            pval.vm = pnorm(abs(teststatVM), mean = 0, sd = sqrt(resvm$vars/resvm$VM))
        }
        else {
            pval.tmp = apply(as.matrix(1:resvm$nmixt), 1, FUN = function(x) pnorm(abs(teststatVM), 
                mean = 0, sd = sqrt(resvm$vars[x]/resvm$VM)) * 
                resvm$tau[, x])
            pval.vm <- rowSums(pval.tmp)
            
        }
        pval.vm <- 2 * (1 - pval.vm)
        resultatVM <- sapply(procs, function(meth) p.adjust(pval.vm, 
            meth))
        PvalueAdj <- data.frame(GenFic, Delta, VarByGene1, VarByGene2, 
            PooledVar, NbObs1, NbObs2, resvm$tau, VarianceVM = resvm$VM, 
            StatOfTestVM = teststatVM, PValueVM = pval.vm, resultatVM)
        names(PvalueAdj) <- sub("X", "tau", names(PvalueAdj))
        write.table(PvalueAdj, file = paste("VM-Complete", fileOUT, 
            sep = ""), row.names = FALSE, sep = sep.write, dec = dec.write)
        diffrtVM <- apply(as.matrix(1:length(procs)), 1, FUN = function(x) PvalueAdj[which(PvalueAdj[, 
            grep(procs[x], names(PvalueAdj))] <= alpha[x]), ])
        for (compteur in 1:length(diffrtVM)) {
            fichierResVM <- diffrtVM[[compteur]]
            cat("Multiple testing procedure : ", procs[compteur], 
                paste("(", alpha[compteur] * 100, "%)", sep = ""), 
                " \n")
            cat("Number of differentially expressed genes - VM : ", 
                nrow(fichierResVM), "\n")
            cat("-------------------------------------------------------------\n")
            if (nrow(fichierResVM) != 0) {
                write.table(fichierResVM, file = paste(procs[compteur], 
                  "-", alpha[compteur], "-VM-", fileOUT, sep = ""), 
                  row.names = FALSE, sep = sep.write, dec = dec.write)
            }
        }
    }
    invisible(PvalueAdj)
}

