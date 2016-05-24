duncan <-
function (y, trt, DFerror, SSerror, alpha = 0.05, group = TRUE,    main = NULL)
{
    MSerror <- SSerror/DFerror
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means <- tapply.stat(junto[, 1], junto[, 2], stat = "mean")
    sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
    nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
    means <- data.frame(means, std.err = sds[, 2]/sqrt(nn[, 2]),
        replication = nn[, 2])
    names(means)[1:2] <- c(name.t, name.y)
    ntr <- nrow(means)

    # criando um contador de medias abrangidas para se calcular as prob de Tukey
    k.snk<-ntr-1
    Tprob <- vector(mode="integer",k.snk)
    kk <- 1
    for (kk in 1:k.snk)  {
        alphap=1-(1-alpha)^((kk+1)-1)
        Tprob[kk] <- qtukey(1 - alphap, kk+1, DFerror)
        }
    # Tprob <- qtukey(1 - alpha, ntr, DFerror)

    nr <- unique(nn[, 2])
    nfila <- c("Alpha", "Error Degrees of Freedom", "Error Mean Square")
    nfila1 <- c("Distances between averages", "Critical Value of Studentized Range")
    # coloquei o vetor Tprop transposto
    nvalor <- c(alpha, DFerror, MSerror)
    nvalor1 <- rbind(t(seq(2,ntr)),t(Tprob))

    #cat("\nStudy:", main)
    #cat("\n\nHSD Test for", name.y, "\n")
    xtabla <- data.frame(...... = nvalor)
    xtabla1 <- data.frame(...... = nvalor1)
    row.names(xtabla) <- nfila
    row.names(xtabla1) <- nfila1
    #print(xtabla1)
    #print(xtabla)
    #cat("\nTreatment Means\n")
    #print(data.frame(row.names = NULL, means))
    HSD <- vector(mode="integer",k.snk)
    if (group) {
        if (length(nr) == 1) {
            kk <- 1
            for (kk in 1:k.snk) {
              HSD[kk] <- Tprob[kk] * sqrt(MSerror/nr)
                 }
            #HSD <- Tprob[1] * sqrt(MSerror/nr)
            #cat("\nHonestly Significant Difference", HSD)
        }
        else {
            nr1 <- 1/mean(1/nn[, 2])
            # criado o for para se ter um vetor de DMS
            kk <- 1
            for (kk in 1:k.snk) {
              HSD[kk] <- Tprob[kk] * sqrt(MSerror/nr1)
                 }
            #cat("\nHonestly Significant Difference", HSD)
            #cat("\nHarmonic Mean of Cell Sizes ", nr1)
            #cat("\n\nDifferent HSD for each comparison")
        }
        cat("\nDuncan's test \n------------------------------------------------------------------------")
        cat("\nGroups  Treatments  Means\n")
        output <- order.stat.SNK(means[, 1], means[, 2], HSD)
        cat('------------------------------------------------------------------------\n')
    }
}
