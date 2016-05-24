plotMoreSamples <-
function (dataframe,
          LMfits,
          sampleLabels = levels(unlist(dataframe["Sample"])),
          indexOfReference = 1,
          pdfName = "SpecifiedModel.pdf",
          joinReplicates = TRUE, 
          showRho = TRUE,
          xlab = "Log(Dosis)",
          ylab = "Respons: Response", 
          pch = 14 + as.numeric(unlist(dataframe["Replicate"])), 
          cex = 2,
          lwd = 4,
          colTst = "black",
          colRef = "grey",
          colRho = "grey80", 
          main = "Parallel Line Model",
          sub = "Specified Model") 
{
    selectS <- !is.na(match(as.character(unlist(dataframe["Sample"])),
                            sampleLabels))
    dataframe <- dataframe[selectS, ]
    alpha <- coef(LMfits)["(Intercept)"]
    namesCoefs <- names(coef(LMfits))
    i <- (nchar(namesCoefs) == 15) &
        (substr(namesCoefs, 1, 14) == "factor(Sample)")
    i <- (substr(namesCoefs, 1, 14) == "factor(Sample)") &
        unlist(lapply(strsplit(namesCoefs, ":"), length)) < 2
    gammas <- coef(LMfits)[i]
    beta <- coef(LMfits)["Z"]
    Treatments <- unlist(dataframe["indexOfTreatment"])
    if (length(gammas) == length(colTst)) 
        colTsts <- colTst
    else colTsts <- rep(colTst, length(gammas))
    COLORS <- c(colRef, colTsts)
    COLORS[indexOfReference]  <- colRef
    COLORS[-indexOfReference] <- colTsts
    colors <- COLORS[Treatments]
    if (pdfName != "None") 
        pdf(pdfName)
    else {
        main <- ""
        cex <- 3/4
        lwd <- 2
        ylab <- paste("Resp.: ", sub)
    }
    plot(Response ~ Zjitter, col = colors, data = dataframe, xlab = xlab, 
        ylab = ylab, pch = pch, cex = cex, cex.lab = 0.75, main = main, 
        sub = sub)
    replicate <- as.numeric(unlist(dataframe["Replicate"]))
    if (joinReplicates) 
        for (j in 1:(0 + length(unique(Treatments)))) {
            for (i in 1:length(unique(replicate))) {
                dataPlate <- dataframe[replicate == i &
                                       unlist(
                                           dataframe["indexOfTreatment"]) == j, ]
                lines(Response ~ Zjitter, lty = i + 1,
                      col = COLORS[j], data = dataPlate)
            }
        }
    response <- unlist(c(dataframe["Response"]))
    Z <- unlist(c(dataframe["Z"]))
    Z <- Z[!is.na(response)]
    Treatments <- Treatments[!is.na(response)]
    minZ <- min(Z)
    maxZ <- max(Z)
    meanZ = (minZ + maxZ)/2
    q1Z = minZ + (maxZ - minZ)/4
    for (j in unique(Treatments)) {
        i <- Treatments == j
        lines(Z[i], fitted(LMfits)[i], col = colors[i], lwd = lwd, lty = 1)
    }
    drawlogRho <- function(gamma, colRho) {
        a1 <- q1Z - gamma/beta/2
        a2 <- q1Z + gamma/beta/2
        f1 <- alpha + gamma + beta * a1
        f1
        lines(c(a1, a2), rep(f1, 2), col = colRho, lwd = lwd/2)
        if (abs(a1 - a2) > 0.1) {
            arrows(a1, f1, a2, col = colRho, lwd = lwd/2)
            arrows(a2, f1, a1, col = colRho, lwd = lwd/2)
        }
        text((a1 + a2)/2, f1 + 0.2, expression(log(rho)), col = colRho, 
            cex = cex)
    }
    if (showRho) 
        for (i in 1:length(gammas)) drawlogRho(gammas[i], colRho = colTsts[i])
    deltas <- (q1Z - minZ) * beta * gammas/abs(gammas)
    a4 <- q1Z
    f4 <- alpha + beta * a4 - 0 * max(deltas)
    text(a4, f4, paste("Standard", sampleLabels[indexOfReference], sep = ": "), 
        col = colRef, cex = cex)
    for (i in 1:length(gammas)) {
        a5 <- q1Z
        f5 <- alpha + gammas[i] * 1 + beta * a5 + 0 * deltas[i]
        text(a5, f5, paste("Test",
                           sampleLabels[-indexOfReference][i], sep = ": "), 
            col = colTsts[i], cex = cex)
    }
}
