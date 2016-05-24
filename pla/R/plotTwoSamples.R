plotTwoSamples <-
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
    gamma <- coef(LMfits)["factor(Sample)T"]
    namesCoefs <- names(coef(LMfits))
    i <- (nchar(namesCoefs) == 15) &
        (substr(namesCoefs, 1, 14) == "factor(Sample)")
    i <- (substr(namesCoefs, 1, 14) == "factor(Sample)") &
           unlist(lapply(strsplit(namesCoefs, ":"), length)) < 2
    gamma <- coef(LMfits)[i]
    beta <- coef(LMfits)["Z"]
    if (pdfName != "None") 
        pdf(pdfName)
    else {
        main <- ""
        cex <- 3/4
        lwd <- 2
        ylab <- paste("Response: ", sub)
    }
    plot(Response ~ Zjitter,
         col = ifelse(unlist(dataframe["indexOfTreatment"]) == indexOfReference,
             colRef, colTst),
         data = dataframe, xlab = xlab, ylab = ylab, 
         pch = pch, cex = cex, cex.lab = 0.75, main = main, sub = sub)
    replicate <- as.numeric(unlist(dataframe["Replicate"]))
    if (joinReplicates) 
        for (j in 1:2) {
            for (i in 1:length(unique(replicate))) {
                dataPlate <- dataframe[replicate == i &
                                       unlist(
                                           dataframe["indexOfTreatment"]) == j, ]
                lines(Response ~ Zjitter, lty = i + 1,
                      col = ifelse(unlist(
                          dataPlate["indexOfTreatment"]) == indexOfReference,
                          colRef, colTst),
                      data = dataPlate)
            }
        }
    response <- unlist(c(dataframe["Response"]))
    Z <- unlist(c(dataframe["Z"]))
    Z <- Z[!is.na(response)]
    minZ <- min(Z)
    maxZ <- max(Z)
    meanZ = (minZ + maxZ)/2
    q1Z = minZ + (maxZ - minZ)/4
    colors <- c(colRef, colTst)
    Treatments <- unlist(dataframe["indexOfTreatment"])
    Treatments <- Treatments[!is.na(response)]
    for (j in unique(Treatments)) {
        i <- Treatments == j
        lines(Z[i], fitted(LMfits)[i], col = colors[j], lwd = lwd, lty = 1)
    }
    drawlogRho <- function() {
        a1 <- q1Z - gamma/beta/2
        a2 <- q1Z + gamma/beta/2
        f1 <- alpha + gamma + beta * a1
        f1
        lines(c(a1, a2), rep(f1, 2), col = colRho, lwd = lwd)
        if (abs(a1 - a2) > 0.1) {
            arrows(a1, f1, a2, col = colRho, lwd = lwd)
            arrows(a2, f1, a1, col = colRho, lwd = lwd)
        }
        text((a1 + a2)/2, f1 + 0.2, expression(log(rho)), col = colRho, 
            cex = cex)
    }
    if (showRho) 
        drawlogRho()
    delta <- (q1Z - minZ) * beta * gamma/abs(gamma)
    a4 <- q1Z
    f4 <- alpha + gamma * 0 + beta * a4 - delta
    text(a4, f4, paste("Standard", sampleLabels[indexOfReference], sep = "-"), 
        col = colRef, cex = cex)
    a5 <- q1Z
    f5 <- alpha + gamma * 1 + beta * a5 + delta
    text(a5, f5, paste("Test", sampleLabels[-indexOfReference], sep = "-"),
         col = colTst, cex = cex)
}
