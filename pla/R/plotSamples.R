plotSamples <-
function (dataframe,
          LMfits,
          sampleLabels = levels(unlist(dataframe["Sample"])),
          indexOfReference = 1,
          pdfName = "SpecifiedModel.pdf",
          joinReplicates = TRUE, 
          showRho = TRUE,
          xlab = "Log(Dosis)",
          ylab = "Response: Response", 
          pch = 14 + as.numeric(unlist(dataframe["Replicate"])), 
          cex = 2,
          lwd = 4,
          colTst = "black",
          colRef = "grey",
          colRho = "grey80", 
          main = "Parallel Line Model",
          sub = "Specified Model") 
{
    if (length(sampleLabels) == 2) 
        plotTwoSamples(dataframe, LMfits,
                       sampleLabels = sampleLabels,
                       indexOfReference = indexOfReference,
                       pdfName = pdfName, 
                       joinReplicates = joinReplicates, showRho = showRho, 
                       xlab = xlab, ylab = ylab, pch = pch, cex = cex,
                       lwd = lwd, colTst = colTst, colRef = colRef,
                       colRho = colRho, main = main, sub = sub)
    else plotMoreSamples(dataframe, LMfits,
                         sampleLabels = sampleLabels,
                         indexOfReference = indexOfReference,
                         pdfName = pdfName, 
                         joinReplicates = joinReplicates, showRho = showRho, 
                         xlab = xlab, ylab = ylab, pch = pch, cex = cex,
                         lwd = lwd, colTst = colTst, colRef = colRef,
                         colRho = colRho, main = main, sub = sub)
}
