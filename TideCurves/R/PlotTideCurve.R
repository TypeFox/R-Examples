#' Plots the tide curve
#' @description Plots the tide curve and residuum to give an overview of the fitting results.
#' @param data The results from TideCurve. Warning: The synthesis period must overlap with the analysis period.
#' @return Generates a plot
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @export

PlotTideCurve <- function(data){

  numm_imm  <- NULL
  i         <- NULL
  k         <- NULL
  numm      <- NULL
  imm       <- NULL
  Diff      <- NULL
  height    <- NULL
  i.height  <- NULL
  splint    <- NULL
  sl <- data.table()
  dm <- data.table()

  sl <- copy(data$synthesis.lunar)
  dm <- copy(data$data.matrix)
  sl[, numm_imm := paste(i, k, sep="_")]
  dm[, numm_imm := paste(numm, imm, sep="_")]

  setkey(sl, numm_imm)
  setkey(dm, numm_imm)
  ttable <- dm[sl]

  sTable <- ttable[order(imm)][order(numm)]
  sTable[, Diff := height - i.height]

  par(mfrow = c(2, 1))
  plot(x = sTable$date_time, y = sTable$i.height, type = "o", pch = 0,
       cex = .1, ylim = c(min(sTable$i.height,sTable$height,na.rm = TRUE),
       max(sTable$height,sTable$i.height,na.rm = TRUE)),
       xlab = "date_time", ylab = "height" )
  lines(x = sTable$date_time, y = sTable$height, type = "p", col = "red", cex = .5)
  legend(x = 'top', legend = c("Synthesis/Prediction","Observation"),
         lty = c(1, 0), pch = c(0, 1), col = c("black", "red"), inset = c(-0, -0.4), xpd = TRUE)

  plot(x = sTable$date_time, y = sTable$Diff, type = "p", col = "green", pch = 17, cex = .7,
       xlab = "date_time", ylab = "observation  minus synthesis")
  legend(x = 'top', legend = "Residuum", pch = 17, col = "green", inset = c(-0, -0.3), xpd = TRUE)

}