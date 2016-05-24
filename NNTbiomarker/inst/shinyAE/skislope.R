skisloplot <- function (RSinput= 30, ytop = 100) {
  plot(RSvector, nnt,   xaxs="i",   yaxs="i", ##  log='y',
       xlim=c(0, 50),
       ylim=c(0, ytop), type="l",  lwd=5,
       main = "Number Needed to Treat by Recurrence Score",
       xlab="Recurrence score", ylab="Number needed to treat")

  AERiskTable = c(.029, .15, .57, .20, .05, .001)   #Data from B20 trial; CMFT
  names(AERiskTable) = c("G0", "G1", "G2", "G3", "G4", "G5")

  Gvec = kronecker(matrix(AERiskTable, nrow=1), matrix(nnt-1))
  #str(Gvec)
  Gcum = t(apply(Gvec, 1, cumsum))+1
  colnames(Gcum) <- c("G0", "G1", "G2", "G3", "G4", "G5")
  #str(Gcum)

  for (i in 1:ncol(Gcum))
    lines(x = RSvector, y = Gcum[, i], col = boxcolors[i])

  rect(0, 0, 50, 1, col = 'green')

  for (i in 1:ncol(Gcum)) {
    if (i == 1)
      polygon(x = c(0:100, 100:0), y = c(rep(1, 101), Gcum[101:1, 1]),
              border = FALSE, col = boxcolors[1])
    else
      polygon(x = c(0:100, 100:0), y = c(Gcum[1:101, i-1], Gcum[101:1, i]),
              border = FALSE, col = boxcolors[i])
  }

  points(x = RSinput, y = nnt[RSinput + 1], type = "h", lwd = 3, col = "CornflowerBlue")
  graphics::text(x = RSinput, y = nnt[RSinput + 1], "RS", col = "CornflowerBlue",
                 cex = 2, adj = c(0,0))
  legendDF = data.frame(legend = c("NNT", "User-selected RS", "Benefitted", "No AE", "Mild", "Moderate", "Severe", "Life-threatening", "Died", "-------"),
                        colors = c("black", "CornflowerBlue", "green", boxcolors, "white"),
                        lwd = 9,
                        stringsAsFactors = FALSE
  )
  legendDF = legendDF[ c(2, 1, 10:3), ]
  legend("topright", legend = legendDF$legend,
         cex=1,
         text.col = legendDF$colors,
         lwd=legendDF$lwd, col = legendDF$colors
  )
  ### Remember to remove the helped patient!  NNT-1.
}
