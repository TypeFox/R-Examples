doubleDiff <- function(dd, file="doubleDiff.txt") {
   psi <- dd@psi.df
   psi[2:dd@nt,] <- psi[2:dd@nt,] - psi[1:(dd@nt-1),]
   dd@psi.df <- psi[-1,]
   dd@x <- dd@x[-1]
   dd@nt <- length(dd@x)
   write.table(t(dd@psi.df), file=file,
               quote=FALSE, row.names=dd@x2,
               col.names=dd@x)
   dd
}
