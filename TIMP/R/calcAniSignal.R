"calcAniSignal" <- function (m, plotoptions)
{
 # need to implement analysis for more than 3 datasets simult.
  perpind <- which(m[[1]]@anispec$angle == "PERP")[1]
  parind <- which(m[[1]]@anispec$angle == "PAR")[1]
  
  anisig <- (m[[parind]]@psi.df - m[[perpind]]@psi.df)/ 
    (m[[parind]]@psi.df + (2 * m[[perpind]]@psi.df))
  
  write.table(anisig, file=paste(plotoptions@makeps,
                        "calculated_ani_signal.txt", sep=""),
              row.names = m[[1]]@x, 
              col.names = m[[1]]@x2, quote=FALSE) 
}
