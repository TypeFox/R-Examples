"boa.plot.par" <-
function(mfdim = c(1, 1), title = TRUE)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   top <- ifelse(title, 3 + mfdim[1], 0)

   val <- boa.par("par")
   val$mfrow <- mfdim
   val$mfcol <- NULL
   val$oma <- c(0, 0, top, 0)
   par(val)
       
   invisible()
}
