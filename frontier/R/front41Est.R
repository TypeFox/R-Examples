front41Est <- function( command = ifelse( .Platform$OS.type == "windows",
      "front41.exe", "front41.bin" ), ... ) {

   # write input file for FRONTIER 4.1
   front41Ins <-front41WriteInput( ... )

   # estimate the model by FRONTIER 4.1
   fullCommand <- paste( command, front41Ins$insFile )
   front41Messages <- system( fullCommand, intern = TRUE )

   # read output of FRONTIER 4.1
   front41Out <- front41ReadOutput( front41Ins$outFile )

   front41Out$input    <- front41Ins
   front41Out$messages <- front41Messages

   return( front41Out )
}

