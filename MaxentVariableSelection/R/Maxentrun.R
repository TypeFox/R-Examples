Maxentrun <-
function(maxent,outdir,gridfolder,occurrencesites,backgroundsites,additionalargs){

     system(paste("java -jar", maxent,
      "-o", outdir,
      "-j", gridfolder,
      "-s", occurrencesites,
      "-e", backgroundsites,
      additionalargs,sep=" "),
            intern=TRUE,
            ignore.stdout = TRUE,
             show.output.on.console = FALSE,
            )
      }


# jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
