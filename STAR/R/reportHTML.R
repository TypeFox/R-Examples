reportHTML <- function(object,
                       filename,
                       extension,
                       directory,
                       Title,
                       ...)
### Generic reportHTML method definition
### object should contain data like a spikeTrain object
### or a repeatedTrain object on whih some standard analysis
### should be carried out first before exporting the results,
### that is, numerical and graphical summaries as well as,
### if necessary actual R objects onto the disk. The summaries
### are formated in html for easy inspection.
### Arguments filename, extension, directory and Title are passed
### to HTMLInitFile function (from R2HTML), directory corresponds
### to argument outdir of the latter.
{

  UseMethod("reportHTML")

}
