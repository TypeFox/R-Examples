readData <- function(filename, markers = "ALL", n.vars = 0, sep = " ", allele.sep = ";", na.strings = "NA", filehash = FALSE,
design = "triad"){
##
## READS RAW DATA IN HAPLIN FORMAT. PREPARES AN info OBJECT FROM CALL
## USES f.read.data TO READ AND PROCESS. 
## TO BE CALLED USING ARGUMENTS IDENTICAL TO HAPLIN'S FILE SPECIFICATIONS
##
## NOTE: use.missing IS ALWAYS SET TO TRUE HERE 
## USER CAN SELECT MARKERS, BUT NOT USE sel.sex OR SET use.missing  = F

# IMPLEMENTER FILEHASH!
# PROBLEM: HVIS MAN GLEMMER n.vars FAAR MAN BESKJED OM FEIL I allele.sep

## CREATE info OBJECT
.info <- f.catch(match.call(), formals())
.info$model$use.missing <- TRUE # OVERRIDE DEFAULT
#
## READ FULL DATA FILE
if(.info$control$verbose)	cat("\nReading data from file...  ")
.data.read <- f.read.data(info = .info) 
if(.info$control$verbose)	cat("Done\n")
#
##
return(.data.read)
}
