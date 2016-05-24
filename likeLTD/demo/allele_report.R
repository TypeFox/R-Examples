# Demo for running allele report.

# Load library!
library(likeLTD) 

# Get path to libray, cos that's where the data is.
likeLTD.path <- path.package('likeLTD')

# Case we are going to be looking at.
caseName = 'hammer'
datapath <- file.path(file.path(likeLTD.path, 'extdata'), caseName)
# Construct input: crime scene profile
mixedFile = file.path(datapath, 'hammer-CSP.csv')
# Construct input: reference profiles
refFile = file.path(datapath, 'hammer-reference.csv')
# Construct input: output path in the R temp directory for now
outputPath = tempdir()

# Construct list of all administrative input
admin = pack.admin.input( caseName=caseName,
                          cspFile=cspFile,
                          refFile=refFile,
                          outputPath=outputPath )
# Finally call allele.report
allele.report(admin)
