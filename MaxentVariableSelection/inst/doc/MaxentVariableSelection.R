## ------------------------------------------------------------------------
occurrencelocations <- system.file("extdata", "Occurrencedata.csv",
                                   package="MaxentVariableSelection")
occurrencelocations <- read.csv(occurrencelocations,header=TRUE)
head(occurrencelocations)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages('raster')

## ---- warning=FALSE, eval=FALSE------------------------------------------
#  # load the raster package
#  library(raster)
#  
#  # Load the occurrence records
#  occurrencelocations <- system.file("extdata", "Occurrencedata.csv",
#                                     package="MaxentVariableSelection")
#  occurrencelocations <- read.csv(occurrencelocations,header=TRUE)
#  LonLatData <- occurrencelocations[,c(2,3)]
#  
#  # Then load the environmental variables into R with the help of the
#  # stack function of the 'raster' package. You can not just copy the
#  # following line but have to adjust the filepath to your own.
#  
#  files <- list.files("/home/alj/Downloads/BioORACLEVariables",pattern='asc',
#  full.names=TRUE)
#  Grids <- raster::stack(files)
#  
#  # Extracting the variables for all occurrencelocations
#  VariablesAtOccurrencelocations <- raster::extract(Grids,LonLatData)
#  
#  # Combining the extracted values with the longitude and latitude values
#  Outfile <- as.data.frame(cbind("Fucusdistichus", LonLatData,
#                                 VariablesAtOccurrencelocations))
#  colnames(Outfile) <- c("species","longitude","latitude",
#                         colnames(VariablesAtOccurrencelocations))
#  
#  
#  #writing this table to a csv file:
#  
#  write.csv(Outfile, file =
#  "VariablesAtOccurrencelocations.csv", append = FALSE,sep = ",", eol =
#  "\n", na = "NA", dec = ".",col.names = TRUE,row.names=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  VariableSelection(maxent, outdir, gridfolder,
#  occurrencelocations, backgroundlocations, additionalargs,
#  contributionthreshold, correlationthreshold, betamultiplier)

## ------------------------------------------------------------------------
maxent <- ("/home/alj/Downloads/maxent.jar")


## ---- eval=FALSE---------------------------------------------------------
#  maxent <- ("C:/.../maxent.jar")

## ------------------------------------------------------------------------
outdir <- ("/home/alj/Downloads/OutputDirectory")

## ------------------------------------------------------------------------
gridfolder <- ("/home/alj/Downloads/BioORACLEVariables")

## ------------------------------------------------------------------------
occurrencelocations <- system.file("extdata", "Occurrencedata.csv",
                                   package="MaxentVariableSelection")

## ---- eval=FALSE---------------------------------------------------------
#  occurrencelocations <- "/home/alj/Downloads/Occurrencedata.csv"

## ------------------------------------------------------------------------
backgroundlocations <- system.file("extdata", "Backgrounddata.csv",
                                   package="MaxentVariableSelection")

## ---- eval=FALSE---------------------------------------------------------
#  backgroundlocations <- "/home/alj/Downloads/Backgrounddata.csv"

## ------------------------------------------------------------------------
additionalargs="nolinear noquadratic noproduct nothreshold noautofeature"

## ------------------------------------------------------------------------
contributionthreshold <- 5 

## ------------------------------------------------------------------------
correlationthreshold <- 0.9

## ------------------------------------------------------------------------
 betamultiplier=seq(2,6,0.5)

## ---- eval=FALSE---------------------------------------------------------
#  library("MaxentVariableSelection")
#  VariableSelection(maxent,
#                    outdir,
#                    gridfolder,
#                    occurrencelocations,
#                    backgroundlocations,
#                    additionalargs,
#                    contributionthreshold,
#                    correlationthreshold,
#                    betamultiplier
#                    )

