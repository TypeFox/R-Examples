#
# MigClim.validate: Validate a genetic clusters migration output file.
#
MigClim.validate <- function (validateFile="Validation.txt", nrPoints=0,
                              simFile="out1.asc", nrClusters=4)
{
  #
  # Call the validation C function.
  #
  validate <- .C("validate", validateFile, as.integer(nrPoints), simFile,
                 as.integer(nrClusters), score=double(2))
  return (validate$score)
}

