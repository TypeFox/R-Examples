# Hack to avoid NOTES in R CMD check
# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c("Height", "Id","Xmin","Xmax","Text")) ## Neded in function generateEPG
}

## Needed to avoid notes when using data.table in functions:
## calibrateLb
if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c("PCR.Amplicon","Sample.Name","Marker","TPH",
                           "LPH","LB","phSum"))
}