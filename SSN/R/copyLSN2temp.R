copyLSN2temp <- function() {
  ## Create temporary .ssn directory to work with
  old.wd <- getwd()
  setwd(system.file("lsndata/MiddleFork04.ssn",package = "SSN"))
  file.list <- list.files()
  dir.create(paste0(tempdir(),'/MiddleFork04.ssn'), showWarnings = FALSE)
  file.copy(file.list, paste0(tempdir(),'/MiddleFork04.ssn'), recursive = TRUE)
  setwd(old.wd)
}

