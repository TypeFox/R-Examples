# test that the fame package is working

#  (at BoC I need environment  export FAME=/apps/fame92r2 ) for fame
# export FAME=/apps/fame92r2

if(identical(as.logical(Sys.getenv("_R_CHECK_HAVE_FAME_")), TRUE)) {

cat("**************   test tis\n")
  require("tis")

  seriesA  <- tis(1:24,	   start = c(2002, 1), freq = 12)
  seriesC <- tis(rnorm(104), start = c(2002, 1), tif = "wfriday")

cat("**************   test fame\n")

  require("fame")

cat("**************        connecting zotTestFame.db\n")

  scratch.db <- paste(getwd(),"/zotTestFame.db", sep="")
  unlink(scratch.db, recursive = TRUE)

  if(!fameRunning()) fameStart(workingDB = FALSE)
  Id <- try(fameDbOpen(scratch.db, accessMode = "create"))
  if(inherits(Id, "try-error") ) 
      stop("Could not establish TSfameConnection to ", scratch.db)

  fameDbClose(Id) 

cat("**************        writing to zotTestFame.db\n")

  Id <- fameDbOpen(scratch.db, accessMode = "update")
  ok <- 0 == fameWriteSeries(Id, "seriesA", seriesA,
  		       update=FALSE, checkBasisAndObserved=FALSE)

  ok <- ok & 0==fameWriteSeries(Id, "seriesC", seriesC,
  		       update=FALSE, checkBasisAndObserved=FALSE)

  fameDbClose(Id) 

  if(!ok) stop("Could not write series.")

cat("**************        reading from zotTestFame.db\n")

   r <- getfame("seriesA", scratch.db, save = FALSE, envir = parent.frame(),
             start = NULL, end = NULL, getDoc = FALSE)[[1]]
   if(1e-10 < max(abs(r - seriesA))) stop("Error reading seriesA.")
   if(any(start(r) != start(seriesA))) stop("Error reading seriesA.")

   r <- getfame("seriesC", scratch.db, save = FALSE, envir = parent.frame(),
             start = NULL, end = NULL, getDoc = FALSE)[[1]]
   if(1e-10 < max(abs(r - seriesC))) stop("Error reading seriesC.")
   if(any(start(r) != start(seriesC))) stop("Error reading seriesC.")

cat("**************        delete\n")
  unlink(scratch.db, recursive = TRUE)

} else  cat("FAME not available. Skipping tests.")
