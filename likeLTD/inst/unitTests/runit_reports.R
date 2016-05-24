## Test unit 'reports'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests//runit_reports.R"
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
  # Sets threads to 2 or less. This is a CRAN requirement.
  if(.Call(likeLTD::.cpp.nbthreads) > 2) {
    nb_threads_in_test = .Call(likeLTD::.cpp.nbthreads)
    .Call(likeLTD::.cpp.set_nbthreads, as.integer(2))
  }
}

.tearDown <-
function () {
	## Specific actions for svUnit: clean up context
	if ("package:svUnit" %in% search()) {
		.Log$..Unit <- ""
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
  # Reset number of threads to what it was.
  if('nb_threads_in_test' %in% ls()) {
    .Call(likeLTD::.cpp.set_nbthreads, as.integer(nb_threads_in_test))
    rm(nb_threads_in_test)
  }
}


##############################################
# Then two functions to help design unit tests
##############################################

temporary.directory <- function(expr) {
  # Creates a temporary working directory.
  #
  # The temp directory is created. The working directory is set to it. The
  # expression is evaluated. Whatever happens, the directory should be
  # removed when we exit this function.
  startdir = getwd()
  tryCatch( { 
              directory = tempfile()
              dir.create(directory, recursive=TRUE)
              setwd(directory)
              ;expr;
            },
            finally = {
              try(unlink(directory, TRUE, TRUE), silent=TRUE)
              setwd(startdir)
            } )
}

checkNoException <- function(expr, msg="") {
  foundError <- FALSE
  tryCatch(expr, error=function(n) foundError <- TRUE)
  checkTrue(!foundError, msg=msg)
}


###################################
# Finally, the unit-test themselves
###################################

"test_pack.admin.input" <- svTest(function() {
  # Tests pack.admin.input interface
  # 
  # Checks it packs the data as expected.
  # Checks it does not raise an exception if files exist
  # Checks it raises an exception if files do not exist

  find.package.files <- function(path_vector) {
    system.file(Reduce(file.path, path_vector), package="likeLTD")
  }
  dataFiles = list.files(system.file("data",package="likeLTD"))
  databaseName = dataFiles[grep("DNA17-db",dataFiles)]
  
  databaseFile = find.package.files(c("data", databaseName))
  cspFile    = find.package.files(c("extdata", "hammer", "hammer-CSP.csv"))
  refFile      = find.package.files(c("extdata", "hammer", "hammer-reference.csv"))
  caseName   = 'hammer'
  outputPath = 'hammer'

  correct <- function() {
    pack.admin.input( caseName=caseName,		
                      databaseFile=databaseFile,
       		    cspFile=cspFile,
                      refFile=refFile, 
                      outputPath=outputPath)
  				}

  incorrect <- function() {
    pack.admin.input( caseName=caseName,		
                      databaseFile=databaseFile,
       		    cspFile='nonexistant-file.csv', # only need one to be wrong
                      refFile=refFile, 
                      outputPath=outputPath)
  				}

  # Checks it packs the data as expected.
  admin <- correct()
  checkEquals(length(admin), 7)
  checkEquals(admin$caseName, caseName)
  checkEquals(admin$databaseFile, databaseFile)
  checkEquals(admin$cspFile, cspFile)
  checkEquals(admin$refFile, refFile)
  checkEquals(admin$outputPath, outputPath)

  # Checks it does not raise an exception if files exist
  checkNoException(correct(), msg="Files do exist.")

  # Checks it raises an exception if files do not exist
  checkException(incorrect(), msg="One file doesn't exist")
})


ref.data <- function() {
  # Reference profile used throughout the tests. 
  ref = list( D3   = c("14,16", "16,16", "15,17"),
              vWA  = c("15,19", "15,16", "16,19"),
              D16  = c("11,14", "13,13", "12,13"),
              D2   = c("24,25", "20,20", "18,25"),
              D8   = c("12,13", "11,15", "11,13"),
              D21  = c("28,31", "29,30", "29,30"),
              D18  = c("14,17", "17,17", "15,17"),
              D19  = c("15.2,17.2", "12,14", "14,14"),
              TH01 = c("9,9.3", "6,8", "6,7"),
              FGA  = c("22,23", "22,25", "20,22") )
  ref = data.frame( ref, row.names=c('Suspect', 'Victim 1', 'Victim 2'),
                    stringsAsFactors=FALSE )
  return(ref)
} 
csp.data <- function() {
  # Crime Scene Profile data used throughout the tests.

  csp = list( D3=c("14,16", "", "14,16", ""), 
               vWA=c("15,16,19", "", "15,16,17,19", ""),
               D16=c("11,13,14", "", "11,13,14", ""),
               D2=c("20,23,24,25", "", "20,24,25", ""),
               D8=c("11,12,13,15", "", "11,12,13,15", ""),
               D21=c("28,31", "", "28,29,30,31,31.2", ""),
               D18=c("", "", "13,14,16,17", ""),
               D19=c("12,14,15.2,17.2", "", "12,13,14,15.2,17.2",""),
               TH01=c("6,8,9,9.3", "", "6,8,9,9.3", ""),
               FGA=c("22", "", "22,23,25", "") )
  csp = data.frame(csp, stringsAsFactors=FALSE)
  return(csp)
}


test_queried.vs.known <- svTest(function() {
  # Tests reading queried vs known profiles
  # Path to data files.                                 
  path = Reduce(file.path, c("extdata", "hammer", "hammer-reference.csv"))
  path = system.file(path, package="likeLTD")
  # Now read data.
  if(! "queried.vs.known" %in% ls(.GlobalEnv))
    queried.vs.known = getFromNamespace( "queried.vs.known", "likeLTD")
  result = queried.vs.known(path)
  checkEquals(result, c(TRUE, FALSE, FALSE))
})




