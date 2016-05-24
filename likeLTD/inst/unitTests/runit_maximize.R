## Test unit 'maximize'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_maximize.R"
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

###############################################################
# Then data functions
###############################################################

ref.data.path <- function() {
  path = Reduce(file.path, c("extdata", "hammer", "hammer-reference.csv"))
  system.file(path, package="likeLTD")
}
csp.data.path <- function() {
  path = Reduce(file.path, c("extdata", "hammer", "hammer-CSP.csv"))
  system.file(path, package="likeLTD")
}

###################################
# Finally, the unit-test themselves
###################################
test_estimates <- svTest(function() {

  cspProfile = read.csp.profile(csp.data.path())
  knownProfiles = read.known.profiles(ref.data.path())

  if(! "estimates" %in% ls(.GlobalEnv))
    estimates <- getFromNamespace("estimates", "likeLTD")

  checkEquals(array(c(0.575, 0.500)),
              estimates(knownProfiles[1, ], cspProfile))
  checkEquals(array(c(0.625, 0.500)),
              estimates(knownProfiles[2, ], cspProfile))
  checkEquals(array(c(0.750, 0.675)),
              estimates(knownProfiles[3, ], cspProfile))
})
