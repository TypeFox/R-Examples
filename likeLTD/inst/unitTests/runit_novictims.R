## Test unit 'novictims'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_novictims.R"
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
  path = Reduce(file.path, c("extdata", "novictim", "hammer-reference.csv"))
  system.file(path, package="likeLTD")
}
csp.data.path <- function() {
  path = Reduce(file.path, c("extdata", "novictim", "hammer-CSP.csv"))
  system.file(path, package="likeLTD")
}


###################################
# Finally, the unit-test themselves
###################################

test_novictim_read.known.profiles <- svTest(function() {
  # Checks reading from profiles without Victim still returns a matrix.
  knownProfiles = read.known.profiles(ref.data.path())
  checkEquals(nrow(knownProfiles), 1)
  checkEquals(ncol(knownProfiles), 11)
  checkTrue(knownProfiles[[1]])
  checkTrue(is.matrix(knownProfiles))
  checkEquals(rownames(knownProfiles), "Suspect")
})
test_novictim_determine.dropout <- svTest(function() {
  # Test that determine.dropout can be fed an empty knownProfiles.
  cspProfile = read.csp.profile(csp.data.path())
  knownProfiles = read.known.profiles(ref.data.path())
  knownProfiles = knownProfiles[!unlist(knownProfiles[, "queried"]), , drop=FALSE]
  
  result = determine.dropout(knownProfiles, cspProfile)
  checkEquals(result, logical(0))
})

test_novictim_regression_prosecution <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "novictim")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'hammer-CSP.csv'),
    refFile = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns = 1,
    doDropin = TRUE,
    ethnic = "EA1",
    adj = 1.0,
    fst = 0.02,
    relatedness = c(0, 0)/4,
    combineRare  = FALSE
  )

  # Create hypothesis for defence and prosecution.
  prosecuHyp = do.call(prosecution.hypothesis, args)

  # Create and call a likelihood function
  prosecuModel <- create.likelihood.vectors(prosecuHyp)

  argsP = list(rcont = c(# 0.57993104483562263, 0.00305578926998365,
                         1.00000000000000000, 0.21751450573279751),
               dropin = 0.13105537723029814,
               degradation = c(# 0.009812112214182461, 0.004750597605281430,
                               0.005091442765999815, 0.000161894965244169),
               locusAdjustment = c(1.021102533224349, 1.049445500333389,
                                   1.061196279429053, 1.099315646424478,
                                   1.070630778996032, 0.926136743765615,
                                   0.871148514301201, 1.092422691938998,
                                   0.906893821819100, 0.901707489767785),
               power = -4.15575811952846,
               dropout = c(0.09167336880621775, 0.00139282368638074))

  newP <- do.call(prosecuModel, argsP)$objectives
  checks = c(3.44999418614564e-01, 1.11714577308687e-03, 6.39008514828798e-03,
             2.38667220196423e-05, 3.34391350299337e-06, 6.61358628290887e-04,
             1.66664947299532e-03, 1.54323319590758e-06, 3.89361671672519e-06,
             7.26597711768314e-03)
  names(checks) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  checkEquals(newP, checks)
})

test_novictim_regression_defence <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "novictim")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'hammer-CSP.csv'),
    refFile = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns = 1,
    doDropin = TRUE,
    ethnic = "EA1",
    adj = 1.0,
    fst = 0.02,
    relatedness = c(0, 0)/4,
    combineRare  = FALSE
  )

  # Create hypothesis for defence and prosecution.
  defenceHyp = do.call(defence.hypothesis, args)

  # Create and call a likelihood function
  defenceModel <- create.likelihood.vectors(defenceHyp)

  argsD = list(rcont = c(# 0.840590259524432, 0.028546165926857,
                         1.000000000000000, 0.946753652037012),
               dropin = 0.158577060475070,
               degradation = c(# 1.38014699093319e-02, 1.76501060235913e-05,
                               2.26110867647914e-03, 2.08110714678324e-04),
               locusAdjustment = c(1.005338141210508, 1.013192796176000,
                                   1.065720552666928, 1.056464566516383,
                                   1.036911990905030, 0.928442923242878,
                                   0.876845227077896, 1.022970368741970,
                                   1.052948710773730, 0.941164722688676),
               power = -4.42426563254477,
               dropout = c(0.415989161376334515, 0.000161235911814838))

  newD <- do.call(defenceModel, argsD)$objectives
  checks = c(1.45452762382673e-02, 1.39059953090237e-03, 4.12069957025512e-03,
             2.23676993312579e-05, 1.39764825435972e-03, 6.55421895552490e-05,
             1.69972376120242e-03, 4.95286954519019e-06, 4.54379488982241e-03,
             1.34217276423583e-03)
  names(checks) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  checkEquals(newD, checks)
})
