## Test unit 'nodropout'
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


###################################
# Finally, the unit-test themselves
###################################
test_dropout_regression_prosecution <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "nodropout")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'CSP.csv'),
    refFile = file.path(datapath, 'reference.csv'),
    nUnknowns = 1,
    doDropin = FALSE,
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

  argsP = list( locusAdjustment=c( 0.980883952583997, 0.980464467636125,
                                   0.980243910764109, 0.980280452794204,
                                   0.980771919430644, 0.980489238622726,
                                   0.981327128481614, 0.979711058314168,
                                   0.980826168761904, 0.980795415358621 ),
                power=-4.35195406059307,
                dropout=c(1e-04, 1e-04),
                degradation=0.00110893573762344,
                dropin=1e-04 )

  newP <- do.call(prosecuModel, argsP)$objectives
  checks = c(0.23340789850456878, 0.11616629711845862, 0.14594965969702192,
             0.11182655329884779, 0.06969260540520721, 0.05045585876831847,
             0.00135846108675573, 0.38503530829085297, 0.12637764623617148,
             0.00805568848771995)
  names(checks) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  checkEquals(newP, checks)
})

test_dropout_regression_defence <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "nodropout")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'CSP.csv'),
    refFile = file.path(datapath, 'reference.csv'),
    nUnknowns = 1,
    doDropin = FALSE,
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

  argsD = list( locusAdjustment=c( 0.979178047455148, 0.981265026384822,
                                   0.981380871718712, 0.980188885321909,
                                   0.979999998544042, 0.979999469147536,
                                   0.979178390673367, 0.979820452384997,
                                   0.979274819180087, 0.979181441869338 ),
                 power=-4.36675942037325,
                 dropout=c(1e-04, 1e-04),
                 degradation=c(0.00146422147787705, 0.00178103375279096),
                 rcont=1.68052278396074, 
                 dropin=1e-4 )

  newD <- do.call(defenceModel, argsD)$objectives
  checks = c(1.14050706665492e-01, 1.08871059314561e-02, 4.99753850183928e-02,
             2.41384344770380e-02, 2.39365424629160e-02, 7.74146362090388e-03,
             7.62049153669593e-05, 1.22837189292898e-01, 6.98334653935746e-02,
             2.64909789599108e-03)
  names(checks) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  checkEquals(newD, checks)
})
