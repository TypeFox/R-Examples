## Test unit 'nodropin'
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

test_dropin_genotype_sizes <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "nodropin")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'hammer-CSP.csv'),
    refFile = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns = 1,
    doDropin = FALSE,
    ethnic = "EA1",
    adj = 1.0,
    fst = 0.02,
    relatedness = c(0, 0)/4,
    combineRare  = FALSE
  )

  # Create hypothesis for defence and prosecution.
  hypothesis = do.call(defence.hypothesis, args)

  # Create and call a likelihood function
  model <- create.likelihood.vectors(hypothesis, addAttr=TRUE)
  functions <- attr(model, "functions")
  
  checkEquals(ncol(attr(functions[["D3S1358"]],   "constructs")$genotypes),  512) 
  checkEquals(ncol(attr(functions[["vWA"]],  "constructs")$genotypes),  512) 
  checkEquals(ncol(attr(functions[["D16S539"]],  "constructs")$genotypes),  169) 
  checkEquals(ncol(attr(functions[["D2S1338"]],   "constructs")$genotypes),  397) 
  checkEquals(ncol(attr(functions[["D8S1179"]],   "constructs")$genotypes), 1331) 
  checkEquals(ncol(attr(functions[["D21S11"]],  "constructs")$genotypes),   72) 
  checkEquals(ncol(attr(functions[["D18S51"]],  "constructs")$genotypes),   78) 
  checkEquals(ncol(attr(functions[["D19S433"]],  "constructs")$genotypes),   72) 
  checkEquals(ncol(attr(functions[["TH01"]], "constructs")$genotypes),  127) 
  checkEquals(ncol(attr(functions[["FGA"]],  "constructs")$genotypes), 2744) 
})

test_dropin_regression_prosecution <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "nodropin")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'hammer-CSP.csv'),
    refFile = file.path(datapath, 'hammer-reference.csv'),
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

  argsP = list(rcont = c(0.57993104483562263, 0.00305578926998365,
                         1.00000000000000000, 0.21751450573279751),
               degradation = c(0.009812112214182461, 0.004750597605281430,
                               0.005091442765999815, 0.000161894965244169),
               locusAdjustment = c(1.021102533224349, 1.049445500333389,
                                   1.061196279429053, 1.099315646424478,
                                   1.070630778996032, 0.926136743765615,
                                   0.871148514301201, 1.092422691938998,
                                   0.906893821819100, 0.901707489767785),
               power = -4.15575811952846,
               dropout = c(0.09167336880621775, 0.00139282368638074))

  newP <- do.call(prosecuModel, argsP)$objectives
  checks = c(3.82454205615816e-01, 1.71125534624295e-01, 2.68681620184101e-01,
             1.32857835412624e-04, 5.23329092369162e-01, 1.40904975902118e-02,
             1.63050400130066e-07, 1.76168459367054e-01, 4.17460002336186e-01,
             6.47510311051041e-02)
  names(checks) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  checkEquals(newP, checks)
})

test_dropin_regression_defence <- svTest(function() {
  # Regression test over cases without victims
  datapath = file.path(system.file("extdata", package="likeLTD"), "nodropin")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile = file.path(datapath, 'hammer-CSP.csv'),
    refFile = file.path(datapath, 'hammer-reference.csv'),
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

  argsD = list(rcont = c(0.840590259524432, 0.028546165926857,
                         1.000000000000000, 0.946753652037012),
               degradation = c(1.38014699093319e-02, 1.76501060235913e-05,
                               2.26110867647914e-03, 2.08110714678324e-04),
               locusAdjustment = c(1.005338141210508, 1.013192796176000,
                                   1.065720552666928, 1.056464566516383,
                                   1.036911990905030, 0.928442923242878,
                                   0.876845227077896, 1.022970368741970,
                                   1.052948710773730, 0.941164722688676),
               power = -4.42426563254477,
               dropout = c(0.415989161376334515, 0.000161235911814838))

  newD <- do.call(defenceModel, argsD)$objectives
  checks = c(1.87064469330605e-02, 1.13464140439553e-02, 8.56137314546463e-03,
             1.64066893120773e-07, 6.89490152498635e-02, 7.27356001552238e-04,
             1.87654288015654e-04, 4.10268638515752e-04, 9.86365555197756e-02,
             2.81428078288624e-03)
  names(checks) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  checkEquals(newD, checks)
})
