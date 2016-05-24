## Test unit 'objectives per locus'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_objectivePerLocus.R"
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


test_empty.alleles = svTest(function() {
  
  if(! "all.genotypes.per.locus" %in% ls(.GlobalEnv))
    all.genotypes.per.locus <- getFromNamespace("all.genotypes.per.locus",
                                                "likeLTD")
  if(! "empty.alleles" %in% ls(.GlobalEnv))
    empty.alleles <- getFromNamespace("empty.alleles", "likeLTD")
  genotypes = matrix(nrow=4, ncol=0)
  dropoutProfs = cbind(c(TRUE, FALSE, FALSE, TRUE),
                       c(FALSE, TRUE, FALSE, FALSE))
  result = empty.alleles(genotypes, dropoutProfs, 0)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 0)
  checkTrue(nrow(result) == 4)
  
  result = empty.alleles(genotypes, dropoutProfs, 1)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 0)
  checkTrue(nrow(result) == 4)

  genotypes = all.genotypes.per.locus(4, 2)
  result = empty.alleles(genotypes, dropoutProfs, 0)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == ncol(genotypes))
  checkTrue(nrow(result) == 4)
  checkTrue(all(result == c(FALSE, FALSE, TRUE, FALSE)))

  result = empty.alleles(genotypes, dropoutProfs, 1)
  indices = c(1, 2, 4, 5, 7, 10, 11, 12, 14, 15, 17, 20, 31, 32, 34, 35, 37,
              40, 41, 42, 44, 45, 47, 50, 61, 62, 64, 65, 67, 70, 91, 92, 94,
              95, 97, 100)
  checkTrue(ncol(result) == ncol(genotypes))
  checkTrue(nrow(result) == 4)
  checkTrue(all(!result[1, ]))
  checkTrue(all(!result[2, ]))
  checkTrue(all(!result[4, ]))
  checkTrue(all(result[3, indices]))
  result[3, indices] = FALSE
  checkTrue(all(!result))
})

test_TH01.regression.with.dropin = svTest(function() {
  # Case we are going to be looking at.
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 1,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(0.0, 0),
    combineRare  = FALSE
  )
  if(! "defence.hypothesis" %in% ls(.GlobalEnv))
    defence.hypothesis <- getFromNamespace("defence.hypothesis", "likeLTD")
  hypothesis = do.call(defence.hypothesis, args)
  if(! "transform.to.locus.centric" %in% ls(.GlobalEnv))
    transform.to.locus.centric <-
      getFromNamespace("transform.to.locus.centric", "likeLTD")
  hypothesisTH01 = transform.to.locus.centric(hypothesis)$TH01
  
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565),
                   dropin = 1e0, #0.108695652173913,
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.175, 0.105) )

  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")
  hypothesisTH01$nUnknowns = 2
  hypothesisTH01$doDropin = TRUE
  objective.function <- create.likelihood.per.locus(hypothesisTH01)

  checkEquals(do.call(objective.function, arguments), 0.0204764693571788)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 0.0017054449886482)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 0.0147496568283615)
})

test_TH01.regression.no.dropin = svTest(function() {
  # Case we are going to be looking at.
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 1,
    doDropin     = FALSE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(0.0, 0),
    combineRare  = FALSE
  )
  if(! "defence.hypothesis" %in% ls(.GlobalEnv))
    defence.hypothesis <- getFromNamespace("defence.hypothesis", "likeLTD")
  hypothesis = do.call(defence.hypothesis, args)
  if(! "transform.to.locus.centric" %in% ls(.GlobalEnv))
    transform.to.locus.centric <-
      getFromNamespace("transform.to.locus.centric", "likeLTD")
  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")
  hypothesisTH01 = transform.to.locus.centric(hypothesis)$TH01

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565), 
                   dropin=1e0,
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.175, 0.105) )

  hypothesisTH01$nUnknowns = 2
  hypothesisTH01$doDropin = FALSE
  objective.function <- create.likelihood.per.locus(hypothesisTH01)

  checkEquals(do.call(objective.function, arguments), 0.0194512081797547)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 0.00166299033316889)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 0.0140009673609186)
})

test_D18.regression.with.dropin = svTest(function() {
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 1,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(0.0, 0),
    combineRare  = FALSE
  )
  if(! "defence.hypothesis" %in% ls(.GlobalEnv))
    defence.hypothesis <- getFromNamespace("defence.hypothesis", "likeLTD")
  hypothesis = do.call(defence.hypothesis, args)
  if(! "transform.to.locus.centric" %in% ls(.GlobalEnv))
    transform.to.locus.centric <-
      getFromNamespace("transform.to.locus.centric", "likeLTD")
  hypothesisD18 = transform.to.locus.centric(hypothesis)$D18
  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")

  
  hypothesisD18$nUnknowns = 0
  hypothesisD18$doDropin = TRUE
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565),
                   dropin=1e0,
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.175, 0.105) )
  objective.function <- create.likelihood.per.locus(hypothesisD18)

  arguments$rcont=c(0.923913043478261, 0.565217391304348)
  arguments$degradation=rep(3e-3, 2)
  checkEquals(do.call(objective.function, arguments), 9.09495530595364e-06)
  arguments$degradation = rep(2e-2, 2)
  checkEquals(do.call(objective.function, arguments), 5.75530620496417e-05)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047)
  checkEquals(do.call(objective.function, arguments), 5.13626195539709e-05)

  hypothesisD18$nUnknowns = 1
  objective.function <- create.likelihood.per.locus(hypothesisD18)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600) 
  arguments$degradation = rep(3e-3, 3)
  arguments$rcont=c(0.923913043478261, 0.565217391304348, 1e0)
  checkEquals(do.call(objective.function, arguments), 9.06669340994184e-06)
  arguments$degradation = rep(2e-2, 3)
  checkEquals(do.call(objective.function, arguments), 4.27972472968122e-05)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600) 
  checkEquals(do.call(objective.function, arguments), 3.77630662967064e-05)

  hypothesisD18$nUnknowns = 2
  objective.function <- create.likelihood.per.locus(hypothesisD18)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  arguments$rcont=c(0.923913043478261, 0.565217391304348, 1e0,
                    0.543478260869565)
  checkEquals(do.call(objective.function, arguments), 1.3537562256385e-05)
})


test_relatedness.factors <- svTest(function() {

  alleleDb   = matrix(c(1, 0), nrow=14, ncol=2)
  row.names(alleleDb) = c("10", "11", "12", "13", "14", "14.2", "15", "16",
                          "17", "18", "19", "20", "21", "22")
  
  queriedAlleles = matrix(list(c("11", "21")), nrow=1, ncol=1)
  alleleDb[2, 1] = 2e0
  alleleDb[13, 1] = 3e0
  
  if(! "all.genotypes.per.locus" %in% ls(.GlobalEnv))
    all.genotypes.per.locus <- getFromNamespace("all.genotypes.per.locus",
                                                "likeLTD")
  if(! "relatedness.factors" %in% ls(.GlobalEnv))
    relatedness.factors <- getFromNamespace("relatedness.factors", "likeLTD")

  genotypes = all.genotypes.per.locus(nrow(alleleDb), 2)
  # No relatedness
  result = array(1.0, ncol(genotypes))
  result = relatedness.factors(result, genotypes, alleleDb, queriedAlleles,
                               c(0, 0))
  checkEquals(result, array(1.0, ncol(genotypes)))
  # First, only one relatedness 
  result = array(1.0, ncol(genotypes))
  result = relatedness.factors(result, genotypes, alleleDb, queriedAlleles,
                               c(0.9, 0))
  checkEquals(length(result), ncol(genotypes))
  na = nrow(alleleDb)
  nComb = na * (na+1) / 2
  nCont = 2
  nTrue = (2*na - 1) * nComb^(nCont-1)
  checkEquals(sum(abs(result - 1.0 * (1.0 - sum(c(0.9, 0)))) > 1e-8), nTrue)

  hasFirst = genotypes[1, ] %in% 2 | genotypes[2, ] %in% 2
  hasSecond = genotypes[1, ] %in% 13 | genotypes[2, ] %in% 13
  het = genotypes[1, ] == genotypes[2, ]

  r = result[hasFirst & (!hasSecond) & (!het)]
  check = 1e0 - 0.9  + 0.9 * 0.5 * 0.5 / alleleDb[2, 1]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  r = result[hasFirst & (!hasSecond) & het]
  check = 1e0 - 0.9 + 0.9 * 0.5 / alleleDb[2, 1]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  r = result[(!hasFirst) & hasSecond & (!het)]
  check = 1e0 - 0.9 + 0.9 * 0.5 * 0.5 / alleleDb[13, 1]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  check = 1e0 - 0.9 + 0.9 * 0.5 / alleleDb[13, 1]
  r = result[(!hasFirst) & hasSecond & het]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)
  
  r = result[hasFirst & hasSecond]
  check = 1e0 - 0.9 + 0.9 * 0.5 * 0.5 / alleleDb[2, 1] +
                      0.9 * 0.5 * 0.5 / alleleDb[13, 1]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  result = array(1.0, ncol(genotypes))
  result = relatedness.factors(result, genotypes, alleleDb, queriedAlleles,
                               c(0, 0.9))
  check = 1e0 - 0.9 + 0.9 * 0.5 / alleleDb[2, 1] / alleleDb[13, 1]
  r = result[hasFirst & hasSecond]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)
})
