## Test unit 'objectives'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_objective.R"
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


test_regression1 <- svTest(function() {
  # Case we are going to be looking at.
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 0,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relationship  = 0,
    combineRare  = FALSE
  )
  if(! "defence.hypothesis" %in% ls(.GlobalEnv))
    defence.hypothesis <- getFromNamespace("defence.hypothesis", "likeLTD")
  hypothesis = do.call(defence.hypothesis, args)
  hypothesis$nUnknowns = 0

  likelihood <- create.likelihood.vectors(hypothesis)
  objectives = c(3.10250372325746e-04, 1.17224578453062e-02,
                 4.76863464507366e-05, 4.29822197384531e-06,
                 1.76350988087800e-03, 1.43444739873622e-08,
                 9.09495530595362e-06, 7.11137407679358e-08,
                 3.01744911691531e-04, 6.89041414610852e-03)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  penalties = c(0.00303192703176332, 0.00303192703176332, 0.00303192703176332,
                0.00303192703176332, 0.00303192703176332, 0.00303192703176332, 
                0.00303192703176332, 0.00303192703176332, 0.00303192703176332,
                0.00303192703176332)

  check = list(objectives=objectives, penalties=penalties)

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348),
                   dropin = 1.0,
                   degradation=rep(3e-3, 2),
                   locusAdjustment=1,
                   power=-4.35, beta=-4.35,
                   dropout=c(0.175, 0.105) )
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, check$objectives)
})

test_regression.zerounknown <- svTest(function() {
  # Case we are going to be looking at.
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    kit = "SGMplus",
    cspFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 0,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(0.0, 0),
    combineRare  = FALSE
  )
  defenceHypothesis = do.call(defence.hypothesis, args)
  prosecutionHypothesis = do.call(prosecution.hypothesis, args)

  likelihood <- create.likelihood.vectors(prosecutionHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0),
                   dropin = 0.543478260869565,
                   degradation=c(3e-3, 3e-3, 3e-3),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.15, 0.01) )
  objectives = c(2.64877312615097e-04, 8.30010229904578e-02,
                 6.57219727208952e-02, 3.55156330421480e-03,
                 5.58990000167620e-01, 2.58579653708631e-05,
                 6.88667211400942e-07, 5.73149656581839e-02,
                 2.70886382744855e-02, 3.64303077108431e-03)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 5.82387768745086e-25)

  likelihood <- create.likelihood.vectors(defenceHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0),
                   dropin = 0.543478260869565,
                   degradation=c(3e-3, 3e-3, 3e-3),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(2.01470581026018e-03, 1.13394042752114e-02, 4.69287546906059e-03,
                 2.47293709587760e-04, 5.90787709164540e-02, 9.97729081224379e-07,
                 6.25814358291258e-06, 1.04230422586896e-04, 1.46838789834590e-02,
                 1.43741838557565e-03)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 6.97617046679607e-32)
})

test_regression.oneunknown <- svTest(function() {
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
  defenceHypothesis = do.call(defence.hypothesis, args)
  prosecutionHypothesis = do.call(prosecution.hypothesis, args)

  likelihood <- create.likelihood.vectors(prosecutionHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0, 
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.15, 0.01) )
  objectives = c(5.24609987450467e-05, 1.28369768074416e-01, 3.19566044368621e-02, 
                 9.37696749226556e-04, 3.54990168257584e-01, 3.61714632221897e-05, 
                 2.78024377866913e-06, 7.42370399489160e-02, 2.15249772203121e-02, 
                 4.27467268573391e-04)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 1.9762409429182e-25)

  likelihood <- create.likelihood.vectors(defenceHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   locusAdjustment=1,
                   power=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(5.40928563410066e-04, 1.24028883413425e-02, 4.32400562747164e-03,
                 3.71369198539134e-04, 7.20410130310739e-02, 1.08896493758788e-06,
                 2.81245481065190e-06, 1.83320557026132e-04, 1.94351600773495e-02,
                 5.75805705402794e-04)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 3.63417836957934e-31)
})

test_regression.relatedness <- svTest(function() {
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
    relationship  = 2,
    combineRare  = FALSE
  )
  defenceHypothesis = do.call(defence.hypothesis, args)
  prosecutionHypothesis = do.call(prosecution.hypothesis, args)

  likelihood <- create.likelihood.vectors(prosecutionHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   locusAdjustment=rep(1,times=10),
                   power=-4.35,
                   dropout=c(0.15, 0.01) )
  objectives = c(5.2460998745046723702e-05, 1.2836976807441610737e-01, 3.1956604436862108554e-02, 
9.3769674922655570342e-04, 3.5499016825758356042e-01, 3.6171463222189653881e-05, 
2.7802437786691249654e-06, 7.4237039948916033749e-02, 2.1524977220312099119e-02, 
4.2746726857339074597e-04)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  penalty = 2.6078646479705969163
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), prod(objectives * penalty), checkNames=FALSE)
  # old result was 1.9762409429182e-25 (different penalties)

  likelihood <- create.likelihood.vectors(defenceHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   locusAdjustment=rep(1,times=10),
                   power=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(2.0324320110552429076e-03, 4.0692306250085193142e-02, 4.4262872552010612548e-02, 
1.5005784145824823013e-03, 1.9188770082595565936e-01, 7.4396398032556112094e-06, 
2.1466961591367787368e-06, 1.3734648258538964885e-02, 5.5349543045990240442e-02, 
2.2848765693775441657e-03)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  penalty = 2.6078646479705969163
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), prod(objectives * penalty), checkNames=FALSE)
  # old result was 5.75589320779353e-25 (different penalties)

  defenceHypothesis$relationship = 2
  likelihood <- create.likelihood.vectors(defenceHypothesis)
  # This also tests that refIndiv works. The objective function should be
  # inserting 1 at position 3 in rcont.
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   locusAdjustment=rep(1,times=10),
                   power=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(2.0324320110552429076e-03, 4.0692306250085193142e-02, 4.4262872552010612548e-02, 
1.5005784145824823013e-03, 1.9188770082595565936e-01, 7.4396398032556112094e-06, 
2.1466961591367787368e-06, 1.3734648258538964885e-02, 5.5349543045990240442e-02, 
2.2848765693775441657e-03)
  names(objectives) = c("D3S1358", "vWA", "D16S539", "D2S1338", "D8S1179", "D21S11", "D18S51", "D19S433",
                        "TH01", "FGA")
  penalty = 2.6078646479705969163
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), prod(objectives * penalty), checkNames=FALSE)
  # old result was 7.07569814196917e-26 (different penalties)
})
