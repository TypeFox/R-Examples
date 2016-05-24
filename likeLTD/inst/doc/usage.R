### R code from vignette source 'usage.Rnw'

###################################################
### code chunk number 1: setthreads
###################################################
  if(.Call(likeLTD::.cpp.nbthreads) > 2) {
    .Call(likeLTD::.cpp.set_nbthreads, as.integer(2))
  }


###################################################
### code chunk number 2: inputs (eval = FALSE)
###################################################
##   require(likeLTD)
## 
##   # Case we are going to be evaluating
##   caseName = "Laboratory"
##   datapath = file.path(system.file("extdata", package="likeLTD"),
##       'laboratory')
## 
##   # File paths and case name for allele report
##   admin = pack.admin.input.peaks(
##       peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
##       refFile = file.path(datapath, 'laboratory-reference.csv'),
##       caseName = caseName,
##       detectionThresh = 20
##       )


###################################################
### code chunk number 3: alleleReport (eval = FALSE)
###################################################
##   # Generate allele report
##   allele.report.peaks(admin)


###################################################
### code chunk number 4: detectThresh (eval = FALSE)
###################################################
## # Specifying locus specific detection thresholds
## admin = pack.admin.input.peaks(
##             peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
##             refFile = file.path(datapath, 'laboratory-reference.csv'),
##             caseName = "Laboratory",
##             detectionThresh = list(D10S1248=20,vWA=20,D16S539=20,D2S1338=20, # blue  
##                                    D8S1179=30,D21S11=30,D18S51=30, # green
##                                    D22S1045=40,D19S433=40,TH01=40,FGA=40, # black  
##                                    D2S441=50,D3S1358=50,D1S1656=50,D12S391=50,SE33=50) # red
##             )


###################################################
### code chunk number 5: hyps (eval = FALSE)
###################################################
##   # Enter arguments
##   args = list(
##         nUnknowns = 1,
##         doDropin = FALSE,
##         ethnic = "NDU1",
##         adj = 1,
##         fst = 0.03,
##         relationship = 0
##         )
## 
##   # Create hypotheses
##   hypP = do.call(prosecution.hypothesis.peaks, append(admin,args))
##   hypD = do.call(defence.hypothesis.peaks, append(admin,args))


###################################################
### code chunk number 6: params (eval = FALSE)
###################################################
##   # Generate likelihood functions and optimisation parameters
##   paramsP = optimisation.params.peaks(hypP,verbose=FALSE)
##   paramsD = optimisation.params.peaks(hypD,verbose=FALSE)
## 
##   # reduce number of iterations for demonstration purposes
##   paramsP$control$itermax=25
##   paramsD$control$itermax=25


###################################################
### code chunk number 7: optimise (eval = FALSE)
###################################################
##   # Run optimisation
##   results = evaluate.peaks(paramsP, paramsD, n.steps=1, 
##       converge=FALSE)


###################################################
### code chunk number 8: output (eval = FALSE)
###################################################
##   # Generate output report
##   output.report.peaks(hypP,hypD,results)


###################################################
### code chunk number 9: optimisedResults (eval = FALSE)
###################################################
##   results = list(Def=list(optim=list(bestmem = c(
##     -3.459508761,     -2.599777111,     -3.385884129,    902.540872861, 
##    171.908479331,    993.204056721,     49.207185386,      0.005254712, 
##      1.016233455,      1.368054382,      0.984844585,      1.371468905, 
##      0.535551035,      1.230762849,      1.358293511,      1.264302595, 
##      0.747843015,      0.775721009,      1.000498903,      0.908954169, 
##      0.822063630,      0.946761244,      0.985732220,      1.135499511, 
##      0.001433327,      0.001742212))))
##   names(results$Def$optim$bestmem) = c(
##     "degradation1",     "degradation2",     "degradation3",         "DNAcont1", 
##         "DNAcont2",         "DNAcont3",            "scale",        "gradientS", 
##  "gradientAdjust1",  "gradientAdjust2",  "gradientAdjust3",  "gradientAdjust4", 
##  "gradientAdjust5",  "gradientAdjust6",  "gradientAdjust7",  "gradientAdjust8", 
##  "gradientAdjust9", "gradientAdjust10", "gradientAdjust11", "gradientAdjust12", 
## "gradientAdjust13", "gradientAdjust14", "gradientAdjust15", "gradientAdjust16", 
##            "meanD",            "meanO") 


###################################################
### code chunk number 10: likely (eval = FALSE)
###################################################
##   # Get the most likely single-contributor genotypes
##   gensMarginal = get.likely.genotypes.peaks(hypD,paramsD,
##       results$Def)
##   # Return joint genotypes and probabilities
##   gensJoint = get.likely.genotypes.peaks(hypD,paramsD,
##       results$Def,joint=TRUE)


###################################################
### code chunk number 11: posterior (eval = FALSE)
###################################################
##   # Get the posterior likelihoods for all genotype combinations
##   gensPosterior = get.likely.genotypes.peaks(hypD,paramsD,
##       results$Def,posterior=TRUE)


###################################################
### code chunk number 12: run
###################################################
  require(likeLTD)

  # Case we are going to be evaluating
  caseName = "Laboratory"
  datapath = file.path(system.file("extdata", package="likeLTD"),
      'laboratory')

  # File paths and case name for allele report
  admin = pack.admin.input.peaks(
      peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
      refFile = file.path(datapath, 'laboratory-reference.csv'),
      caseName = caseName,
      detectionThresh = 20
      )
  # Enter arguments
  args = list(
        nUnknowns = 1,
        doDropin = FALSE,
        ethnic = "NDU1",
        adj = 1,
        fst = 0.03,
        relationship = 0
        )

  # Create hypotheses
  hypP = do.call(prosecution.hypothesis.peaks, append(admin,args))
  hypD = do.call(defence.hypothesis.peaks, append(admin,args))
  # Generate likelihood functions and optimisation parameters
  paramsP = optimisation.params.peaks(hypP,verbose=FALSE)
  paramsD = optimisation.params.peaks(hypD,verbose=FALSE)

  # reduce number of iterations for demonstration purposes
  paramsP$control$itermax=25
  paramsD$control$itermax=25
  results = list(Def=list(optim=list(bestmem = c(
    -3.459508761,     -2.599777111,     -3.385884129,    902.540872861, 
   171.908479331,    993.204056721,     49.207185386,      0.005254712, 
     1.016233455,      1.368054382,      0.984844585,      1.371468905, 
     0.535551035,      1.230762849,      1.358293511,      1.264302595, 
     0.747843015,      0.775721009,      1.000498903,      0.908954169, 
     0.822063630,      0.946761244,      0.985732220,      1.135499511, 
     0.001433327,      0.001742212))))
  names(results$Def$optim$bestmem) = c(
    "degradation1",     "degradation2",     "degradation3",         "DNAcont1", 
        "DNAcont2",         "DNAcont3",            "scale",        "gradientS", 
 "gradientAdjust1",  "gradientAdjust2",  "gradientAdjust3",  "gradientAdjust4", 
 "gradientAdjust5",  "gradientAdjust6",  "gradientAdjust7",  "gradientAdjust8", 
 "gradientAdjust9", "gradientAdjust10", "gradientAdjust11", "gradientAdjust12", 
"gradientAdjust13", "gradientAdjust14", "gradientAdjust15", "gradientAdjust16", 
           "meanD",            "meanO") 


###################################################
### code chunk number 13: diagnose (eval = FALSE)
###################################################
##   # Plot CSP with most likely genotypes
##   peaks.results.plot(hypD,results$Def,replicate=1)


###################################################
### code chunk number 14: plot
###################################################
  # Plot CSP with most likely genotypes
  peaks.results.plot(hypD,results$Def,replicate=1)


