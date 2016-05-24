### --- Test setup ---
`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA = function(x) all(is.na(x))
 
if(FALSE) {
  library("RUnit")
  library("roxygen2")
  library("devtools")
  setwd(".."); setwd(".."); getwd()
  document()
  load_all("./", create = FALSE) # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # tmlenet:::debug_set() # SET TO DEBUG MODE

  setwd("..");
  install("tmlenet", build_vignettes = FALSE) # INSTALL W/ devtools:

  # system("echo $PATH") # see the current path env var
  # system("R CMD Rd2pdf tmlenet")  # just create the pdf manual from help files

  # CHECK AND BUILD PACKAGE:
  getwd()
  # setwd("./tmlenet"); setwd(".."); getwd()
  devtools::check(args = "--as-cran")
  # devtools::check() # runs check with devtools
  # devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs check with devtools
  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  devtools::build()
  devtools::build_win(args = "--as-cran") # build package on CRAN servers (windows os?)
  # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  # devtools::build(args = "--compact-vignettes") # build package tarball compacting vignettes
  # devtools::build(args = "--no-build-vignettes") # build package tarball compacting vignettes
  # devtools::build() # build package tarball
  setwd("..")
  system("R CMD check --as-cran tmlenet_0.0.9.tar.gz") # check R package tar ball prior to CRAN submission
      ## system("R CMD check --no-manual --no-vignettes tmlenet") # check without building the pdf manual and not building vignettes
      ## system("R CMD build tmlenet --no-build-vignettes")
      ## system("R CMD build tmlenet")  
  # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # INSTALLING FROM SOURCE:
  # install.packages("./tmlenet_0.2.0.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # library(tmlenet)
  # tmlenet:::debug_set() # SET TO DEBUG MODE
  # tmlenet:::debug_off() # SET DEBUG MODE OFF

  # To install a specific branch:
  # devtools::install_github('osofr/simcausal', ref = "simnet", build_vignettes = FALSE)
  # devtools::install_github('osofr/tmlenet', ref = "master", build_vignettes = FALSE)

  # TEST COVERATE:
  # if your working directory is in the packages base directory
  # package_coverage()
  # or a package in another directory
  # cov <- package_coverage("tmlenet")
  # view results as a data.frame
  # as.data.frame(cov)
  # zero_coverage() can be used to filter only uncovered lines.
  # zero_coverage(cov)
}

sample_checks <- function() {   # doesn't run, this is just to show what test functions can be used
  print("Starting tests...")
    checkTrue(1 < 2, "check1")     ## passes fine
     ## checkTrue(1 > 2, "check2")  ## appears as failure in the test protocol
     v <- 1:3
     w <- 1:3
     checkEquals(v, w)               ## passes fine
     names(v) <- c("A", "B", "C")
     ## checkEquals(v, w)            ## fails because v and w have different names
     checkEqualsNumeric(v, w)        ## passes fine because names are ignored
     x <- rep(1:12, 2)
     y <- rep(0:1, 12)
     res <- list(a=1:3, b=letters, LM=lm(y ~ x))
     res2 <- list(a=seq(1,3,by=1), b=letters, LM=lm(y ~ x))
     checkEquals( res, res2)        ## passes fine
     checkIdentical( res, res)
     checkIdentical( res2, res2)
     ## checkIdentical( res, res2)  ## fails because element 'a' differs in type
     fun <- function(x) {
       if(x)
       {
        stop("stop conditions signaled")
       }
       return()
     }
     checkException(fun(TRUE))      ## passes fine
     ## checkException(fun(FALSE))  ## failure, because fun raises no error
     checkException(fun(TRUE), silent=TRUE)
     ##  special constants
     ##  same behaviour as for underlying base functions
     checkEquals(NA, NA)
     checkEquals(NaN, NaN)
     checkEquals(Inf, Inf)
     checkIdentical(NA, NA)
     checkIdentical(NaN, NaN)
     checkIdentical(-Inf, -Inf)
}


test.bugfixes <- function() {

}

# Add a bug with automatic interval/bin detection on binary variable (all vals get placed in one bin)
test.bin01bug <- function() {

}

test.nullnodesbug <- function() {
  # When is.null(DatNet$nodes) no checks are made and regression is attempted anyways
  # results in uninterpretable errror
  # ....................................
  # Fixed by setting DatNet.sWsA$get.outvar to throw an error when outvar is not found
  # ....................................
}

test.opts.misfun.chkpkgs <- function() {
  checkException(old_opts <- tmlenet_options(bin.method = "blah"))

  funmiss <- tmlenet:::testmisfun()
  checkTrue(funmiss(NA))
  checkTrue(is.na(tmlenet:::get.misval()))
  
  tmlenet:::set.misval(tmlenet:::gvars, NaN)
  checkTrue(is.nan(tmlenet:::gvars$misval))
  tmlenet:::set.misval(tmlenet:::gvars, NA)
  checkTrue(is.na(tmlenet:::gvars$misval))

  checkException(tmlenet:::checkpkgs("blahblah"))

  warns <- tmlenet:::GetWarningsToSuppress()

  testdat <- data.frame(a = rnorm(5), b = rep("str", 5), stringsAsFactors=TRUE)
  checkTrue(tmlenet:::CheckExistFactors(data = testdat)%in%"b")
}

# making bin indicator matrix from ordinal sVar:
test.makebins.fromord <- function() {
  # testing bin indicator matrix for ordinal variable that's in range 1:nbins
  ord_vec_1 <- sample(c(1:5), 200, replace=TRUE)
  nbins <- length(unique(ord_vec_1))
  # levels <- sort(unique(ord_vec_1))
  bin.nms <- "B_"%+%(1:nbins)
  bin_ord_vec_1 <- tmlenet:::make.bins_mtx_1(ord_vec_1, nbins = nbins, bin.nms = bin.nms)
  out_mat_1 <- cbind(ord_vec_1 = ord_vec_1, bin_ord_vec_1)

  # testing bin indicator matrix for ordinal variable that's not in standard range 1:nbins
  ord_vec_2 <- sample(c(-2,1:5,8,9), 200, replace=TRUE)
  nbins <- length(unique(ord_vec_2))
  levels <- sort(unique(ord_vec_2))
  bin.nms <- "B_"%+%levels
  bin_ord_vec_2 <- tmlenet:::make.bins_mtx_1(ord_vec_2, nbins = nbins, bin.nms = bin.nms, levels = levels)
  out_mat_2 <- cbind(ord_vec_2 = ord_vec_2, bin_ord_vec_2)
}


# helper function for generating some data
get.testDat <- function(nsamp = 100000) {
  `%+%` <- function(a, b) paste0(a, b)
  require(simcausal)
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
  D <- set.DAG(D)
  datO <- sim(D, n = nsamp, rndseed = 12345)
}

# helper function for generating network and tmlenet data storage objects
get.testDatNet <- function(datO) {
  Kmax <- 1
  nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"), nFnode = "nF")
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA")
  netind_cl <- simcausal::NetIndClass$new(nobs = nrow(datO))
  # Define datNetObs:
  datnetW <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sW)
  datnetA <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sA)
  datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
  return(list(datNetObs = datNetObs, netind_cl = netind_cl, def_sA = def_sA, def_sW = def_sW, nodes = nodes))
}

test.RegressionClass <- function() {
  # Tests for RegressionClass:
  reg_test1 <- RegressionClass$new(outvar.class = c(tmlenet:::gvars$sVartypes$bin, tmlenet:::gvars$sVartypes$bin),
                                  outvar = c("A1", "A2"),
                                  predvars = c("W1", "W2"))
  class(reg_test1)
  reg_test1$subset
  model1 <- SummariesModel$new(reg = reg_test1)
  # [1] "Init BinOutModel:"
  # [1] "P(A1|W1,W2)"
  # [1] "Init BinOutModel:"
  # [1] "P(A2|A1,W1,W2)"
  checkTrue("BinOutModel" %in% class(model1$getPsAsW.models()$`P(sA|sW).1`))

  reg_test2 <- RegressionClass$new(outvar.class = c(tmlenet:::gvars$sVartypes$bin, tmlenet:::gvars$sVartypes$bin),
                                  outvar = c("A1", "A2"),
                                  predvars = c("W1", "W2"),
                                  subset = list("A1"))
  class(reg_test2)
  reg_test2$subset
  model2 <- SummariesModel$new(reg = reg_test2)

  checkException(
    reg_test3 <- RegressionClass$new(outvar.class = c(tmlenet:::gvars$sVartypes$cont, tmlenet:::gvars$sVartypes$cont),
                                    outvar = c("sA"), predvars = c("W1", "W2", "W3"),
                                    subset = list(quote(sA==0)))
    )
  reg_test3 <- RegressionClass$new(outvar.class = tmlenet:::gvars$sVartypes$cont,
                                  outvar = c("sA"), predvars = c("W1", "W2", "W3"),
                                  subset = list("sA"))
                                  # subset = list(quote(sA==0)))

  reg_test_new <- reg_test3$clone()
  checkTrue(all.equal(reg_test3, reg_test_new))
  checkTrue("RegressionClass" %in% class(reg_test3))
  nsamp <- 1000
  datO <- get.testDat(nsamp)
  head(datO)
  nodeobjs <- get.testDatNet(datO)
  
  # ContinSummaryModel$new
  # class(reg_test3)
  # Detect intervals with a minimal RegressionClass inputs:
  intervls <- nodeobjs$datNetObs$detect.sVar.intrvls(reg_test3$outvar,
                                                      nbins = reg_test3$nbins,
                                                      bin_bymass = reg_test3$bin_bymass,
                                                      bin_bydhist = reg_test3$bin_bydhist,
                                                      max_nperbin = reg_test3$max_nperbin)
  
  class(reg_test3$subset[[1]])
  model3 <- SummariesModel$new(reg = reg_test3, DatNet.sWsA.g0 = nodeobjs$datNetObs)
}


# ------------------------------------------------------------------------------------------------
# BUG in $fit WHEN USING SUBSET EXPRESSIONS (INSTEAD OF VAR NAMES) FOR CONTINUOUS OUTVAR:
# Tries to add "TRUE" to the current bin indicator name: self$subset_expr = c("sA_B.1", "TRUE")
# ------------------------------------------------------------------------------------------------
test.butregquote <- function() {
  gvars <- tmlenet:::gvars
  reg_test <- RegressionClass$new(outvar.class = c(gvars$sVartypes$cont),
                                  outvar = c("sA"),
                                  predvars = c("W1", "W2", "W3"),
                                  subset = list(quote(TRUE)))
  datO <- get.testDat(nsamp = 1000)
  nodeobjs <- get.testDatNet(datO)
  datNetObs <- nodeobjs$datNetObs
  model3 <- SummariesModel$new(reg = reg_test, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  
  model3$getPsAsW.models()[[1]]$getPsAsW.models()
  model3$getPsAsW.models()[[1]]$reg$subset
  # model3$fit(data = nodeobjs$datNetObs)
}


# ------------------------------------------------------------------------------------------------
# Test for pooled fitting of the bin indicators
# AN IDEA FOR TESTING pooled regression: 
# USE IT TO ESTIMATE POOLED IPTW FOR LONGITUDINAL DATA WITH SEVERAL TIME POINTS (RN's simulation)
# ------------------------------------------------------------------------------------------------
test.PoolContRegression <- function() {
  require(data.table)
  gvars <- tmlenet:::gvars
  reg_test <- RegressionClass$new(outvar.class = c(gvars$sVartypes$cont),
                                  outvar = c("sA"),
                                  predvars = c("W1", "W2", "W3"),
                                  # subset = list(quote(TRUE)))
                                  subset = list("sA"),
                                  pool_cont = TRUE, nbins = 10, bin_bymass = FALSE)
  datO <- get.testDat(nsamp = 1000)
  nodeobjs <- get.testDatNet(datO)
  datNetObs <- nodeobjs$datNetObs
  class(datNetObs) # [1] "DatNet.sWsA" "DatNet"      "R6"
  model3 <- SummariesModel$new(reg = reg_test, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  # Matrix of all summary measures: (sW,sA)
  head(nodeobjs$datNetObs$mat.sVar); class(nodeobjs$datNetObs$mat.sVar)
  head(datNetObs$mat.bin.sVar)
  binfit_time <- system.time(
    model3$fit(data = nodeobjs$datNetObs)
    # Error in 1L:self$nbins : argument of length 0
  )
  binfit_time
  binpredict_time <- system.time(
    probAeqa <- model3$predictAeqa(newdata = nodeobjs$datNetObs)
  )
  binpredict_time
  # [1] "fit (10K)"
  # $coef
  #  Intercept     bin_ID         W1         W2         W3 
  # -2.7756215  0.1553186 -1.0014477 -0.5720651 -0.3339728 
  # [1] "res_DT: "
  #            ID ProbAeqa_long
  #      1:     1    0.97396036
  #      2:     1    0.96971742
  #      3:     1    0.96480811
  #      4:     1    0.95913645
  #      5:     1    0.95259565
  #     ---                    
  # 104496:  9998    0.07668105
  # 104497:  9999    0.93215687
  # 104498:  9999    0.92165035
  # 104499:  9999    0.09032560
  # 104500: 10000    0.06784313
  # [1] "res_DT_short: "
  #           ID    cumprob
  #     1:     1 0.06099655
  #     2:     2 0.06145986
  #     3:     3 0.03836225
  #     4:     4 0.05821479
  #     5:     5 0.07303417
  #    ---                 
  #  9996:  9996 0.05119563
  #  9997:  9997 0.05896735
  #  9998:  9998 0.06414013
  #  9999:  9999 0.07760077
  # 10000: 10000 0.06784313
  # [1] "head(ProbAeqa, 50)"
  #  [1] 0.060996548 0.061459862 0.038362248 0.058214786 0.073034166 0.064140127 0.060658764 0.050023002 0.026039639 0.075325033 0.029168620
  # [12] 0.054538219 0.054031618 0.083549453 0.008653412 0.029594466 0.077600772 0.081220201 0.068319822 0.061459862 0.071357407 0.039453938
  # [23] 0.075325033 0.039007914 0.057871503 0.077600772 0.057871503 0.058967354 0.064140127 0.043973691 0.046655735 0.079794387 0.074434114
  # [34] 0.058967354 0.067843133 0.063492979 0.033237556 0.064138704 0.056974041 0.065426910 0.037236039 0.029168620 0.056974041 0.047226347
  # [45] 0.043973691 0.084256432 0.060173071 0.073034166 0.029168620 0.060183301
}

## ---------------------------------------------------------------------
# Test set for detecting vector types for continuous sVar
## ---------------------------------------------------------------------
test.intervals <- function() {
  test_mat <- as.matrix(data.frame(a = c(0,1,0,0,1), b = rep(5,5), c = c(1,2,3,4,5), d = rnorm(5)))
  correct.types <- list(a = tmlenet:::gvars$sVartypes$bin, 
                        b = tmlenet:::gvars$sVartypes$bin, 
                        c = tmlenet:::gvars$sVartypes$cat, 
                        d = tmlenet:::gvars$sVartypes$cont)
  out.types <- tmlenet:::detect.col.types(test_mat)
  checkTrue(all.equal(correct.types, out.types))

  # ------------------------------------------------------------------------------
  # Simulate some network data
  # ------------------------------------------------------------------------------
  n <- 100
  # NetInd_k <- matrix(NA, nrow=n, ncol=3)
  # NetInd_k[1,] <- c(2,3,NA)
  # NetInd_k[4,] <- c(2,NA,NA)
  # nF <- rep.int(0,n)
  # nF[1] <- 2; nF[4] <- 1
  Kmax <- 3
  NetInd_k <- t(replicate(n, sample(1:n, Kmax, replace = FALSE)))
  nF <- rep.int(Kmax, n)
  ## ---------------------------------------------------------------------
  # TEST SET FOR NetIndClass class
  ## ---------------------------------------------------------------------
  testdf <- data.frame(a = rnorm(100), b = rnorm(100))
  nettest <- simcausal::NetIndClass$new(nobs = nrow(testdf), Kmax = 3)
  nettest$getNetInd

  NetInd_k_toobig <- t(replicate(1000, sample(1:1000, Kmax, replace = FALSE)))
  checkException(nettest$NetInd <- NetInd_k_toobig)
  nettest$NetInd <- NetInd_k
  checkEquals(nettest$NetInd, NetInd_k)

  ## ---------------------------------------------------------------------
  # TESTING ContinSummaryModel class and newsummarymodel.contin constructor
  ## ---------------------------------------------------------------------
  # TEST 1: Binary outvar (treated as continuous). Testing that results match with binary class prediction.
  ## ---------------------------------------------------------------------
  nbins <- 10L
  binvar <- rbinom(n=100, size = 1, prob = 0.3)

  # binning by equal length
  int_bylen <- c(0, 0.1, 0.9, 1)
  ord1 <- findInterval(x = binvar, vec = int_bylen, rightmost.closed = TRUE)
  bins_1 <- tmlenet:::make.bins_mtx_1(x.ordinal = ord1, nbins = nbins, bin.nms = "B_"%+%1:nbins)
  cbind(binvar = binvar, ord1 = ord1, bins_1)

  # binning by equal mass
  int_bymass <- quantile(binvar, prob = c(0, 0.1, 0.9, 1))
  ord2 <- findInterval(x = binvar, vec = int_bymass, rightmost.closed = TRUE)
  bins_2 <- tmlenet:::make.bins_mtx_1(x.ordinal = ord1, nbins = nbins, bin.nms = "B_"%+%1:nbins)
  cbind(binvar = binvar, ord2 = ord2, bins_2)
}


# ------------------------------------------------------------------------------
# Detecting sVar types and intervals for continuous and categorical sVar
# ------------------------------------------------------------------------------
test.detect.int.sA <- function() {
  nsamp <- 1000
  gvars <- tmlenet:::gvars
  # reg_test <- RegressionClass$new(outvar.class = c(gvars$sVartypes$cont),
  #                                 outvar = c("sA"),
  #                                 predvars = c("W1", "W2", "W3"),
  #                                 # subset = list(quote(TRUE)))
  #                                 subset = list("sA"),
  #                                 pool_cont = TRUE, nbins = 10, bin_bymass = FALSE)

  makedat <- function(nsamp, Kmax = 3) {
    datO <- get.testDat(nsamp = nsamp)
    # ------------------------------------------------------------------------------
    # Simulate a network mat
    # ------------------------------------------------------------------------------
    n <- nsamp
    NetInd_k <- t(replicate(n, sample(1:n, Kmax, replace = FALSE)))
    nF <- rep.int(Kmax, n)
    NetInd_k[1,3] <- NA
    NetInd_k[2,3] <- NA
    NetInd_k[3,c(2,3)] <- NA
    NetInd_k[4,c(1:Kmax)] <- NA
    NetInd_k[5,c(3)] <- NA
    NetInd_k[6,c(Kmax:(Kmax-1))] <- NA
    NetInd_k[7,c(Kmax:(Kmax-2))] <- NA
    NetInd_k[8,c(Kmax:(Kmax-3))] <- NA

    head(NetInd_k)
    netind_cl <- simcausal::NetIndClass$new(nobs = nrow(datO), Kmax = Kmax)
    netind_cl$NetInd <- NetInd_k
    netind_cl$make.nF()
    nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"), nFnode = "nF")
    def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
    def_sA <- def.sA(sA = "sA")
    # Define datNetObs:
    datnetW <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sW)
    datnetA <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sA)
    datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
    return(datNetObs)
  }
  
  print("Binning by mass")
  nbins <- 10
  oldopts <- tmlenet_options(maxncats = 5, nbins = nbins)

  # ----------------------------------------------------------------------------------------
  # Continuous
  # ---------------------------------------------------------------------------------------- 
  datNetObs <- makedat(nsamp=nsamp, Kmax=3)
  obsdat.sW <- datNetObs$datnetW$dat.sVar
  print("head(obsdat.sW)"); print(head(obsdat.sW,10))
  head(datNetObs$netind_cl$NetInd, 10)

  print("Var types: "); str(datNetObs$datnetW$type.sVar)
  print("Contin covars names: "); print(datNetObs$datnetW$names.c.sVar)
  defints1 <- datNetObs$detect.sVar.intrvls("sA",
                                            nbins = nbins,
                                            bin_bymass = TRUE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  print("No normalization bin intervals by mass: ");
  print(defints1)

  # ----------------------------------------------------------------------------------------
  # Categorical w ncats < maxncats
  # ----------------------------------------------------------------------------------------
  defints2 <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = nbins,
                                            bin_bymass = TRUE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)

  # Testing ordinals with ncats < nbins get nbins = ncats:
  checkTrue((length(defints2)-1) < nbins)
  # Testing all categories in ordinal are represented for categorical with 3 values (nF):
  checkTrue(all(unique(datNetObs$dat.sWsA[,"nF"])%in%defints2))


  # ----------------------------------------------------------------------------------------
  # ******** THIS IS PRETTY BAD. THE RESULT COLLAPSES nF TO 3/4 categories only ***********
  # Categorical w maxncats < ncats, so gets detected as contin and categories get collapsed
  # ----------------------------------------------------------------------------------------
  datNetObs <- makedat(nsamp=nsamp, Kmax=10)
  tmlenet_options(maxncats = 5, nbins = 10)
  obsdat.sW <- datNetObs$datnetW$dat.sVar
  head(obsdat.sW,10)
  head(datNetObs$netind_cl$NetInd, 10)
  table(obsdat.sW[,"nF"])
  # 0   6   7   8   9  10
  # 1   1   1   2   3 992
  str(datNetObs$datnetW$type.sVar)
 # $ W1: chr "binary"
 # $ W2: chr "binary"
 # $ W3: chr "binary"
 # $ nF: chr "contin"
  defints3a <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = nbins,
                                            bin_bymass = TRUE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  defints3a
  # [1] -1000.000000     0.000000     1.833333     3.666667     5.500000     7.333333     9.166667    11.000000  1011.000000
  defints3b <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = nbins,
                                            bin_bymass = FALSE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  defints3b
  # [1] -1000.000000     0.000000     1.833333     3.666667     5.500000     7.333333     9.166667    11.000000  1011.000000

  defints3c <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = length(unique(obsdat.sW[,"nF"])),
                                            bin_bymass = FALSE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  defints3c
  # [1] -1000.000000     0.000000     1.833333     3.666667     5.500000     7.333333     9.166667    11.000000  1011.000000
  defints3d <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = 2,
                                            bin_bymass = FALSE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  defints3d
  # [1] -1000.0     0.0     5.5    11.0  1011.0
  defints3e <- datNetObs$detect.sVar.intrvls("nF",
                                            nbins = 30,
                                            bin_bymass = FALSE,
                                            bin_bydhist = FALSE,
                                            max_nperbin = 100)
  defints3e
  # [1] -1000.000000     0.000000     1.833333     3.666667     5.500000     7.333333     9.166667    11.000000  1011.000000

 do.call(tmlenet_options, oldopts)
}

test.NetIndClassFromString <- function() {
	# ----------------------------------------------------------------------------------------
	# TESTING NetIndClass CLASS
	# ----------------------------------------------------------------------------------------
	k <- 2
	dftestW <- data.frame(W = as.integer(c(6,7,8,9,10))) # W_netF1 = rep(6,5), W_netF2 = rep(8,5)
	dftestA <- data.frame(A = as.integer(c(1,2,3,4,5))) # , A_netF1 = rep(1,5), A_netF2 = rep(3,5)
	class(dftestW$W)
	class(dftestA$A)

	NET_id <- c(rep("1 3", nrow(dftestW)-1), "1")
	class(dftestW$A)
	dftest1 <- data.frame(dftestW, dftestA, NETID = NET_id, stringsAsFactors = FALSE)
	is.factor(dftest1$NETID)
  netindcl <- simcausal::NetIndClass$new(nobs = nrow(dftest1), Kmax = k)
  netindcl$makeNetInd.fromIDs(Net_str = dftest1[,"NETID"])
	netindcl$NetInd_k
	checkTrue("matrix" %in% class(netindcl$NetInd_k))
	checkTrue("integer" %in% class(netindcl$NetInd_k[,1]))
}

# TESTING sVar expressions parser:
test.DefineSummariesClass <- function() {
  # ----------------------------------------------------------------------------------------
  # TEST DATA:
  # ----------------------------------------------------------------------------------------
  k <- 2
  Kmax <- 2
  dftestW <- data.frame(W = as.integer(c(6,7,8,9,10)))
  dftestA <- data.frame(A = as.integer(c(1,2,3,4,5)))
  NET_id <- c(rep("1 3", nrow(dftestW)-1), "1 ")
  # NET_id <- c(rep("1 3", nrow(dftestW)-1), "1")
  class(dftestW$W)
  class(dftestA$A)

  dfnet <- data.frame(dftestW, W_netF1 = rep(6,5), W_netF2 = c(rep(8,4), NA), dftestA, A_netF1 = rep(1,5), A_netF2 = c(rep(3,4), NA))
  dfnet

  # ----------------------------------------------------------------------------------------
  # TESTING sVar expressions parser
  # ----------------------------------------------------------------------------------------
  dftest <- data.frame(dftestW, dftestA, NETID = NET_id, stringsAsFactors = FALSE)
  netind_cl <- simcausal::NetIndClass$new(nobs = nrow(dftest), Kmax = k)
  netind_cl$makeNetInd.fromIDs(Net_str = dftest[,"NETID"])
  netind_cl$Kmax
  netind_cl$NetInd_k
  netind_cl$nF

  # **** Example TESTING Kmax substitute ****
  defsVar.expr0 <- def.sW(sA.1 = A[[Kmax]])
  evaled.sVar.expr0 <- defsVar.expr0$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(all(as.vector(evaled.sVar.expr0[,"A_netF2"]) %in% c(rep(3,4),NA)))

  defsVar.expr0 <- def.sW(A)
  evaled.sVar.expr0 <- defsVar.expr0$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)

  defsVar.expr0 <- def.sW(A[[0:Kmax]])
  evaled.sVar.expr0 <- defsVar.expr0$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(is.matrix(evaled.sVar.expr0))

  # Example 1.
  defsVar.expr1a <- def.sW(sA.1 = rowSums(A[[0:Kmax]]))
  evaled.sVar.expr1a <- defsVar.expr1a$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  defsVar.expr1b <- def.sW(sA.1 = sum(A[[0:Kmax]]))
  evaled.sVar.expr1b <- defsVar.expr1b$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(all.equal(evaled.sVar.expr1a, evaled.sVar.expr1b))
  # w/ NA for missing vars:
  checkTrue(is.na(evaled.sVar.expr1a[5,1]))

  # Example 2. Using a variable to pass sVar expression.
  # ************************************************************************************************************************
  # ANY FUNCTION INSIDE QUOTE WILL NOT BE MODIFIED BY THE PARSER
  # => Using quote(sum(A[[0:Kmax]])) WILL ACTUALLY apply sum to the entire matrix A[[0:Kmax]] and will produce ONE NUMBER!!!
  # ************************************************************************************************************************
  testexpr_call <- quote(rowSums(A[[0:Kmax]]))
  # testexpr_call <- quote(A[[0]])
  defsVar.expr2 <- def.sW(sA.1 = eval(testexpr_call)) # this works
  # defsVar.expr2 <- def.sW(W = testexpr_call) # doesn't work
  defsVar.expr2$exprs_list
  evaled.sVar.expr2 <- defsVar.expr2$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  res1 <- as.integer(c(5, 6, 7, 8, NA))
  checkTrue(all.equal(as.vector(evaled.sVar.expr2[,1]), res1))

  defsVar.expr1 <- def.sW(sA.1 = sum(A[[0:Kmax]]))
  evaled.sVar.expr1 <- defsVar.expr1$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(all.equal(evaled.sVar.expr1, evaled.sVar.expr2))
  checkTrue(all(res1 %in% as.vector(evaled.sVar.expr1)))

  # Example 3. Generate a matrix of sVar[1], ..., sVar[j] from one sVar expression.
  defsVar.expr1 <- def.sW(W = W[[0:Kmax]])
  evaled.sVar.expr1 <- defsVar.expr1$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)

  # k is not a known constant (Kmax is), will parse this as having two parents: W and k => will throw exception
  defsVar.expr1 <- def.sW(W[[0:k]])
  checkException(evaled.sVar.expr1 <- defsVar.expr1$eval.nodeforms(data.df = dftest, netind_cl = netind_cl))
  # correct way to do above (k is defined in this environment):
  defsVar.expr1a <- def.sW(netW = W[[0:k]], replaceNAw0=TRUE)
  evaled.sVar.expr1a <- defsVar.expr1a$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  # this will parse as having only one parent (Kmax is a known constant) => no exception
  defsVar.expr1b <- def.sW(W[[0:Kmax]], replaceNAw0=TRUE)
  evaled.sVar.expr1b <- defsVar.expr1b$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(!any((evaled.sVar.expr1a-evaled.sVar.expr1b)!=0))

  checkTrue(is.matrix(evaled.sVar.expr1))
  checkTrue(all(evaled.sVar.expr1[,1] == dftest$W))

  # Example 4a. Generate a matrix of sVar[1], ..., sVar[j] from one sVar expression that is a combination of different Vars in Odata.
  defsVar.expr <- def.sA(sA.1 = W[[0:Kmax]] + sum(A[[1:Kmax]]), replaceNAw0 = TRUE)
  evaled.sVar.expr <- defsVar.expr$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  class(evaled.sVar.expr)
  colnames(evaled.sVar.expr)

  testres1_cl <- def.sA(netW = W[[0:Kmax]], replaceNAw0 = TRUE)
  evaled.testres1 <- testres1_cl$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  testres2_cl <- def.sA(sA.1 = rowSums(A[[1:Kmax]]), replaceNAw0 = TRUE)
  evaled.testres2 <- testres2_cl$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(all((evaled.testres1 + as.vector(evaled.testres2)) == evaled.sVar.expr))

  # Example 4b. Generate a matrix of sVar[1], ..., sVar[j] from one sVar expression that is a combination of different Vars in Odata.
  defsVar.expr <- def.sW(W = "W[[0:Kmax]] + rowSums(A[[1:Kmax]])", replaceNAw0 = TRUE)
  class(defsVar.expr$sVar.exprs[["W"]])
  defsVar.expr$sVar.exprs[["W"]]
  evaled.sVar.expr2 <- defsVar.expr$eval.nodeforms(data.df = dftest,  netind_cl = netind_cl)
  checkTrue(all(evaled.sVar.expr2[,c(1,2,3)]==evaled.sVar.expr))

  # Example 5. sum of prod of netA and netW:
  defsVar.expr <- def.sA(sumAWnets = sum(A[[1:Kmax]] * W[[1:Kmax]]), replaceNAw0 = TRUE)
  evaled.sVar.expr <- defsVar.expr$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)
  checkTrue(all(as.integer(as.vector(evaled.sVar.expr)) == c(30,30,30,30,6)))

  # Example 6. More than one summary measure
  defsVar.expr <- def.sA(A = A, sumAnets = sum(A[[1:Kmax]]), sumAWnets = sum(A[[1:Kmax]] * W[[1:Kmax]]), replaceNAw0 = TRUE)
  evaled.sVar.expr <- defsVar.expr$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)

  # Example 7. No names
  defsVar.expr <- def.sA(A, sumAnets = sum(A[[1:Kmax]]), sumAWnets = sum(A[[1:Kmax]] * W[[1:Kmax]]))
  evaled.sVar.expr <- defsVar.expr$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)

  defsVar.expr <- def.sA(A[[0:Kmax]], sumAnets = sum(A[[1:Kmax]])) + def.sA(sumAWnets = sum(A[[1:Kmax]] * W[[1:Kmax]]))
  evaled.sVar.expr <- defsVar.expr$eval.nodeforms(data.df = dftest, netind_cl = netind_cl)

  defsVar.expr <- def.sA(A, sumAnets = sum(A[[1:Kmax]])) + def.sA(sum(A[[1:Kmax]] * W[[1:Kmax]]), replaceNAw0 = TRUE)
  checkException(evaled.sVar.expr <- defsVar.expr$eval.nodeforms(data.df = dftest, netind_cl = netind_cl))
}