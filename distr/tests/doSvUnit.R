# we only run the tests, if svUnit is available
if (require(svUnit, quietly=TRUE)) {
  pkg <- "distr"
  require("distr")
  
  # we must investigate whether R CMD check is running or not
  #   and if the check is running, whether a time limit exists
  RCMDCHECK <- FALSE
  RCMDCHECKCRAN <- FALSE
  
  for (actual.name in names(Sys.getenv())) {
    if (substr(actual.name, 1, 9) == "_R_CHECK_") {
      RCMDCHECK <- TRUE
      
      if (actual.name == "_R_CHECK_TIMINGS_") {
        RCMDCHECKCRAN <- (as.numeric(Sys.getenv("_R_CHECK_TIMINGS_")) > 0)
      }
    }
  }
  
  # we must determine the path for tests in the installation and outside installation
  if (RCMDCHECK) {
    ## Path to unit tests for R CMD check
    ## PKG.Rcheck/tests/../PKG/unitTests
    ## PKG.Rcheck/tests/unitTests
    
    # we determine the two paths
    pathTestsInInstallation <- system.file(package=pkg, "unitTests")
    pathTestsOutsideInstallation <- file.path(getwd(), "unitTests")
  } else {
    ## Path to unit tests for standalone running as script with "PKG/tests" as working directory
    ## PKG/tests/../inst/unitTests
    ## PKG/tests/unitTests
    
    # we determine the two paths
    pathTestsInInstallation <- file.path(getwd(), "..", "inst", "unitTests")
    pathTestsOutsideInstallation <- file.path(getwd(), "unitTests")
  }
  
  print(pathTestsInInstallation)
  print(pathTestsOutsideInstallation)
  
  # it depends whether we want to skip the long running tests or not
  if (RCMDCHECKCRAN) {
    mypkgSuite <- svSuiteList(packages=pkg, dirs=pathTestsInInstallation)
  } else {
    mypkgSuite <- svSuiteList(packages=pkg, dirs=c(pathTestsInInstallation, pathTestsOutsideInstallation))
  }
  
  unlink("report.txt")  # Make sure we generate a new report
  
  print(svSuiteList(packages=FALSE, dirs=c(pathTestsInInstallation, pathTestsOutsideInstallation)))
  
  runTest(mypkgSuite, name = pkg)  # Run them...
  
  ## makeTestListFromExamples is in svUnit 0.7.8 or more
  #doRunExamples <- TRUE
  #svUnitVersion = as.integer(strsplit(installed.packages()[which(installed.packages()[, 'Package'] == "svUnit"), "Version"], "[\\.-]")[[1]])
  #if (svUnitVersion[1] == 0) {
  #  if (svUnitVersion[2] < 7) {
  #    doRunExamples <- FALSE
  #  } else {
  #    if (svUnitVersion[2] == 7)
  #      doRunExamples <- svUnitVersion[3] >= 8
  #  }
  #}
  #if(doRunExamples)
  #  runTest(tryCatch(makeTestListFromExamples(pkg, "../../pkg/man/"), error=function(e) NULL))
  
  
  protocol(Log(), type = "text", file = "report.txt")  # ... and write report
}

