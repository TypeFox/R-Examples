## Test unit 'hypothesis'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_hypothesis.R"
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
test_determine.dropout <- svTest(function() {

  cspProfile = read.csp.profile(csp.data.path())
  knownProfiles = read.known.profiles(ref.data.path())
  
  result = determine.dropout(knownProfiles, cspProfile)
  check = list("Suspect"=TRUE, "Victim 1"=TRUE, "Victim 2"=TRUE)
  checkTrue(result == check)
  checkTrue(names(result) == names(check))

  # construct fake matrix with one replicate containing alleles from 1 and 2.
  other = mapply(union, knownProfiles[2, colnames(cspProfile), drop=FALSE],
                 knownProfiles[1, colnames(cspProfile), drop=FALSE])
  other = matrix(other, nrow=1)
  colnames(other) = colnames(cspProfile)
  result = determine.dropout(knownProfiles, other)
  check = list("Suspect"=FALSE, "Victim 1"=FALSE, "Victim 2"=TRUE)
  checkTrue(result == check)
  checkTrue(names(result) == names(check))

  # construct fake matrix with one replicate containing alleles from 1 and 2
  # and another replicate containing alleles from 1 and 3.
  other2 = mapply(union, knownProfiles[3, colnames(cspProfile), drop=FALSE],
                  knownProfiles[1, colnames(cspProfile), drop=FALSE])
  other = rbind(other, other2)
  result = determine.dropout(knownProfiles, other)
  check = list("Suspect"=FALSE, "Victim 1"=TRUE, "Victim 2"=TRUE)
  checkTrue(result == check)
  checkTrue(names(result) == names(check))
})

test_missing.alleles <- svTest(function() {

# Load missing.alleles from likeLTD namespace
  if(! "missing.alleles" %in% ls(.GlobalEnv))
    missing.alleles <- getFromNamespace("missing.alleles", "likeLTD")
  
  alleleDb   = list(D2=matrix(c(2, 2), nrow=14, ncol=2))
  row.names(alleleDb$D2) = c("10", "11", "12", "13", "14", "14.2", "15", "16",
                             "17", "18", "19", "20", "21", "22")
  cspProfile = matrix(list("D2"=list(c("10", "16"), c("14", "16", "17"))))
  colnames(cspProfile) <- c("D2")
  dropoutProfiles = matrix(list("D2"=list(c("10", "16"), c("14", "16"))))
  colnames(dropoutProfiles) <- c("D2")
  queriedProfile = matrix(list("D2"=list(c("10","16"))))
  colnames(queriedProfile) <- c("D2")
  
# Check cspProfile is a matrix
  checkException(missing.alleles(alleleDb,list("D2"=list(c("10", "16"), c("14", "16", "17"))),queriedProfile,dropoutProfiles))

# Check dropoutProfiles are a matrix
  checkException(missing.alleles(alleleDb,cspProfile,queriedProfile,list("D2"=list(c("10", "16"), c("14", "16")))))

# Check queriedProfile is a matrix
  checkException(missing.alleles(alleleDb,cspProfile,list("D2"=list(c("10","16"))),dropoutProfiles))
 
# Check colnames are present
  tmpcspProfile = cspProfile; colnames(tmpcspProfile) = NULL
  checkException(missing.alleles(alleleDb,tmpcspProfile,queriedProfile,dropoutProfiles))
		
# Does not change if no missing alleles
  checkEquals(alleleDb,missing.alleles(alleleDb,cspProfile,queriedProfile,dropoutProfiles))

# 
  cspProfile = matrix(list("D2"=list(c("10", "16"), c("14", "16", "999"))))
  colnames(cspProfile) <- c("D2")
  tmp = missing.alleles(alleleDb,cspProfile,queriedProfile,dropoutProfiles)
# Check 999 is in new alleleDb
  checkTrue("999" %in% rownames(tmp$D2))
# Check have only added 999
  checkEquals(alleleDb$D2, tmp$D2[-which(rownames(tmp$D2)=="999"),])
# Check values of 999 are correct
  checkEquals(tmp$D2["999",][1], 1) 
  checkEquals(tmp$D2["999",][2], 0) 

# Check behaviours of NA
  cspProfile = matrix(list("D2"=list(NA, NA)))
  colnames(cspProfile) <- c("D2")
  checkEquals(alleleDb,missing.alleles(alleleDb,cspProfile,queriedProfile,dropoutProfiles))

  cspProfile = matrix(list("D2"=list(NA)))
  colnames(cspProfile) <- c("D2")
  checkEquals(alleleDb,missing.alleles(alleleDb,cspProfile,queriedProfile,dropoutProfiles))

# Check behaviour of NA and exceptional
  queriedProfile = matrix(list("D2"=list(c("10","999"))))
  colnames(queriedProfile) <- c("D2")
  tmp = missing.alleles(alleleDb,cspProfile,queriedProfile,dropoutProfiles)
# Check 999 is in new alleleDb
  checkTrue("999" %in% rownames(tmp$D2))
# Check have only added 999
  checkEquals(alleleDb$D2, tmp$D2[-which(rownames(tmp$D2)=="999"),])
# Check values of 999 are correct
  checkEquals(tmp$D2["999",][1], 1) 
  checkEquals(tmp$D2["999",][2], 0) 

  dropoutProfiles = matrix(list("D2"=list(c("10", "16"), c("14", "998"))))
  colnames(dropoutProfiles) <- c("D2")
  tmp = missing.alleles(alleleDb,cspProfile,queriedProfile,dropoutProfiles)
# Check 999 and 998 is in new alleleDb
  checkTrue("999" %in% rownames(tmp$D2))
  checkTrue("998" %in% rownames(tmp$D2))
# Check have only added 999 and 998
  checkEquals(alleleDb$D2, tmp$D2[-which(rownames(tmp$D2)%in%c("999","998")),])
# Check values of 999 and 998 are correct
  checkEquals(tmp$D2["999",][1], 1) 
  checkEquals(tmp$D2["999",][2], 0) 
  checkEquals(tmp$D2["998",][1], 1) 
  checkEquals(tmp$D2["998",][2], 0) 
})


test_add.args.to.hypothesis <- svTest(function() {

  start = list(a="hello", b="world")
  if(! "add.args.to.hypothesis" %in% ls(.GlobalEnv))
    add.args.to.hypothesis <- getFromNamespace("add.args.to.hypothesis", "likeLTD")
  result = add.args.to.hypothesis(start, c=5, a=2)
  
  checkEquals(result, list(a=2, b="world", c=5))
})

test_read.csp.profile <- svTest(function() {
  
  data = list( D3S1358=list(c("14", "16"), c("14", "16")),
               vWA=list(c("15", "16", "19"), c("15", "16", "17", "19")),
               D16S539=list(c("11", "13", "14"), c("11", "13", "14")),
               D2S1338=list(c("20", "23", "24", "25"), c("20", "24", "25")),
               D8S1179=list(c("11", "12", "13", "15"), c("11", "12", "13", "15")),
               D21S11=list(c("28", "31"), c("28", "29", "30", "31", "31.2")),
               D18S51=list(character(0), c("13", "14", "16", "17")),
               D19S433=list(c("12", "14", "15.2", "17.2"), c("12", "13", "14", "15.2", "17.2")),
               TH01=list(c("6", "8", "9", "9.3"), c("6", "8", "9", "9.3")),
               FGA=list(c("22"), c("22", "23", "25")) )
  if(! "read.csp.profile" %in% ls(.GlobalEnv))
    read.csp.profile = getFromNamespace("read.csp.profile", "likeLTD")
  # Path to data files.                                 
  for(filename in c("hammer", "space")) {
    path = Reduce(file.path, c("extdata", "hammer", paste(filename, "-CSP.csv", sep="")))
    path = system.file(path, package="likeLTD")
    # Now read data.
    result = read.csp.profile(path)
  
    # And check it
    checkTrue(is.matrix(result))
    checkEquals(nrow(result), 2)
    checkEquals(ncol(result), length(data))
    checkTrue(setequal(names(data), colnames(result)))
    for(col in names(data)) {
      checkEquals(data[[col]][[1]], result[[1, col]])
      checkEquals(data[[col]][[2]], result[[2, col]])
    }
  }
})

test_read.known.profiles <- svTest(function() {
  
  data = list( queried=c(TRUE, FALSE, FALSE),
               D3S1358  =list(c("14", "16"), c("16", "16"), c("15", "17")),
               vWA =list(c("15", "19"), c("15", "16"), c("16", "19")),
               D16S539 =list(c("11", "14"), c("13", "13"), c("12", "13")),
               D2S1338  =list(c("24", "25"), c("20", "20"), c("18", "25")),
               D8S1179  =list(c("12", "13"), c("11", "15"), c("11", "13")),
               D21S11 =list(c("28", "31"), c("29", "30"), c("29", "30")),
               D18S51 =list(c("14", "17"), c("17", "17"), c("15", "17")),
               D19S433 =list(c("15.2", "17.2"), c("12", "14"), c("14", "14")),
               TH01=list(c("9", "9.3"), c("6", "8"), c("6", "7")),
               FGA =list(c("22", "23"), c("22", "25"), c("20", "22")) )
  if(! "read.known.profiles" %in% ls(.GlobalEnv))
    read.known.profiles = getFromNamespace("read.known.profiles", "likeLTD")
  # Path to data files.                                 
  for(filename in c("hammer", "space")) {
    path = Reduce(file.path, c("extdata", "hammer", paste(filename, "-reference.csv", sep="")))
    path = system.file(path, package="likeLTD")
    # Now read data.
    result = read.known.profiles(path)
  
    # And check it
    checkTrue(is.matrix(result))
    checkEquals(nrow(result), 3)
    checkEquals(ncol(result), length(data))
    checkTrue(setequal(names(data), colnames(result)))
    checkTrue(setequal(c("Suspect", "Victim 1", "Victim 2"), rownames(result)))
    for(col in names(data)) {
      checkEquals(data[[col]][[1]], result[["Suspect", col]])
      checkEquals(data[[col]][[2]], result[["Victim 1", col]])
      checkEquals(data[[col]][[3]], result[["Victim 2", col]])
    }
  }
})
