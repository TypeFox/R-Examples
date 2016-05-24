## ---- echo = FALSE-------------------------------------------------------
pkgUrl <- gsub("\n", "", packageDescription("seroincidence")$URL)

## ------------------------------------------------------------------------
# Load package "seroincidence"
library(seroincidence)

## ---- eval=FALSE---------------------------------------------------------
#  # List all objects (functions and data) exposed by package "seroincidence"
#  ls("package:seroincidence")

## ---- echo=FALSE---------------------------------------------------------
ls("package:seroincidence")

## ------------------------------------------------------------------------
# Show first rows of "salmonellaSerologyData" data.frame
head(salmonellaSerologyData)

# Show first rows of "campylobacterSerologyData" data.frame
head(campylobacterSerologyData)

## ------------------------------------------------------------------------
# Assign data.frame "salmonellaSerologyData" to object named "serologyData"
serologyData <- salmonellaSerologyData

## ------------------------------------------------------------------------
# Assign output of function "simulateSerologyData" to object named "serologyData"
serologyData <- simulateSerologyData(n = 300)

# Show first rows of object "serologyData"
head(serologyData)

## ---- eval=FALSE---------------------------------------------------------
#  # Read content of file "c:\\cross-sectional-data.csv" into object named "serologyData"
#  serologyData <- read.csv(file = "c:\\cross-sectional-data.csv")

## ------------------------------------------------------------------------
# Show first rows of data.frame "A" in list "campylobacterResponseParams"
head(campylobacterResponseParams$A)

# Show first rows of data.frame "k" in list "campylobacterResponseParams"
head(campylobacterResponseParams$k)

## ------------------------------------------------------------------------
responseParams <- campylobacterResponseParams

## ------------------------------------------------------------------------
responseParams <- simulateSalmonellaResponseParams()

## ----eval=FALSE----------------------------------------------------------
#  AData <- read.csv(file = "A.csv")
#  kData <- read.csv(file = "k.csv")
#  
#  # Create a list named "responseData" containing objects named "A" and "k"
#  responseParams <- list(A = AData, k = kData)

## ---- eval=FALSE---------------------------------------------------------
#  # Assign output of function "estimateSeroincidence" to object named "seroincidenceData"
#  serologyData <- salmonellaSerologyData
#  responseParams <- simulateSalmonellaResponseParams()
#  seroincidenceData <- estimateSeroincidence(data = serologyData,
#                                  antibodies = c("IgG", "IgM", "IgA"),
#                                  strata = "sex",
#                                  Ak = responseParams,
#                                  censorLimits = list(IgG = 0.25, IgM = 0.25, IgA = 0.25))
#  
#  # Show content of the output variable
#  print(seroincidenceData)
#  # or simply type in the console: 'seroincidenceData' (without "'") and press ENTER

## ---- eval=FALSE---------------------------------------------------------
#  censorLimits <- list(IgG = 0, IgM = 0, IgA = 0)

## ---- eval=FALSE---------------------------------------------------------
#  summary(seroincidenceData)

## ---- eval=FALSE---------------------------------------------------------
#  # Compute seroincidence summary and assign to object "seroincidenceSummary"
#  seroincidenceSummary <- summary(seroincidenceData)
#  
#  # Show the results
#  seroincidenceSummary$Results

## ---- eval=FALSE---------------------------------------------------------
#  # Calculate seroincidence rates for Salmonella
#  
#  # 1. Define cross-sectional data
#  serologyData <- salmonellaSerologyData
#  
#  # 2. Define longitudinal response data (simulate 1000 observations)
#  responseParams <- simulateSalmonellaResponseParams(n = 1000)
#  
#  # 3. Define cut-offs
#  cutoffs <- list(IgG = 0.25, IgM = 0.25, IgA = 0.25)
#  
#  # 4a. Calculate seroincidence rates by age (triplet of titres)...
#  seroincidenceData <- estimateSeroincidence(data = serologyData,
#                          antibodies = c("IgG", "IgM", "IgA"),
#                          strata = "age",
#                          Ak = responseParams,
#                          censorLimits = cutoffs,
#                          showProgress = TRUE)
#  
#  # 4b. ...or calculate a single seroincidence rate for all serum samples (triplet of titres)...
#  seroincidenceData <- estimateSeroincidence(data = serologyData,
#                          antibodies = c("IgG", "IgM", "IgA"),
#                          strata = "",
#                          Ak = responseParams,
#                          censorLimits = cutoffs)
#  
#  # 4c. ...or calculate a single seroincidence rate for a single serum sample (triplet of titres)...
#  seroincidenceData <- estimateSeroincidence(data = serologyData[1, ],
#                          antibodies = c("IgG", "IgM", "IgA"),
#                          strata = "",
#                          Ak = responseParams,
#                          censorLimits = cutoffs)
#  
#  # 4d. ...or calculate a single seroincidence rate for all serum samples (only IgG)
#  seroincidenceData <- estimateSeroincidence(data = serologyData,
#                          antibodies = c("IgG"),
#                          strata = "",
#                          Ak = responseParams,
#                          censorLimits = cutoffs)
#  
#  # 5a. Produce summary of the results with 2.5% and 97.5% bounds...
#  summary(seroincidenceData)
#  
#  # 5b. ...or produce summary of the results with 5% and 95% bounds, do not show convergence...
#  summary(seroincidenceData, quantiles = c(0.05, 0.95), showConvergence = FALSE)
#  
#  # 5c. ...or produce summary and assign to an object...
#  seroincidenceSummary <- summary(seroincidenceData)
#  # ...and work with the results object from now on (here: display the results).
#  seroincidenceSummary$Results

