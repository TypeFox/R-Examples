## ---- results = "hide", echo=FALSE, message = FALSE----------------------
library(matchingR)

## ------------------------------------------------------------------------
uM = matrix(c(1.0, 0.5, 0.0,
              0.5, 0.0, 0.5,
              0.0, 1.0, 1.0), nrow = 3, ncol = 3, byrow = TRUE)

## ---- echo=FALSE---------------------------------------------------------
dimnames(uM) = list(rows = c('Woman 1', 'Woman 2', 'Woman 3'), cols = c('Man 1', 'Man 2', 'Man 3'))
uM

## ------------------------------------------------------------------------
uW = matrix(c(0.0, 1.0, 0.0,
              0.5, 0.0, 0.5,
              1.0, 0.5, 1.0), nrow = 3, ncol = 3, byrow = TRUE)

## ---- echo=FALSE---------------------------------------------------------
dimnames(uW) = list(rows = c('Man 1', 'Man 2', 'Man 3'), cols = c('Woman 1', 'Woman 2', 'Woman 3'))
uW

## ------------------------------------------------------------------------
prefM = matrix(c(1, 3, 3,
                 2, 1, 2,
                 3, 2, 1), nrow = 3, ncol = 3, byrow = TRUE)

## ---- echo=FALSE---------------------------------------------------------
dimnames(prefM) = list(rows = c('Woman 1', 'Woman 2', 'Woman 3'), cols = c('Man 1', 'Man 2', 'Man 3'))
prefM

## ------------------------------------------------------------------------
prefW = matrix(c(3, 1, 3,
                 2, 3, 2,
                 1, 2, 1), nrow = 3, ncol = 3, byrow = TRUE)

## ---- echo=FALSE---------------------------------------------------------
dimnames(prefW) = list(rows = c('Man 1', 'Man 2', 'Man 3'), cols = c('Woman 1', 'Woman 2', 'Woman 3'))
prefW

## ------------------------------------------------------------------------
matching = galeShapley.marriageMarket(uM, uW)

## ------------------------------------------------------------------------
matching = galeShapley.marriageMarket(proposerPref = prefM, reviewerPref = prefW)

## ---- echo=FALSE---------------------------------------------------------
dimnames(matching$proposals) = list(rows = c("Man 1", "Man 2", "Man 3"),
                                    cols = c("Proposed to Woman"))

## ------------------------------------------------------------------------
matching$proposals

## ---- echo=FALSE---------------------------------------------------------
dimnames(matching$engagements) = list(rows = c("Woman 1", "Woman 2", "Woman 3"),
                                      cols = c("Engaged to Man"))

## ------------------------------------------------------------------------
matching$engagements

## ------------------------------------------------------------------------
galeShapley.checkStability(uM, uW, matching$proposals, matching$engagements)

## ------------------------------------------------------------------------
# set seed
set.seed(1)
# set number of men
nmen = 2500
# set number of women
nwomen = 2000
# generate preferences
uM = matrix(runif(nmen*nwomen), nrow = nwomen, ncol = nmen)
uW = matrix(runif(nmen*nwomen), nrow = nmen, ncol = nwomen)
# male-optimal matching
resultsM = galeShapley.marriageMarket(uM, uW)
str(resultsM)
galeShapley.checkStability(uM, uW, resultsM$proposals, resultsM$engagements)
# female-optimal matching
resultsW = galeShapley.marriageMarket(uW, uM)
str(resultsW)
galeShapley.checkStability(uW, uM, resultsW$proposals, resultsW$engagements)

## ------------------------------------------------------------------------
# set seed
set.seed(1)
# set number of students
nstudents = 1000
# set number of colleges
ncolleges = 400
# generate preferences
uStudents = matrix(runif(ncolleges*nstudents), nrow = ncolleges, ncol = nstudents)
uColleges = matrix(runif(nstudents*ncolleges), nrow = nstudents, ncol = ncolleges)
# student-optimal matching
results = galeShapley.collegeAdmissions(studentUtils =  uStudents, collegeUtils =  uColleges, slots = 2)
str(results)
# check if matching is stable
galeShapley.checkStability(uStudents, uColleges, results$matched.students, results$matched.colleges)

## ------------------------------------------------------------------------
pref = matrix(c(3, 6, 2, 5, 3, 5,
                4, 5, 4, 2, 1, 1,
                2, 4, 5, 3, 2, 3,
                6, 1, 1, 6, 4, 4,
                5, 3, 6, 1, 6, 2), nrow = 5, ncol = 6, byrow = TRUE)

## ------------------------------------------------------------------------
roommate.checkPreferences(pref)

## ------------------------------------------------------------------------
results = roommate(pref = pref)
results

## ------------------------------------------------------------------------
roommate.checkStability(pref = pref, matching = results)

## ------------------------------------------------------------------------
# generate preferences
N = 10
u = matrix(runif(N^2),  nrow = N, ncol = N)
results = roommate(utils = u)
results
roommate.checkStability(utils = u, matching = results)

## ------------------------------------------------------------------------
set.seed(1)
N = 512
u = matrix(runif(N^2),  nrow = N, ncol = N)
results = roommate(utils = u)
print(results)

## ------------------------------------------------------------------------
pref = matrix(c(4, 4, 2, 4,
                2, 1, 1, 1,
                1, 2, 3, 3,
                3, 3, 4, 2), nrow = 4, ncol = 4, byrow = TRUE)

## ------------------------------------------------------------------------
results = toptrading(pref = pref)
results

## ------------------------------------------------------------------------
toptrading.checkStability(pref = pref, matchings = results)

