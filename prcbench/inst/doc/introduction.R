## ------------------------------------------------------------------------
library(prcbench)

## A single tool
toolsetA <- create_toolset("ROCR")

## Multiple tools
toolsetB <- create_toolset(c("PerfMeas", "PRROC"))
                           
## Tool sets can be manually combined to a single set
toolsetAB <- c(toolsetA, toolsetB)

## ------------------------------------------------------------------------
library(prcbench)

## A single tool - lower case
toolsetA2 <- create_toolset("rocr")

## Multiple tools - lower case and partially matched 
toolsetB2 <- create_toolset(c("perf", "prr"))

## ------------------------------------------------------------------------
## Use 'set_names'
toolsetC <- create_toolset(set_names = "auc5")

## Multiple sets are automatically combined to a single set
toolsetD <- create_toolset(set_names = c("auc5", "crv4"))

## ------------------------------------------------------------------------
## A balanced data set with 50 positives and 50 negatives
testset1A <- create_testset("bench", "b100")

## An imbalanced data set with 2500 positives and 7500 negatives
testset1B <- create_testset("bench", "i10k")

## Test data sets can be manually combined to a single set
testset1AB <- c(testset1A, testset1B)

## Multiple sets are automatically combined to a single set
testset1C <- create_testset("bench", c("i10", "b10"))

## ------------------------------------------------------------------------
## C1 test set
testset2A <- create_testset("curve", "c1")

## C2 test set
testset2B <- create_testset("curve", "c2")

## Test data sets can be manually combined to a single set
testset2AB <- c(testset2A, testset2B)

## Multiple sets are automatically combined to a single set
testset2C <- create_testset("curve", c("c1", "c2"))

## ------------------------------------------------------------------------
## Run microbenchmark for aut5 on b10
testset <- create_testset("bench", "b10")
toolset <- create_toolset(set_names = "auc5")
res <- run_benchmark(testset, toolset)
res

## ------------------------------------------------------------------------
## Evaluate Precision-Recall curves for ROCR and precrec with c1 test set
testset <- create_testset("curve", "c1")
toolset <- create_toolset(c("ROCR", "precrec"))
scores <- run_evalcurve(testset, toolset)
scores

## ------------------------------------------------------------------------
## Print all results
print(scores, data_type = "all")

## ---- fig.width=7, warning=FALSE, fig.show='hold'------------------------
## ggplot2 is necessary to use autoplot
library(ggplot2)

## Plot base points and the result of precrec on c1, c2, and c3 test sets
testset <- create_testset("curve", c("c1", "c2", "c3"))
toolset <- create_toolset("precrec")
scores1 <- run_evalcurve(testset, toolset)
autoplot(scores1)

## Plot the results of PerfMeas and PRROC on c1, c2, and c3 test sets
toolset <- create_toolset(c("PerfMeas", "PRROC"))
scores2 <- run_evalcurve(testset, toolset)
autoplot(scores2, base_plot = FALSE)

## ------------------------------------------------------------------------
## Create a new tool set for 'xyz' 
toolname <- "xyz"
calcfunc <- create_example_func()
toolsetU <- create_usrtool(toolname, calcfunc)

## User-defined tools can be combined with predefined tools
toolsetA <- create_toolset("ROCR")
toolsetU2 <- c(toolsetA, toolsetU)

## ---- fig.width=7, warning=FALSE, fig.show='hold'------------------------
## Curve evaluation
testset3 <- create_testset("curve", "c2")
scores3 <- run_evalcurve(testset3, toolsetU2)
autoplot(scores3, base_plot = FALSE)

## ------------------------------------------------------------------------
## Show an example of the second argument
calcfunc <- create_example_func()
print(calcfunc)

## ------------------------------------------------------------------------
## Create a test dataset 'b5' for benchmarking
testsetB <- create_usrdata("bench", scores = c(0.1, 0.2), labels = c(1, 0),
                           tsname = "b5")

## ------------------------------------------------------------------------
## Run microbenchmark for ROCR and precrec on a predefined test dataset
toolset <- create_toolset(c("ROCR", "precrec"))
res <- run_benchmark(testsetB, toolset)
res

## ------------------------------------------------------------------------
## Create a test dataset 'c5' for benchmarking
testsetC <- create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                           tsname = "c5", base_x = c(0.0, 1.0), 
                           base_y = c(0.0, 0.5))

## ---- fig.width=7, warning=FALSE, fig.show='hold'------------------------
## Run curve evaluation for ROCR and precrec on a predefined test dataset
toolset2 <- create_toolset(c("ROCR", "precrec"))
scores2 <- run_evalcurve(testsetC, toolset2)
autoplot(scores2, base_plot = FALSE)

