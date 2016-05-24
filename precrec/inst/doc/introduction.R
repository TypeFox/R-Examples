## ------------------------------------------------------------------------
library(precrec)

# Load a test dataset
data(P10N10)

# Calculate ROC and Precision-Recall curves
sscurves <- evalmod(scores = P10N10$scores, labels = P10N10$labels)

## ---- fig.width=7, fig.show='hold'---------------------------------------
# Show ROC and Precision-Recall plots
plot(sscurves)

# Show a Precision-Recall plot
plot(sscurves, "PRC")

## ---- fig.width=7, fig.show='hold'---------------------------------------
# The ggplot2 package is required 
library(ggplot2)

# Show ROC and Precision-Recall plots
autoplot(sscurves)

# Show a Precision-Recall plot
autoplot(sscurves, "PRC")

## ------------------------------------------------------------------------
# Get a data frame with AUC scores
aucs <- auc(sscurves)

# Use knitr::kable to display the result in a table format
knitr::kable(aucs)

# Get AUCs of Precision-Recall
aucs_prc <- subset(aucs, curvetypes == "PRC")
knitr::kable(aucs_prc)


## ------------------------------------------------------------------------
s1 <- c(1, 2, 3, 4)
s2 <- c(5, 6, 7, 8)
s3 <- matrix(1:8, 4, 2)

# Join two score vectors
scores1 <- join_scores(s1, s2)

# Join two vectors and a matrix
scores2 <- join_scores(s1, s2, s3)

## ------------------------------------------------------------------------
l1 <- c(1, 0, 1, 1)
l2 <- c(1, 0, 1, 1)
l3 <- c(1, 0, 1, 0)

# Join two label vectors
labels1 <- join_labels(l1, l2)
labels2 <- join_labels(l1, l3)

## ------------------------------------------------------------------------
# Create an input dataset with two score vectors and one label vector
msmdat <- mmdata(scores1, labels1)

# Specify dataset IDs
smmdat <- mmdata(scores1, labels2, dsids = c(1, 2))

# Specify model names and dataset IDs
mmmdat <- mmdata(scores1, labels2, modnames = c("mod1", "mod2"), dsids = c(1, 2))

## ------------------------------------------------------------------------
# A dataset with 10 positives and 10 negatives for the random performance level
samps1 <- create_sim_samples(1, 10, 10, "random")

#  A dataset for five different performance levels
samps2 <- create_sim_samples(1, 10, 10, "all")

# A dataset with 20 samples for the good early retrieval performance level
samps3 <- create_sim_samples(20, 10, 10, "good_er")

# A dataset with 20 samples for five different performance levels
samps4 <- create_sim_samples(20, 10, 10, "all")

## ------------------------------------------------------------------------
# Use a list with multiple score vectors and a list with a single label vector
msmdat1 <- mmdata(scores1, labels1)

# Explicitly specify model names
msmdat2 <- mmdata(scores1, labels1, modnames = c("mod1", "mod2"))

# Use a sample dataset created by the create_sim_samples function
msmdat3 <- mmdata(samps2[["scores"]], samps2[["labels"]], modnames = samps2[["modnames"]])

## ------------------------------------------------------------------------
# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(msmdat3)

## ---- fig.width=7, fig.show='hold'---------------------------------------
# Show ROC and Precision-Recall curves with the ggplot2 package
autoplot(mscurves)

## ------------------------------------------------------------------------
# Specify test dataset IDs names
smmdat1 <- mmdata(scores1, labels2, dsids = c(1,2))

# Use a sample dataset created by the create_sim_samples function
smmdat2 <- mmdata(samps3[["scores"]], samps3[["labels"]], dsids = samps3[["dsids"]])

## ------------------------------------------------------------------------
# Calculate curves for multiple test datasets and keep all the curves
smcurves <- evalmod(smmdat2, raw_curves = TRUE)

## ---- fig.width=7, fig.show='hold'---------------------------------------
# Show an average Precision-Recall curve with the 95% confidence bounds
autoplot(smcurves, "PRC")

# Show raw Precision-Recall curves
autoplot(smcurves, "PRC", raw_curves = TRUE)

## ------------------------------------------------------------------------
# Specify model names and test dataset IDs names
mmmdat1 <- mmdata(scores1, labels2, modnames= c("mod1", "mod2"), dsids = c(1, 2))

# Use a sample dataset created by the create_sim_samples function
mmmdat2 <- mmdata(samps4[["scores"]], samps4[["labels"]], 
                  modnames = samps4[["modnames"]], dsids = samps4[["dsids"]])

## ------------------------------------------------------------------------
# Calculate curves for multiple models and multiple test datasets
mmcurves <- evalmod(mmmdat2)

## ---- fig.width=7, fig.show='hold'---------------------------------------
# Show average Precision-Recall curves
autoplot(mmcurves, "PRC")

# Show average Precision-Recall curves with the 95% confidence bounds
autoplot(mmcurves, "PRC", show_cb = TRUE)

## ------------------------------------------------------------------------
# Balanced dataset
samps5 <- create_sim_samples(100, 100, 100, "all")
simmdat1 <- mmdata(samps5[["scores"]], samps5[["labels"]], 
                   modnames = samps5[["modnames"]], dsids = samps5[["dsids"]])

# Imbalanced dataset
samps6 <- create_sim_samples(100, 25, 100, "all")
simmdat2 <- mmdata(samps6[["scores"]], samps6[["labels"]], 
                   modnames = samps6[["modnames"]], dsids = samps6[["dsids"]])


## ------------------------------------------------------------------------
# Balanced dataset
simcurves1 <- evalmod(simmdat1)

# Imbalanced dataset
simcurves2 <- evalmod(simmdat2)

## ---- fig.width=7, fig.show='hold'---------------------------------------
# Balanced dataset
autoplot(simcurves1)

# Imbalanced dataset
autoplot(simcurves2)

## ------------------------------------------------------------------------
# Calculate basic evaluation measures
mmpoins <- evalmod(mmmdat2, mode = "basic")

## ---- fig.width=7, fig.show='hold'---------------------------------------
# Show normalized threshold values vs. error rate and accuracy
autoplot(mmpoins, c("error", "accuracy"))

# Show normalized threshold values vs. specificity and sensitivity
autoplot(mmpoins, c("specificity", "sensitivity"))

# Show normalized threshold values vs. precision
autoplot(mmpoins, "precision")

