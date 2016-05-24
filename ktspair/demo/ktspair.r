# --------------------------------------------------------------------------
## Demo file for ktspair; start with 'demo(ktspair)'
# --------------------------------------------------------------------------


data(ktspdata)
dim(dat)

ktsp1 <- ktspcalc(dat, grp, 3)
ktsp1

## Can also use Expression Sets

readline(paste("Hit return to continue.\n"))

ktsp2 <- ktspcalc(eSet, grp, 3)
ktsp2
ktsp3 <- ktspcalc(eSet, 1,3)
ktsp3

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------

## Control of the length of the list in the C code.

ktsp <- ktspcalc(dat, grp, 5, length=8)
ktsp

## Here length=8 is too short, need to use a bigger number. 

readline(paste("Hit return to continue.\n"))

ktsp <- ktspcalc(dat, grp, 5, length=100)
ktsp

## length=100 seems to be a big enough number. Nevertheless using high
## value of length will reduce the computation speed of the programm.
## Using length=20 would be a good compromise.

readline(paste("Hit return to continue.\n"))

ktsp <- ktspcalc(dat, grp, 5, length=20)
ktsp

## The parameter length is available in almost all function of this package.
## It has to be adapted by the user on the analysed dataset.

readline(paste("Hit return to continue.\n"))

## An alternative method:  Compute for each gene the average of median within each group and subtract this value from gene's expression value.

ktsp_median <- ktspcalc(dat, grp, 5, length=20, med=TRUE)
ktsp_median

## All the functions in the package are able to deal this substraction, and can also compute it.
## This is set via the variable med.

# --------------------------------------------------------------------------

## Obtain a summary of the model based on the dataset.
## Either for the all list of TSPs using the majority voting algorithm

summary(ktsp)

## Or for every single TSP in the list
readline(paste("Hit return to continue.\n"))

summary(ktsp, printall=TRUE)

## Or for a selected pair.
readline(paste("Hit return to continue.\n"))

summary(ktsp, select=3)

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------

## Display a graphical summary of the results
## Either plot every TSP

dev.new()
plot(ktsp)

readline(paste("Hit return to continue.\n"))

## Or select one

plot(ktsp, select=1)

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------

## Predict invidual based on the estimated model
## By default, predict the individual used to construct the model

predict(ktsp)


## It is also possible to predict classification for new data
## using the variable dat.
readline(paste("Hit return to continue.\n"))

predict(ktsp, dat=dat[,20:30])

## Finally, it is also possible to select one pair of genes and to control 
## the warnings displayed in the shell (useful for multiple uses of the method,
## for example in bootstrap, see later)
## with the parameter select and display respectively.

readline(paste("Hit return to continue.\n"))

predict(ktsp, select=1)
predict(ktsp, dat=dat[,20:30], display=FALSE)

## The parameter display can be chosen in almost all function of this package.
## It can be useful to check if the parameter length was correctly set. 

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------

## The parameter k can be chosen by crossvalidation to minimize
## the error rate prediction
## By default the number of fold is set to 5, but can be chosen
## by the user through the variable cross.

cv <- cv(dat, grp, seed=1)
ktsp <- ktspcalc(dat, grp, cv$k)
ktsp

## This function returns which k should be used as well as the accuracy
## reached by the k-TSP with different values of k (1,3,5,7,9).

readline(paste("Hit return to continue.\n"))

## It is also possible to call the function cv directly in the
## function ktspcalc() if the value of k is not specified.

ktsp <- ktspcalc(dat, grp, seed=1, display=FALSE)
ktsp

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------

## In order to study the robustness of the k-TSPon the original dataset,
## a bootstrap procedure may be useful.

bootstrap1 <- bootstrap.ktsp(dat, grp, k=3, n=5, seed=1)
bootstrap2 <- bootstrap.ktsp(dat, grp, n=5, seed=1)
 
## Here the number of iteration were chosen to be quite small. This is to speed up the computation.
## It is obvious that, for more accurate results, the number of iterations has to be bigger.


bootstrap.graphic.ktsp(bootstrap1, title="Without adjustment of the graphic parameters")

readline(paste("Hit return to continue.\n"))


bootstrap.graphic.ktsp(bootstrap1, para1=0.3, para2=0.3, title="With adjustment of the graphic parameters")

readline(paste("Hit return to continue.\n"))


bootstrap.graphic.ktsp(bootstrap2, para1=0.3, para2=0.3, title="With k chosen by crossvalidation at every step", 
mtext="With adjustment of the graphic parameters")

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------


## Based on the results of the bootstrap, one may want to choose some pairs of genes to construct a k-TSP object
## This is possible with the function ktspcalc2().

ktsp2 <- ktspcalc2(c(74, 704, 298, 587, 34, 266), dat, grp, performance = TRUE, healthy="healthy")
ktsp2

ktsp2$accuracy
ktsp2$sensitivity
ktsp2$specificity

readline(paste("Hit return to continue.\n"))

# --------------------------------------------------------------------------

## In order to study the accuracy of the classifier k-TSP, a ROC curve may be useful.
## It is possible to transform the majority voting algorithm by using different cutoff.
## The variable mult.cutoff allows to use the cutoff 0.25, 0.5 and 0.75. 

roc1 <- ROC.offset(dat, grp, n=5, healthy="healthy", mult.cutoff=FALSE)
roc2 <- ROC.offset(dat, grp, n=5, healthy="healthy", mult.cutoff=TRUE, para1=50, para2=2)

dev.off()
ROC.graphic.ktsp(roc1, boxplot = TRUE, AUC = TRUE, maintitle="ROC curve")

readline(paste("Hit return to continue.\n"))

dev.off()
ROC.graphic.ktsp(roc1, boxplot = FALSE, AUC = TRUE, maintitle="ROC curve", mtext="Without boxplot")

readline(paste("Hit return to continue.\n"))

dev.off()
ROC.graphic.ktsp(roc2, maintitle = "ROC curve for the k-TSP with several cutoffs")

readline(paste("Hit return to continue.\n"))

dev.off()
roc3 <- ROC.voting(dat, grp, n=5, healthy="healthy")
ROC.graphic.ktsp(roc3, maintitle = "ROC curve with the voting method")
# --------------------------------------------------------------------------

################################END#########################################
