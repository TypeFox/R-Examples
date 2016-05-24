### R code from vignette source 'spls-example.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ")


###################################################
### code chunk number 2: spls-prelim
###################################################
library("spls")


###################################################
### code chunk number 3: spls-data
###################################################
data(yeast)
yeast$x[1:5,1:5]
yeast$y[1:5,1:5]


###################################################
### code chunk number 4: plot1 (eval = FALSE)
###################################################
## set.seed(1)
## cv <- cv.spls( yeast$x, yeast$y, eta = seq(0.1,0.9,0.1), K = c(5:10) )


###################################################
### code chunk number 5: cv
###################################################
set.seed(1)
cv <- cv.spls( yeast$x, yeast$y, eta = seq(0.1,0.9,0.1), K = c(5:10) )


###################################################
### code chunk number 6: spls-fn
###################################################
f <- spls( yeast$x, yeast$y, eta = cv$eta.opt, K = cv$K.opt )
print(f)
coef.f <- coef(f)
coef.f[1:5,1:5]


###################################################
### code chunk number 7: plot2 (eval = FALSE)
###################################################
## plot.spls( f, yvar=1 )


###################################################
### code chunk number 8: plot
###################################################
plot.spls( f, yvar=1 )


###################################################
### code chunk number 9: plot3 (eval = FALSE)
###################################################
## coefplot.spls( f, nwin=c(2,2), xvar=c(1:4) )


###################################################
### code chunk number 10: coef
###################################################
coefplot.spls( f, nwin=c(2,2), xvar=c(1:4) )


###################################################
### code chunk number 11: eqtl-data
###################################################
data(mice)
mice$x[1:5,1:5]
mice$y[1:5,1:5]


###################################################
### code chunk number 12: eqtl-cv (eval = FALSE)
###################################################
## set.seed(1)
## cv <- cv.spls( mice$x, mice$y, eta = seq(0.1,0.9,0.1), K = c(1:5) )


###################################################
### code chunk number 13: eqtl-fn (eval = FALSE)
###################################################
## f <- spls( mice$x, mice$y, eta = cv$eta.opt, K = cv$K.opt )
## print(f)


###################################################
### code chunk number 14: eqtl-fn
###################################################
f <- spls( mice$x, mice$y, eta = 0.6, K = 1 )
print(f)


###################################################
### code chunk number 15: plot5 (eval = FALSE)
###################################################
## set.seed(1)
## ci.f <- ci.spls( f, plot.it=TRUE, plot.fix='x', plot.var=20 )


###################################################
### code chunk number 16: ciplot
###################################################
set.seed(1)
ci.f <- ci.spls( f, plot.it=TRUE, plot.fix='x', plot.var=20 )


###################################################
### code chunk number 17: eqtl-ci
###################################################
cis <- ci.f$cibeta
cis[[20]][1:5,]


###################################################
### code chunk number 18: cor1 (eval = FALSE)
###################################################
## cf <- correct.spls( ci.f )


###################################################
### code chunk number 19: cor2
###################################################
cf <- correct.spls( ci.f, plot.it=FALSE )


###################################################
### code chunk number 20: cor-out
###################################################
cf[15:20,1:5]


###################################################
### code chunk number 21: plot6 (eval = FALSE)
###################################################
## heatmap.spls( mat=f$betahat, xlab='Predictors', ylab='Responses',
##             main='Original Coefficient Estimates', coln=16, as='i' )


###################################################
### code chunk number 22: plot7 (eval = FALSE)
###################################################
## heatmap.spls( mat=cf, xlab='Predictors', ylab='Responses',
##             main='Corrected Coefficient Estimates', coln=16, as='i' )        


###################################################
### code chunk number 23: coef1
###################################################
heatmap.spls( mat=f$betahat, xlab='Predictors', ylab='Responses',
            main='Original Coefficient Estimates', coln=16, as='i' )


###################################################
### code chunk number 24: coef2
###################################################
heatmap.spls( mat=cf, xlab='Predictors', ylab='Responses',
            main='Corrected Coefficient Estimates', coln=16, as='i' )        


