## ----knitropts,echo=FALSE,message=FALSE----------------------------------
if (require('knitr')) opts_chunk$set(fig.width = 5, fig.height = 5, fig.align = 'center', tidy = FALSE, warning = FALSE, cache = TRUE)

## ----prelim,echo=FALSE---------------------------------------------------
enpls.version = '1.1'

## ----load-package--------------------------------------------------------
require(enpls)
data(alkanes)
x = alkanes$x
y = alkanes$y

## ----enpls.fs,fig.cap='Top ten important variables of the \\texttt{alkanes} dataset'----
set.seed(42)
varimp = enpls.fs(x, y, MCtimes = 100)
print(varimp, nvar = 10L)
plot(varimp, nvar = 10L)

## ----enpls.od,fig.cap='Outlier detection result of the \\texttt{alkanes} dataset'----
od = enpls.od(x, y, MCtimes = 100)
plot(od, criterion = 'sd')

## ----enpls.en------------------------------------------------------------
enpls.fit = enpls.en(x, y, MCtimes = 100)

## ----predict.enpls.en,fig.cap='Experimental values vs. predicted values'----
y.pred = predict(enpls.fit, newx = x)
plot(y, y.pred, xlim = range(y), ylim = range(y),
     xlab = 'Experimental', ylab = 'Predicted')
abline(a = 0L, b = 1L)

## ----cv.enpls,fig.cap='Cross validation result: experimental values vs. predicted values'----
cv.enpls.fit = cv.enpls(x, y, MCtimes = 20)
print(cv.enpls.fit)
plot(cv.enpls.fit)

