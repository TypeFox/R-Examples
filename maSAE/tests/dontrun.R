#!/usr/bin/Rscript --vanilla
warning('tests in tests/dootrun.R a skipped for the sake of cpu time on CRAN.')
if (FALSE){
## from system.file("test-tools-1.R",      package = "Matrix"):
assertError <- function(expr) {
  d.expr <- deparse(substitute(expr))
  t.res <- tryCatch(expr, error = function(e) e)
  print(t.res)
  if(!inherits(t.res, "error"))
    stop(d.expr, "\n\t did not give an error", call. = FALSE)
  invisible(t.res)
}
assertWarning <- function(expr) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, warning = function(w)w)
    if(!inherits(t.res, "warning"))
        stop(d.expr, "\n\t did not give a warning", call. = FALSE)
    invisible(t.res)
}
options(warn = 2)
library('maSAE')
library('methods')
## ## ## ## ## ## real data
data('s0')
phase0 <- s0; phase0$x1 <- NULL
data('s1')
data('s2')

s2$s1 <- TRUE; s2$s2 <- TRUE
s1$clustid <- s1$y <- NA
s1$s1 <- TRUE; s1$s2 <- FALSE
s12 <- rbind(s1, s2)
phase0$clustid <- phase0$y <-phase0$x1 <- NA
phase0$s1 <- phase0$s2 <- FALSE
s012 <- rbind(phase0, s12)
s0$clustid <- s0$y <- NA
s0$s1 <- s0$s2 <- FALSE
s0f12 <- rbind(s0, s12)

preds <- paste('x',1:3, sep='')

trueMeans.s0 <- as.data.frame(
    rbind(
        colMeans(subset(s0f12, g =='a')[, preds])
        , colMeans(subset(s0f12, g =='b')[, preds])
        )
    ); trueMeans.s0$g=c('a', 'b')
partialMeans.s0 <- trueMeans.s0[,-2]
trueMeans.s1 <- as.data.frame(
    rbind(
        colMeans(subset(s12, g =='a')[, preds])
        , colMeans(subset(s12, g =='b')[, preds])
        )
    ); trueMeans.s1$g=c('a', 'b')
partialMeans.s1 <- trueMeans.s1[,-2]

str(predict(saObj(data = s2, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1)))
str(predict(saObj(data = s2, f = y ~ NULL | g)))
## ## ## ## ## missing data
## ## ## ## missing predictand in s2
missing <- s2
missing[missing$clustid <= 4000, ]$y <- NA
## ## ## design
assertError(
predict(saObj(data = missing, f = y ~ NULL| g))
  )
assertError(
predict(saObj(data = missing, f = y ~ NULL| g, cluster = 'clustid'))
  )
## ## ## exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1))
  )
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1, s2 = 's2', cluster = 'clustid'))
  )

## ## ## partially exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', cluster = 'clustid'))
  )

## ## ## non-exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, s2 = 's2', cluster = 'clustid'))
  )
## ## ## three-phase
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', s1 = 's1'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', s1 = 's1', cluster = 'clustid'))
  )


## ## ## ## missing predictor in s2
missing <- s2
missing[missing$clustid <= 4000, ]$x2 <- NA
## ## ## design, should work
predict(saObj(data = missing, f = y ~ NULL| g))
predict(saObj(data = missing, f = y ~ NULL| g, cluster = 'clustid', include = 'inclusion'))
## ## ## exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1))
  )
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1, s2 = 's2', cluster = 'clustid'))
  )

## ## ## partially exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', cluster = 'clustid'))
  )

## ## ## non-exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, s2 = 's2', cluster = 'clustid'))
  )
## ## ## three-phase
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', s1 = 's1'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', s1 = 's1', cluster = 'clustid'))
  )


## ## ## ## missing predictor in s1
missing <- s1
missing[1:floor(nrow(missing)/10), ]$x2 <- NA
missing <- rbind(missing, s2)

## ## ## partially exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', cluster = 'clustid'))
  )

## ## ## non-exhaustive
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, s2 = 's2'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, s2 = 's2', cluster = 'clustid'))
  )
## ## ## three-phase
## ## noncluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', s1 = 's1'))
  )
## ## cluster
assertError(
  predict(saObj(data = missing, f = y ~x1 + x2 + x3 | g, smallAreaMeans = partialMeans.s1, s2 = 's2', s1 = 's1', cluster = 'clustid'))
  )

## ## ## ## missing predictor in s0
missing <- phase0
missing[1:floor(nrow(missing)/10), ]$x2 <- NA
missing <- rbind(missing, s12)
## ## ## three-phase, should work
## ## noncluster
a <- predict(saObj(data =   missing, f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1'))
b <- predict(saObj(data = s0f12, f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1'))
assertError( identical(a,b)|| stop('NA should give different result'))
## ## cluster
a <- predict(saObj(data =   missing, f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1', cluster = 'clustid'))
b <- predict(saObj(data = s0f12, f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1', cluster = 'clustid'))
assertError( identical(a,b)|| stop('NA should give different result'))


## ## ## ## ## plausibility
## ## unclustered
## EXHAUSTIVE
(ex <- predict(saObj(data = s2, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1)))
ex.0 <- predict(saObj(data = s0f12, f = y ~x1 + x2 + x3 | g, s2 = 's2', smallAreaMeans = trueMeans.s1))
all.equal(ex[,3], ex.0[,3]) ||  stop('s2-option should guarantee identical results')
## PARTIALLY EXHAUSTIVE
(part <- predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g,  s2 = 's2', smallAreaMeans = partialMeans.s1)))
## THREE-PHASE
(tp <- predict(saObj(data = s012,  f = y ~x1 + x2 + x3 | g, s1 = 's1', s2 = 's2')))
## NON-EXHAUSTIVE
(nonex <- predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g, s2 = 's2')))

if (! all.equal(as.numeric(ex[, 2]), as.numeric(part[, 2]))) stop('smallAreaMeans were set to sample means in s1, all but THREE-PHASE predictions should be identical!')
if (! all(ex[,3] < part[,3]) )  warning('at least one EXHAUSTIVE variance not smaller than PARTIALLY EXHAUSTIVE variance')
if (! all(ex[,3] < nonex[,3]) )  stop('at least one EXHAUSTIVE variance not smaller than NON EXHAUSTIVE variance')
if (! all(ex[,3] < tp[,3]) )  stop('at least one EXHAUSTIVE variance not smaller than THREE-PHASE variance')
if (! all(part[,3] < nonex[,3]) ) warning('at least one PARTIALLY EXHAUSTIVE variance not smaller than NON-EXHAUSTIVE variance')
if (! all(tp[,3] < nonex[,3]) ) stop('at least one THREE-PHASE variance not smaller than NON-EXHAUSTIVE variance')
if (! all(part[,3] < tp[,3])) warning('at least one PARTIALLY EXHAUSTIVE variance not smaller than THREE-PHASE variance. This seems o.k., since they use different data')

## THREE-PHASE and NON-EXHAUSTIVE
(a <- predict(saObj(data = s2,  f = y ~x1 + x2 + x3 | g, s2 = 's1')))
(b <- predict(saObj(data = s2,  f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1')))
if (! all.equal(as.numeric(unlist(a)),as.numeric(unlist(b)) )  )  stop('NON-EXHAUSTIVE and THREE-PHASE should be the same given identical data')
## THREE-PHASE and EXHAUSTIVE
(a <- predict(saObj(data = s2, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s0)))
(b <- predict(saObj(data = s0f12,  f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1')))
if (! all.equal(as.numeric(a[,2]), as.numeric(b[,2])))  stop('EXHAUSTIVE and THREE-PHASE predictions should be the same given identical data')
if (! all(a[,3] < b[,3]) )  stop('at least one EXHAUSTIVE variance not smaller THREE-PHASE than variance')
## THREE-PHASE and PARTIALLY EXHAUSTIVE
(a <- predict(saObj(data = s0f12, f = y ~x1 + x2 + x3 |  g,  s2 = 's2', smallAreaMeans = partialMeans.s0)))
(b <- predict(saObj(data = s0f12,  f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1')))
if (! all.equal(as.numeric(a[,2]), as.numeric(b[,2])))  stop('PARTIALLY EXHAUSTIVE and THREE-PHASE predictions should be the same given identical data')
if (! all(a[,3] < b[,3]) )  stop('at least one PARTIALLY EXHAUSTIVE variance not smaller THREE-PHASE than variance')


## ## clustered
## EXHAUSTIVE
(ex.cl <- predict(saObj(data = s2, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s1, cluster = 'clustid', s2 = 's2')))
ex.0.cl <- predict(saObj(data = s0f12, f = y ~x1 + x2 + x3 | g, s2 = 's2', smallAreaMeans = trueMeans.s1, cluster = 'clustid'))
all.equal(ex.cl[,3], ex.0.cl[,3]) ||  stop('s2-option should guarantee identical results')
## PARTIALLY EXHAUSTIVE
(part.cl <- predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g,  s2 = 's2', smallAreaMeans = partialMeans.s1, cluster = 'clustid')))
## THREE-PHASE
(tp.cl <- predict(saObj(data = s012,  f = y ~x1 + x2 + x3 | g, s1 = 's1', s2 = 's2', cluster = 'clustid')))
## NON-EXHAUSTIVE
(nonex.cl <- predict(saObj(data = s12, f = y ~x1 + x2 + x3 | g, s2 = 's2', cluster = 'clustid')))

if (! all.equal(as.numeric(ex.cl[, 2]), as.numeric(part.cl[, 2])) ) stop('smallAreaMeans were set to sample means in s1, all but THREE-PHASE predictions should be identical!')
if (! all(ex.cl[,3] < part.cl[,3]) )  warning('at least one EXHAUSTIVE variance not smaller than PARTIALLY EXHAUSTIVE variance')
if (! all(ex.cl[,3] < nonex.cl[,3]) )  stop('at least one EXHAUSTIVE variance not smaller than NON EXHAUSTIVE variance')
if (! all(ex.cl[,3] < tp.cl[,3]) )  stop('at least one EXHAUSTIVE variance not smaller than THREE-PHASE variance')
if (! all(part.cl[,3] < nonex.cl[,3]) ) warning('at least one PARTIALLY EXHAUSTIVE variance not smaller than NON-EXHAUSTIVE variance')
if (! all(tp.cl[,3] < nonex.cl[,3]) ) stop('at least one THREE-PHASE variance not smaller than NON-EXHAUSTIVE variance')
if (! all(part.cl[,3] < tp.cl[,3])) warning('at least one PARTIALLY EXHAUSTIVE variance not smaller than THREE-PHASE variance. This seems o.k., since they use different data')

## THREE-PHASE and NON-EXHAUSTIVE
(a <- predict(saObj(data = s2,  f = y ~x1 + x2 + x3 | g, s2 = 's1', cluster = 'clustid')))
(b <- predict(saObj(data = s2,  f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1', cluster = 'clustid')))
all.equal(as.numeric(unlist(a)),as.numeric(unlist(b)) )  ||  stop('NON-EXHAUSTIVE and THREE-PHASE should be the same given identical data')
## THREE-PHASE and EXHAUSTIVE
(a <- predict(saObj(data = s2, f = y ~x1 + x2 + x3 | g, smallAreaMeans = trueMeans.s0, cluster = 'clustid', s2 = 's2')))
(b <- predict(saObj(data = s0f12,  f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1', cluster = 'clustid')))
if (! all.equal(as.numeric(a[,2]), as.numeric(b[,2])))  stop('EXHAUSTIVE and THREE-PHASE predictions should be the same given identical data')
if (! all(a[,3] < b[,3]) )  stop('at least one EXHAUSTIVE variance not smaller THREE-PHASE than variance')
## THREE-PHASE and PARTIALLY EXHAUSTIVE
(a <- predict(saObj(data = s0f12, f = y ~x1 + x2 + x3 |  g,  s2 = 's2', smallAreaMeans = partialMeans.s0, cluster = 'clustid')))
(b <- predict(saObj(data = s0f12,  f = y ~x1 + x2 + x3 | g, s2 = 's2', s1 = 's1', cluster = 'clustid')))
if (! all.equal(as.numeric(a[,2]), as.numeric(b[,2])))  stop('PARTIALLY EXHAUSTIVE and THREE-PHASE predictions should be the same given identical data')
if (! all(a[,3] < b[,3]) )  stop('at least one PARTIALLY EXHAUSTIVE variance not smaller THREE-PHASE than variance')

## ## unclustered vs clustered
if (! all(ex[,3] < ex.cl[,3]) ) stop('at least one UNCLUSTERED variance not smaller than CLUSTERED variance')
if (! all(part[,3] < part.cl[,3]) ) stop('at least one UNCLUSTERED variance not smaller than CLUSTERED variance')
if (! all(tp[,3] < tp.cl[,3]) ) stop('at least one UNCLUSTERED variance not smaller than CLUSTERED variance')
if (! all(nonex[,3] < nonex.cl[,3]) ) stop('at least one UNCLUSTERED variance not smaller than CLUSTERED variance')

if (! all( abs(1 - (ex[,2] / ex.cl[,2] )) < 0.01) ) warning('at least one area with UNCLUSTERED and CLUSTERED predictions differing by more than .01 ')
if (! all( abs(1 - (part[,2] / part.cl[,2] )) < 0.01) ) warning('at least one area with UNCLUSTERED and CLUSTERED predictions differing by more than .01 ')
if (! all( abs(1 - (tp[,2] / tp.cl[,2] )) < 0.01) ) warning('at least one area with UNCLUSTERED and CLUSTERED predictions differing by more than .01 ')
if (! all( abs(1 - (nonex[,2] / nonex.cl[,2] )) < 0.01) ) warning('at least one area with UNCLUSTERED and CLUSTERED predictions differing by more than .01 ')
}
