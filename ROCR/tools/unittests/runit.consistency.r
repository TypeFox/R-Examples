##
##
##


library(RUnit)
library(ROCR)

# source("tests/runit.simple.r")    # needed for .get.performance.measures

## predict performance measures on random data and check their consistency

testConsistency <- function() {

    for (i in 1:100) {
        n.folds <- sample(1:10,1)
        fold.sizes <- sample(10:100, n.folds, replace=T)
        error.rates <- runif( n.folds )
        pp <- .mock.prediction( fold.sizes, error.rates )
        pred <- prediction( pp$predictions, pp$labels )
        .check.prediction.object(pred)
        a <- .get.performance.measures( pred, c('acc','err','fpr','tpr','fnr','tnr','prec','pcfall','npv','pcmiss','rpp','rnp'))
        .check.consistency( a)
    }
}

testCombining <- function() {
    measures <- c('tpr','fpr','acc','err','rec','sens','fnr','tnr','spec',
                  'ppv','prec','npv','fall','miss','pcfall','pcmiss','rpp','rnp',
                  'phi','mat','mi','chisq','odds','lift')
#                  'auc','prbe','rch','mxe','rmse','phi','mat','mi','chisq',
#                  'odds','lift','f','sar','ecost','cost')

    for (measure1 in measures) {
        print(measure1)
        for (measure2 in measures) {
            n.folds <- sample(1:2,1)
            fold.sizes <- sample(10:20, n.folds, replace=T)
            error.rates <- runif( n.folds )
            pp <- .mock.prediction( fold.sizes, error.rates )
            pred <- prediction( pp$predictions, pp$labels )
            .check.prediction.object(pred)
            perf1 <- performance( pred, measure1 )
            perf2 <- performance( pred, measure2 )
            perf3 <- performance( pred, measure2, measure1 )
            .check.performance.object(perf1)
            .check.performance.object(perf2)
            .check.performance.object(perf3)

            for (i in 1:n.folds) {
                #check elements
                checkEquals(setequal( c( perf1@x.values[[i]], perf2@x.values[[i]]), perf3@alpha.values[[i]] ),T)
                checkEquals(setequal( perf1@y.values[[i]], perf3@x.values[[i]] ),T)
                checkEquals(setequal( perf2@y.values[[i]], perf3@y.values[[i]] ),T)

                #check order
                ind <- sapply( perf1@x.values[[i]], function(x) { min(which(x==perf3@alpha.values[[i]]))})
                checkEquals( perf1@y.values[[i]], perf3@x.values[[i]][ind] )
                checkEquals( perf2@y.values[[i]], perf3@y.values[[i]][ind] )
                
            }
        }
    }
}

.get.performance.measures <- function(pred, measures) {

    ans <- list()
    for (measure in measures) {
        ## need to enclose y.values into a list to avoid flattening
        perf <- performance(pred, measure)
        .check.performance.object( perf )
        ans <- c(ans, list(perf@y.values))
    }
    names(ans) <- measures
    ans
}

.check.consistency <- function(measures) {
    ## check entries of contingency table for consistency

    for (measure in c("acc", "err", "fnr", "tpr", "fpr", "tnr", "pcfall", "prec", "npv", "pcmiss",'rpp','rnp')) {
        if (!measure %in% names(measures)) {
            stop(paste("Performance measure", measure, "not in argument list."))
        }
    }
    for (i in 1:length(measures$acc)) {
        finite.bool <- is.finite(measures$acc[[i]]) & is.finite(measures$err[[i]])
        checkEquals(measures$acc[[i]][finite.bool] + measures$err[[i]][finite.bool],
                    rep(1,length(measures$acc[[i]]))[finite.bool])
        
        finite.bool <- is.finite(measures$fnr[[i]]) & is.finite(measures$tpr[[i]])
        checkEquals(measures$fnr[[i]][finite.bool] + measures$tpr[[i]][finite.bool],
                    rep(1,length(measures$fnr[[i]]))[finite.bool])

        finite.bool <- is.finite(measures$fpr[[i]]) & is.finite(measures$tnr[[i]])
        checkEquals(measures$fpr[[i]][finite.bool] + measures$tnr[[i]][finite.bool],
                    rep(1,length(measures$fpr[[i]]))[finite.bool])

        finite.bool <- is.finite(measures$prec[[i]]) & is.finite(measures$pcfall[[i]])
        checkEquals(measures$prec[[i]][finite.bool] + measures$pcfall[[i]][finite.bool],
                    rep(1,length(measures$acc[[i]]))[finite.bool])

        finite.bool <- is.finite(measures$npv[[i]]) & is.finite(measures$pcmiss[[i]])
        checkEquals(measures$npv[[i]][finite.bool] + measures$pcmiss[[i]][finite.bool],
                    rep(1,length(measures$acc[[i]]))[finite.bool])

        checkEquals(measures$rpp[[i]] + measures$rnp[[i]], rep(1, length(measures$rpp[[i]])))
    }
    
}

  
## use consistency checks to validate results on pathological input cases


performance.measures <- c('tpr','fpr','acc','err','rec','sens','fnr','tnr','spec',
                          'ppv','prec','npv','fall','miss','pcfall','pcmiss','rpp','rnp',
                          'auc','prbe','rch','mxe','rmse','phi','mat','mi','chisq',
                          'odds','lift','f','sar','ecost','cost')

testPathological <- function() {

    # mxe needs 0,1 labels (warning otherwise),
    # rmse needs numeric labels (warning otherwise), sar as well
    pred <- prediction( c(0.1, 0.2, 0.5), c("a", "a", "b"))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'mxe' &
                                                  performance.measures != 'rmse' &
                                                  performance.measures != 'sar']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(0.1, 0.2, 0.5), c(F, F, T))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'mxe' &
                                                  performance.measures != 'rmse' &
                                                  performance.measures != 'sar']

    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(0.1, 0.2, 0.5), c("1", "1", "0"))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(0.1, 0.2, 0.5), c(T, F, F))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'mxe' &
                                                  performance.measures != 'rmse' &
                                                  performance.measures != 'sar' ]
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    # prbe cannot be computed, because only one prec/rec pair available.
    pred <- prediction( c(0,0,0), c(0,1,1))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'prbe' ]
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(0,0,0), ordered(c(0,0,0), levels=c(0,1)))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'auc' &
                                                  performance.measures != 'prbe' &
                                                  performance.measures != 'rch' &
                                                  performance.measures != 'sar' &
                                                  performance.measures != 'ecost']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(-1,-0.2,-0.6), ordered(c(1,0,1), levels=c(0,1)))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'mxe' ]
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(-1,-0.2,-0.6), c(-1,1,-1))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'mxe']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(-1,-0.2,-0.6), c(3,2,3))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'mxe']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction( c(1), ordered(c("a"),levels=c('a','b')))
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[ performance.measures != 'auc' &
                                                  performance.measures != 'prbe' &
                                                  performance.measures != 'rch' &
                                                  performance.measures != 'mxe' &
                                                  performance.measures != 'rmse' &
                                                  performance.measures != 'sar' &
                                                  performance.measures != 'ecost']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    
    .check.consistency( measures)

}



############################################################
# test length of performance measures

.check.performance.object <- function(perf) {
    ylen <- length(perf@y.values)
    xlen <- length(perf@x.values)
    alphalen <- length(perf@alpha.values)

    checkEquals( (xlen==0 || xlen==ylen) && (alphalen==0 || (alphalen==xlen && alphalen==ylen)), T )

    if (xlen==ylen) {
        for (i in 1:ylen) checkEquals( length(perf@x.values[[i]]), length(perf@y.values[[i]])  )
    }
    if (alphalen==ylen) {
        for (i in 1:ylen) checkEquals( length(perf@alpha.values[[i]]), length(perf@y.values[[i]])  )
    }
}

.check.prediction.object <- function( pred) {
    # 1. all entries in prediction object must have equals number of cross-validation runs
    lenvec <- c(length(pred@predictions), length(pred@labels), length(pred@cutoffs), length(pred@fp),
                length(pred@tp), length(pred@fn), length(pred@tn), length(pred@n.pos),
                length(pred@n.neg), length(pred@n.pos.pred), length(pred@n.neg.pred))
    checkEquals( length(unique(lenvec)), 1)

    # 2. inside: xval runs:
    for (i in 1:length(pred@predictions)) {
        checkEquals( length(pred@predictions[[i]]), length(pred@labels[[i]]))
        lenvec <- c(length(pred@cutoffs[[i]]), length(pred@fp[[i]]),
                    length(pred@tp[[i]]), length(pred@fn[[i]]),
                    length(pred@tn[[i]]), length(pred@n.pos.pred[[i]]), length(pred@n.neg.pred[[i]]))
        checkEquals( length(unique(lenvec)), 1)
        checkEquals( unique(lenvec), length(unique(pred@predictions[[i]]))+1 )
    }
    # 3. cutoffs sorted in descending order?
    for (i in 1:length(pred@predictions)) {
        checkEquals( sort(pred@cutoffs[[i]], decreasing=TRUE ), pred@cutoffs[[i]] ) 
    }

    # 4. check 2x2 table for consistency with marginal sums
    for (i in 1:length(pred@predictions)) {
        checkEquals( pred@tp[[i]] + pred@fp[[i]], pred@n.pos.pred[[i]] )
        checkEquals( pred@fn[[i]] + pred@tn[[i]], pred@n.neg.pred[[i]] )
        checkEquals( pred@tp[[i]] + pred@fn[[i]], rep( pred@n.pos[[i]], length(pred@tp[[i]]))  )
        checkEquals( pred@fp[[i]] + pred@tn[[i]], rep( pred@n.neg[[i]], length(pred@tp[[i]]))  )
        checkEquals(pred@n.pos.pred[[i]] + pred@n.neg.pred[[i]],
                    rep( pred@n.pos[[i]] + pred@n.neg[[i]], length(pred@n.pos.pred[[i]])) )
        checkEquals(pred@n.pos[[i]] + pred@n.neg[[i]], length(pred@labels[[i]]))
    }
        
}
 




############################################################
# test measures for consistency on supplied data sets

testDatabase <- function() {
    data(ROCR.simple)
    pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    data(ROCR.xval)
    pred <- prediction(ROCR.xval$predictions, ROCR.xval$labels)
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    data(ROCR.hiv)
    pred <- prediction(ROCR.hiv$hiv.nn$predictions, ROCR.hiv$hiv.nn$labels)
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[performance.measures != 'mxe' &
                                                 performance.measures != 'cal']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)

    pred <- prediction(ROCR.hiv$hiv.svm$predictions, ROCR.hiv$hiv.svm$labels)
    .check.prediction.object(pred)
    measures.to.evaluate <- performance.measures[performance.measures != 'mxe' &
                                                 performance.measures != 'cal']
    measures <- .get.performance.measures(pred, measures.to.evaluate)
    .check.consistency( measures)
}


## remove XXX below to include testDatabaseCombine (currently disabled for speed
## reasons) in the test suite.
tXXXestDatabasesCombine <- function() {
    measures <- c('tpr','fpr','acc','err','rec','sens','fnr','tnr','spec',
                  'ppv','prec','npv','fall','miss','pcfall','pcmiss','rpp','rnp',
                  'phi','mat','mi','chisq','odds','lift')
                  #'auc','prbe','rch','mxe','rmse','phi','mat','mi','chisq',
                  #'odds','lift','f','sar','ecost','cost')
    print("Database combine test deactivated.")
    data(ROCR.simple)
    data(ROCR.xval)
    data(ROCR.hiv)
    all.pred <- list(prediction(ROCR.simple$predictions, ROCR.simple$labels),
                     prediction(ROCR.xval$predictions, ROCR.xval$labels),
                     prediction(ROCR.hiv$hiv.nn$predictions, ROCR.hiv$hiv.nn$labels),
                     prediction(ROCR.hiv$hiv.svm$predictions, ROCR.hiv$hiv.svm$labels))
    lapply(all.pred, .check.prediction.object)
    for (pred in all.pred) {
        for (measure1 in measures) {
            print(measure1)
            for (measure2 in measures) {
                perf1 <- performance( pred, measure1 )
                perf2 <- performance( pred, measure2 )
                perf3 <- performance( pred, measure2, measure1 )
                .check.performance.object(perf1)
                .check.performance.object(perf2)
                .check.performance.object(perf3)
                
                for (i in 1:length(pred@labels)) {
                    #check elements
                    checkEquals(setequal( c( perf1@x.values[[i]], perf2@x.values[[i]]), perf3@alpha.values[[i]] ),T)
                    checkEquals(setequal( perf1@y.values[[i]], perf3@x.values[[i]] ),T)
                    checkEquals(setequal( perf2@y.values[[i]], perf3@y.values[[i]] ),T)

                    # check order
                    ind <- sapply( perf1@x.values[[i]], function(x) { min(which(x==perf3@alpha.values[[i]]))})
                    checkEquals( perf1@y.values[[i]], perf3@x.values[[i]][ind] )
                    checkEquals( perf2@y.values[[i]], perf3@y.values[[i]][ind] )
                    
                }
            }
        }
    }
}
############################################################

crashCases <- list(       ## cases that are ok to crash:
                          list(pred= c(0), lab= c(0)), #-> Number of classes is not equal to 2.
                          list(pred= c(1), lab= c(1)), #-> Number of classes is not equal to 2.
                          list(pred= c(0.1, 0.2, 0.5), lab= c(1,1,1)), #-> Number of classes is not equal to 2.
                          list(pred= c(0.1, 0.2, 0.5), lab= c(0,0,0)), #-> Number of classes is not equal to 2.
                          list(pred= c(0.1, 0.2, 0.5), lab= c("a", "a", "a")),  #-> Number of classes is not equal to 2.
                          list(pred= c(0.1, 0.2, 0.5), lab= c(T, T, T)), #-> Number of classes is not equal to 2.
                          list(pred= c(0.1, 0.2, 0.5), lab= c(F, F, F))  #-> Number of classes is not equal to 2.
)


#                          list(pred= c(), lab= c()),   #-> Number of classes is not equal to 2.

testCrash <- function() {
    for (case in crashCases) {
        cat(case$pred, " ", case$lab, "\n")
        checkException(pred <- prediction(case$pred, case$lab))
        #checkException(measures <- .get.performance.measures(pred))
    }
}


#  
.mock.prediction <- function( n.predictions, error.rate ) {
    if ( length(n.predictions) > 1 && length(error.rate)==1) {
        error.rate <- rep(error.rate, length(n.predictions) )
    }

    if (length(n.predictions)>1) {
        predictions <- list()
        labels <- list()
    } else {
        predictions <- c()
        labels <- c()
    }

    for (i in 1:length(n.predictions)) {
        current.predictions <- runif( n.predictions[i] )
        current.labels <- as.numeric( current.predictions >= 0.5)
        flip.indices <- sample( n.predictions[i], round( error.rate[i] * n.predictions[i] ))
        current.labels[ flip.indices ] <- !current.labels[ flip.indices ]
        # current.labels[ current.labels=="1" ] <- "+"
        # current.labels[ current.labels=="0" ] <- "-"

        if (length(n.predictions)>1) {
            predictions <- c( predictions, list( current.predictions ))
            labels <- c( labels, list( current.labels ))
        }
    }

    if (length( n.predictions)==1) {
        predictions <- list(current.predictions)
        labels <- list(current.labels)
    }

    ans <- list(predictions= predictions, labels= labels)

    # ensure, that random labels have exactly two levels
    if (any(   sapply(labels, function(run) {length(unique(run))})  !=  rep(2, length(labels)) )) {
        print(paste("XXX", labels, str(n.predictions), str(error.rate)))
        return(.mock.prediction(n.predictions, error.rate))
    }
    else return( ans )
}
