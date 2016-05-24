perfScatterPlot <-
function(metric, modelx, modely, datax, datay=datax, addCostsx=NULL, addCostsy=NULL, pargs=NULL, ...) {
    if(is.null(metric)) {
        stop("Need evaluation metric for plotting!")
    }
    if(is.null(modelx) || is.null(modely)) {
        stop("Need models to plot performances!")
    }
    if(is.null(datax)) {
        stop("Need data to plot performances!")
    }

    idCols = intersect(c(datax$ids, datax$extra), c(datay$ids, datay$extra))
    if(length(idCols) == 0) {
        stop("Cannot match up the two data frames!")
    }

    if(length(datax$test) > 0) {
        edatax = rbind.fill(lapply(datax$test, function(x) {
            datax$data[x,idCols,drop=F]
        }))
    } else {
        edatax = datax$data[idCols]
    }
    if(!is.null(datax$extra)) {
        edatax[datax$extra] = datax$data[datax$extra]
    }
    if(length(datay$test) > 0) {
        edatay = rbind.fill(lapply(datay$test, function(x) {
            datay$data[x,idCols,drop=F]
        }))
    } else {
        edatay = datay$data[idCols]
    }
    if(!is.null(datay$extra)) {
        edatay[datay$extra] = datay$data[datay$extra]
    }

    scoresx = cbind(edatax, data.frame(scorex=metric(datax, modelx, addCosts=addCostsx, ...)))
    scoresy = cbind(edatay, data.frame(scorey=metric(datay, modely, addCosts=addCostsy, ...)))

    d = merge(scoresx, scoresy, by=idCols)
    ggplot(d, aes_string(x="scorex", y="scorey")) +
        geom_abline(slope = 1) +
        geom_point(pargs)
}
