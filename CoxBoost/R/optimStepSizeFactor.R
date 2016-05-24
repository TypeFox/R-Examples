optimStepSizeFactor <- function(time,status,x,direction=c("down","up","both"),start.stepsize=0.1,iter.max=10,
                                constant.cv.res=NULL,parallel=FALSE,trace=FALSE,...) 
{
    if (parallel) {
        if (!require(snowfall)) {
            parallel <- FALSE
            warning("package 'snowfall' not found, i.e., parallelization cannot be performed")
        } else {
            snowfall::sfExport("x")
        }
    }

    direction <- match.arg(direction)

    factor.list <- switch(direction,down=c(1,1-start.stepsize,1-2*start.stepsize),
                                    up=c(1,1+start.stepsize,1+2*start.stepsize),
                                    both=c(1-start.stepsize,1,1+start.stepsize))

    critmat <- NULL

    folds <- NULL
    if (!is.null(constant.cv.res)) folds <- constant.cv.res$folds

    i <- 1
    step.size <- start.stepsize
    reduction.done <- FALSE
    iter.count <- 0

    while (i <= length(factor.list) && iter.count < iter.max) {
        iter.count <- iter.count + 1
        if (trace) cat("iteration ",iter.count,": evaluating factor ",factor.list[i],"\n",sep="")
        
        if (factor.list[i] == 1 && !is.null(constant.cv.res)) {
            cv.res.act.path <- constant.cv.res
        } else {
            if (is.null(folds)) {
                cv.res.act.path <- cv.CoxBoost(time=time,status=status,x=x,stepsize.factor=factor.list[i],parallel=parallel,upload.x=FALSE,trace=trace,...)
            } else {
                cv.res.act.path <- cv.CoxBoost(time=time,status=status,x=x,stepsize.factor=factor.list[i],parallel=parallel,upload.x=FALSE,folds=folds,trace=trace,...)
            }
        }

        critmat <- rbind(critmat,cv.res.act.path$mean.logplik)

        if (is.null(folds)) folds <- cv.res.act.path$folds
    
        i <- i + 1
        if (i > length(factor.list)) {
            if (reduction.done) break

            actual.max.val <- apply(critmat,1,max)
            actual.max <- which.max(actual.max.val)

            do.reduction <- TRUE

            if (direction %in% c("down","both") && factor.list[actual.max] == min(factor.list)) {
                if (min(factor.list) < step.size*2) break
                factor.list <- c(factor.list,min(factor.list)-step.size)
                do.reduction <- FALSE
            }

            if (direction %in% c("up","both") && factor.list[actual.max] == max(factor.list)) {
                factor.list <- c(factor.list,max(factor.list)+step.size)
                do.reduction <- FALSE
            }
            
            if (do.reduction) {
                sort.factor.list <- sort(factor.list)
                max.pos <- which(sort.factor.list == factor.list[actual.max])
                max.factor <- sort.factor.list[max.pos]
                
                if (max.pos == 1) {
                    second.max.pos <- 2
                } else {
                    if (max.pos == length(sort.factor.list)) {
                        second.max.pos <- length(sort.factor.list) - 1
                    } else {
                        candidate1.pos <- max.pos - 1
                        candidate2.pos <- max.pos + 1

                        if (actual.max.val[which(sort.factor.list[candidate1.pos] == factor.list)] > actual.max.val[which(sort.factor.list[candidate2.pos] == factor.list)]) {
                            second.max.pos <- candidate1.pos
                        } else {
                            second.max.pos <- candidate2.pos
                        }
                    }
                }
                
                second.max.factor <- sort.factor.list[second.max.pos]
                
                factor.list <- c(factor.list,mean(c(max.factor,second.max.factor)))
        
                reduction.done <- TRUE
            }
        }
    }

    optimal.factor.index <- which.max(apply(critmat,1,max))
    
    list(factor.list=factor.list,
         critmat=critmat,
         optimal.factor.index=optimal.factor.index,
         optimal.factor=factor.list[optimal.factor.index],
         optimal.step=which.max(critmat[optimal.factor.index,])-1)
}
