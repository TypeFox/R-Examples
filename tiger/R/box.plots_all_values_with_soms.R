box.plots.false <- function(result, solution, show.measures =
1:num.measures, new.order=1:solution, show.synthetic.peaks=FALSE, 
synthetic.peaks.col = c(2:8,3), show.timestep=NA, show.cell=NA,  ref=NULL,
ref.new.order=new.order, ref.solutions=solution
){

    num.measures <- length(result$names)
    n.errors <- dim(result$measures.synthetic.peaks)[1]
    palette("default")

    cluster.rer <- change.order.clusters(result$cluster.assignment[[solution]], new.order)
    do.ref <- FALSE
    if(!is.null(ref)){
        do.ref<-TRUE
        cluster.rer.ref <- change.order.clusters(ref$cluster.assignment[[ref.solutions]], ref.new.order)
    }
    
    centers.box <- result$best.value.location$all.values.reranged
    centers.rer <- result$best.value.location$central.best.reranged

    the.best <- rep("", solution)
    the.worst <- rep("", solution)
    the.low <- rep("", solution)
    the.high <- rep("", solution)
    bp.medians <- matrix(nrow=num.measures, ncol=solution)

    for(i in show.measures){
        main.expr <- substitute(expression(' '*b),list(b=result$names[[i]]))
        if(do.ref){
            if(result$use.som){
                data <- ref$som$code[,i] 
            } else {
                data <- ref$measures.uniform[!ref$na.rows,i]
            }
            bp <- boxplot(data ~ cluster.rer.ref$cluster, col="lightgray", xlab="region", main=eval(main.expr), names=LETTERS[1:ref.solutions], ylim=c(0,1))
        }
        if(result$use.som){
            clusters <- merge(result$som$visual,
                              cbind(result$som$code.sum,
                                    cluster=cluster.rer$cluster))$cluster
        } else {
            clusters <- cluster.rer$cluster
        }
        data <- result$measures.uniform[!result$na.rows,i]
        bp <- boxplot(data ~ clusters , col=rainbow(solution), xlab="region", main=eval(main.expr), names=LETTERS[1:solution], ylim=c(0,1), add=do.ref , boxwex=0.6 )
        bp.medians[i,]<-bp$stats[3,]
        lines(c(0,solution+1), c(centers.box[i],centers.box[i]), lwd=2, col="yellow")
        if(show.synthetic.peaks){
            plot.window(c(0,7),c(0,1))
            data <- result$measures.uniform.synthetic.peaks[,,i]
            for(error in 1:n.errors){
                points(data[error,], pch=error, col=synthetic.peaks.col[error])
            }
        }
        if(!any(is.na(show.cell))){
             l.dat <- result$som$code[1+show.cell$x+show.cell$y * result$som.dim[1],i]
             lines(c(0,solution+1), c(l.dat,l.dat), lwd=2, col="blue")
        }
        if(!is.na(show.timestep)){
             l.dat <- result$measures.uniform[show.timestep,i]
             lines(c(0,solution+1), c(l.dat,l.dat), lwd=2, col="blue")
        }
        dist.from.center <- centers.box[i] - bp$stats[3,]
        if(is.na(centers.box[i])) stop("Please define full centers.box")

        best <- which.min(abs(dist.from.center))  #center - median
        best <- include.others(best, centers.box[i], bp$stats, best=TRUE)
        the.best[best] <- paste(the.best[best] , result$names[i], sep=", ")

        #do we have errors to both sides?
        two.sides <- FALSE
        if(!is.na(centers.rer[i]) & 
           length(unique(sign(dist.from.center[-best])))==2){
               two.sides <-TRUE
        }

        if(two.sides){
            low <- which.min(bp$stats[3,])
            low <- include.others(low, centers.box[i], bp$stats)
            low <- low[! low %in% best]
            the.low[low] <- paste(the.low[low] , result$names[i], sep=", ")
            high <- which.max(bp$stats[3,])
            high <- include.others(high, centers.box[i], bp$stats)
            high <- high[! high %in% best]
            the.high[high] <- paste(the.high[high] , result$names[i], sep=", ")

        } else {
            worst <- which.max(abs(dist.from.center))
            worst <- include.others(worst, centers.box[i], bp$stats)
            worst <- worst[! worst %in% best]
            the.worst[worst] <- paste(the.worst[worst] , result$names[i], sep=", ")
        }
    }


  return(paste(LETTERS[1:solution], " & ", 
    "{\bf best:} ", substring(the.best, 3), 
    "; {\bf worst:} ", substring(the.worst, 3) ,
    "; {\bf low:} ", substring(the.low , 3), 
    "; {\bf high:} ", substring(the.high, 3), "\\", sep=""))
 }
