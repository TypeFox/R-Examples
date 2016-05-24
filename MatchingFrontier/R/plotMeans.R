plotMeans <-
function(frontier.object,
         xlab = 'Number of Observations Pruned',
         main = 'Means Plot',
         xlim = c(1,max(frontier.object$frontier$Xs)),
         ylim = c(0, 1),
         cols = rainbow(length(frontier.object$match.on)),
         diff.in.means = FALSE,
         ...){

    # NOTE: This borrows the basic structure of
    # estimateEffects. 
    
    covs.mat <- matrix(nrow = 0, ncol = length(frontier.object$match.on), byrow = FALSE)
    colnames(covs.mat) <- frontier.object$match.on

    for(col in frontier.object$match.on){
        frontier.object$dataset[,colnames(frontier.object$dataset) == col] <- range01(
                             frontier.object$dataset[,colnames(frontier.object$dataset) == col])
    }
    
    # Calculate means

    # For variable ratio
    for(i in 1:length(frontier.object$frontier$Xs)){        
        this.dat.inds <- unlist(frontier.object$frontier$drop.order[i:length(frontier.object$frontier$drop.order)])
        dataset <- frontier.object$dataset[this.dat.inds,]       
    
        if(frontier.object$ratio == 'variable'){
            w <- makeWeights(dataset, frontier.object$treatment)
            dataset$w <- w
            
            new.row <- c()
            for(col in colnames(covs.mat)){
                if(diff.in.means){
                    treated.inds <- dataset[[frontier.object$treatment]] == 1
                    treated.mean <- weighted.mean(dataset[treated.inds, col], w = dataset$w[treated.inds])

                    control.inds <- dataset[[frontier.object$treatment]] == 0
                    control.mean <- weighted.mean(dataset[control.inds, col], w = dataset$w[control.inds])

                    new.row <- c(new.row, abs(treated.mean - control.mean))
                }else{
                    new.row <- c(new.row, weighted.mean(dataset[,col], w = dataset$w))
                }
            }
            covs.mat <- rbind(covs.mat, new.row)
        }else{ # For fixed ratio
            new.row <- c()
            for(col in colnames(covs.mat)){
                if(diff.in.means){
                    treated.inds <- dataset[[frontier.object$treatment]] == 1
                    treated.mean <- mean(dataset[treated.inds, col])
                    
                    control.inds <- dataset[[frontier.object$treatment]] == 0
                    control.mean <- mean(dataset[control.inds, col])

                    new.row <- c(new.row, abs(treated.mean - control.mean))
                }else{
                    new.row <- c(new.row, mean(dataset[,col]))
                }

            }
            covs.mat <- rbind(covs.mat, new.row)            
        }
    }

    if(!(exists('ylab'))){
        if(diff.in.means){
            ylab <- 'Difference in scaled means'
        }else{
            ylab <- 'Scaled means'
        }
    }
    
    # make plot
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    x <- frontier.object$frontier$Xs
    plot(1, type='n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
    for(i in 1:ncol(covs.mat)){
        this.y <- covs.mat[, colnames(covs.mat)[i]]
        points(x, this.y, type = 'l', col = cols[i])
    }

    legend("topright", inset=c(-0.325,0), lty=c(1,1), legend=colnames(covs.mat), col = cols, title="Group")
}    

