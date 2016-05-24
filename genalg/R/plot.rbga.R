plot.rbga <- function(x, type="default", breaks=10, ...) {
    # should do object type checking here
    rbga.object = x

    if ((type == "hist") & (rbga.object$type == "floats chromosome")) {
        vars = length(rbga.object$stringMin)
        par(mfcol=c(vars,1));
        for (var in 1:vars) {
            hist(rbga.object$population[,var], main=paste("Variable", var), 
                 xlab="variable value", breaks=breaks, ...);
        }
        par(mfcol=c(1,1));
    } else if ((type == "hist") & (rbga.object$type == "binary chromosome")) {
        plot(x=c(1:ncol(rbga.object$population)), y=apply(rbga.object$population, 2, sum), 
             main="Sum of Selected Genes", xlab="gene", ylab="times selected", 
             type="h", ...);
    } else if ((type == "vars") & (rbga.object$type == "floats chromosome")) {
        vars = length(rbga.object$stringMin)
        par(mfcol=c(vars,1));
        for (var in 1:vars) {
            plot(x=rbga.object$population[,var], y=rbga.object$evaluations,
                 main=paste("Variable", var), 
                 xlab="variable value", ylab="evaluation value", ...);
        }
        par(mfcol=c(1,1));
    } else {
        if (type != "default") {
            warning(paste("Plot type", type, "not supported for a RBGA object of type",
                    rbga.object$type));
        }
        max = max(rbga.object$best, rbga.object$mean);
        min = min(rbga.object$best, rbga.object$mean);
        plot(rbga.object$best, type="l", main="Best and mean evaluation value",
             ylim=c(min, max), xlab="generation", ylab="evaluation value", ...);
        lines(rbga.object$mean, col="blue", ...);
    }
}

