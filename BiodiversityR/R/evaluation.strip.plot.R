`evaluation.strip.plot` <- function(
    data, TrainData=NULL,
    modelnames=c("MAXENT", "GBM", "GBMSTEP", "RF", "GLM", "GLMSTEP", "GAM", "GAMSTEP", "MGCV", "MGCVFIX",
        "EARTH", "RPART", "NNET", "FDA", "SVM", "SVME", "BIOCLIM", "DOMAIN", "MAHAL"),
    variable=NULL, model=NULL, 
    dev.new.width=7, dev.new.height=7, ...
) 
{
    if (is.null(TrainData) == F) {
        TrainData <- TrainData[TrainData[, "pb"]==1, ]
        TrainData[, "pb"] <- as.factor(TrainData[, "pb"])
    }
    if(is.null(variable)==F) {
        v <- (which(names(data) == variable))
        v <- v-2
        f <- data[,1]==v
        vars <- max(data[,1])
# plot for all models
        if (is.null(model) ==T) {
            modelnames <- c(modelnames, "ENSEMBLE")
            nmodels <- length(modelnames)
# models with data
            models <- 0
            for (j in 1:nmodels) {
                if (any(is.na(data[f, 2+vars+j])==F)) {models <- models + 1}
            }
            dim1 <- max(1, ceiling(sqrt(models)))
            dim2 <- max(1, ceiling(models/dim1))
            par.old <- graphics::par(no.readonly=T)
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            graphics::par(mfrow=c(dim1,dim2))
            for (j in 1:models) {
                if (any(is.na(data[v, 2+vars+j])==F)) {
                    if (is.null(TrainData)==T  || is.factor(TrainData[, which(names(TrainData) == variable)])==T) {
                        graphics::plot(data[f,v+2], data[f, 2+vars+j], main=variable, xlab="", ylab=names(data)[2+vars+j],...)
                    }else{
                        graphics::plot(data[f,v+2], data[f, 2+vars+j], main=variable, xlab="", ylab=names(data)[2+vars+j], ylim=c(0, 1.25), ...)
                        graphics::boxplot(TrainData[, which(names(TrainData) == variable)] ~ TrainData[,"pb"], add=T, horizontal=T)
                    }
                }
            }
            graphics::par(par.old)
        }else{
            m <- names(data) == model
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            if (any(is.na(data[v, m])==F)) {
                    if (is.null(TrainData)==T  || is.factor(TrainData[, which(names(TrainData) == variable)])==T) {
                        graphics::plot(data[f,v+2], data[f, m], main=variable, xlab="", ylab=names(data)[2+vars+j],...)
                    }else{
                        graphics::plot(data[f,v+2], data[f, m], main=variable, xlab="", ylab=names(data)[2+vars+j], ylim=c(0, 1.25), ...)
                        graphics::boxplot(TrainData[, which(names(TrainData) == variable)] ~ TrainData[,"pb"], add=T, horizontal=T)
                    }
            }
        }
    }
    if(is.null(model)==F && is.null(variable)==T) {
        m <- names(data) == model
# models with data
        if(is.na(sum(data[, m]))) { 
            cat(paste("NOTE: No data for model: ",  model, "\n", sep = ""))
        }else{
            vars <- max(data[,1])
            dim1 <- max(1, ceiling(sqrt(vars)))
            dim2 <- max(1, ceiling(vars/dim1))
            if (dev.new.width > 0 && dev.new.height > 0) {grDevices::dev.new(width=dev.new.width, height=dev.new.height)}
            par.old <- graphics::par(no.readonly=T)
            graphics::par(mfrow=c(dim1,dim2))
            for (i in 1:vars) {
                f <- data[,1]==i
                if (is.null(TrainData)==T  || is.factor(TrainData[, which(names(TrainData) == names(data)[i+2])])==T) {
                    graphics::plot(data[f,i+2], data[f, m], main=names(data)[i+2], xlab="", ylab=model,...)
                }else{
                    graphics::plot(data[f,i+2], data[f, m], main=names(data)[i+2], xlab="", ylab=model, ylim=c(0, 1.25), ...)
                    varfocal <- names(data)[i+2]
                    graphics::boxplot(TrainData[, which(names(TrainData) == varfocal)] ~ TrainData[,"pb"], add=T, horizontal=T)
                }
            }
            graphics::par(par.old)
        }
    }
}

