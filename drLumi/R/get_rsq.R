# Estimate of Rsquared
# summary(model)$sigma = sqrt(sum(residuals(model)^2)/df)
get_rsq <- function(model, adjusted=TRUE){
    lhs <- as.character(model$m$formula()[[2]])
    parameters <- names(model$m$getPars())
    allobj <- ls(model$m$getEnv())
    rhs <- allobj[-match(c(parameters,lhs),allobj)]
    ndf <- data.frame(get(lhs, model$m$getEnv()), get(rhs, model$m$getEnv()))
    names(ndf) <- c(lhs, rhs)

    tss.fit <- var(ndf[,lhs])

    if(adjusted==TRUE){
        rss.df <- summary(model)$df[2]-1
        rss.fit <- sum(residuals(model)^2)/rss.df    
    } else {
        rss.fit <- sum(residuals(model)^2)/(length(residuals(model))-1)
    }

    rsquare.fit <- 1 - (rss.fit/tss.fit)
    return(rsquare.fit)  
}









