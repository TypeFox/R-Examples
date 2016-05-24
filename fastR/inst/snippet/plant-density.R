# plant densities (lambda) for simulations
density <- c(0.01,0.1,0.25,0.5,1,2,4,10,100)       
sizes <- c(10,20,50)
simulate <- function(lambda=1,size=10){
    x <- rpois(size, pi*lambda)
    plants.x <- sum(x)
    area.x <- pi * length(x)
    y <- rweibull( size, shape=2, scale=1/sqrt(pi * lambda) )
    plants.y <- length(y)
    area.y <- pi * sum(y^2)
    mle.x <- plants.x / area.x  
    mle.y <- plants.y / area.y  
    return( c(size, lambda, mle.x,mle.y, 
              plants.x, area.x, plants.y,area.y) )
}

results <- c()
for (lambda in density) {
    for (size in sizes) {
        results <- rbind(results, 
                         t(replicate(1000, simulate(lambda,size))))
    }
}

results <- data.frame(results)
names(results) <- c("size","lambda","mle.count", "mle.measure",
    "plants.count", "area.count", "plants.measure", "area.measure")

results2 <- data.frame(
    size=c(results$size,results$size),
    lambda=c(results$lambda,results$lambda),
    mle=c(results$mle.count,results$mle.measure),
    plants=c(results$plants.count,results$plants.measure),
    area=c(results$area.count,results$area.measure),
    method=rep(c("count","measure"),each=nrow(results))
    )

require(Hmisc)
# have a bake-off:
summary( 
    as.numeric(abs(mle.count-lambda) > abs(mle.measure - lambda)) ~
        size+lambda,
    data=results,
    method="cross",
    fun = function(x) { round(mean(x),3) } 
    )
plot1 <- bwplot( (mle - lambda)/lambda ~ method | factor(lambda),
    subset=size==20,
    data=results2 
    )
# compare standard deviations
summary( 
    mle.count ~ size+lambda,
    data=results,
    method="cross",
    fun = function(x) {round(sd(x),3)}
    )
summary( 
    mle.measure ~ size+lambda,
    data=results,
    method="cross",
    fun = function(x) {round(sd(x),3)}
    )
