mymclust <- function (formula = .~., diagonal = TRUE) 
{    
    require("mvtnorm")
    retval <- new("FLXMC", weighted = TRUE,
                  formula = formula, dist = "mvnorm",
                  name = "my model-based clustering")
    retval@defineComponent <- expression({
        logLik <- function(x, y) {
            dmvnorm(y, mean = center, sigma = cov, log = TRUE)
          }
        predict <- function(x) {
            matrix(center, nrow = nrow(x),
                   ncol = length(center), byrow = TRUE)
        }
        new("FLXcomponent",
            parameters = list(center = center, cov = cov),
            df = df, logLik = logLik, predict = predict)
    })
    retval@fit <- function(x, y, w, ...) {
        
        para <- cov.wt(y, wt = w)[c("center", "cov")]
        df <- (3 * ncol(y) + ncol(y)^2)/2
        
        if (diagonal) {
            para$cov <- diag(diag(para$cov))
            df <- 2 * ncol(y)
          }
        with(para, eval(retval@defineComponent))
    }
    retval
}
