library(gnm)
set.seed(1)
DNase1 <- subset(DNase, Run == 1)

test <- gnm(density ~ -1 +
            Mult(1, Inv(Const(1) + Exp(Mult(1 + offset(-log(conc)),
                                                Inv(1))))),
            start = c(NA, 0, 1), data = DNase1, trace = TRUE)
coef(test)

Logistic <- function(x, inst = NULL){
    list(predictors = list(Asym = 1, xmid = 1, scal = 1),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste(predLabels[1], "/(1 + exp((", predLabels[2], "-",
                   varLabels[1], ")/", predLabels[3], "))")
         },
         start = function(theta){
             theta[3] <- 1
             theta
         }
         )
}
class(Logistic) <- "nonlin"

test <- gnm(density ~ -1 + Logistic(log(conc)),
            data = DNase1, trace = TRUE)
coef(test)
