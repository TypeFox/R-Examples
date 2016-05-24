eta2 <- function(x,y) {
return(summary(lm(as.formula(x~y)))$r.squared)
}
