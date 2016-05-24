f <- function(x) { 
        coef(ut.lm3)[4] * x + coef(ut.lm3)[3] - coef(ut.lm2)[3] 
}
uniroot( f, c(20,50) )$root
