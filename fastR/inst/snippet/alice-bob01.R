f <- function(prob) { 
        prob  + prob * (1-prob)^3/(1-(1-prob)^2) 
}
g <- function(prob) { f(prob) - 0.50 }
uniroot(g,c(0.20,0.5))$root          # when g=0, f=0.5
