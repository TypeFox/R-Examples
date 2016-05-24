n <- nrow(football)
x <- sum(football$Helium > football$Air)
binom.test(x,n)
