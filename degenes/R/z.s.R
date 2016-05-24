z.s <-
function(a, b, m1, m2) {
cat("Calculation of z", fill=TRUE)
z <- ((apply(a, 1, sum)/m1) - (apply(b, 1, sum)/m2))  / sqrt ((apply(a, 1, var)/m1) + (apply(b, 1, var)/m1))
return(z)
}

