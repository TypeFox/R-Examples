"rounding" <-
function(y, m){
    a <- y / m - 0.5 / (m + 1)
    a <- ifelse(a < 0, 0, a)
    a <- ifelse(y == m, 1 - 0.5 / (m + 1), a)
    b <- y / m + 0.5 / (m + 1)
    b <- ifelse(b > 1, 1, b)
    b <- ifelse(y == 0, 0.5 / (m + 1), b)
    cbind(a, b)
}

