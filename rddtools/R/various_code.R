### MISC
is_even <- function(a) {
    a%%2 == 0
}


Kernel_tri <- function(X, center, bw) {
    ifelse(abs(X - center) > bw, 0, 1 - (abs(X - center)/bw))
}

Kernel_uni <- function(X, center, bw) {
    ifelse(abs(X - center) > bw, 0, 1)
}
