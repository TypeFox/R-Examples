## Ackley
tfAckley <- function(x) {
    n <- length(x)
    20 + exp(1) - 20 * exp(-0.2*sqrt(x %*% x /n)) - exp(sum(cos(2*pi*x))/n)
}

## tfGriewank
tfGriewank <- function(x) {
    n <- length(x)
    x %*% x / 4000 - prod(cos(x / sqrt(1L:n))) + 1
}

##tfRastrigin
tfRastrigin <- function(x) {
    n <- length(x)
    10*n + sum(x^2 - 10 * cos(2*pi*x))
}

## Rosenbrock
tfRosenbrock <- function(x) {
    n <- length(x)
    xi <- x[seq_len(n-1L)]
    sum(100 * (x[2L:n] - xi * xi)^2 + (1 - xi)^2)
}

## tfSchwefel
tfSchwefel <- function(x)
    sum(-x * sin(sqrt(abs(x))))

## Trefethen
tfTrefethen <- function(x) {
    y <- x[2L]
    x <- x[1L]
    exp(sin(50*x)) + sin(60*exp(y)) + sin(70*sin(x)) +
        sin(sin(80*y)) - sin(10*(x+y))  + (x^2+y^2)/4
}
