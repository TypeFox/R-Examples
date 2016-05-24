# tria.R -- version 2010-12-16
# variant 1
Min <- 0; Max <- 3; Mode <- 0.75
rtria1 <- function(u1, u2, Min, Max, Mode){
    Mode + (Min + u1 * (Max - Min) - Mode) * sqrt(u2)
}
trials <- 1000L
T <- rtria1(runif(trials), runif(trials), Min,Max,Mode)
hist(T)


# variant 2
Min <- 0; Max <- 3; Mode <- 0.75
rtria2 <- function(u1,u2,u3,Min,Max,Mode) {
    Mode + (Min + u1 * (Max - Min) - Mode) * pmax(u2, u3)
}
trials <- 1000L
T <- rtria2(runif(trials), runif(trials), runif(trials), 
            Min, Max, Mode)
hist(T, breaks = 100L)


# inverse
triaInverse <- function(u,min,max,mode) {
    range <- max-min
    r <- ifelse(u <= (mode-min)/(range),
            min + sqrt(   u  * range * (mode - min)),
            max - sqrt((1-u) * range * (max  - mode)))
    return(r)
}
trials <- 1000L
u <- runif(trials)
Min <- 0; Max <- 3; Mode <- 0.75
system.time(T <- triaInverse(u,Min,Max,Mode))
hist(T)


# correlation
N <- 1000L; p <- 4L
rho <- 0.7; rho <- 2 * sin(rho * pi / 6) # spearman
X <- array(rnorm(N * p), dim = c(N, p))
C <- matrix(rho, nrow = p, ncol = p); 
diag(C) <- 1; C <- chol(C)

X <- X %*% C
cor(X, method = "spearman")

U <- pnorm(X)
T <- triaInverse(U,Min,Max,Mode)

# graphic (see ?pairs)
panel.hist <- function(x) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1L:2L], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col=grey(.5))
}
par(las = 1, mar = c(2,2,0.5,0.5), ps = 10, tck = 0.01, 
    mgp = c(3, 0.2, 0), pch = '.')
pairs(T, diag.panel = panel.hist, gap = 0)
