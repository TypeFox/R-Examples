
set.seed(123)
     
n.seq <- 250
p <- 5
K <- 2
mix.prop <- c(0.3, 0.7)

TP <- array(rep(NA, p * p * K), c(p, p, K))   
TP[,,1] <- matrix(c(0.20, 0.10, 0.15, 0.15, 0.40,
                    0.10, 0.10, 0.20, 0.20, 0.40,
                    0.15, 0.10, 0.20, 0.20, 0.35,
                    0.15, 0.10, 0.20, 0.20, 0.35,
                    0.30, 0.30, 0.10, 0.10, 0.20), byrow = TRUE, ncol = p)
     
TP[,,2] <- matrix(c(0.15, 0.35, 0.20, 0.20, 0.10,
                    0.40, 0.10, 0.20, 0.20, 0.10,
                    0.25, 0.20, 0.15, 0.15, 0.25,
                    0.25, 0.20, 0.15, 0.15, 0.25,
                    0.10, 0.20, 0.20, 0.20, 0.30), byrow = TRUE, ncol = p)
     
# DATA SIMULATION
     
A <- click.sim(n = n.seq, int = c(10, 50), alpha = mix.prop, gamma = TP)

C <- click.read(A$S)

N2 <- click.EM(X = C$X, K = 2); N2$BIC
M2 <- click.EM(X = C$X, y = C$y, K = 2); M2$BIC
F2 <- click.forward(X = C$X, K = 2); F2$BIC
B2 <- click.backward(X = C$X, K = 2); B2$BIC
