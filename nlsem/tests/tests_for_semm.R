# Testing functions for semm
# ===========================

library(nlsem)

# SEMM model for structural equation models
# ==========================================

# create model
mod <- specify_sem(num.x=4, num.y=4, num.xi=2, num.eta=2,
                     xi="x1-x2,x3-x4", eta="y1-y2,y3-y4",
                     constraints="direct1", num.classes=2,
                     interaction="none",
                     rel.lat="eta1~xi1,eta2~xi2,eta2~eta1,eta1~eta2")

dat <- as.data.frame(mod)

dat[c(2,8,10,16),2:3]   <- 1    # Lambda.x.21 & 42 and Lambda.y.21 & 42
dat[c(58,62),2:3]       <- 0    # psi.21 and phi.21
dat[65:72,2:3]          <- 1    # nu.x and nu.y

model <- create_sem(dat)

# simulate data
parameters <- c(
                # class 1
                rep(0.5, 2),    # Gamma
                c(-0.3, 0.7),   # Beta
                rep(0.5, 10),   # Theta.d, Theta.e, Psi
                rep(1, 2),      # Phi
                rep(1, 2),      # alpha
                rep(1, 2),      # tau
                # class 2
                rep(-0.5, 2),   # Gamma
                c(0.7, -0.3),   # Beta
                rep(0.5, 10),   # Theta.d, Theta.e, Psi
                rep(1, 2),      # Phi
                rep(1, 2),      # alpha
                rep(4, 2)       # tau
)
data <- simulate(model, seed=7, parameters=parameters)

# estimate model

set.seed(8)
parameters <- runif(count_free_parameters(model), 0.1, 1.5)

# system.time(
#     res <- em(model, data, parameters, verbose=TRUE)
# )

