
## ----chgt_package, warning=FALSE, message=FALSE--------------------------
library(ESGtoolkit)


## ----example_simshocks_1-------------------------------------------------
# Number of simulations
nb <- 1000

# Number of risk factors
d <- 2

# Number of possible combinations of the risk factors (here : 1)
dd <- d*(d-1)/2

# Family : Gaussian copula 
fam1 <- rep(1,dd)
# Correlation coefficients between the risk factors (d*(d-1)/2)
par0.1 <- 0.1
par0.2 <- -0.9


## ----example_simshocks_simul---------------------------------------------
set.seed(2)
# Simulation of shocks for the d risk factors
s0.par1 <- simshocks(n = nb, horizon = 4, 
family = fam1, par = par0.1)

s0.par2 <- simshocks(n = nb, horizon = 4, 
family = fam1, par = par0.2)


## ----example_simshocks_2-------------------------------------------------
# Correlation test
esgcortest(s0.par1)


## ----example_simshocks_3, fig.align='center', fig.height=4.5-------------
test <- esgcortest(s0.par2)
par(mfrow=c(1, 2))
esgplotbands(esgcortest(s0.par1))
esgplotbands(test)


## ----example_simshocks_4-------------------------------------------------
# Family : Rotated Clayton (180 degrees)
fam2 <- 13
par0.3 <- 2

# Family : Rotated Clayton (90 degrees)
fam3 <- 23
par0.4 <- -2

# number of simulations
nb <- 200

# Simulation of shocks for the d risk factors
s0.par3 <- simshocks(n = nb, horizon = 4, 
family = fam2, par = par0.3)

s0.par4 <- simshocks(n = nb, horizon = 4, 
family = fam3, par = par0.4)


## ----example_simshocks_5, fig.align='center', fig.height=6---------------
esgplotshocks(s0.par3, s0.par4)


## ----chgtfOptions, message=FALSE, warning=FALSE--------------------------
library(fOptions)


## ----<example_SVJD_1-----------------------------------------------------
# Spot variance
V0 <- 0.1372
# mean-reversion speed
kappa <- 9.5110/100
# long-term variance
theta <- 0.0285
# volatility of volatility
volvol <- 0.8010/100
# Correlation between stoch. vol and prices
rho <- -0.5483
# Intensity of the Poisson process
lambda <- 0.3635
# mean and vol of the merton jumps diffusion
mu.J <- -0.2459
sigma.J <- 0.2547/100
m <- exp(mu.J + 0.5*(sigma.J^2)) - 1
# Initial stock price
S0 <- 4468.17
# Initial short rate
r0 <- 0.0357


## ----<example_SVJD_2-----------------------------------------------------
n <- 300
horizon <- 1
freq <- "weekly" 

# Simulation of shocks, with antithetic variates
shocks <- simshocks(n = n, horizon = horizon, 
          frequency = freq, 
          method = "anti", 
          family = 1, par = rho)

# Vol simulation
sim.vol <- simdiff(n = n, horizon = horizon,
                   frequency = freq, model = "CIR", x0 = V0,
                   theta1 = kappa*theta, theta2 = kappa, 
                   theta3 = volvol,
                    eps = shocks[[1]])

# Plotting the volatility (only for a low number of simulations)
esgplotts(sim.vol)


## ----<example_SVJD_price-------------------------------------------------
# prices simulation
sim.price <- simdiff(n = n, horizon = horizon,
                     frequency = freq, model = "GBM", x0 = S0,
                     theta1 = r0 - lambda*m, theta2 = sim.vol,
                     lambda = lambda, mu.z = mu.J, 
                     sigma.z = sigma.J, 
                     eps = shocks[[2]])


## ----<example_SVJD_4-----------------------------------------------------
par(mfrow=c(2,1))
matplot(time(sim.price), sim.price, type = 'l', 
        main = "with matplot")
esgplotbands(sim.price, main = "with esgplotbands", xlab = "time", 
             ylab = "values")


## ----<example_SVJD_3-----------------------------------------------------
# Discounted Monte Carlo price
as.numeric(esgmcprices(r0, sim.price, 2/52))
# Inital price
S0
# pct. difference
as.numeric((esgmcprices(r0, sim.price, 2/52)/S0 - 1)*100)


## ----<example_SVJD_cvS0, fig.height=4------------------------------------
# convergence of the discounted price
esgmccv(r0, sim.price, 2/52, 
        main = "Convergence towards the initial \n asset price")


## ----<example_SVJD_martingale_1, results='hide'--------------------------
martingaletest.sim.price <- esgmartingaletest(r = r0, 
                                              X = sim.price, 
                                              p0 = S0)


## ----<example_SVJD_martingale_2------------------------------------------
esgplotbands(martingaletest.sim.price)


## ----<example_pricing----------------------------------------------------
# Option pricing

# Strike
K <- 3400
Kts <- ts(matrix(K, nrow(sim.price), ncol(sim.price)), 
               start = start(sim.price), 
          deltat = deltat(sim.price),
          end = end(sim.price))

# Implied volatility
sigma.imp <- 0.6625

#Maturity
maturity <- 2/52

# payoff at maturity
payoff <- (sim.price - Kts)*(sim.price > Kts)
payoff <- window(payoff, 
             start = deltat(sim.price), 
             deltat = deltat(sim.price),
             names = paste0("Series ", 1:n))

# True price
c0 <- GBSOption("c", S = S0, X = K, Time = maturity, r = r0, 
                b = 0, sigma = sigma.imp)
c0@price

# Monte Carlo price
as.numeric(esgmcprices(r = r0, X = payoff, maturity))

# pct. difference
as.numeric((esgmcprices(r = r0, X = payoff, 
             maturity = maturity)/c0@price - 1)*100)

# Convergence towards the option price
esgmccv(r = r0, X = payoff, maturity = maturity, 
        main = "Convergence towards the call \n option price")



