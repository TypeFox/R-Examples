context("Count data regression charts")

## Simulation parameters
S <- 1 ; t <- 1:120 ; m <- length(t)
beta <- c(1.5,0.6,0.6)
omega <- 2*pi/52
#log mu_{0,t}
alpha <- 0.2
base <- beta[1] + beta[2] * cos(omega*t) + beta[3] * sin(omega*t) 
#Generate example data with changepoint and tau=tau
tau <- 100
kappa <- 0.4
mu0 <- exp(base)
mu1 <- exp(base  + kappa) 

## Generate counts
set.seed(42)
x <- rnbinom(length(t),mu=mu0*(exp(kappa)^(t>=tau)),size=1/alpha)
s.ts <- create.disProg(week=t, observed=x, state=(t>=tau))

## Define control object
cntrl1 <- list(range=t,c.ARL=5, mu0=mu0, alpha=alpha,
               change="intercept", ret="value", dir="inc")

## Run algorithm
glr.ts1 <- algo.glrnb(s.ts, control=cntrl1)

## Correct upperbound (rounded)
## dput(signif(c(glr.ts1$upperbound), 7))
correctUpperbound <- c(
    0.0933664, 0, 0.001387989, 0.4392282, 1.239898, 2.983766, 1.954988, 
    1.722341, 1.586777, 0.7331938, 0.9337575, 0.7903225, 1.104522, 
    1.425098, 1.24129, 1.633672, 2.033343, 1.788079, 1.397671, 0.9081794, 
    0.797097, 0.7270934, 0.5248943, 0.3093548, 0.2622768, 0.2301054, 
    0.1595651, 0.1484989, 0.06889605, 0.1504776, 0.04138495, 0.02219845, 
    0.0231524, 0.009575689, 0.1504776, 0.5827537, 0.0357062, 0.005011513, 
    0, 1.390972, 0.3167743, 0.5717088, 0.1053871, 0.003442552, 0.0005934715, 
    0, 0, 0.05509335, 0.1375619, 0.2449853, 0.6840703, 0.5427538, 
    0.05675776, 0.06656547, 0.09036596, 0.209314, 0.1392091, 0.03494786, 
    0.026216, 0.277202, 0.01762547, 0, 0, 0, 3.564077, 1.41019, 0.290548, 
    0.3740241, 0.4269062, 0.1296794, 0.1298662, 0.6322042, 0.2115204, 
    0.107457, 0.9366399, 0.1379007, 0.1509654, 0.03392803, 0.005775552, 
    0, 0, 0, 0, 0, 0.001143512, 0.001637927, 1.021689, 1.965804, 
    1.83044, 1.017412, 0.3033473, 0.1689957, 0.4051742, 0.1247774, 
    0.1460143, 0.03590031, 0.9459381, 0.4189531, 0.2637725, 0.03925406, 
    0.01374443, 0.2283519, 2.535301, 1.406133, 1.692899, 2.021258, 
    2.951635, 4.25683, 4.77543, 3.90064, 3.646361, 3.680106, 4.236502, 
    5.522696, 0.1221651, 0.4054735, 0.6761779, 0.8039129, 0.3913383, 
    0.1261521)

test_that("upperbound equals pre-computed value",
          expect_that(c(glr.ts1$upperbound),
                      equals(correctUpperbound, tolerance=1e-6)))
