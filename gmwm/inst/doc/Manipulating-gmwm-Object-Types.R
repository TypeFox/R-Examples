## ----load_pkg, echo = F, message=F---------------------------------------
library(gmwm)

## ------------------------------------------------------------------------
vals = rnorm(100)

# Basic
ex1.basic = gts(vals)

# Advanced
ex1.adv = gts(vals, freq = 2, unit = "sec", name = "Example 1 Data")

## ------------------------------------------------------------------------
N = 1000
model = WN(sigma2=1) + AR1(phi = 0.1, sigma2 = .7)

# Basic
ex2.basic = gen.gts(model, N)

# Advanced
ex2.adv = gen.gts(model, freq = 10, unit = "sec", name = "Example 2 Data")

## ----eval = F------------------------------------------------------------
#  is.gts(x)

## ----eval = F------------------------------------------------------------
#  length(x)
#  nrow(x)

## ------------------------------------------------------------------------
x = gen.gts(WN(sigma2=1), 50)
plot(x)

## ------------------------------------------------------------------------
plot(x, title = "Example of `gts` plot", axis.x.label = "Time", axis.y.label = "Value")

## ----eval = F------------------------------------------------------------
#  N = 10000
#  vals = cbind(rnorm(N),
#               rnorm(N, mean = 3, sd = 1),
#               rnorm(N, mean = 2, sd = 3),
#               rnorm(N, mean = 1, sd = 2))
#  
#  # Basic
#  x = imu(vals,
#          gyros = 1:2,
#          accels = 3:4)
#  
#  # Advanced
#  x = imu(vals,
#          gyros = 1:2,
#          accels = 3:4,
#          axis = c("A","B"),
#          freq = 5, unit = "sec", name = "Example 3 Data")

## ----eval = F------------------------------------------------------------
#  # Option 1: Relative Directory
#  setwd("C:/Users/James")
#  x = read.imu("Documents/imu_data.imu", type = "IXSEA")
#  
#  # Option 2: Fixed Directory
#  x = read.imu("C:/Users/James/Documents/imu_data.imu", type = "IXSEA")

## ----eval = F------------------------------------------------------------
#  is.imu(x)

## ----eval = F------------------------------------------------------------
#  nrow(x)

## ----eval = F------------------------------------------------------------
#  ncol(x)
#  length(x)

## ----eval = F------------------------------------------------------------
#  # Numeric Return
#  value(x,"accel")   # Number of Accelerometers
#  value(x,"gyro")    # Number of Gyroscopes
#  value(x,"sensors") # Total Number of Sensors
#  
#  # Logical Return (T/F)
#  has(x, "accel")
#  has(x, "gyro")
#  has(x, "sensors")

## ----eval = F------------------------------------------------------------
#  # Start of data set
#  head(x, 5)
#  
#  # End of data set
#  tail(x, 5)

## ----eval = F------------------------------------------------------------
#  # Creates an AR1 modeling component with the program set to guess initial values.
#  model.guided = AR1()
#  
#  # Makes an AR1 modeling component with user supplied initial values.
#  model.adv = AR1(phi = .3, sigma2 = 1)

## ----eval = F------------------------------------------------------------
#  # Creates an AR1 + WN model with the program set to guess initial values.
#  model.guided = AR1() + WN()
#  
#  # Builds an AR1 + WN model with user supplied initial values.
#  model.adv = AR1(phi = .9, sigma2 = 1) + WN(sigma2 = .1)
#  
#  # Creates an ARMA(2,2) Process
#  arma.guided = ARMA(2,2)
#  
#  # Specifying parameters for an arma 2,2
#  arma.adv = ARMA(ar = c(0.3,.7), ma = c(0.5,.1), sigma2 = 1)

## ----eval = F------------------------------------------------------------
#  # Creates an 3AR1 model with the program set to guess initial values.
#  model.rep = 3*AR1()
#  
#  # Equivalent to:
#  model.add = AR1() + AR1() + AR1()

