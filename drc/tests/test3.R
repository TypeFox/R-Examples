## Test provided by Nicholas Lewin-Koh 2007-03-01

library(drc)

# First dataset
dat1 <-
structure(list(conc = c(500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 
3.90625, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 
500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625), 
response =
c(2.756, 2.167, 1.38, 0.873, 0.571, 0.43, 0.361, 0.326, 2.82, 2.174, 1.402, 
0.911, 0.593, 0.458, 0.387, 0.348, 2.732, 2.143, 1.419, 0.874, 
0.582, 0.442, 0.366, 0.331)), 
.Names = c("conc", "response"), 
row.names = as.integer(c(1, 2, 3, 4, 5, 6, 7, 8, 17, 18, 19, 20, 21, 22, 23, 24, 33, 34, 
35, 36, 37, 38, 39, 40)), 
class = "data.frame")

m1 <- drm(response~conc, data=dat1, fct=LL.4())
#m2 <- drm(response~conc, data=dat1, fct=LL.4(), adjust="vp")  # huge standard errors


## Second dataset
dat2 <-
structure(list(conc = c(500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 
3.90625, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 
500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625), 
response = c(2.943, 2.337, 1.521, 0.989, 0.669, 0.481, 0.413, 0.36, 2.952, 2.272, 
1.518, 0.974, 0.648, 0.493, 0.413, 0.36, 2.943, 2.309, 1.505, 
0.979, 0.649, 0.478, 0.387, 0.34)), 
.Names = c("conc", "response" ), 
row.names = as.integer(c(49, 50, 51, 52, 53, 54, 55, 56, 65, 
66, 67, 68, 69, 70, 71, 72, 81, 82, 83, 84, 85, 86, 87, 88)), 
class = "data.frame")

m3 <- drm(response~conc, data=dat2, fct=LL.4())
#m4 <- drm(response~conc, data=dat2, fct=LL.4(), adjust="vp")  # huge standard errors
