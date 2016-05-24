
 library(rcdd)

 ##### H-representation #####

 hrep <- rbind(
 c(0,  0,  1,  1,  0),
 c(0,  0, -1,  0,  0),
 c(0,  0,  0, -1,  0),
 c(0,  0,  0,  0, -1),
 c(0,  0, -1, -1, -1))

 redundant(d2q(hrep), representation = "H")

 redundant(hrep, representation = "H")

 ##### V-representation #####

 foo <- c(1, 0, -1)
 hrep <- cbind(0, 1, rep(foo, each = 9), rep(foo, each = 3), foo)
 print(hrep)

 redundant(d2q(hrep), representation = "V")

 redundant(hrep, representation = "V")

 ##### another V-representation #####

 hrep <- rbind(
 c(0,  0,  1,  0,  0),
 c(0,  0,  0,  1,  0),
 c(0,  0,  0,  0,  1),
 c(0,  0, -1, -1, -1))

 redundant(d2q(hrep), representation = "V")

 redundant(hrep, representation = "V")

 ##### negative new position #####

 hrep <- rbind(
 c(1, 0, 1, 0, 0),
 c(1, 0, 0, 1, 0),
 c(1, 0, 0, 0, 1),
 c(1, 0, 0, 1, 0),
 c(1, 0, 0, 0, 1),
 c(1, 0, 0, 0, 1))

 redundant(hrep, representation = "V")

