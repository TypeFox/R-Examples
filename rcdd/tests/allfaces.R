
 library(rcdd)

 hrep <- rbind(c(0, 1,  1,  1, -1),
               c(0, 1,  1, -1, -1),
               c(0, 1, -1, -1, -1),
               c(0, 1, -1,  1, -1),
               c(0, 0,  0,  0,  1))

 qout <- allfaces(d2q(hrep))
 print(qout)

 dout <- allfaces(hrep)

 identical(qout$dimension, dout$dimension)
 identical(qout$active.set, dout$active.set)
 all.equal(lapply(qout$relative.interior.point, FUN = q2d),
     dout$relative.interior.point)

