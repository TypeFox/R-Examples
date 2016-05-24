
##
library("slam")

##
x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 2,
            dimnames = list(A = 1:2, B = 1:3))
x

a <- as.simple_sparse_array(x)
a

##
z <- rollup(x, 2L, c(1,2,1), na.rm = TRUE)
z
identical(as.array(z),
	  as.array(rollup(a, 2L, c(1,2,1), na.rm = TRUE)))
identical(as.array(z),
	  as.array(rollup(a, 2L, c(1,2,1), na.rm = TRUE, EXPAND = "dense")))
identical(as.array(z),
	  as.array(rollup(a, 2L, c(1,2,1), na.rm = TRUE, EXPAND = "all")))

##
z <- rollup(x, 2L, c(1,NA,1), na.rm = TRUE)
z
identical(as.array(z),
          as.array(rollup(a, 2L, c(1,NA,1), na.rm = TRUE)))
identical(as.array(z),
          as.array(rollup(a, 2L, c(1,NA,1), na.rm = TRUE, EXPAND = "dense")))
identical(as.array(z),
          as.array(rollup(a, 2L, c(1,NA,1), na.rm = TRUE, EXPAND = "all")))

##
z <- rollup(x, 2L, c(1,NA,1), na.rm = TRUE, DROP = TRUE)
identical(as.array(z),
          as.array(rollup(a, 2L, c(1,NA,1), na.rm = TRUE, DROP = TRUE)))


##
z <- rollup(x, 1:2, list(1:2, c(1,2,1)), na.rm = TRUE)
identical(as.array(z),
	  as.array(rollup(a, 1:2, list(1:2, c(1,2,1)), na.rm = TRUE)))

##
s <- as.simple_triplet_matrix(a)
z <- rollup(x, 2L, FUN = min, na.rm = TRUE)
identical(as.matrix(z),
	  as.matrix(rollup(s, 2L, FUN = min, na.rm = TRUE, EXPAND = "dense")))

###
