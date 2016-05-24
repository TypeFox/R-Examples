
library("slam")

##
.Call("_part_index", factor(1:4))
.Call("_part_index", factor(c(1L,2L,2L,1L)))
.Call("_part_index", factor(c(1L,2L,NA,1L)))

##
i <- 1:27
x <- arrayInd(i, .dim = c(3L,3L,3L))
.Call("_vector_index", c(3L,3L,3L), x)
x[14L, 2L] <- NA
.Call("_vector_index", c(3L,3L,3L), x)

##
v <- c(1L,1L)
p <- matrix(c(1L,2L,3L, 2L,2L,2L), nrow = 2L, byrow = TRUE)
.Call("_ini_array", c(3L,3L,3L), p, v, 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.logical(v), 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.double(v), 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.raw(v), 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.complex(v), 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.character(v), 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.list(v), 2L)
.Call("_ini_array", c(3L,3L,3L), p, as.expression(v), 2L)

.Call("_ini_array", 3L, c(1L,2L), c(1L,1L), 2L)

.Call("_split_col", array(c(1L,2L), c(2L, 2L)))

##
x <- matrix(c(1L,1L,1L,1L,1L,2L,1L,3L,1L,2L), 
	    ncol = 2, byrow = TRUE)
x
.Call("_match_matrix", x, NULL, NULL)
.Call("_match_matrix", x, x[1:3,], 0L)


##
x <- matrix(c(1L,2L,2L,2L,NA,1L,NA,2L,NA,NA), 
	    ncol = 2, byrow = TRUE)
x
.Call("_all_row", x > 1L, FALSE)
.Call("_all_row", x > 1L, TRUE)

###
