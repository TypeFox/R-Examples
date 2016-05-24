
library("slam")

x <- simple_triplet_diag_matrix(1, nrow = 3L)
rownames(x) <- letters[1:3]

identical(as.matrix(cbind(x, x)),
	  cbind(as.matrix(x), as.matrix(x)))
identical(as.matrix(rbind(t(x), t(x))),
	  rbind(as.matrix(t(x)), as.matrix(t(x))))

identical(as.matrix(cbind(x, t(x))),
	  cbind(as.matrix(x), as.matrix(t(x))))
identical(as.matrix(rbind(t(x), x)),
	  rbind(as.matrix(t(x)), as.matrix(x)))

###
