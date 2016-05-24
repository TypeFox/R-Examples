
 library(aster)

 set.seed(42)

 m <- 10
 n <- 5
 a <- matrix(rnorm(m * n), nrow = m)
 b <- rnorm(n)

 out <- .C("aster_mat_vec_mult",
     nrow = as.integer(m),
     ncol = as.integer(n),
     a = as.double(a),
     b = as.double(b),
     c = double(m))

 all.equal(out$c, as.numeric(a %*% b))

 ##########

 b <- rnorm(m)

 out <- .C("aster_vec_mat_mult",
     nrow = as.integer(m),
     ncol = as.integer(n),
     a = as.double(a),
     b = as.double(b),
     c = double(n))

 all.equal(out$c, as.numeric(b %*% a))

 ##########

 out <- .C("aster_mat_vec_mat_mult",
     nrow = as.integer(m),
     ncol = as.integer(n),
     a = as.double(a),
     b = as.double(b),
     c = matrix(as.double(0), n, n))

 all.equal(out$c, t(a) %*% diag(b) %*% a)

 ##########

 b <- matrix(rnorm(n * n), n)

 out <- .C("aster_diag_mat_mat_mat_mult",
     nrow = as.integer(m),
     ncol = as.integer(n),
     a = as.double(a),
     b = as.double(b),
     c = double(m))

 all.equal(out$c, diag(a %*% b %*% t(a)))

