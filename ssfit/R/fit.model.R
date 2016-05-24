fit.model <-
function (p, q, n, r, starting_values, h_vector, data_true, sim_data, 
    features, n_iter, print_results = TRUE, variances = TRUE) 
{
    h <- list()
    for (s in 1:p) {
        h[[s]] <- h_vector[s]
    	}

    theta_0 <- list()
    for (s in 1:p) {
        theta_0[[s]] <- starting_values[s]
    	}

    k <- 2^ceiling(log2(p + 1))

    # library(survey)
    if (p == 1) {
        D <- matrix(c(1, -1), nrow = 2, ncol = 1)
    	}
    if (p > 1) {
        D <- 2 * hadamard(k - 1) - 1
        D <- D[, 2:(p + 1)]
    	}

    sim <- list()
    z <- list()

    results <- matrix(NA, nrow = n_iter, ncol = p)
    k_ones <- matrix(1, nrow = k, ncol = 1)
    r_ones <- matrix(1, nrow = r, ncol = 1)
    design.points <- list()
    z_data <- matrix(NA, nrow = 1, ncol = q)

    for (j in 1:q) {
        z_data[1, j] <- features(data_true)[j]
    }
    z_data <- t(z_data)

    if (variances == TRUE) {
	var_theta <- list()
        }

    for (iter in 1:n_iter) {

        for (d in 1:k) {
            sim[[d]] <- list()
            z[[d]] <- matrix(NA, nrow = r, ncol = q)
	    }

        for (d in 1:k) {
            design.points[[d]] <- list()
            for (s in 1:p) {
                design.points[[d]][[s]] <- theta_0[[s]] + (D[d, s] * h[[s]])
            	}
            for (i in 1:r) {
                sim[[d]][[i]] <- sim_data(n, design.points[[d]])
                }
            }

        for (d in 1:k) {
            for (i in 1:r) {
                z[[d]][i, ] <- features(sim[[d]][[i]])
            	}
            }

        zbar <- list()
        for (d in 1:k) {
            zbar[[d]] <- (t(r_ones) %*% z[[d]])/r
	    }

	Zbar <- matrix(unlist(zbar), nrow = q, byrow = F)

        ztilde <- list()
        sigma <- list()
        for (d in 1:k) {
		ztilde[[d]] <- z[[d]] - (r_ones %*% zbar[[d]])
 	        sigma[[d]] <- (t(ztilde[[d]]) %*% ztilde[[d]])/(r - 1)
	        }

	V <- Zbar %*% D

  	sigma_vv <- (1/k) * Reduce("+", sigma)
        sigma_inv <- solve(sigma_vv)

	lambda <- rep(NA, times = p)
        for (s in 1:p) {
            lambda[s] <- sqrt((t(V[, s]) %*% sigma_inv %*% V[, s]))
	    }

	L <- (1/lambda) * (sigma_inv %*% V)
	Ybar <- t(L) %*% Zbar
 
   	Y_data <- t(L) %*% z_data

        A <- (2/k) * (t(V) %*% L)
        for (u in 1:dim(A)[1]) {
            for (v in 1:dim(A)[2]) {
                A[u, v][u != v] <- 0
            	}
            }

        H <- matrix(0, nrow = p, ncol = p)
        for (s in 1:p) {
            H[s, s] <- h[[s]]
            }

        theta_hat <- matrix(NA, nrow = p, ncol = 1)
        th_0 <- matrix(unlist(theta_0))

        theta_hat <- th_0 + (2 * H %*% solve(A) %*% (Y_data - ((Ybar %*% k_ones)/k)))

        for (s in 1:p) {
 	       theta_0[[s]] <- theta_hat[s]
 	       }

        results[iter, ] <- t(theta_hat)

        if (variances == TRUE) {
		var_theta[[iter]] <- 4 * H %*% solve(A) %*% t(L) %*% sigma_vv %*% L %*% H %*% solve(A)
        	}

        if (print_results == TRUE) {
            	print(results[iter, ])
        	}
    }

        if (variances == FALSE) {
	    var_theta <- 4 * H %*% solve(A) %*% t(L) %*% sigma_vv %*% L %*% H %*% solve(A)
            }

    return(list(estimates = results, var_estimates = var_theta, L = L, sigma = sigma_vv, zbar = Zbar, z_D = z_data, ybar = Ybar, y_D = Y_data))

}

