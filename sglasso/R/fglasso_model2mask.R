fglasso_model2mask <- function(fglasso.model, tp, p){
	mask <- matrix(".", p * tp, p * tp)
	s_lbl <- rep("s", length = p)
	s_u_lbl <- paste(s_lbl, "_", 1:p, sep = "")
	n_lbl <- rep("n", length = p * (p - 1) / 2)
	lbl <- outer(1:p, 1:p, function(i,j) paste(i, j, sep = "."))
	U <- upper.tri(lbl)
	U <- outer(1:p, 1:p, "<")
	L <- t(U)
	n_u_lbl <- paste(n_lbl, "_", lbl[U], sep = "")
	n_l_lbl <- paste(n_lbl, "_", lbl[L], sep = "")
	for(rb in 0:(tp - 1)){
		k <- rb * p + 1:p
		diag(mask[k, k]) <- switch(fglasso.model[1, 1],
								   c = paste(s_lbl, "-h0", sep = ""),
								   u = paste(s_u_lbl, "-h0", sep = ""),
								   t = paste(s_lbl, "-t", rb + 1, "-h0", sep = ""),
								   ut = paste(s_u_lbl, "-t", rb + 1, "-h0", sep = ""),
								   . = ".")
		mask[k, k][U] <- switch(fglasso.model[1, 2],
								 c = paste(n_lbl, "-h0", sep = ""),
								 u = paste(n_u_lbl, "-h0", sep = ""),
								 t = paste(n_lbl, "-t", rb + 1, "-h0", sep = ""),
								 ut = paste(n_u_lbl, "-t", rb + 1, "-h0", sep = ""),
								 . = ".")
	}
	nlag <- dim(fglasso.model)[1]
	if(nlag > 1){
		for(lag in 1:(nlag - 1)){
			for(rb in 0:(tp - lag - 1)){
				s <- rb * p + 1:p
				e <- (rb + lag) * p + 1:p
				diag(mask[s, e]) <- switch(fglasso.model[lag + 1, 1],
										   c = paste(s_lbl, "-h", lag, sep = ""),
										   u = paste(s_u_lbl, "-h", lag, sep = ""),
										   t = paste(s_lbl, "-t", rb + 1, "-h", lag, sep = ""),
										   ut = paste(s_u_lbl, "-t", rb + 1, "-h", lag, sep = ""),
										   . = ".")
				mask[s, e][U] <- switch(fglasso.model[lag + 1, 2],
										c = paste(n_lbl, "-h", lag, sep = ""),
										u = paste(n_u_lbl, "-h", lag, sep = ""),
										t = paste(n_lbl, "-t", rb + 1, "-h", lag, sep = ""),
										ut = paste(n_u_lbl, "-t", rb + 1, "-h", lag, sep = ""),
										. = ".")
				mask[s, e][L] <- switch(fglasso.model[lag + 1, 2],
										c = paste(n_lbl, "-h", lag, sep = ""),
										u = paste(n_l_lbl, "-h", lag, sep = ""),
										t = paste(n_lbl, "-t", rb + 1, "-h", lag, sep = ""),
										ut = paste(n_l_lbl, "-t", rb + 1, "-h", lag, sep = ""),
										. = ".")
			}
		}
	}
	mask[lower.tri(mask)] <- t(mask)[lower.tri(mask)]
	mask
}
	
	
	
