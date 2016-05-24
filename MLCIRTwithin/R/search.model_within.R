search.model_within <- function (S, yv = rep(1, ns), kv1, kv2, X = NULL, link = c("global","local"),
                                  disc = FALSE, difl = FALSE, multi1, multi2, fort = FALSE, 
                                  tol1 = 10^-6, tol2 = 10^-10, glob = FALSE, disp = FALSE,
									output = FALSE, out_se = FALSE, nrep = 2){
									
# preliminaries
	link = match.arg(link)
	ns = dim(S)[1]
	J = dim(S)[2]
# store results with tol1
	k = 0
	if(min(kv1) > 1) kv1 = c(1, kv1)
	if(min(kv2) > 1) kv2 = c(1, kv2)
    k1m = length(kv1); k2m = length(kv2)
	out = vector("list", k1m*k2m)
	lkv = aicv = bicv = entv = necv = errv = rep(NA, k1m*k2m)
	for(k1 in kv1) for(k2 in kv2){
	    	cat("***************************************************************************\n")
		cat(c(k1,k2), "\n")
		k = k+1
		if(k1==1 & k2==1){
			out[[k]] = try(est_multi_poly(S = S, yv = yv, k = 1,tol=tol2))
       	}else if(k1==1 & k2>1){
			out[[k]] = try(est_multi_poly_between(S = S, yv = yv, k = k2,
							X = X, start = "deterministic", link = link, disc = disc, difl = difl,
        	               multi = multi2, fort = fort, tol = tol1, glob = glob, disp = disp))
       	}else if(k1>1 & k2==1){
			out[[k]] = try(est_multi_poly_between(S = S, yv = yv, k = k1,
							X = X, start = "deterministic", link = link, disc = disc, difl = difl,
        	               multi = multi1, fort = fort, tol = tol1, glob = glob, disp = disp))
       	}else{ 
			out[[k]] = try(est_multi_poly_within(S = S, yv = yv, k1 = k1, k2 = k2,
							X = X, start = "deterministic", link = link, disc = disc, difl = difl,
        	               multi1 = multi1, multi2 = multi2, fort = fort, tol = tol1, glob = glob, disp = disp))
		}
		if (!inherits(out[[k]], "try-error")){
			errv[[k]] = FALSE
		}else{
			errv[[k]] = TRUE
			if (k > 1) out[[k]] = out[[k - 1]]
		}
		lktrace = out[[k]]$lk
		lkv[k] = out[[k]]$lk
		aicv[k] = out[[k]]$aic
		bicv[k] = out[[k]]$bic
		entv[k] = out[[k]]$ent
		if (k == 1) necv[1] = 1
		else necv[k] = entv[k]/(lkv[k] - lkv[1])
		cat("lktrace = ", sort(lktrace), "\n")
		cat("lk = ", lkv, "\n")
		cat("aic = ", aicv, "\n")
		cat("bic = ", bicv, "\n")
		cat("ent = ", entv, "\n")
		cat("nec = ", necv, "\n")
		save(file = "search.model_within.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
		if(k1>1 | k2>1){
			if(nrep>0) for (h in 1:(nrep * (k1*k2 - 1))){
				cat("***************************************************************************\n")
				cat(c(k1,k2,h), "\n")
				if(k1==1 & k2>1){
					outh = try(est_multi_poly_between(S = S, yv = yv, k = k2,
							X = X, start = "random", link = link, disc = disc, difl = difl,
        	               multi = multi2, fort = fort, tol = tol1, glob = glob, disp = disp))
	       		}else if(k1>1 & k2==1){
					outh = try(est_multi_poly_between(S = S, yv = yv, k = k1,
							X = X, start = "random", link = link, disc = disc, difl = difl,
        	               multi = multi1, fort = fort, tol = tol1, glob = glob, disp = disp))
       			}else{ 
					outh = try(est_multi_poly_within(S = S, yv = yv, k1 = k1, k2 = k2,
							X = X, start = "random", link = link, disc = disc, difl = difl,
        	               multi1 = multi1, multi2 = multi2, fort = fort, tol = tol1, glob = glob, disp = disp))
				}
           	    if(!inherits(outh, "try-error")){
					lktrace = c(lktrace, outh$lk)
					if (outh$lk > out[[k]]$lk) out[[k]] = outh
				}
				lkv[k] = out[[k]]$lk
				aicv[k] = out[[k]]$aic
				bicv[k] = out[[k]]$bic
				entv[k] = out[[k]]$ent
				necv[k] = entv[k]/(lkv[k] - lkv[1])
				cat("lktrace = ", sort(lktrace), "\n")
				cat("lk = ", lkv, "\n")
				cat("aic = ", aicv, "\n")
				cat("bic = ", bicv, "\n")
				cat("ent = ", entv, "\n")
				cat("nec = ", necv, "\n")
				save(file = "search.model_within.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
          	}
          	if(tol2<tol1 || output || out_se){	
          		cat("***************************************************************************\n")
          		cat(k1,k2,"\n")
				if(k1==1 & k2>1){
					outh = try(est_multi_poly_between(S = S, yv = yv, k = k2,
							X = X, start = "external", link = link, disc = disc, difl = difl,
        	               multi = multi2, fort = fort, tol = tol2, glob = glob, disp = disp,
							output=output,out_se=out_se,piv=out[[k]]$piv,Phi=out[[k]]$Phi,
							gac=out[[k]]$gac,De=out[[k]]$De))
	       		}else if(k1>1 & k2==1){
					outh = try(est_multi_poly_between(S = S, yv = yv, k = k1,
							X = X, start = "external", link = link, disc = disc, difl = difl,
        	               multi = multi1, fort = fort, tol = tol2, glob = glob, disp = disp,
        	               output=output,out_se=out_se,piv=out[[k]]$piv,Phi=out[[k]]$Phi,
							gac=out[[k]]$gac,De=out[[k]]$De))
       			}else{ 
					outh = try(est_multi_poly_within(S = S, yv = yv, k1 = k1, k2 = k2,
							X = X, start = "external", link = link, disc = disc, difl = difl,
        	               multi1 = multi1, multi2 = multi2, fort = fort, tol = tol2, glob = glob, disp = disp,
        	               output=output,out_se=out_se,piv1=out[[k]]$piv1,piv2=out[[k]]$piv2,Phi=out[[k]]$Phi,
							ga1c=out[[k]]$ga1c,ga2c=out[[k]]$ga2c,De1=out[[k]]$De1,De2=out[[k]]$De2))
				}
           	    if(!inherits(outh, "try-error")){
					lktrace = c(lktrace, outh$lk)
					if (outh$lk > out[[k]]$lk) out[[k]] = outh
				}
				lkv[k] = out[[k]]$lk
				aicv[k] = out[[k]]$aic
				bicv[k] = out[[k]]$bic
				entv[k] = out[[k]]$ent
          	 	necv[k] = entv[k]/(lkv[k] - lkv[1])
          	 	cat("lktrace = ", sort(lktrace), "\n")
          	 	cat("lk = ", lkv, "\n")
          	 	cat("aic = ", aicv, "\n")
          	 	cat("bic = ", bicv, "\n")
          	 	cat("ent = ", entv, "\n")
          	 	cat("nec = ", necv, "\n")
          	 	save(file = "search.model_within.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
          	}
		}
        out[[k]]$lktrace = lktrace
        out[[k]]$k1 = k1
        out[[k]]$k2 = k2
        save(file = "search.model_between.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
    }
    out = list(out.single = out, aicv = aicv, bicv = bicv, entv = entv, necv = necv, lkv = lkv, errv = errv)
}