search.model_between <- function (S, yv = rep(1, ns), kv, X = NULL, link = c("global","local"),
                                  disc = FALSE, difl = FALSE, multi = 1:J, fort = FALSE, 
                                  tol1 = 10^-6, tol2 = 10^-10, glob = FALSE, disp = FALSE,
									output = FALSE, out_se = FALSE, nrep = 2){
									
# preliminaries
	link = match.arg(link)
	ns = dim(S)[1]
	J = dim(S)[2]
# store results with tol1
	out = vector("list", max(kv))
	if(min(kv) > 1) kv = c(1, kv)
	lkv = aicv = bicv = entv = necv = errv = rep(NA, max(kv))
	for(k in kv){
    	cat("***************************************************************************\n")
		cat(k, "\n")
		if(k==1){
			out[[k]] = try(est_multi_poly(S = S, yv = yv, k = 1,tol=tol2))
       	}else{
			out[[k]] = try(est_multi_poly_between(S = S, yv = yv, k = k,
							X = X, start = "deterministic", link = link, disc = disc, difl = difl,
        	               multi = multi, fort = fort, tol = tol1, glob = glob, disp = disp))
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
		else if (1 %in% kv) necv[k] = entv[k]/(lkv[k] - lkv[1])
		cat("lktrace = ", sort(lktrace), "\n")
		cat("lk = ", lkv, "\n")
		cat("aic = ", aicv, "\n")
		cat("bic = ", bicv, "\n")
		cat("ent = ", entv, "\n")
		cat("nec = ", necv, "\n")
		save(file = "search.model_bewteen.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
		if(k > 1){
			if(nrep>0) for (h in 1:(nrep * (k - 1))){
				cat("***************************************************************************\n")
				cat(c(k, h), "\n")
				outh = try(est_multi_poly_between(S = S, yv = yv, k = k, 
						    X = X, start = "random", link = link, disc = disc, difl = difl, 
           	               multi = multi, fort = fort, tol = tol1, glob = glob, disp = disp))
           	    if(!inherits(outh, "try-error")){
					lktrace = c(lktrace, outh$lk)
					if (outh$lk > out[[k]]$lk) out[[k]] = outh
				}
				lkv[k] = out[[k]]$lk
				aicv[k] = out[[k]]$aic
				bicv[k] = out[[k]]$bic
				entv[k] = out[[k]]$ent
				if (1 %in% kv) necv[k] = entv[k]/(lkv[k] - lkv[1])
				cat("lktrace = ", sort(lktrace), "\n")
				cat("lk = ", lkv, "\n")
				cat("aic = ", aicv, "\n")
				cat("bic = ", bicv, "\n")
				cat("ent = ", entv, "\n")
				cat("nec = ", necv, "\n")
				save(file = "search.model_between.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
          	}
          	if(tol2<tol1 || output || out_se){	
          		cat("***************************************************************************\n")
          		cat(k, "\n")
          		outh = try(est_multi_poly_between(S = S, yv = yv, k = k, 
           	            X = X, start = "external", link = link, disc = disc, difl = difl, 
           	            multi = multi, fort = fort, tol = tol2, glob = glob, disp = disp,
           	            output=output,out_se=out_se,piv=out[[k]]$piv,Phi=out[[k]]$Phi,
           	            gac=out[[k]]$gac,De=out[[k]]$De))
           	    if (!inherits(outh, "try-error")){
					lktrace = c(lktrace, outh$lk)
          	     	if (outh$lk > out[[k]]$lk) out[[k]] = outh
				}
				lkv[k] = out[[k]]$lk
				aicv[k] = out[[k]]$aic
				bicv[k] = out[[k]]$bic
				entv[k] = out[[k]]$ent
          	 	if (1 %in% kv) necv[k] = entv[k]/(lkv[k] - lkv[1])
          	 	cat("lktrace = ", sort(lktrace), "\n")
          	 	cat("lk = ", lkv, "\n")
          	 	cat("aic = ", aicv, "\n")
          	 	cat("bic = ", bicv, "\n")
          	 	cat("ent = ", entv, "\n")
          	 	cat("nec = ", necv, "\n")
          	 	save(file = "search.model_between.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
          	}
		}
        out[[k]]$lktrace = lktrace
        save(file = "search.model_between.temp.RData", out, aicv, bicv, entv, necv, errv, lkv)
    }
    out = list(out.single = out, aicv = aicv, bicv = bicv, entv = entv, necv = necv, lkv = lkv, errv = errv)
}