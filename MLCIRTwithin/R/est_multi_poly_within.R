est_multi_poly_within <- function (S, yv = rep(1, ns), k1, k2, X = NULL, 
    start = c("deterministic","random","external"), link = c("global","local"), 
    disc = FALSE, difl = FALSE, multi1, multi2, piv1 = NULL, piv2 = NULL, Phi = NULL, ga1c = NULL, 
    ga2c = NULL, De1 = NULL, De2 = NULL, fort = FALSE, tol = 10^-10, disp = FALSE, 
    output = FALSE, out_se = FALSE, glob = FALSE){
    	
# check input
    if (k1 == 1 & k2 == 1) stop("--> for LC models use est_multi_poly (package MultiLCIRT)")
    if (k1 == 1 | k2 == 1) stop("--> for k1=1 or k2=1 use est_multi_poly_between")
    if (max(S, na.rm = TRUE) == 1 & difl != 0) {
        warning("with binary data put difl=0\n")
        difl = FALSE
    }
    link0 = match.arg(link)
	if(link0=="global") link=1
	if(link0=="local") link=2    
    start0 = match.arg(start)
	if(start0=="deterministic") start=0
	if(start0=="random") start=1    
	if(start0=="external") start=2    
    multi1v = as.vector(multi1)
    if(any(multi1v==0)) multi1v = sort(multi1v[-which(multi1v==0)])
    multi2v = as.vector(multi2)
    if(any(multi2v==0)) multi2v = sort(multi2v[-which(multi2v==0)])
    multiv = union(multi1v,multi2v)
	if(length(multi1v)==length(multi2v)) if(all(multi1v==multi2v)) warning("multi1 and multi2 cannot be equal")
	if(is.vector(multi1)) cons1 = multi1[1] else cons1 = multi1[,1]
	if(is.vector(multi2)) cons2 = multi2[1] else cons2 = multi2[,1]
	if(length(intersect(cons1,cons2))>0) warning("There are repetitions in the constrained items")
# adjust covariates	
    cov = !is.null(X)
    if (cov) {
        X = as.matrix(X)
        namesX = colnames(X)
        if (glob) logit_cov = "g" else logit_cov = "m"
    }else{
        logit_cov = "m"
    }
    miss = any(is.na(S))
    ns = nrow(S)
    J = ncol(S)
    if (miss) {
        cat("Missing data in the dataset, units and items without responses are removed\n")
        ind = which(apply(is.na(S), 1, all))
        if (length(ind) > 0) {
            S = S[-ind, ]
            yv = yv[-ind]
            if (!is.null(X)) 
                X = as.matrix(X[-ind, ])
            ind = which(apply(is.na(S), 2, all))
            if (length(ind) > 0) {
                S = S[, -ind]
                miss = any(is.na(S))
            }
        }
    }
    if (miss) {
        R = 1 * (!is.na(S))
        S[is.na(S)] = 0
    }
    lv = apply(S, 2, max) + 1
    if(difl){
    	flag = FALSE
    	if(is.vector(multi1)){
    		tmp = lv[multi1]
    		if(max(tmp)>min(tmp)) flag=TRUE
    	}else{
    		for(h in 1:nrow(multi1)){
    			tmp = multi1[h,]; tmp = tmp[tmp>0]
	    		tmp = lv[multi1[h,]]
    			if(max(tmp)>min(tmp)) flag=TRUE    			
    		}
    	} 
    	if(is.vector(multi2)){
    		tmp = lv[multi2]
    		if(max(tmp)>min(tmp)) flag=TRUE
    	}else{
    		for(h in 1:nrow(multi2)){
    			tmp = multi2[h,]; tmp = tmp[tmp>0]
	    		tmp = lv[multi2[h,]]
    			if(max(tmp)>min(tmp)) flag=TRUE    			
    		}
    	}    	 
    	if(flag) stop("Option difl=TRUE not allowed in the presence of items with a different number of response categories within the same dimension")
    }
    lm = max(lv)
    ns = nrow(S)
    J = ncol(S)
    n = sum(yv)
    if (cov) {
        ncov = ncol(X)
        out = aggr_data(X, fort = fort)
        Xdis = as.matrix(out$data_dis)
        Xlabel = out$label
        Xndis = max(out$label)
        if (glob) {
            XX1dis = array(0, c(k1 - 1, k1 - 1 + ncov, Xndis))
            for (i in 1:Xndis) XX1dis[, , i] = cbind(diag(k1 - 
                1), rep(1, k1 - 1) %o% Xdis[i, ])
        }
        else {
            XX1dis = array(0, c(k1 - 1, (k1 - 1) * (ncov + 1), 
                Xndis))
            if (k1 == 2) 
                II = 1
            else II = diag(k1 - 1)
            for (i in 1:Xndis) XX1dis[, , i] = II %x% t(c(1, 
                Xdis[i, ]))
        }
        if (glob) {
            XX2dis = array(0, c(k2 - 1, k2 - 1 + ncov, Xndis))
            for (i in 1:Xndis) XX2dis[, , i] = cbind(diag(k2 - 
                1), rep(1, k2 - 1) %o% Xdis[i, ])
        }
        else {
            XX2dis = array(0, c(k2 - 1, (k2 - 1) * (ncov + 1), 
                Xndis))
            if (k2 == 2) 
                II = 1
            else II = diag(k2 - 1)
            for (i in 1:Xndis) XX2dis[, , i] = II %x% t(c(1, 
                Xdis[i, ]))
        }
    }
    else {
        ncov = 0
        XX1dis = array(diag(k1 - 1), c(k1 - 1, k1 - 1, 1))
        XX2dis = array(diag(k2 - 1), c(k2 - 1, k2 - 1, 1))
        Xlabel = rep(1, ns)
    }
    Aggr1 = diag(k1) %x% matrix(1, 1, k2)
    Aggr2 = matrix(1, 1, k1) %x% diag(k2)
    if (link == 1) 
        ltype = "g"
    else if (link == 2) 
        ltype = "l"
    if (link == 1 || link == 2) {
        items1 = sort(unique(as.vector(multi1)))
        if (any(items1 == 0)) 
            items1 = items1[items1 > 0]
        J1 = length(items1)
        if (is.vector(multi1)) 
            rm1 = 1
        else rm1 = nrow(multi1)
        items2 = sort(unique(as.vector(multi2)))
        if (any(items2 == 0)) 
            items2 = items2[items2 > 0]
        J2 = length(items2)
        if (is.vector(multi2)) 
            rm2 = 1
        else rm2 = nrow(multi2)
        Dem1 = matrix(0, J, rm1)
        if (rm1 == 1) {
            Dem1 = 1
            fv1 = multi1[1]
        }
        else {
            for (r in 1:rm1) {
                ind = multi1[r, ]
                ind = ind[ind > 0]
                Dem1[ind, r] = 1
            }
            fv1 = multi1[, 1]
        }
        fv1e = NULL
        count = 0
        for (j in 1:J) {
            if (j %in% fv1) 
                fv1e = c(fv1e, count + 1)
            count = count + lv[j] - 1
        }
        Dem2 = matrix(0, J, rm2)
        if (rm2 == 1) {
            Dem2 = 1
            fv2 = multi2[1]
        }
        else {
            for (r in 1:rm2) {
                ind = multi2[r, ]
                ind = ind[ind > 0]
                Dem2[ind, r] = 1
            }
            fv2 = multi2[, 1]
        }
        fv2e = NULL
        count = 0
        for (j in 1:J) {
            if (j %in% fv2) 
                fv2e = c(fv2e, count + 1)
            count = count + lv[j] - 1
        }
        fv = union(fv1, fv2)
        fve = union(fv1e, fv2e)
        rm = length(fve)
        indga1 = 1:J
        indga1 = indga1[setdiff(items1, fv1)]
        indga2 = 1:J
        indga2 = indga2[setdiff(items2, fv2)]
        indth1 = 1:(k1 * rm1)
        indth2 = (k1 * rm1 + 1):(k1 * rm1 + k2 * rm2)
        if(!difl){
            indbe = k1 * rm1 + k2 * rm2 + (1:(sum(lv - 1) - rm1 - rm2))
            indbec = 1:sum(lv - 1)
            indbec = indbec[-fve]
        }
        else {
            indbe = k1 * rm1 + k2 * rm2 + (1:(J - rm + sum(lv[fv] - 
                2)))
            indbec = 1:J
            indbec = indbec[-fv]
        }
        abils1 = rep(0, J)
        if (rm1 == 1) {
            abils1[multi1] = 1
        }
        else {
            for (h in 1:rm1) {
                ind = multi1[h, ]
                ind = ind[ind > 0]
                abils1[ind] = h
            }
        }
        abils2 = rep(0, J)
        if (rm2 == 1) {
            abils2[multi2] = 1
        }
        else {
            for (h in 1:rm2) {
                ind = multi2[h, ]
                ind = ind[ind > 0]
                abils2[ind] = h
            }
        }
        abils = rep(0, J)
        for (j in 1:J) for (h in 1:rm) if (abils1[j] == abils1[fv[h]] & 
            abils2[j] == abils2[fv[h]]) 
            abils[j] = h
        if (!difl) 
            ZZ = array(NA, c(lm - 1, k1 * rm1 + k2 * rm2 + sum(lv - 1) - rm, J * k1 * k2))
        if (difl) 
            ZZ = array(NA, c(lm - 1, k1 * rm1 + k2 * rm2 + J - rm + sum(lv[fv] - 2), J * k1 * k2))
        cont = 0
        refitem = matrix(0, J * k1 * k2, 1)
        for (c1 in 1:k1) {
            u11 = matrix(0, 1, k1)
            u11[c1] = 1
            for (c2 in 1:k2) {
                u12 = matrix(0, 1, k2)
                u12[c2] = 1
                for (j in 1:J) {
                  u21 = matrix(0, 1, rm1)
                  u21[abils1[j]] = 1
                  u22 = matrix(0, 1, rm2)
                  u22[abils2[j]] = 1
                  v = matrix(0, 1, J)
                  v[j] = 1
                  cont = cont + 1
                  if (!difl) {
                    Te = matrix(0, lv[j] - 1, sum(lv - 1))
                    if (j == 1) 
                      Te[, 1:(lv[j] - 1)] = diag(lv[j] - 1)
                    else Te[, sum(lv[1:(j - 1)] - 1) + (1:(lv[j] - 
                      1))] = diag(lv[j] - 1)
                    Te = matrix(Te[, -fve], lv[j] - 1, dim(Te)[2] - 
                      length(fve))
                  }
                  else if (difl) {
                    Te = matrix(0, lv[j] - 1, sum(lv[fv] - 2))
                    if (lv[j] > 2) {
                      if (abils[j] == 1) Te[, 1:(lv[j] - 2)] = diag(lv[j] - 1)[, -1]
                      else Te[, sum(lv[1:(abils[j] - 1)] - 2) + (1:(lv[j] - 2))] = diag(lv[j] - 1)[, -1]
                    }
                    Te = cbind(v %x% rep(1, lv[j] - 1), Te)
                    Te = matrix(Te[, -fv], lv[j] - 1)
                  }
                  ZZ[1:(lv[j] - 1), , cont] = cbind(rep(1, lv[j] - 
                    1) %*% (u11 %x% u21), rep(1, lv[j] - 1) %*% 
                    (u12 %x% u22), -Te)
                  refitem[cont] = j
                }
            }
        }
    }
    ZZ0 = ZZ
    if (glob) {
        II1 = diag(k1 - 1)
        II1 = cbind(0, II1) - cbind(II1, 0)
        if (rm1 > 1) 
            II1 = II1 %x% matrix(c(1, rep(0, rm1 - 1)), 1, rm1)
        II2 = diag(k2 - 1)
        II2 = cbind(0, II2) - cbind(II2, 0)
        if (rm2 > 1) 
            II2 = II2 %x% matrix(c(1, rep(0, rm1 - 1)), 1, rm1)
        II = rbind(cbind(II1, matrix(0, k1 - 1, dim(II2)[2])), 
            cbind(matrix(0, k2 - 1, dim(II1)[2]), II2))
        Dis = cbind(II, matrix(0, k1 + k2 - 2, dim(ZZ)[2] - (k1 * 
            rm1 + k2 * rm2)))
    }
    else Dis = NULL
    if (start == 0) {
        if(cov){
            de1 = de2 = NULL
            Piv = matrix(1/(k1 * k2), ns, k1 * k2)
            piv = NULL
        }else{
            be = de1 = de2 = NULL
            piv1 = rep(1, k1)/k1
            piv2 = rep(1, k2)/k2
        }
        if (k1 == 1) 
            grid1 = 0
        else grid1 = seq(-k1, k1, 2 * k1/(k1 - 1))
        if (k2 == 1) 
            grid2 = 0
        else grid2 = seq(-k2, k2, 2 * k2/(k2 - 1))
        Phi = array(NA, c(lm, J, k1 * k2))
        for (j in 1:J) {
            dist = rep(0, lv[j])
            for (y in 0:(lv[j] - 1)) dist[y + 1] = (sum(yv[S[,j] == y]) + 0.5)/n
            out = matr_glob(lv[j])
            Co = out$Co
            Ma = out$Ma
            eta = Co %*% log(Ma %*% dist)
            count = 0
            for (c1 in 1:k1) for (c2 in 1:k2) {
                count = count + 1
                Phi[1:lv[j], j, count] = inv_glob(eta + grid1[c1] + 
                  grid2[c2])$p
            }
        }
    }
    if (start == 1) {
        if (cov) {
            if (glob){
                de1 = de2 = NULL
                Piv = matrix(runif(ns * k1 * k2), ns, k1 * k2)
                Piv = Piv * (1/rowSums(Piv))
                piv = NULL
            }else {
                de1 = rnorm((k1 - 1) * (ncov + 1))/rep(c(1, apply(X, 
                  2, sd)), (k1 - 1))
                de2 = rnorm((k2 - 1) * (ncov + 1))/rep(c(1, apply(X, 
                  2, sd)), (k2 - 1))
                if (k1 > 1) 
                  Piv1 = prob_multi_glob(XX1dis, logit_cov, de1, 
                    Xlabel)$P
                if (k2 > 1) 
                  Piv2 = prob_multi_glob(XX2dis, logit_cov, de2, 
                    Xlabel)$P
                if (k1 == 1) 
                  Piv = Piv2
                if (k2 == 1) 
                  Piv = Piv1
                if (k1 > 1 & k2 > 1) {
                  Piv = matrix(0, ns, k1 * k2)
                  for (i in 1:ns) Piv[i, ] = Piv1[i, ] %x% Piv2[i, 
                    ]
                }
                piv = NULL
            }
        } else {
        	de1 = de2 = NULL
            piv1 = runif(k1)
            piv1 = piv1/sum(piv1)
            piv2 = runif(k2)
            piv2 = piv2/sum(piv2)
        }
        Phi = array(runif(lm * J * k1 * k2), c(lm, J, k1 * k2))
        for (c in 1:(k1 * k2)) for (j in 1:J) {
            Phi[1:lv[j], j, c] = Phi[1:lv[j], j, c]/sum(Phi[1:lv[j], 
                j, c])
            if (lv[j] < lm) 
                Phi[(lv[j] + 1):lm, j, c] = NA
        }
        if (glob) {
            if (runif(1) > 0.5) 
                for (j in 1:J) {
                  mPhi = (0:(lv[j] - 1)) %*% Phi[1:lv[j], j, 
                    ]
                  ind = order(mPhi)
                  Phi[, j, ] = Phi[, j, ind]
                }
        }
    }
    if (start == 2){ 
        de1 = as.vector(De1)
        de2 = as.vector(De2)
    }
    if (start == 0 || start == 1){
        ga1 = rep(1, J1 - rm1)
        ga2 = rep(1, J2 - rm2)
    }else {
        ga1 = ga1c[indga1]
		ga2 = ga2c[indga2]
    }
    Psi = matrix(1, ns, k1 * k2)
    if (miss) {
        for (j in 1:J) for (c in 1:(k1 * k2)) Psi[, c] = Psi[, 
            c] * (Phi[S[, j] + 1, j, c] * R[, j] + (1 - R[, j]))
    }
    else {
        for (j in 1:J) for (c in 1:(k1 * k2)) Psi[, c] = Psi[, 
            c] * Phi[S[, j] + 1, j, c]
    }
    if(cov){
        if (start == 2) {
	        if (k1 > 1) Piv1 = prob_multi_glob(XX1dis, logit_cov, de1, Xlabel)$P
   			if (k2 > 1) Piv2 = prob_multi_glob(XX2dis, logit_cov, de2, Xlabel)$P
        	if (k1 == 1) Piv = Piv2
       		if (k2 == 1) Piv = Piv1
        	if (k1 > 1 & k2 > 1) {
        		Piv = matrix(0, ns, k1 * k2)
           		for (i in 1:ns) Piv[i, ] = Piv1[i, ] %x% Piv2[i,]
        	}
        }
    }
    else{
        Piv = rep(1, ns) %o% (piv1 %x% piv2)
    }
    if (k1 * k2 == 1) Pj = Psi
    else Pj = Psi * Piv
    pm = rowSums(Pj)
    lk = sum(yv * log(pm))
    cat("*-------------------------------------------------------------------------------*\n")
    if (link == 1 || link == 2) {
        cat(c("Model with multidimensional structure\n\n"))
        cat("* First latent variable\n")
        names11 = NULL
        for (j in 1:rm1) {
            names11 = c(names11, paste("Dimension", j))
        }
        multi1_out = as.matrix(multi1)
        if (rm1 == 1) 
            multi1_out = t(multi1_out)
        rownames(multi1_out) = names11
        print(multi1_out)
        cat("\n * Second latent variable\n")
        names12 = NULL
        for (j in 1:rm2) {
            names12 = c(names12, paste("Dimension", j))
        }
        multi2_out = as.matrix(multi2)
        if (rm2 == 1) 
            multi2_out = t(multi2_out)
        rownames(multi2_out) = names12
        print(multi2_out)
    }
    cat("\n")
    cat(c("Link of type =                 ", link0, "\n"))
    cat(c("Discrimination index =         ", disc, "\n"))
    cat(c("Constraints on the difficulty =", difl, "\n"))
    cat(c("Type of initialization =       ", start0, "\n"))
    cat("*-------------------------------------------------------------------------------*\n")
    if (disp) {
        if (!disc || (length(ga1) == 0 & length(ga2) == 0)) {
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
            cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |     dis     |   min(par)  |   max(par)  |\n")
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
        }
        if (disc) {
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
            cat("  iteration |   classes   |    model    |    lk       |    lk-lko   |      dis    |   min(ga)   |   max(ga)   |   min(par)  |   max(par)  |\n")
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
        }
        cat(sprintf("%11g", c(0, k1+ k2/10, link, lk)), "\n", sep = " | ")
    }
    it = 0
    lko = lk - 10^10
    dis = 0
    par = NULL
    dga = NULL
    lkv = NULL
    while (((abs(lk - lko)/abs(lko) > tol) && it < 10^4) || it < 
        2) {
        it = it + 1
        paro = par
        ga1o = ga1
        ga2o = ga2
        piv1o = piv1
        piv2o = piv2
        de1o = de1
        de2o = de2
        lko = lk
        V = ((yv/pm) %o% rep(1, k1 * k2)) * Piv * Psi
        sV = colSums(V)
        YY = matrix(NA, J * k1 * k2, lm)
        count = 0
        for (c in 1:(k1 * k2)) for (j in 1:J) {
            count = count + 1
            for (y in 1:lv[j]) {
                ind = (S[, j] == (y - 1))
                if (miss) 
                  YY[count, y] = sum(V[ind, c] * R[ind, j])
                else YY[count, y] = sum(V[ind, c])
            }
        }
        if (disc) {
            if (it > 1 & rm < J) {
                ZZ1 = ZZ2 = array(NA, c(lm - 1, J, J * k1 * k2))
                count = 0
                for (c1 in 1:k1) for (c2 in 1:k2) for (j in 1:J) {
                  count = count + 1
                  ZZ1[1:(lv[j] - 1), , count] = 0
                  ZZ1[1:(lv[j] - 1), j, count] = ZZ0[1:(lv[j] - 
                    1), 1:(k1 * rm1), count] %*% par[1:(k1 * 
                    rm1)]
                  ZZ2[1:(lv[j] - 1), , count] = 0
                  ZZ2[1:(lv[j] - 1), j, count] = ZZ0[1:(lv[j] - 
                    1), k1 * rm1 + 1:(k2 * rm2), count] %*% par[k1 * 
                    rm1 + 1:(k2 * rm2)]
                }
                ZZ1 = array(ZZ1[, indga1, ], c(lm - 1, length(ga1), 
                  J * k1 * k2))
                ZZ2 = array(ZZ2[, indga2, ], c(lm - 1, length(ga2), 
                  J * k1 * k2))
                if (!difl) ind = k1 * rm1 + k2 * rm2 + 1:(sum(lv - 1) - rm)
                if (difl) ind = k1 * rm1 + k2 * rm2 + 1:(J - rm + sum(lv[fv] - 2))
                ZZ1Int = ZZ2Int = array(NA, c(lm - 1, J * k1 * 
                  k2))
                count = 0
                for (c1 in 1:k1) for (c2 in 1:k2) for (j in 1:J) {
                  count = count + 1
                  ZZ1Int[1:(lv[j] - 1), count] = ga2c[j] * ZZ0[1:(lv[j] - 
                    1), k1 * rm1 + 1:(k2 * rm2), count] %*% par[k1 * 
                    rm1 + 1:(k2 * rm2)] + ZZ0[1:(lv[j] - 1), 
                    ind, count] %*% par[ind]
                  ZZ2Int[1:(lv[j] - 1), count] = ga1c[j] * ZZ0[1:(lv[j] - 
                    1), 1:(k1 * rm1), count] %*% par[1:(k1 * 
                    rm1)] + ZZ0[1:(lv[j] - 1), ind, count] %*% 
                    par[ind]
                }
                ga1 = est_multi_glob_gen(YY, ZZ1, ltype, be = ga1, Int = ZZ1Int)$be
                ga2 = est_multi_glob_gen(YY, ZZ2, ltype, be = ga2, Int = ZZ2Int)$be
            }
            ga1c = rep(1, J)
            ga1c[indga1] = ga1
            ga2c = rep(1, J)
            ga2c[indga2] = ga2
            ZZ = ZZ0
            for (j in 1:J) {
                ind = (refitem == j)
                ZZ[, 1:(k1 * rm1), ind] = ZZ[, 1:(k1 * rm1), 
                  ind] * ga1c[j]
                ZZ[, k1 * rm1 + 1:(k2 * rm2), ind] = ZZ[, k1 * 
                  rm1 + 1:(k2 * rm2), ind] * ga2c[j]
            }
        }
        if(start==2 & it==1) maxit = 250 else maxit = 10
        out = est_multi_glob_gen(YY, ZZ, ltype, be = par, Dis = Dis, maxit=maxit)
        par = out$be
        P = out$P
        Phi = array(t(P), c(lm, J, k1 * k2))
        if (cov) {
            if (k1 > 1) {
                out = est_multi_glob(V %*% t(Aggr1), XX1dis, 
                  logit_cov, Xlabel, de1)
                de1 = out$be
                P1dis = out$Pdis
                Piv1 = out$P
            }
            if (k2 > 1) {
                out = est_multi_glob(V %*% t(Aggr2), XX2dis, 
                  logit_cov, Xlabel, de2)
                de2 = out$be
                P2dis = out$Pdis
                Piv2 = out$P
            }
            if (k1 == 1) 
                Piv = Piv2
            if (k2 == 1) 
                Piv = Piv1
            if (k1 > 1 & k2 > 1) 
                for (i in 1:ns) Piv[i, ] = Piv1[i, ] %x% Piv2[i,]
        }
        else {
            piv1 = as.vector(Aggr1 %*% sV)/n
            piv2 = as.vector(Aggr2 %*% sV)/n
            Piv = rep(1, ns) %o% (piv1 %x% piv2)
        }
        Psi = matrix(1, ns, k1 * k2)
        if (miss) {
            for (j in 1:J) for (c in 1:(k1 * k2)) Psi[, c] = Psi[, 
                c] * (Phi[S[, j] + 1, j, c] * R[, j] + (1 - R[, 
                j]))
        }
        else {
            for (j in 1:J) for (c in 1:(k1 * k2)) Psi[, c] = Psi[, 
                c] * Phi[S[, j] + 1, j, c]
        }
        if (k1 * k2 == 1) 
            Pj = Psi
        else Pj = Psi * Piv
        pm = rowSums(Pj)
        lk = sum(yv * log(pm))
        if (it > 1 & link > 0) 
            dis = max(c(abs(par - paro), abs(ga1 - ga1o), abs(ga2 - 
                ga2o), abs(piv1 - piv1o), abs(piv2 - piv2o)))
        if (disp) {
            if (it/10 == floor(it/10)) {
                if (!disc || (length(ga1) == 0 & length(ga2) == 0)) 
                  cat(sprintf("%11g", c(it, k1 + k2/10, link, 
                    lk, lk - lko, dis, min(par), max(par))), "\n", sep = " | ")
                if (disc) cat(sprintf("%11g", c(it, k1 + k2/10, link, 
                    lk, lk - lko, dis, min(c(ga1, ga2)), max(c(ga1, 
                      ga2)), min(par), max(par))), "\n", sep = " | ")
            }
        }
        lkv = c(lkv, lk)
    }
    if (disp) {
        if (it/10 > floor(it/10)) {
            if (!disc || (length(ga1) == 0 & length(ga2) == 0)) 
                cat(sprintf("%11g", c(it, k1+ k2/10, link, 
                  lk, lk - lko, dis, min(par), max(par))), "\n", 
                  sep = " | ")
            if (disc) cat(sprintf("%11g", c(it, k1+k2/10, link, 
                  lk, lk - lko, dis, min(c(ga1, ga2)), max(c(ga1, 
                    ga2)), min(par), max(par))), "\n", sep = " | ")
        }
        if (!disc || (length(ga1) == 0 & length(ga2) == 0)) 
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
        if (disc) cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
    }
    np = length(par)
    if (disc) np=np+length(ga1)+length(ga2)
    if (cov) {
        if (glob) {
            np = np + k1 + k2 - 2 + 2 * ncov
        }
        else {
            np = np + (k1 + k2 - 2) * (ncov + 1)
        }
    }
    else {
        np = np + k1 - 1 + k2 - 1
    }
    aic = -2 * lk + 2 * np
    bic = -2 * lk + np * log(n)
    th1 = par[indth1]
    th2 = par[indth2]
    be = par[indbe]
    if (!difl) {
        bec = rep(0, sum(lv - 1))
        bec[indbec] = be
        Bec = matrix(NA, J, lm - 1)
        count = 0
        for (j in 1:J) {
            Bec[j, (1:lv[j] - 1)] = bec[count + (1:(lv[j] - 1))]
            count = count + lv[j] - 1
        }
        dimnames(Bec) = list(item=1:J,cutoff=1:(lm-1))
    }
    else if (difl) {
        bec1 = rep(0, J)
        bec1[indbec] = be[1:(J - rm)]
        if (rm == 1) {
            bec2 = rep(0, lm - 1)
            bec2[2:(lm - 1)] = be[J - rm + (1:(lm - 2))]
        }
        else {
            bec2 = matrix(NA, lm - 1, rm)
            bec2[1, ] = 0
            count = 0
            for (h in 1:rm) if (lv[fv[h]] > 2) {
                bec2[2:(lv[h] - 1), h] = be[J - rm + count + 
                  (1:(lv[h] - 2))]
                count = count + lv[h] - 2
            }
            dimnames(bec2) = list(level = 1:(lm - 1), dim = 1:rm)
        }
        Bec = list(difficulties = bec1, cutoffs = bec2)
    }
    ga1c = rep(1, J)
    ga1c[indga1] = ga1
    ga1c[setdiff(1:J,multi1v)] = NA
    ga2c = rep(1, J)
    ga2c[indga2] = ga2
    ga2c[setdiff(1:J,multi2v)] = NA
    Th1 = matrix(th1, rm1, k1)
    Th2 = matrix(th2, rm2, k2)
    dimnames(Th1) = list(dimension=1:rm1,class=1:k1)
    dimnames(Th2) = list(dimension=1:rm2,class=1:k2)
    Pp = ((1/pm) %o% rep(1, k1 * k2)) * Piv * Psi
    Pp1 = Pp %*% t(Aggr1)
    Pp2 = Pp %*% t(Aggr2)
    ent = -sum(V * log(pmax(Pp, 10^-100)))
    if (cov) {
        if (glob) {
            if (k1 == 1) {
                De1 = NULL
            }
            else {
                De1 = matrix(de1, ncov + k1 - 1, 1)
                names_cutoff = paste("cutoff", 1:(k1 - 1), sep = "")
                if (is.null(namesX)) {
                  namesX1 = c(names_cutoff, paste("X", 1:ncov, 
                    sep = ""))
                }
                else {
                  namesX1 = c(names_cutoff, namesX)
                }
                rownames(De1) = namesX1
            }
            if (k2 == 1) {
                De2 = NULL
            }
            else {
                De2 = matrix(de2, ncov + k2 - 1, 1)
                names_cutoff = paste("cutoff", 1:(k2 - 1), sep = "")
                if (is.null(namesX)) {
                  namesX2 = c(names_cutoff, paste("X", 1:ncov, 
                    sep = ""))
                }
                else {
                  namesX2 = c(names_cutoff, namesX)
                }
                rownames(De2) = namesX2
            }
        }
        else {
            if (is.null(namesX)) {
                namesX = c("intercept", paste("X", 1:ncov, sep = ""))
            }
            else {
                namesX = c("intercept", namesX)
            }
            if (k1 == 1) 
                De1 = NULL
            else {
                De1 = matrix(de1, ncov + 1, k1 - 1)
                dimnames(De1) = list(namesX, logit = 2:k1)
            }
            if (k2 == 1) 
                De2 = NULL
            else {
                De2 = matrix(de2, ncov + 1, k2 - 1)
                dimnames(De2) = list(namesX, logit = 2:k2)
            }
        }
        piv1 = colMeans(Piv1)
        piv2 = colMeans(Piv2)
    }
    else {
        De1 = De2 = NULL
    }
    dimnames(Phi) = list(category = 0:(lm - 1), item = 1:J, class = 1:(k1*k2))
    if (!cov) {
        de1 = De1 = log(piv1[-1]/piv1[1])
        de2 = De2 = log(piv2[-1]/piv2[1])
    }
    if (out_se) {
        lde1 = length(de1)
        lde2 = length(de2)
        lpar = length(par)
        lga = 0
        par_comp = c(de1, de2, par)
        if (disc) {
            lga1 = length(ga1)
            lga2 = length(ga2)
            par_comp = c(par_comp, ga1, ga2)
        }
        if (disp) {
            cat("computation of derivatives\n")
            cat(length(par_comp), "parameters\n")
        }
        out = lk_obs_score_within(par_comp, lde1, lde2, lpar, 
            lga1, lga2, S, R, yv, k1, k2, rm1, rm2, lv, J, fv, 
            link, disc, indga1, indga2, glob, refitem, miss, 
            ltype, XX1dis, XX2dis, Xlabel, ZZ0, fort)
        scn = rep(0, length(par_comp))
        Jn = NULL
        for (j in 1:length(par_comp)) {
            par_comp1 = par_comp
            par_comp1[j] = par_comp1[j] + 10^-6
            out1 = lk_obs_score_within(par_comp1, lde1, lde2, 
                lpar, lga1, lga2, S, R, yv, k1, k2, rm1, rm2, 
                lv, J, fv, link, disc, indga1, indga2, glob, 
                refitem, miss, ltype, XX1dis, XX2dis, Xlabel, 
                ZZ0, fort)
            scn[j] = (out1$lk - lk) * 10^6
            Jn = cbind(Jn, (out1$sc - out$sc) * 10^6)
            if (disp){
            	if (j/10 > floor(j/10)) cat(".") else cat(round(j/10))
                if (j/100 == floor(j/100)) cat("\n")
            }
        }
        if (disp) 
            cat("\n")
        Jn = -(Jn + t(Jn))/2
        Vn = ginv(Jn)
        se = sqrt(abs(diag(Vn)))
        if (k1 > 1) 
            sede1 = se[1:lde1]
        if (k2 > 1) 
            sede2 = se[lde1 + 1:lde2]
        	separ = se[lde1 + lde2 + (1:lpar)]
        if (disc) {
            sega1 = se[lpar + lde1 + lde2 + (1:lga1)]
            sega2 = se[lpar + lde1 + lde2 + lga1 + (1:lga2)]
    		sega1c = rep(0, J)
    		sega1c[indga1] = sega1
			sega1c[setdiff(1:J,multi1v)] = NA
			sega2c = rep(0, J)
			sega2c[indga2] = sega2
    		sega2c[setdiff(1:J,multi2v)] = NA   
        }else{
            sega1c = NULL
            sega2c = NULL
        }
        seth1 = separ[indth1]
 		seth2 = separ[indth2]
		seTh1 = matrix(seth1, rm1, k1)
    	seTh2 = matrix(seth2, rm2, k2)
		sebe = separ[indbe]
		if (!difl) {
        	sebec = rep(0, sum(lv - 1))
        	sebec[indbec] = sebe
			seBec = matrix(NA, J, lm - 1)
        	count = 0
        	for (j in 1:J) {
            	seBec[j, (1:lv[j] - 1)] = sebec[count + (1:(lv[j] - 1))]
				count = count + lv[j] - 1
        	}
        	dimnames(seBec) = list(item=1:J,cutoff=1:(lm-1))
    	}else if (difl){
			sebec1 = rep(0, J)
        	sebec1[indbec] = sebe[1:(J - rm)]
        	if (rm == 1) {
            	sebec2 = rep(0, lm - 1)
            	sebec2[2:(lm - 1)] = sebe[J - rm + (1:(lm - 2))]
        	}else {
            	sebec2 = matrix(NA, lm - 1, rm)
            	sebec2[1, ] = 0
            	count = 0
            	for (h in 1:rm) if (lv[fv[h]] > 2) {
                	sebec2[2:(lv[h] - 1), h] = sebe[J - rm + count + (1:(lv[h] - 2))]
                	count = count + lv[h] - 2
				}
            	dimnames(bec2) = list(level = 1:(lm - 1), dim = 1:rm)
        	}
        	seBec = list(difficulties = sebec1, cutoffs = sebec2)
    	}
    	dimnames(seTh1) = list(dimension=1:rm1,class=1:k1)
    	dimnames(seTh2) = list(dimension=1:rm2,class=1:k2)
        if (glob) {
            if (k1 == 1) 
                seDe1 = NULL
            else seDe1 = matrix(sede1, ncov + k1 - 1, 1)
            if (k2 == 1) 
                seDe2 = NULL
            else seDe2 = matrix(sede2, ncov + k2 - 1, 1)
        }
        else {
            if (k1 == 1) 
                seDe1 = NULL
            else seDe1 = matrix(sede1, ncov + 1, k1 - 1)
            if (k2 == 1) 
                seDe2 = NULL
            else seDe2 = matrix(sede2, ncov + 1, k2 - 1)
        }
    }
    out = list(piv1 = piv1, piv2 = piv2, Th1 = Th1, Th2 = Th2, 
        Bec = Bec, ga1c = ga1c, ga2c = ga2c, fv1 = fv1, fv2 = fv2, 
        De1 = De1, De2 = De2, Phi = Phi, lk = lk, np = np, aic = aic, 
        bic = bic, ent = ent, call = match.call())
    if (output) {
        out$Pp = Pp
        out$Pp = Pp1
        out$Pp = Pp2
        out$lkv = lkv
        if(cov){
        	out$Xlabel = Xlabel
        	out$XX1dis = XX1dis
        	out$XX2dis = XX2dis
	        out$Piv1 = Piv1
    	    out$Piv2 = Piv2
        }
    }
    if (out_se) {
        out$seDe1 = seDe1
        out$seDe2 = seDe2
        out$seTh1 = seTh1
        out$seTh2 = seTh2
        out$seBec = seBec
        out$sega1c = sega1c
        out$sega2c = sega2c
        out$Vn = Vn
    }
    class(out) = "est_multi_poly_within"
    return(out)
}
