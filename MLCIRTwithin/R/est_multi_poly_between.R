est_multi_poly_between <- function (S, yv = rep(1, ns), k, X = NULL, start = c("deterministic","random","external"),
    link = c("global","local"), disc = FALSE, difl = FALSE, multi = 1:J, piv = NULL, Phi = NULL, gac = NULL, De = NULL,
    fort = FALSE, tol = 10^-10, disp = FALSE, output = FALSE, out_se = FALSE, glob = FALSE){
    	
    if (k == 1) stop("--> use est_multi_poly")
    if (max(S, na.rm = TRUE) == 1 & difl) {
        warning("with binary data put difl=FALSE\n")
        difl = FALSE
    }
    link0 = match.arg(link)
	if(link0=="global") link=1
	if(link0=="local") link=2   
	start0 = match.arg(start)
	if(start0=="deterministic") start=0
	if(start0=="random") start=1
	if(start0=="external") start=2 
    cov = !is.null(X)
    if (cov) {
        X = as.matrix(X)
        namesX = colnames(X)
        if (glob) 
            logit_cov = "g"
        else logit_cov = "m"
    }
    else {
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
            XXdis = array(0, c(k - 1, k - 1 + ncov, Xndis))
            for (i in 1:Xndis) XXdis[, , i] = cbind(diag(k - 
                1), rep(1, k - 1) %o% Xdis[i, ])
        }
        else {
            XXdis = array(0, c(k - 1, (k - 1) * (ncov + 1),Xndis))
            if (k == 2) II = 1
            else II = diag(k - 1)
            for (i in 1:Xndis) XXdis[, , i] = II %x% t(c(1,Xdis[i, ]))
        }
    }
    else {
        ncov = 0
        XXdis = array(diag(k - 1), c(k - 1, k - 1, 1))
        Xlabel = rep(1, ns)
    }
    if (link == 1) ltype = "g"
    else if (link == 2) ltype = "l"
    if (link == 1 || link == 2) {
        items = sort(unique(as.vector(multi)))
        if (any(items == 0)) items = items[items > 0]
        Jitems = length(items)
        if (is.vector(multi)) rm = 1
        else rm = nrow(multi)
        Dem = matrix(0, J, rm)
        # index of the discriminant parameters that are constrained
        if (rm == 1) {
            Dem = 1
            fv = multi[1]
        }
        else {
            for (r in 1:rm) {
                ind = multi[r, ]
                ind = ind[ind > 0]
                #Dem1[ind, r] = 1
                Dem[ind, r] = 1
            }
            fv = multi[, 1]
        }
        # index of the difficulty parameters that are constrained
        fve = NULL
        count = 0
        for (j in 1:J) {
            if (j %in% fv) fve = c(fve, count + 1)
            count = count + lv[j] - 1
        }
        # index of free discriminant prameters
        indga = 1:J
        indga = indga[setdiff(items, fv)]
        # index of free ability parameters
        indth = 1:(k * rm)
        # index of free difficulty parameters
        if (!difl) {
            indbe = k * rm + (1:(sum(lv - 1) - rm))
            indbec = 1:sum(lv - 1)
            indbec = indbec[-fve]
        }
        else {
            indbe = k * rm + (1:(J - rm + sum(lv[fv] - 2)))
            indbec = 1:J
            indbec = indbec[-fv]
        }
        abils = rep(0, J)
        if (rm == 1) {
            abils[multi] = 1
        }
        else {
            for (h in 1:rm) {
                ind = multi[h, ]
                ind = ind[ind > 0]
                abils[ind] = h
            }
        }
        if (!difl) ZZ = array(NA, c(lm - 1, k * rm + sum(lv - 1) - rm, J * k))
        if (difl) ZZ = array(NA, c(lm - 1, k * rm + J - rm + sum(lv[fv] - 2), J * k))
        cont = 0
        refitem = matrix(0, J * k, 1)
		for (c in 1:k){
			u1 = matrix(0, 1, k)
			u1[c] = 1
            for(j in 1:J){
				u2 = matrix(0, 1, rm)
				u2[abils[j]] = 1
				v = matrix(0, 1, J)
				v[j] = 1
				cont = cont + 1
				if (!difl){
					Te = matrix(0, lv[j] - 1, sum(lv - 1))
					if (j == 1) Te[, 1:(lv[j] - 1)] = diag(lv[j] - 1)
					else Te[, sum(lv[1:(j - 1)] - 1) + (1:(lv[j] - 1))] = diag(lv[j] - 1)
					Te = matrix(Te[, -fve], lv[j] - 1, dim(Te)[2] - length(fve))
				}
                else if (difl) {
                		Te = matrix(0, lv[j] - 1, sum(lv[fv] - 2))
                		if (lv[j] > 2) {
                			if (abils[j] == 1) Te[, 1:(lv[j] - 2)] = diag(lv[j] - 1)[, -1]
                			else Te[, sum(lv[1:(abils[j] - 1)] - 2) + (1:(lv[j] - 2))] = diag(lv[j] - 1)[,-1]
                		}
                		Te = cbind(v %x% rep(1, lv[j] - 1), Te)
                		Te = matrix(Te[, -fv], lv[j] - 1)
                	}
                	ZZ[1:(lv[j] - 1), , cont] = cbind(rep(1, lv[j] - 1) %*% (u1 %x% u2), -Te)
                	refitem[cont] = j
            }
        }
    }
    ZZ0 = ZZ
    if (glob) {
        II = diag(k - 1)
        II = cbind(0, II) - cbind(II, 0)
        if (rm > 1) II = II %x% matrix(c(1, rep(0, rm - 1)), 1, rm)
        Dis = cbind(II, matrix(0, k - 1, dim(ZZ)[2] - (k * rm)))
    }
    else Dis = NULL
    if (start == 0) {
        if (cov) {
            de = NULL
            Piv = matrix(1/k, ns, k)
            piv = NULL
        }
        else {
            be = de = NULL
            piv = rep(1, k)/k
        }
        if (k == 1) grid = 0
        else grid = seq(-k, k, length.out = k)
        Phi = array(NA, c(lm, J, k))
        for (j in 1:J) {
            dist = rep(0, lv[j])
            for (y in 0:(lv[j] - 1)) dist[y + 1] = (sum(yv[S[, 
                j] == y]) + 0.5)/n
            out = matr_glob(lv[j])
            Co = out$Co
            Ma = out$Ma
            eta = Co %*% log(Ma %*% dist)
            count = 0
            for (c in 1:k){
                count = count + 1
                Phi[1:lv[j], j, count] = inv_glob(eta + grid[c])$p
            }
        }
    }
    if (start == 1) {
        if (cov) {
            if (glob){
                de = NULL
                Piv = matrix(runif(ns * k), ns, k)
                Piv = Piv * (1/rowSums(Piv))
                piv = NULL
            }else {
                de = rnorm((k - 1) * (ncov + 1))/rep(c(1, apply(X, 2, sd)), (k - 1))
                if (k > 1) Piv = prob_multi_glob(XXdis, logit_cov, de, Xlabel)$P
                piv = NULL
            }
        }
        else {
            piv = runif(k)
            piv = piv/sum(piv)
        }
        Phi = array(runif(lm * J * k), c(lm, J, k))
        for (c in 1:k) for (j in 1:J) {
            Phi[1:lv[j], j, c] = Phi[1:lv[j], j, c]/sum(Phi[1:lv[j], j, c])
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
        de = as.vector(De)
    }
    if (start == 0 || start == 1){
        ga = rep(1, Jitems - rm)
    }else {
        ga = gac[indga]
    }
    Psi = matrix(1, ns, k)
    if (miss) {
        for (j in 1:J) for (c in 1:k) Psi[, c] = Psi[, c] * (Phi[S[, j] + 1, j, c] * R[, j] + (1 - R[, j]))
    }
    else {
        for (j in 1:J) for (c in 1:k) Psi[, c] = Psi[,c] * Phi[S[, j] + 1, j, c]
    }
    if(cov){
        if (start == 2) {
	        if (k > 1) Piv = prob_multi_glob(XXdis, logit_cov, de, Xlabel)$P
        }
    }
    else{
        Piv = rep(1, ns) %o% piv
    }
    if (k == 1) Pj = Psi
    else Pj = Psi * Piv
    pm = rowSums(Pj)
    lk = sum(yv * log(pm))
    cat("*-------------------------------------------------------------------------------*\n")
    if (link == 1 || link == 2) {
        cat(c("Model with multidimensional structure\n"))
        names1 = NULL
        for (j in 1:rm) {
            names1 = c(names1, paste("Dimension", j))
        }
        multi_out = as.matrix(multi)
        if (rm == 1) 
            multi_out = t(multi_out)
        rownames(multi_out) = names1
        print(multi_out)
    }
    cat(c("Link of type =                 ", link0, "\n"))
    cat(c("Discrimination index =         ", disc, "\n"))
    cat(c("Constraints on the difficulty =", difl, "\n"))
    cat(c("Type of initialization =       ", start0, "\n"))
    cat("*-------------------------------------------------------------------------------*\n")
    if (disp) {
        if (!disc || length(ga) == 0 ) {
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
            cat("  iteration |   classes   |    model    |      lk     |    lk-lko   |     dis     |   min(par)  |   max(par)  |\n")
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
        }
        if (disc) {
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
            cat("  iteration |   classes   |    model    |    lk       |    lk-lko   |      dis    |   min(ga)   |   max(ga)   |   min(par)  |   max(par)  |\n")
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
        }
        cat(sprintf("%11g", c(0, k, link, lk)), "\n", sep = " | ")
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
        gao = ga
        pivo = piv
        deo = de
        lko = lk
        V = ((yv/pm) %o% rep(1, k)) * Piv * Psi
        sV = colSums(V)
        YY = matrix(NA, J * k, lm)
        count = 0
        for (c in 1:k) for (j in 1:J) {
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
                ZZ = array(NA, c(lm - 1, J, J * k))
                count = 0
                for(c in 1:k) for(j in 1:J) {
                  count = count + 1
                  ZZ[1:(lv[j] - 1), , count] = 0
                  ZZ[1:(lv[j] - 1), j, count] = ZZ0[1:(lv[j] - 1), 1:(k * rm), count] %*% par[1:(k*rm)]
                }
                ZZ = array(ZZ[, indga, ], c(lm - 1, length(ga), J * k))
                if (!difl) ind = k * rm + 1:(sum(lv - 1) - rm)
                if (difl) ind = k * rm + 1:(J - rm + sum(lv[fv] - 2))
                ZZInt = array(NA, c(lm - 1, J * k))
                count = 0
                for (c in 1:k) for (j in 1:J) {
                  count = count + 1
                  ZZInt[1:(lv[j] - 1), count] = ZZ0[1:(lv[j] - 1), ind, count] %*% par[ind]
                  ZZ0[1:(lv[j] - 1), ind, count] %*% par[ind]
                }
                ga = est_multi_glob_gen(YY, ZZ, ltype, be = ga, Int = ZZInt)$be
            }
            gac = rep(1, J)
            gac[indga] = ga
            ZZ = ZZ0
            for (j in 1:J) {
                ind = (refitem == j)
                ZZ[, 1:(k * rm), ind] = ZZ[, 1:(k * rm), ind] * gac[j]
            }
        }
        if(start==2 & it==1) maxit = 250 else maxit = 10
        out = est_multi_glob_gen(YY, ZZ, ltype, be = par, Dis = Dis, maxit=maxit)
        par = out$be
        P = out$P
        Phi = array(t(P), c(lm, J, k))
        if (cov) {
            if (k > 1) {
                out = est_multi_glob(V, XXdis, logit_cov, Xlabel, de)
                de = out$be
                Pdis = out$Pdis
                Piv = out$P
            }
        }
        else {
            piv = sV/n
            Piv = rep(1, ns) %o% piv
        }
        Psi = matrix(1, ns, k)
        if (miss) {
            for (j in 1:J) for (c in 1:k) Psi[, c] = Psi[,c] * (Phi[S[, j] + 1, j, c] * R[, j] + (1 - R[,j]))
        }
        else {
            for (j in 1:J) for (c in 1:k) Psi[, c] = Psi[,c] * Phi[S[, j] + 1, j, c]
        }
        if (k == 1) Pj = Psi
        else Pj = Psi * Piv
        pm = rowSums(Pj)
        lk = sum(yv * log(pm))
        if (it > 1 & link > 0) dis = max(c(abs(par - paro), abs(ga - ga), abs(piv - pivo)))
        if (disp) {
            if (it/10 == floor(it/10)) {
                if (!disc || length(ga) == 0 ) 
                  cat(sprintf("%11g", c(it, k, link, lk, lk - lko, dis, min(par), max(par))), "\n", sep = " | ")
                if (disc) 
                  cat(sprintf("%11g", c(it, k, link, lk, lk - lko, dis, min(ga), max(ga), min(par), max(par))), "\n", sep = " | ")
            }
        }
        lkv = c(lkv, lk)
    }
    if (disp) {
        if (it/10 > floor(it/10)) {
            if (!disc || length(ga) == 0) 
                cat(sprintf("%11g", c(it, k, link, lk, lk - lko, dis, min(par), max(par))), "\n", sep = " | ")
            if (disc) 
                cat(sprintf("%11g", c(it, k, link, lk, lk - lko, dis, min(ga), max(ga), min(par), max(par))), "\n", sep = " | ")
        }
        if (!disc || length(ga) == 0) 
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
        if (disc) 
            cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|-------------|\n")
    }
    np = length(par)
    if (disc==1) np=np+length(ga)
    if (cov) {
        if (glob) {
            np = np + k - 1 + ncov
        }
        else {
            np = np + (k-1) * (ncov + 1)
        }
    }
    else {
        np = np + k - 1
    }
    aic = -2 * lk + 2 * np
    bic = -2 * lk + np * log(n)
    th = par[indth]
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
        dimnames(Bec) = list(item = 1:J, cutoff = 1:(lm-1))
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
    Th = matrix(th, rm, k)
    dimnames(Th) = list(dimension = 1:rm,class=1:k)
    gac = rep(1, J)
    gac[indga] = ga
    Pp = ((1/pm) %o% rep(1, k)) * Piv * Psi
    ent = -sum(V * log(pmax(Pp, 10^-100)))
    if (cov) {
        if (glob) {
            if (k == 1) {
                De = NULL
            }
            else {
                De = matrix(de, ncov + k - 1, 1)
                names_cutoff = paste("cutoff", 1:(k - 1), sep = "")
                if (is.null(namesX)) {
                  namesX1 = c(names_cutoff, paste("X", 1:ncov, 
                    sep = ""))
                }
                else {
                  namesX1 = c(names_cutoff, namesX)
                }
                rownames(De) = namesX
            }
        }
        else {
            if (is.null(namesX)) {
                namesX = c("intercept", paste("X", 1:ncov, sep = ""))
            }
            else {
                namesX = c("intercept", namesX)
            }
            if (k == 1) 
                De = NULL
            else {
                De = matrix(de, ncov + 1, k - 1)
                dimnames(De) = list(namesX, logit = 2:k)
            }
        }
        piv = colMeans(Piv)
    }
    else {
        De = NULL
    }
    dimnames(Phi) = list(category = 0:(lm - 1), item = 1:J, class = 1:k)
    if (!cov) {
        de = De = log(piv[-1]/piv[1])
    }
    if (out_se) {
        lde = length(de)
        lpar = length(par)
        lga = 0
        par_comp = c(de, par)
        if (disc) {
            lga = length(ga)
            par_comp = c(par_comp, ga)
        }
        if (disp) {
            cat("computation of derivatives\n")
            cat(length(par_comp), "parameters\n")
        }
        out = lk_obs_score_between(par_comp, lde, lpar, 
            lga, S, R, yv, k, rm, lv, J, fv, 
            link, disc, indga, glob, refitem, miss, 
            ltype, XXdis, Xlabel, ZZ0, fort)
        scn = rep(0, length(par_comp))
        Jn = NULL
        for (j in 1:length(par_comp)) {
            par_comp1 = par_comp
            par_comp1[j] = par_comp1[j] + 10^-6
            out1 = lk_obs_score_between(par_comp1, lde,  
                lpar, lga, S, R, yv, k, rm,  
                lv, J, fv, link, disc, indga, glob, 
                refitem, miss, ltype, XXdis, Xlabel, ZZ0, fort)
            scn[j] = (out1$lk - lk) * 10^6
            Jn = cbind(Jn, (out1$sc - out$sc) * 10^6)
            if (disp) 
                if (j/10 > floor(j/10)) 
                  cat(".")
                else cat(round(j/10))
            if (disp) 
                if (j/100 == floor(j/100)) 
                  cat("\n")
        }
        if (disp) 
            cat("\n")
        Jn = -(Jn + t(Jn))/2
        Vn = ginv(Jn)
        se = sqrt(abs(diag(Vn)))
        if (k > 1) sede = se[1:lde]
       	separ = se[lde + (1:lpar)]
	    seth = separ[indth]
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
        	dimnames(seBec) = list(item = 1:J, cutoff = 1:(lm-1))
    	}else if (difl) {
        	sebec1 = rep(0, J)
        	sebec1[indbec] = sebe[1:(J - rm)]
			if (rm == 1) {
            	sebec2 = rep(0, lm - 1)
            	sebec2[2:(lm - 1)] = sebe[J - rm + (1:(lm - 2))]
        	}else{
            	sebec2 = matrix(NA, lm - 1, rm)
            	sebec2[1, ] = 0
            	count = 0
            	for (h in 1:rm) if (lv[fv[h]] > 2) {
                	sebec2[2:(lv[h] - 1), h] = sebe[J - rm + count + (1:(lv[h] - 2))]
                	count = count + lv[h] - 2
            	}
            	dimnames(sebec2) = list(level = 1:(lm - 1), dim = 1:rm)
        	}
        	seBec = list(difficulties = sebec1, cutoffs = sebec2)
    	}
    	seTh = matrix(seth, rm, k)
    	dimnames(seTh) = list(dimension = 1:rm,class=1:k)       	
        if (disc) sega = se[lpar + lde + (1:lga)] else sega = NULL
        if (glob){
            if (k == 1) seDe = NULL
            else seDe = matrix(sede, ncov + k - 1, 1)
        }else{
            if (k == 1) seDe = NULL
            else seDe = matrix(sede, ncov + 1, k - 1)
        }
    }
    out = list(piv = piv, Th = Th, Bec = Bec, gac = gac, fv = fv, 
        De = De, Phi = Phi, lk = lk, np = np, aic = aic, 
        bic = bic, ent = ent, call = match.call())
    if (output) {
        out$Piv = Piv
        out$Pp = Pp
        out$lkv = lkv
        if(cov){
          	out$Xlabel = Xlabel
          	out$XXdis = XXdis
        }
    }
    if (out_se) {
        out$seDe = seDe
        out$seBec = seBec
        out$seTh = seTh
        out$sega = sega
        out$Vn = Vn
    }
    class(out) = "est_multi_poly_between"
    return(out)
}
