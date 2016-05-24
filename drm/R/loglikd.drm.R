"loglikd.drm" <-
  function (p, y, X, dep, asscov = NULL, npar, nclass, nrep, tm, 
            w, tlen, lab, equal, inv, ass, link, constant, Ncond, Lclass, offset, 
            wt, Nf, misn, dropx, droplab, drop.cov) 
{
  assp <- matrix(p[(npar + 1):length(p)], nrow = 1)
  alab <- asscov[[4]]
  if (!is.null(equal)) {
    if (length(constant[[1]] > 0)) {
      assc <- matrix(c(equal, constant[[1]]), nrow = 1)
      assc[, constant[[1]]] <- constant[[2]]
      assc[, -constant[[1]]] <- assp[, equal, drop = FALSE]
      assp <- assc
    }
    else (assp <- assp[, equal, drop = FALSE])
  }
  nass <- nrow(assp)
  eta <- offset
  if (nclass == 2) 
    mu <- cbind(1, matrix(inv(eta + X %*% p[seq(npar)]), 
                          ncol = (nclass - 1)))
  else {
    if (link == "bcl") {
      e.eta <- exp(sapply(0:(nclass - 2), function(i, x, 
                                                   p, nclass) {
        x %*% p[c(1 + i, (ncol(x) - 1) * i + (nclass:(ncol(x) + 
                                                      nclass - 2)))]
      }, x = X, p = p, nclass = nclass))
      mu <- cbind(e.eta, 1)/(1 + c(e.eta %*% rep(1, nclass - 
                                                 1)))
      mu <- cbind(1, mu[, -1])
    }
    else {
      theta <- (p[seq(nclass - 1)])
      if (npar > (nclass - 1)) 
        eta <- eta + (X[, -1, drop = FALSE] %*% p[nclass:npar])
      if (link == "cum") {
        eta <- sapply(-theta, function(i, eta) i + eta, 
                      eta = eta)
        mu <- matrix(inv(-eta), ncol = (nclass - 1))
        mu <- cbind(1, cbind(mu[, 2:(nclass - 1)], 1) - 
                    mu[, seq(nclass - 1)])
      }
      if (link == "acl") {
        eta <- sapply(theta, function(i, eta) i + eta, 
                      eta = eta)
        tri <- lower.tri(matrix(1, ncol = (nclass - 1), 
                                nrow = (nclass - 1)), TRUE)
        mu <- exp(cbind(eta %*% tri, 0))/c(1 + exp(eta %*% 
                                                   tri) %*% rep(1, (nclass - 1)))
        mu <- cbind(1, mu[, -1])
      }
    }
  }
  tau <- matrix(1, nrow = nrow(tm), ncol = nass)
  nr <- nrep
  tord <- apply(tm, 1, sum)
  if (length(grep("N", dep)) > 0) {
    nf <- assp[, 1]
    assp <- assp[, -1, drop = FALSE]
    if (!Ncond) 
      mu[, -1] <- mu[, -1]/rep(nf[alab[!duplicated(lab)]], 
                               rep(nrep, length(unique(lab))))
  }
  else (nf <- 1)
  if (length(grep("M", dep)) > 0) {
    tau <- matrix(tau[, 1], nrow = (nclass)^2, ncol = (nrep - 
                                                       1))
    if (is.null(ass)) {
      tau[tord > 1, ] <- assp[, 1:(nclass - 1)^2] * tau[tord > 
                                                        1, ]
      assp <- assp[, -(1:(nclass - 1)^2), drop = FALSE]
    }
    else {
      tau[tord > 1, ] <- do.call("rbind",
                                 lapply(1:length(tord[tord >  1]),
                                        function(i, tlen, assp) {
                                          eval(parse(text = paste("ass[[", i, "]](", paste("assp[,", 
                                                       tlen[[i]], "]", collapse = ","), ")", sep = "")))
                                        }, tlen = tlen, assp = assp))
      assp <- assp[, -(unlist(tlen)), drop = FALSE]
    }
    nr <- 2
  }
  if (dep == "L" || dep == "NL") {
    if(Lclass==2)
      tau <- (sapply(1:nass, function(i, assp, tm, nclass, 
                                      nr) {
        (sapply(1:(nclass^nr), function(l, assp, tm, i) {
          dn <- prod(sapply(1:ncol(tm), function(m, l, 
                                                 assp, tm, i) {
            (assp[i, 1] + (1 - assp[i, 1]) * assp[i, (m + 
                                                      1)])^tm[l, m]
          }, assp = assp, tm = tm, l = l, i = i))
          num <- (assp[i, 1] + (1 - assp[i, 1]) *
                  prod(sapply(1:ncol(tm), 
                              function(m, l, assp, tm, i) {
                                (assp[i, (m + 1)]^tm[l, m])
                  }, assp = assp, tm = tm, l = l, i = i)))
          num/dn
        }, assp = assp, tm = tm, i = i))
        }, assp = assp, tm = tm, nclass = nclass, nr = nr))
    if(Lclass>2)
      tau <- sapply(1:nass, function(i, assp, tm, nr, Lclass) {
        sapply(1:(2^nr), function(l, assp, tm, Lclass, i) {
          sum(assp[i,Lclass-1],c(1 - sum(assp[i, 1:(Lclass-1)]),assp[i,1:(Lclass-2)])*
              assp[i,Lclass:(2*(Lclass-1))]^tm[l,1])/
                sum(assp[i,Lclass-1],c(1 - sum(assp[i, 1:(Lclass-1)]),assp[i,1:(Lclass-2)])*
                    assp[i,Lclass:(2*(Lclass-1))])^tm[l,1]  
        }, assp = assp, tm = tm, Lclass=Lclass, i = i)
      }, assp = assp, tm = tm, Lclass = Lclass, nr = nr)
  }
    if (length(grep("B", dep)) > 0) {
        tm[tm == 0] <- 1
        tau <- (tau * sapply(1:nass, function(i, assp, tm, nclass, 
            nr) {
            sapply(1:(nclass^nr), function(l, assp, tm, nclass, 
                i) {
                prod(sapply(2 * (1:(nclass - 1)), function(m, 
                  l, assp, tm, i) {
                  prod(c(assp[i, m - 1] + (1:tm[l, m/2] - 1), 
                    1/c(assp[i, m - 1] + assp[i, m] + (1:tm[l, 
                      m/2] - 1)))) * (((assp[i, m - 1] + assp[i, 
                    m])/assp[i, m - 1])^tm[l, m/2])
                }, assp = assp, tm = tm, l = l, i = i))
            }, assp = assp, tm = tm, nclass = nclass, i = i)
        }, assp = assp, tm = tm, nclass = nclass, nr = nr))
    }
    if (length(grep("D", dep)) > 0) {
        tau <- (sapply(1:nass, function(i, assp, tm, nclass, 
            nr) {
            sapply(1:(nclass^nr), function(l, assp, tm, nclass, 
                i) {
                prod(sapply(2:nclass, function(m, l, assp, tm, 
                  i) {
                  (prod((assp[i, m] + (0:tm[l, m - 1]) - 1))/(assp[i, 
                    m] - 1)) * ((sum(assp[i, ])/assp[i, m])^tm[l, 
                    m - 1])
                }, assp = assp[, 1:nclass, drop = FALSE], tm = tm, 
                  l = l, i = i)) * ((sum(assp[i, ]) - 1)/prod(sum(assp[i, 
                  ]) + (0:sum(tm[l, ])) - 1))
            }, assp = assp[, 1:nclass, drop = FALSE], tm = tm, nclass = nclass, 
                i = i)
        }, assp = assp[, 1:nclass, drop = FALSE], tm = tm, nclass = nclass, 
            nr = nr))
    }
    dropinit <- assp[, (length(assp) - (ncol(dropx) - 1)):length(assp)] %*% 
        t(dropx)
    dropinit <- 1/(1 + exp(-dropinit))
    dropm <- 1 - sapply(1:(nrep - 1), function(i, dropinit, droplab, 
        drop.cov, y) dropinit[match(paste(y[, i + 1], y[, i], 
        drop.cov[, i], sep = ""), droplab)], dropinit = dropinit, 
        droplab = droplab, drop.cov = drop.cov, y = y)
    dropm <- apply(dropm, 1, prod, na.rm = TRUE)
    dropm <- lapply(seq(dropm), function(i, misn, dropm, dropinit, 
        drop.cov, nclass) {
        if (any(is.na(misn[[i]]))) 
            dropm[i]
        else (dropm[i] * dropinit[match(paste(1:nclass, paste(misn[[i]], 
            collapse = ""), sep = ""), droplab)])
    }, misn = misn, dropm = dropm, dropinit = dropinit, nclass = nclass)
    if (dep == "M" || dep == "NM") {
        Mv <- lapply(unique(lab), function(i, mu, nrep, tau) {
            u <- sapply(1:(nrep - 1), function(j, i, mu, nrep, 
                tau) tau[, j] * (mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])), mu = mu, 
                i = i, nrep = nrep, tau = tau)
            cbind(u, u[, -1])
        }, mu = mu, nrep = nrep, tau = tau)
        lik <- sapply(seq(along = lab), function(i, Mv, lab, 
            w, dropm, Nf, nf, ind, le) {
            l <- apply(array(Mv[[lab[i]]], dim = dim(w[[i]])) * 
                w[[i]], 3, function(i, ind) ind %*% i, ind = ind)
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), , drop = FALSE], 
                  1/l[(1 + (le + 1)/2):le, , drop = FALSE]), 2, prod)
            sum(dropm[[i]] * (Nf[[i]] * (1 - nf) + nf * l))
        }, Mv = Mv, w = w, lab = lab, dropm = dropm, Nf = Nf, 
            le = 2 * nrep - 3, nf = nf, ind = rep(1, nclass^2))
    }
    if (length(grep("LM", dep)) > 0) {
        mul <- cbind(1, sapply(2:nclass, function(i, mu, assp) mu[, 
            i]/(assp[1, 1] + (1 - assp[1, 1]) * assp[1, i]), 
            mu = mu, assp = assp))
        muu <- cbind(1, sapply(2:nclass, function(i, mul, assp) (mul[, 
            i] * assp[1, i]), mul = mul, assp = assp))
        Mvl <- lapply(unique(lab), function(i, mu, nrep, tau) {
            u <- sapply(1:(nrep - 1), function(j, i, mu, nrep, 
                tau) tau[, j] * (mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])), mu = mul, 
                i = i, nrep = nrep, tau = tau)
            cbind(u, u[, -1])
        }, mu = mu, nrep = nrep, tau = tau)
        likl <- sapply(seq(along = lab), function(i, Mv, lab, 
            w, dropm, Nf, nf, ind, le) {
            l <- apply(array(Mv[[lab[i]]], dim = dim(w[[i]])) * 
                w[[i]], 3, function(i, ind) ind %*% i, ind = ind)
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), , drop = FALSE], 
                  1/l[(1 + (le + 1)/2):le, , drop = FALSE]), 2, prod)
            sum(dropm[[i]] * (Nf[[i]] * (1 - nf) + nf * l))
        }, Mv = Mvl, w = w, lab = lab, dropm = dropm, Nf = Nf, 
            le = 2 * nrep - 3, nf = nf, ind = rep(1, nclass^2))
        Mvu <- lapply(unique(lab), function(i, mu, nrep, tau) {
            u <- sapply(1:(nrep - 1), function(j, i, mu, nrep, 
                tau) tau[, j] * (mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])), mu = muu, 
                i = i, nrep = nrep, tau = tau)
            cbind(u, u[, -1])
        }, mu = mu, nrep = nrep, tau = tau)
        lik <- sapply(seq(along = lab), function(i, Mv, lab, 
            w, dropm, Nf, nf, ind, le) {
            l <- apply(array(Mv[[lab[i]]], dim = dim(w[[i]])) * 
                w[[i]], 3, function(i, ind) ind %*% i, ind = ind)
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), , drop = FALSE], 
                  1/l[(1 + (le + 1)/2):le, , drop = FALSE]), 2, prod)
            sum(dropm[[i]] * (Nf[[i]] * (1 - nf) + nf * l))
        }, Mv = Mvu, w = w, lab = lab, dropm = dropm, Nf = Nf, 
            le = 2 * nrep - 3, nf = nf, ind = rep(1, nclass^2))
        lik <- likl + lik
    }
    if (dep == "M2" || dep == "NM2") {
        Mv <- lapply(unique(lab), function(i, mu, nrep, tau) {
            u <- sapply(1:(nrep - 2), function(j, i, mu, nrep, 
                tau) tau[, j] * (c(mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])) %*% t(mu[j + 
                2 + (nrep * (i - 1)), ])), mu = mu, i = i, nrep = nrep, 
                tau = tau)
            cbind(u, u[, -1])
        }, mu = mu, nrep = nrep, tau = tau)
        lik <- sapply(seq(along = lab), function(i, Mv, lab, 
            w, dropm, Nf, nf, ind, le) {
            l <- apply(array(Mv[[lab[i]]], dim = dim(w[[i]])) * 
                w[[i]], 3, function(i, ind) ind %*% i, ind = ind)
            if (le > 1) 
                l <- apply(rbind(dropm[[i]], l[1:((le + 1)/2), 
                  , drop = FALSE], 1/l[(1 + (le + 1)/2):le, , drop = FALSE]), 
                  2, prod)
            sum(dropm[[i]] * (Nf[[i]] * (1 - nf) + nf * l))
        }, Mv = Mv, w = w, lab = lab, dropm = dropm, Nf = Nf, 
            le = 2 * nrep - 5, nf = nf, ind = rep(1, nclass^3))
    }
    if (length(grep("M", dep)) == 0) {
        Mv <- lapply(unique(lab), function(i, mu, nrep) eval(parse(text = paste(paste(rep("c(", 
            nrep - 1), collapse = ""), paste("mu[(1+(", nrep, 
            " * (i - 1))),]", paste("%*% \nmu[(", 2:nrep, "+ (", 
                nrep, " * (i - 1))),drop=FALSE,])", collapse = ""))))), 
            mu = mu, nrep = nrep)
        lik <- sapply(seq(along = lab), function(i, Mv, w, tau, 
            Nf, nf, lab, dropm, le) sum(dropm[[i]] * (Nf[[i]] * 
            (1 - nf) + nf * rep(1, le) %*% (c(tau) * Mv[[lab[i]]] * 
            w[[i]]))), Mv = Mv, w = w, tau = tau, Nf = Nf, nf = nf, 
            dropm = dropm, lab = lab, le = nclass^nrep)
    }
    if (any(is.na(lik)) || any(lik <= 0) || any(lik > 1)) 
        .Machine$double.xmax
    else (-(sum(wt * log(lik))))
}
