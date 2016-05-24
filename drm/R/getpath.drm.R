"getpath.drm" <-
  function (p, v, X, dep, asscov, npar, nclass, nrep, tm, w, lab, 
            equal, inv, ass, link, constant, Ncond, Lclass, offset, wt, tlen, 
            save.profiles) 
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
  if (length(asscov[[1]]) > 0) {
    assm <- sapply(seq(along = asscov[[1]]), function(i, 
                         assp, asscov) asscov[[1]][[i]] %*% assp[asscov[[3]][[i]]], 
                   asscov = asscov, assp = assp)
    assp <- matrix(assp, ncol = length(assp), nrow = nrow(assm), 
                   byrow = TRUE)
    assp[, asscov[[2]]] <- assm
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
    assp <- assp[1, -1, drop = FALSE]
    if (!Ncond) {
      mu[, -1] <- mu[, -1]/rep(nf[alab[!duplicated(lab)]], 
                               rep(nrep, length(unique(lab))))
    }
  }
  if (length(grep("M", dep)) > 0) {
    if (length(grep("M2", dep)) > 0) {
      tau <- matrix(tau[, 1], nrow = (nclass)^3, ncol = 2 * 
                    nrep - 5)
      if (is.null(ass)) {
        tau[tord > 1, ] <- assp[, c(1, 2, 1, 3)] * tau[tord > 
                                                       1, ]
        if (nrep > 3) 
          tau[7, -(nrep - 2)] <- tau[4, -1]
        assp <- assp[, -(1:3), drop = FALSE]
      }
      else {
        tmp <- lapply(1:3, function(i, tlen, assp) {
          eval(parse(text = paste("ass[[", i, "]](", 
                       paste("assp[,", tlen[[i]], "]", collapse = ","), 
                       ")", sep = "")))
        }, tlen = tlen, assp = assp)
        tmp <- rbind(tmp[[1]][-length(tmp[[1]])], tmp[[2]], 
                     tmp[[1]][-1], tmp[[3]])
        if (ncol(tau) > (nrep - 2)) 
          tau[tord > 1, ] <- cbind(tmp, tmp[, -(ncol(tmp))])
        else (tau[tord > 1, ] <- tmp)
        assp <- assp[, -(unlist(tlen)), drop = FALSE]
      }
      nr <- 3
    }
    else {
      tau <- matrix(tau[1, ], nrow = nclass^2, ncol = (nrep - 
                                1))
      if (is.null(ass)) {
        tau[tord > 1, ] <- assp[, 1:(nclass - 1)^2] * 
          tau[tord > 1, ]
        assp <- assp[, -(1:(nclass - 1)^2), drop = FALSE]
      }
      else {
        tau[tord > 1, ] <- do.call("rbind",
                                   lapply(1:length(tord[tord >1]),
                                          function(i, tlen, assp) {
                                            eval(parse(text = paste("ass[[", i, "]](", 
                                                         paste("assp[,", tlen[[i]], "]", collapse = ","), 
                                                         ")", sep = "")))
                                          }, tlen = tlen, assp = assp))[, 1:(nrep - 1)]
        assp <- assp[, -(unlist(tlen)), drop = FALSE]
      }
      nr <- 2
    }
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
          num <- (assp[i, 1] + (1 - assp[i, 1]) * prod(sapply(1:ncol(tm), 
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
        tau <- (sapply(1:nass, function(i, assp, tm, nclass, 
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
    if (dep == "M" || dep == "NM") {
        lik <- sapply(unique(lab), function(i, mu, w, nrep, tau, 
            le) {
            Mv <- sapply(1:(nrep - 1), function(j, i, mu, nrep, 
                tau) tau[, j] * (mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])), mu = mu, 
                i = i, nrep = nrep, tau = tau)
            l <- apply(array(cbind(Mv, Mv[, -1]), dim = dim(w)) * 
                w, 3, function(i, ind) ind %*% i, ind = rep(1, 
                nclass^2))
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), ], 1/l[(1 + 
                  (le + 1)/2):le, ]), 2, prod)
            if (save.profiles) 
                l
            else (ifelse(min(l) < 0, -1, 1))
        }, mu = mu, w = w, nrep = nrep, tau = tau, le = 2 * nrep - 
            3)
    }
    if (length(grep("LM", dep)) > 0) {
        mul <- cbind(rep(1, nrow(mu)), sapply(2:nclass, function(i, 
            mu, assp) mu[, i]/(assp[1, 1] + (1 - assp[1, 1]) * 
            assp[1, i]), mu = mu, assp = assp))
        likl <- sapply(unique(lab), function(i, mu, w, nrep, 
            tau, le) {
            Mv <- sapply(1:(nrep - 1), function(j, i, mu, nrep, 
                tau) tau[, j] * (mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])), mu = mu, 
                i = i, nrep = nrep, tau = tau)
            Mv <- cbind(Mv, Mv[, -1])
            l <- apply(array(Mv, dim = dim(w)) * w, 3, function(i, 
                ind) ind %*% i, ind = rep(1, nclass^2))
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), ], 1/l[(1 + 
                  (le + 1)/2):le, ]), 2, prod)
            l[is.na(l)] <- 0
            l
        }, mu = mul, w = w, nrep = nrep, tau = tau, le = 2 * 
            nrep - 3)
        muu <- cbind(rep(1, nrow(mu)), sapply(2:nclass, function(i, 
            mu, assp) (mu[, i] * assp[1, i])/(assp[1, 1] + (1 - 
            assp[1, 1]) * assp[1, i]), mu = mu, assp = assp))
        liku <- sapply(unique(lab), function(i, mu, w, nrep, 
            tau, le) {
            Mv <- sapply(1:(nrep - 1), function(j, i, mu, nrep, 
                tau) tau[, j] * (mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])), mu = mu, 
                i = i, nrep = nrep, tau = tau)
            Mv <- cbind(Mv, Mv[, -1])
            l <- apply(array(Mv, dim = dim(w)) * w, 3, function(i, 
                ind) ind %*% i, ind = rep(1, nclass^2))
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), ], 1/l[(1 + 
                  (le + 1)/2):le, ]), 2, prod)
            l[is.na(l)] <- 0
            l
        }, mu = muu, w = w, nrep = nrep, tau = tau, le = 2 * 
            nrep - 3)
        lik <- assp[1, 1] * likl + (1 - assp[1, 1]) * liku
        if (!save.profiles) 
            lik <- ifelse(apply(lik, 2, min) < 0, -1, 1)
    }
    if (dep == "M2" || dep == "NM2") {
        lik <- sapply(unique(lab), function(i, mu, w, nrep, tau, 
            le) {
            Mv <- sapply(1:(nrep - 2), function(j, i, mu, nrep, 
                tau) tau[, j] * (c(mu[j + (nrep * (i - 1)), ] %*% 
                t(mu[j + 1 + (nrep * (i - 1)), ])) %*% t(mu[j + 
                2 + (nrep * (i - 1)), ])), mu = mu, i = i, nrep = nrep, 
                tau = tau)
            cbind(Mv, Mv[, -1])
            l <- apply(array(Mv, dim = dim(w)) * w, 3, function(i, 
                ind) ind %*% i, ind = rep(1, nclass^3))
            if (le > 1) 
                l <- apply(rbind(l[1:((le + 1)/2), ], 1/l[(1 + 
                  (le + 1)/2):le, ]), 2, prod)
            if (save.profiles) 
                l
            else (ifelse(min(l) < 0, -1, 1))
        }, mu = mu, w = w, nrep = nrep, tau = tau, le = 2 * nrep - 
            5)
    }
    if (length(grep("M", dep)) == 0) {
        lik <- sapply(unique(lab), function(i, mu, w, tau, alab, 
            nclass, nrep) {
            Mv <- eval(parse(text = paste(paste(rep("c(", nrep - 
                1), collapse = ""), paste("mu[(1+(", nrep, " * (i - 1))),]", 
                paste("%*% \nmu[(", 2:nrep, "+ (", nrep, " * (i - 1))),drop=FALSE,])", 
                  collapse = "")))))
            p <- rep(1, nclass^nrep) %*% (tau[, alab[i]] * Mv * 
                w)
            if (save.profiles) 
                p
            else (ifelse(min(p) < 0, -1, 1))
        }, mu = mu, nclass = nclass, nrep = nrep, w = w, tau = tau, 
            alab = alab[!duplicated(lab)])
    }
    if (save.profiles) {
        if (length(grep("N", dep)) > 0) {
            lik <- nf[alab[!duplicated(lab)]] * lik
            lik[1, ] <- (1 - nf[alab[!duplicated(lab)]]) + lik[1, 
                ]
        }
    }
    if (!save.profiles) 
        lik <- matrix(lik, nrow = 1)
    lik[, lab]
}

