
truncationRegion <- function(obj, ydecomp, type, tol = 1e-15) {

  n <- nrow(obj$x)
  Z <- ydecomp$Z
  if (type == "TC") {
      eta <- ydecomp$eta
  } else {
      Vd <- ydecomp$Vd
      V2 <- ydecomp$V2
      C <- ydecomp$C
      R <- ydecomp$R
  }
  L <- lapply(1:length(obj$action), function(s) {

    Ug <- obj$maxprojs[[s]]
    peng <- obj$maxpens[[s]]
    if (s > 1) {
        Zs <- obj$cumprojs[[s-1]] %*% Z
        if (type == "TC") {
            etas <- obj$cumprojs[[s-1]] %*% eta
        } else {
            Vds <- obj$cumprojs[[s-1]] %*% Vd
            V2s <- obj$cumprojs[[s-1]] %*% V2
        }
    } else {
        Zs <- Z
        if (type == "TC") {
            etas <- eta
        } else {
            Vds <- Vd
            V2s <- V2
        }
    }

    num.projs <- length(obj$projections[[s]])
    if (num.projs == 0) {
        return(list(Intervals(c(-Inf,0))))
    } else {
      lapply(1:num.projs, function(l) {

          Uh <- obj$projections[[s]][[l]]
          penh <- obj$aicpens[[s]][[l]]
          # The quadratic form corresponding to
          # (t*U + Z)^T %*% Q %*% (t*U + Z) \geq 0
          # we find the roots in t, if there are any
          # and return the interval of potential t
          if (type == "TC") {
              coeffs <- quadratic_coefficients(obj$sigma, Ug, Uh, peng, penh, etas, etas, Zs, Zs)
              quadratic_roots(coeffs$A, coeffs$B, coeffs$C, tol)
          } else {
              coeffs <- TF_coefficients(R, Ug, Uh, peng, penh, Zs, Zs, Vds, Vds, V2s, V2s)
              roots <- TF_roots(R, C, coeffs)
              return(roots)
          }
      })
    }
    # LL is a list of intervals
  })
  # L is now a list of lists of intervals
  return(unlist(L, recursive = FALSE, use.names = FALSE))
}

quadratic_coefficients <- function(sigma, Ug, Uh, peng, penh, etag, etah, Zg, Zh) {
    # g indexes minimizer, h the comparison
    Uheta <- t(Uh) %*% etah
    Ugeta <- t(Ug) %*% etag
    UhZ <- t(Uh) %*% Zh
    UgZ <- t(Ug) %*% Zg
    etaZh <- t(etah) %*% Zh
    etaZg <- t(etag) %*% Zg
    if (is.null(sigma)) {
        A <- penh * (sum(etah^2) - sum(Uheta^2)) - peng * (sum(etag^2) - sum(Ugeta^2))
        B <- 2 * penh * (etaZh - t(Uheta) %*% UhZ) - 2 * peng * (etaZg - t(Ugeta) %*% UgZ)
        C <- penh * (sum(Zh^2) - sum(UhZ^2)) - peng * (sum(Zg^2) - sum(UgZ^2))
    } else {
        A <- (sum(etah^2) - sum(Uheta^2)) - (sum(etag^2) - sum(Ugeta^2))
        B <- 2 * (etaZh - t(Uheta) %*% UhZ) - 2 * (etaZg - t(Ugeta) %*% UgZ)
        C <- (sum(Zh^2) - sum(UhZ^2) + penh) - (sum(Zg^2) - sum(UgZ^2) + peng)
    }
    return(list(A = A, B = B, C = C))
}

quadratic_roots <- function(A, B, C, tol) {
    disc <- B^2 - 4*A*C
    b2a <- -B/(2*A)

    if (disc > tol) {
        # Real roots
        pm <- sqrt(disc)/(2*A)
        endpoints <- sort(c(b2a - pm, b2a + pm))

    } else {
        # No real roots
        if (A > -tol) {
          # Quadratic form always positive
            return(Intervals(c(-Inf,0)))
        } else {
          # Quadratic form always negative
            stop("Empty TC support is infeasible")
        }
    }

    if (A > tol) {
        # Parabola opens upward
        if (min(endpoints) > 0) {
          # Both roots positive, union of intervals
            return(Intervals(rbind(c(-Inf,0), endpoints)))
        } else {
          # At least one negative root
            return(Intervals(c(-Inf, max(0, endpoints[2]))))
        }
    } else {
        if (A < -tol) {
          # Parabola opens downward
            if (endpoints[2] < 0) {
            # Positive quadratic form only when t negative
                stop("Negative TC support is infeasible")
            } else {
            # Part which is positive
                if (endpoints[1] > 0) {
                    return(Intervals(rbind(c(-Inf, endpoints[1]), c(endpoints[2], Inf))))
                } else {
                    return(Intervals(c(endpoints[2], Inf)))
                }
            }
        } else {
          # a is too close to 0, quadratic is actually linear
            if (abs(B) > tol) {
                if (B > 0) {
                    return(Intervals(c(-Inf, max(0, -C/B))))
                } else {
                    if (-C/B < 0) stop("Infeasible linear equation")
                    return(Intervals(rbind(c(-Inf, 0), c(-C/B, Inf))))
                }
            } else {
                warning("Ill-conditioned quadratic")
                return(Intervals(c(-Inf,0)))
            }
        }
    }
}

# Helper functions for TF roots
roots_to_checkpoints <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(c(0, (checkpoints + c(checkpoints[-1], 200 + checkpoints[length(checkpoints)]))/2))
}
roots_to_partition <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(list(endpoints = c(checkpoints, Inf), midpoints = (checkpoints + c(checkpoints[-1], 200 + checkpoints[length(checkpoints)]))/2))
}

# Efficiently compute coefficients of one-dimensional TF slice function
TF_coefficients <- function(R, Ug, Uh, peng, penh, Zg, Zh, Vdg, Vdh, V2g, V2h) {

    UhZ <- t(Uh) %*% Zh
    UgZ <- t(Ug) %*% Zg
    UhVd <- t(Uh) %*% Vdh
    UgVd <- t(Ug) %*% Vdg
    UhV2 <- t(Uh) %*% V2h
    UgV2 <- t(Ug) %*% V2g
    VdZh <- sum(Vdh*Zh)
    VdZg <- sum(Vdg*Zg)
    V2Zh <- sum(V2h*Zh)
    V2Zg <- sum(V2g*Zg)

    x0 <- penh * (sum(Zh^2) - sum(UhZ^2)) - peng * (sum(Zg^2) - sum(UgZ^2))
    x1 <- 2*R*(penh * (VdZh - sum(UhZ*UhVd)) - peng * (VdZg - sum(UgZ*UgVd)))
    x2 <- 2*R*(penh * (V2Zh - sum(UhZ*UhV2)) - peng * (V2Zg - sum(UgZ*UgV2)))
    x12 <- 2*R^2*(penh * (sum(Vdh*V2h) - sum(UhVd*UhV2)) - peng * (sum(Vdg*V2g) - sum(UgVd*UgV2)))
    x11 <- R^2*(penh * (sum(Vdh^2) - sum(UhVd^2)) - peng * (sum(Vdg^2) - sum(UgVd^2)))
    x22 <- R^2*(penh * (sum(V2h^2) - sum(UhV2^2)) - peng * (sum(V2g^2) - sum(UgV2^2)))

    return(list(x11=x11, x22=x22, x12=x12, x1=x1, x2=x2, x0=x0))
}

# Numerically solve for roots of TF slice using
# hybrid polyroot/uniroot approach
TF_roots <- function(R, C, coeffs, tol = 1e-8, tol2 = 1e-6) {

    x11 <- coeffs$x11
    x22 <- coeffs$x22
    x12 <- coeffs$x12
    x1 <- coeffs$x1
    x2 <- coeffs$x2
    x0 <- coeffs$x0

    g1 <- function(t) sqrt(C*t/(1+C*t))
    g2 <- function(t) 1/sqrt(1+C*t)
    I <- function(t) x11*g1(t)^2 + x12*g1(t)*g2(t) + x22*g2(t)^2 + x1*g1(t) + x2*g2(t) + x0

    z4 <- complex(real = -x11 + x22, imaginary = -x12)/4
    z3 <- complex(real = x2, imaginary = -x1)/2
    z2 <- complex(real = x11/2+x22/2+x0)
    z1 <- Conj(z3)
    z0 <- Conj(z4)

    zcoefs <- c(z0, z1, z2, z3, z4)
    croots <- polyroot(zcoefs)
    thetas <- Arg(croots)
    # Can't specify polyroot precision :(
    modinds <- Mod(croots) <= 1 + tol2 & Mod(croots) >= 1 - tol2
    angleinds <- thetas >=0 & thetas <= pi/2
    roots <- unique(thetas[which(modinds & angleinds)])
    troots <- tan(roots)^2/C

    checkpoints <- c()
    if (length(troots) > 0) checkpoints <- roots_to_checkpoints(troots)
    checkpoints <- sort(
        c(checkpoints, 0, tol, tol2,
                        seq(from = sqrt(tol2), to = 1, length.out = 50),
                        seq(from = 1.2, to=50, length.out = 20),
                        100, 1000, 10000))
    ## if (length(troots) == 0) {
    ##     # Polyroot didn't catch any roots
    ##     # ad-hoc check:
    ##     checkpoints <- c(0, tol, tol2,
    ##                      seq(from = sqrt(tol2), to = 1, length.out = 50),
    ##                      seq(from = 1.2, to=50, length.out = 20),
    ##                      100, 1000, 10000)
    ## } else {
    ##     checkpoints <- roots_to_checkpoints(troots)
    ## }

    signs <- sign(I(checkpoints))
    diffs <- c(0, diff(signs))
    changeinds <- which(diffs != 0)

    if (length(changeinds) > 0) {

        roots <- unlist(lapply(changeinds, function(ind) {
            uniroot(I, lower = checkpoints[ind-1], upper = checkpoints[ind], tol = tol)$root
        }))

        partition <- roots_to_partition(roots)
        negative <- which(I(partition$midpoints) < 0)

        intervals <- matrix(NA, ncol=2)
        for (i in 1:length(negative)) {
            ind <- negative[i]
            if ((i > 1) && (ind == negative[i-1] + 1)) {
                # There was not a sign change at end of previous interval
                intervals[nrow(intervals), 2] <- partition$endpoints[ind+1]
            } else {
                intervals <- rbind(intervals, c(partition$endpoints[ind], partition$endpoints[ind+1]))
            }
        }

        return(Intervals(intervals[-1,]))
    }

    # Apparently no roots, always positive
    if (I(0) < 0) stop("Infeasible constraint!")
    return(Intervals(c(-Inf,0)))
}
