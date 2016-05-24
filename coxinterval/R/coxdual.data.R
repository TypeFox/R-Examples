### format data for coxdual.c
coxdual.data <- function(id, start, stop, from, to, contrib, z, states, sieve,
                         eps)
{
  p <- max(1, ncol(z))
  z <- data.frame(id, from, to, contrib, z)
  names(z) <- c("id", "from", "to", "contrib", paste("z", 1:p, sep = ""))
  ## type-specific covariates (nb: ? -> 2 presumed to hold values for 1 -> 2)
  z <- merge(merge(subset(z, from %in% states[1] & to == states[2]),
                   subset(z, from %in% states[1] & to == states[3]), by = "id"),
             subset(z, from %in% c(states[2], NA) & to == states[3]),
             by = "id", all = TRUE)
  z <- z[, substr(names(z), 1, 1) == "z"]
  z[is.na(z)] <- 0
  names(z) <- paste(paste("z", 1:p, ".", sep= ""),
                    rep(c("01", "02", "12"), each = p), sep = "")
  rownames(z) <- NULL
  uid <- unique(id)
  n <- length(uid)
  ## NA action permits missing 'start'/'stop' when 'start' = 'stop'
  start[is.na(start) & !is.na(from)] <- stop[is.na(start) & !is.na(from)]
  stop[is.na(stop)] <- start[is.na(stop)]
  ## largest and smallest observation times
  u <- as.vector(by(start, id, min, na.rm = TRUE))
  v <- as.vector(by(stop, id, max))
  vmax <- max(v)
  if (vmax == vmax - eps)
    stop("Observations large relative to epsilon. Use a smaller time scale.")
  eps <- min(eps, 1 / vmax)
  ## T observed?
  absorb <- is.element(uid, id[to == states[3] & contrib == 1])
  ## wlog largest observation is right-censored
  if (!sieve) absorb[v == vmax] <- FALSE
  ## contribution via 0 -> 1 (1), 0 -> 2 (2), both (0)?
  contrib <- rep(2, n)
  contrib[uid %in% id[from %in% states[2]]] <- 1
  contrib[uid %in% id[is.na(from)]] <- 0
  ## (possible) censoring intervals (L, R] with L = R if T01 observed exactly
  left <- stop[to == states[2]]
  left[!absorb & contrib == 2] <- v[!absorb & contrib == 2]
  right <- rep(Inf, n)
  right[contrib == 0] <- stop[is.na(from)]
  right[contrib == 1] <- start[from %in% states[2]]
  right[absorb & contrib == 0] <- v[absorb & contrib == 0]
  right[absorb & right == v] <- right[absorb & right == v] - eps
  if (sieve) {
    t01 <- maximalint(cbind(left, right)[contrib == 1, ], eps)$int[, 2]
    if (!any(contrib == 2 & absorb))
      stop("Support for progression-free survival is ambiguous.")
    if(!any(contrib == 1 & absorb))
      stop("Support for survival following progression is ambiguous.")
    t02 <- v[absorb & contrib == 2]
    t12 <- v[absorb & contrib == 1]
  }
  else {
    t01 <- maximalint(cbind(left, right)[contrib != 2, ], eps)$int[, 2]
    t02 <- v[absorb & contrib != 1]
    t12 <- v[absorb & contrib != 2]
  }
  names(t01) <- NULL
  t02 <- sort(unique(t02))
  t12 <- sort(unique(t12))
  smax <- pmin(right, v)
  sobs <- !sieve | contrib != 0
  r01 <- apply(sapply(t01, function(x) sobs & smax >= x), 2, sum)
  r02 <- apply(sapply(t02, function(x) sobs & smax >= x), 2, sum)
  r12 <- apply(sapply(t12, function(x) sobs & right <= x & x <= v), 2, sum)
  list(supp = list(t01 = t01, t02 = t02, t12 = t12),
       risk = list(r01 = r01, r02 = r02, r12 = r12),
       left = left, right = right, u = u, v = v, contrib = contrib,
       absorb = absorb, z = z)
}
