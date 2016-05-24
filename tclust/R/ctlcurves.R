
ctlcurves <-
function (x, k = 1:4, alpha = seq (0, 0.2, len = 6), 
          restr.fact = 50, trace = 1, ...)
{
## disabled parameter:
  mah.alpha <- 0.05


  stopifnot (length (restr.fact) == 1)

  if (trace >= 1)
  {
    cat ("Depending on arguments x, k and alpha, ",
         "this function needs some time to compute.\n",
         "(Remove this message by setting \"trace = 0\")\n", sep = "")
     flush.console()
  }

  rs <- array (NA, c (length (k), length (alpha)))
  dimnames (rs) <- list (k = k, alpha = round (alpha, 2))
  ndoubt <- drop <- out <- obj <- min.weights <- unrestr.fact <- rs

#  stab <- 
#  stab[,] <- TRUE

  p <- ncol (x)
  n.k <- length (k)
  n.alpha <- length (alpha)
  n <- nrow (x)

  for (i in 1:n.k)
  {
    for (j in 1:n.alpha)
    {
      cur.alpha <- alpha[j]

      if (trace >= 2)
        cat ("k =", k[i], "; alpha =", alpha[j])

      clus <- tclust (x, k = k[i], alpha = alpha[j], restr.fact = restr.fact,
                    warnings = FALSE, trace = .tr (trace), ...)

      out[i, j] <- sum (clus$mah > qchisq (1-mah.alpha, p), na.rm = TRUE)

      drop[i, j] <- clus$warnings$drop | clus$warnings$size |
                  clus$warnings$sizep


      dsc <- DiscrFact (clus)

      ndoubt[i, j] <- sum (dsc$assignfact > dsc$threshold)

      unrestr.fact[i, j] <- clus$unrestr.fact

      obj[i,j] <- clus$obj            ##  the objective function criterion
                                      ##  weights of the smallest group
      min.weights[i,j] <- min(clus$weights)



#     if (stab [i, j])
#    {
#      srf <- 1e7
#      opt <- tclust (x, k = k[i], alpha = alpha[j], restr.fact = srf,
#                   warnings = FALSE, trace = .tr (trace), ...)
#
#    if (opt$unrestr.fact >= 1e5 ||
#        opt$warn$size || opt$warn$sizep)
#      stab[i:n.k, j:n.alpha] <- FALSE
#      }

      if (trace >= 2)
        cat ("; obj = ", clus$obj, ";min (weights) =", min(clus$weights), "\n")
#     else if (trace == 1)
#      cat (".")
    }
  }

#  if (trace == 1)
#    cat ("\n"
#
#  unrestr.fact [!stab] <- NA
#
#  warn <- list (restr.lo = max (unrestr.fact [stab]) > restr.fact)
#
#  if (warn$restr.lo)
#  warning (paste ("The selected restriction factor (", 
#    ceiling (restr.fact), ") appears to be too small.\n",
#    "ctlcurves suggests to increase the argument restr.fact to ", 
#    ceiling (1.5 * max (unrestr.fact [stab])),".", sep = ""))

  par <- list (x = x, k = k, alpha = alpha, mah.alpha = mah.alpha,
               restr.fact = restr.fact)

  ret <- list (obj = obj, min.weights = min.weights, par = par, 
    unrestr.fact = unrestr.fact, ndoubt = ndoubt, out = out, drop = drop)
  #, stable = stab

##
##  calculating a suggestion
##
#
#  dif.obj <- apply (obj, 2, diff)
#  dif.obj [dif.obj < 0] <- 0
#
#  m.dif <- apply (dif.obj, 1, mean)
#
#  sug.alpha <- sug.k <- NULL
#  if (any (dif.obj [,1] / m.dif > 2))                          ## trimming seems to make sense
#    sug.k <- which.max (dif.obj [,1] / m.dif)  
#  else                      ## seems like trimming is not really needed
#     sug.k <- which (rowMeans(dif.obj)/diff(range(obj)) < 0.05)[1]
#
#  if (!is.null(sug.k))
#    sug.alpha <- which (out[sug.k,] / (n * (1 - alpha)) < mah.alpha)[1]
#
#  if (!warn$restr.lo &&
#      length (sug.alpha) &&
#      !is.na (unrestr.fact[sug.k, sug.alpha]))
#    ret$suggestion <- list (k = k[sug.k], 
#                            alpha = alpha[sug.alpha],
#                            restr.fact = ceiling (1.5 * unrestr.fact[sug.k, sug.alpha]))

  class (ret) <- "ctlcurves"
  ret
}

.tr <- function (trace, trace.red = 2)  ##  trace reduction
{
  max (trace - trace.red, 0)
}

plot.ctlcurves <-
function(x, what = c ("obj", "min.weights", "doubtful"),#, "out"
         main, xlab, ylab, xlim, ylim, col, lty = 1, ...)
{
  mark.art = FALSE

  what <- match.arg (what[1], eval (formals ()$what))

  n.k <- length (x$par$k)
  idxb.use <- rep (TRUE, n.k)
  if (what == "obj")
  {
    dat <- x$obj
    set.ylab <- "Objective Function Value"
  set.main <- "CTL-Curves"
  }
  else if (what == "out")
  {
    dat <- x$out
    set.ylab <- "Number of Outliers among Clusters"
  set.main <- "Outlier-Curves"
  }
  else if (what == "doubtful")
  {
    dat <- x$ndoubt
    set.ylab <- "Number of Doubtful Decision"
  set.main <- "Doubtfull Decision-Curves"
  }
  else if (what == "min.weights")
  {
  idxb.use <- x$par$k != 1       # k == 1 -> min.weights == 1 -> we're not interested in that. 
    dat <- x$min.weights
    set.ylab <- "Minimum Weigths"
  set.main <- "Minimum Weight-Curves"
  }

  if (missing (ylim)) ylim <- range (dat[idxb.use,])
  if (missing (xlim)) xlim <- range (x$par$alpha)
  if (missing (xlab)) xlab <- expression (alpha)
  if (missing (ylab)) ylab <- set.ylab
  if (missing (main)) main <- set.main

  n.alpha <- length (x$par$alpha)
  lty <- rep (lty, length (x$par$k))

  n.use.k <- sum (idxb.use)

  if (missing (col))
    col <- 1 + (1:n.k)
  else
    col <- rep (col, len = n.k)

  plot (0, type = "n", ylim = ylim, xlim = xlim ,
        main = main, xlab = xlab, ylab = ylab)
  mtext (paste ("Restriction Factor =", x$par$restr.fact), line = 0.25)

  pch <- rownames (dat)

#  idx.grey <- !x$stable
  if (mark.art)
    idx.grey <- x$par$restr.fact < x$unrestr.fact
  else
    idx.grey <- FALSE

  for (j in n.k:1)
  {
    if (!idxb.use[j])
    next

  cur.col <- rep (col[j], len = n.alpha)
  if (any (idx.grey))
      cur.col [idx.grey[j,]] <- 8

#    lines (x$par$alpha, dat[j,], type="b", col = cur.col, lty = lty[j],
#           pch = as.character (pch[j]))
    lines (x$par$alpha, dat[j,], type = "b", col = col[j], lty = lty[j],
           pch = "  ")

    text (x = x$par$alpha, y = dat[j,], col = cur.col, labels = pch[j])

  }

  if (what == "out")
    lines (x$par$alpha, nrow (x$par$x) * (1-x$par$alpha) * x$par$mah.alpha, lty = 2)
}

print.ctlcurves <-
function (x, ...)
{
  cat ("Computed ", length (x$par$k) * length (x$par$alpha), 
       " solutions (chosen restr.fact = ", x$par$restr.fact, ").\n\n",
     sep = "")

  idx.ar <- x$par$restr.fact < x$unrestr.fact

  if (any (idx.ar | x$drop))
  {
    uf <- array ("", dim = dim (x$unrestr.fact))
    uf[idx.ar] <- "*"
    uf[x$drop] <- paste (uf[x$drop], "k", sep = "")

    attributes (uf) <- attributes (x$unrestr.fact)

  print (uf, justify = "right", quote = FALSE)

    if (any (idx.ar))
      cat ("\n(*) Identified ", sum (idx.ar), " artificially restricted solutions.", sep = "")
    if (any (x$drop))
      cat ("\n(k) Identified ", sum (x$drop), " solutions with very small/dropped clusters.", sep = "")
  }
  else
    cat ("\nNo artificially restricted solutions or dropped clusters found.")

 cat ("\n")
}

              ##  discarded versions of print.ctlcurves.
              ##  Some features might "return in the future..

#{
#print.ctlcurves <-
#function (x, ...)
#{
#  cat ("The models' individual restriction factors:\n")
#  print (x$unrestr.fact)
#  mrf <- max (x$unrestr.fact, na.rm = TRUE)
#  cat ("\n")
#  if (mrf > x$par$restr.fact)
#    cat ("ctlcurves suggests to increase restr.fact to", ceiling (mrf * 1.5), "\n")
#
#  cat ("The number of doubtful decisions for each models:\n\n")
#  print (x$ndoubt)
#  cat ("\n\n")
#
#  cat ("The fraction of outliers found within the clusters:\n\n")
#
#   out <- x$out %*% diag ((1 / (nrow (x$par$x) * (1 - x$par$alpha))))
#  dimnames (out) <- dimnames (x$out)
#  print (round (out, 2))
#  cat ("\n\n")
#}

#print.ctlcurves<-
#function (x, trace = 1, ...)
#{
#  cat ("Analyzed ", length (x$par$k) * length (x$par$alpha), " models.\n",
#       "Found ", sum (x$stable), " stable models.\n"
#     , sep = "")
#  if (!is.null (x$suggestion))
#  {
#    cat ("\nA heuristic analysis of the ctl-curves suggests the following parameters:\n\n",
#       "  k:          ", x$suggestion$k, "\n",
#     "  alpha:      ", x$suggestion$alpha, "\n",
#     "  restr.fact: ", x$suggestion$restr.fact, "\n\n", sep = "")
#    if (trace >= 1)
#     cat ("NOTE: this suggestion is based on the assumption of normality\n",
#     "and moderate contamination.\n",
#     "DON'T trust this suggestion blindly and check the resulting model manually\n",
#     " (e.g. with \"DiscrFact\" or \"plot (ctl)\").\n\n",
#     sep = "")
#  }
#  else
#    cat ("No parameter constellation could be found automatically\n",
#       "Try increasing alpha and k\n\n", sep = "")
#}

#print.ctlcurves <-
#function (x, ...)
#{
##  n.stable <- sum (x$stable)
##       , "Found ", n.stable, " stable models.\n"
#
#  cat ("Computed ", length (x$par$k) * length (x$par$alpha), 
#       " solutions (chosen restr.fact = ", x$par$restr.fact, ").\n", sep = "")
##  cat ("Analyzed ", length (x$par$k) * length (x$par$alpha), " models for values of\n", sep = "")
##  cat ("  k:          ", paste (x$par$k, collapse = ", "), "\n")
##  cat ("  alpha:      ", paste (round (x$par$alpha, 2), collapse = ", "), "\n")
##  cat ("  restr.fact: ", x$par$restr.fact, "\n\n")
#
#  n.restr.art <- sum (x$par$restr.fact < x$unrestr.fact)
#  cat ("The solutions' individual restriction factors:\n\n")
#
#  uf <- array (dim = dim (x$unrestr.fact))
#  for (i in 1:ncol (uf))
#   for (j in 1:nrow (uf))
#    uf[j, i] <- format (x$unrestr.fact[j, i], scientific = x$unrestr.fact[j, i] >= 1e5, digits = 2)
#
#  ufs <- array ("", dim = dim (uf))
#  ufs[x$par$restr.fact < x$unrestr.fact] <- "*"
#
#  ufpr <- paste (uf, ufs, sep = "")
##  dim (ufpr) <- dim (uf)
##  dimnames (ufpr) <- dimnames (uf)
#  attributes (ufpr) <- attributes (x$unrestr.fact)
#  print (noquote (ufpr))
#
#  if (n.restr.art)
#    cat ("\n(*) Identified ", n.restr.art, " artificially restricted solutions.", sep = "")
# cat ("\n")
#}
#
