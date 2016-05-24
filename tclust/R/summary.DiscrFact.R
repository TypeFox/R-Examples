summary.DiscrFact <-
function (object, hide.emtpy = TRUE, show.clust, show.alt, ...)
{
  idx <- object$assignfact > object$threshold
  k <- object$x$k
  
  if (missing (show.clust))
    show.clust <- k <= 20
  if (missing (show.alt))
    show.alt <- k <= 5

  k1 <- k + 1

  dn <- c ("O", 1:k)

  if (show.clust)
  {
    cat ("\nNumber of doubtful assignments in clusters:\n")
    facts1 <- object$ind[idx] + 1
    cha <- as.array (tabulate (facts1, nbins = k1))
    dimnames (cha) <- list (dn)
    if (hide.emtpy)
      cha <- cha[cha != 0]
    print (cha)
  }
  
  if (show.alt)
  {
    facts2 <- object$ind[idx] * k1 + object$ind2[idx]  + 1

    #chma <- chma[apply (chma, 1, sum) != 0, , drop = FALSE]

    cat ("\nObservations with doubtful decision are alternatively assigned to clusters:\n")  
    chma <- matrix (tabulate (facts2, nbins = k1 * k1), ncol =k1, nrow = k1)
    dimnames (chma) <- list (dn, dn)
    if (hide.emtpy)
    {
      chma <- chma[apply (chma, 1, sum) != 0, , drop = FALSE]
      chma <- chma[, apply (chma, 2, sum) != 0, drop = FALSE]
    }
    print (chma)
  }  
}

