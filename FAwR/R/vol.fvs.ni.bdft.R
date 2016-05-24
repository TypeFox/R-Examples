vol.fvs.ni.bdft <-
function(spp, dbh.in, ht.ft) {
  # Parameters from Wykoff et al. (1982)
  bf.params <-
    data.frame(species = c("WP", "WL", "DF", "GF", "WH", "WC", "LP", "ES",
                 "SF", "PP", "MH"),
               b0.small = c(26.729, 29.790, 25.332, 34.127, 37.314, 10.472,
                 8.059, 11.851, 11.403, 50.340, 37.314),
               b1.small = c(0.01189, 0.00997, 0.01003, 0.01293, 0.01203,
                 0.00878, 0.01208, 0.01149, 0.01011, 0.01201, 0.01203),
               b0.large = c(32.516, 85.150, 9.522, 10.603, 50.680, 4.064,
                 14.111, 1.620, 124.425, 298.784, 50.680),
               b1.large = c(0.01181, 0.00841, 0.01011, 0.01218, 0.01306,
                 0.00799, 0.01103, 0.01158, 0.00694, 0.01595, 0.01306))
  b0 <- b1 <- rep(NA, length(spp))
  dimensions <- data.frame(dbh.in,
                           ht.ft,
                           this.order = 1:length(spp),
                           species = spp)
  small <- !(dbh.in > 20.5)
  small[is.na(small)] <- TRUE # Arbitrary; dbh values are missing
  dimensions <- merge(x=dimensions, y=bf.params, all.x=TRUE, all.y=FALSE)
  dimensions <- dimensions[order(dimensions$this.order, decreasing=FALSE),]
  b0[small] <- dimensions$b0.small[small]
  b1[small] <- dimensions$b1.small[small]
  b0[!small] <- dimensions$b0.large[!small]
  b1[!small] <- dimensions$b1.large[!small]
  b0 + b1 * dimensions$dbh.in^2 * dimensions$ht.ft
}

