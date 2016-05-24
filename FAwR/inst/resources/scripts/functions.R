## These functions are used in Chapter 3

ht.fvs.ni.ft <- function(spp, dbh.in) {
  # Parameters from Wykoff et al. (1982)
  ht.params <-
    data.frame(species = c("WP", "WL", "DF", "GF", "WH", "WC", "LP", "ES",
                 "SF", "PP", "MH"),
               c0 = c(5.19988, 4.97407, 4.81519, 5.00233, 4.97331, 4.89564,
                 4.62171, 4.92190, 4.76537, 4.92880, 4.77951),
               c1 = c(-9.26718, -6.78347, -7.29306, -8.19365, -8.19370,
                 -8.39057, -5.32481, -8.30289, -7.61062, -9.32795, -9.31743))
  dimensions <- as.data.frame(dbh.in)
  dimensions$this.order <- c(1:length(spp))
  dimensions$species <- spp
  dimensions <- merge(x=dimensions, y=ht.params, all.x=T, all.y=F)
  dimensions <- dimensions[order(dimensions$this.order, decreasing=F),]
  4.5 + exp(dimensions$c0 + dimensions$c1 / (dimensions$dbh.in + 1))
}


ht.fvs.ni.m <- function(spp, dbh.cm) 
  ht.fvs.ni.ft(spp, dbh.cm/2.54) * 0.3048

vol.fvs.ni.bdft <- function(spp, dbh.in, ht.ft) {
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
  dimensions <- as.data.frame(cbind(dbh.in, ht.ft))
  dimensions$this.order <- c(1:length(spp))
  dimensions$species <- spp
  small <- !(dbh.in > 20.5)
  small[is.na(small)] <- TRUE # Arbitrary; dbh values are missing
  dimensions <- merge(x=dimensions, y=bf.params, all.x=T, all.y=F)
  dimensions <- dimensions[order(dimensions$this.order, decreasing=F),]
  b0[small] <- dimensions$b0.small[small]
  b1[small] <- dimensions$b1.small[small]
  b0[!small] <- dimensions$b0.large[!small]
  b1[!small] <- dimensions$b1.large[!small]
  b0 + b1 * dimensions$dbh.in^2 * dimensions$ht.ft
}

vol.fvs.ni.m3 <- function(spp, dbh.cm, ht.m) 
  vol.fvs.ni.bdft(spp, dbh.cm/2.54, ht.m/0.3048) *
  144 * # Board feet to cubic inches
  2.54^3 / # Cubic inches to cubic centimetres
  100^3  # Cubic centimetres to cubic metres

BAF.imp.2.met <- function(BAF.imp)
  BAF.imp *
  0.3048^2 * # Square feet to square metres
  2.47105381 # Acres to hectares
  
BAF.met.2.imp <- function(BAF.met)
  BAF.met /
  0.3048^2 / # Square metres to square feet
  2.47105381 # Hectares to acres

show.cols.with.na <- function(x) {
  if (class(x) != "data.frame") 
    stop("x must be a dataframe.\n")
  missing.by.column <- apply(is.na(x), 2, sum)
  if (sum(missing.by.column) == 0) {
    cat("No missing values.\n")
  } else {
    missing <- which(missing.by.column > 0)
    return(missing.by.column[missing])
  }
}
