ht.fvs.ni.ft <-
function(spp, dbh.in) {
  # Parameters from Wykoff et al. (1982)
  ht.params <-
    data.frame(species = c("WP", "WL", "DF", "GF", "WH", "WC", "LP", "ES",
                 "SF", "PP", "MH"),
               c0 = c(5.19988, 4.97407, 4.81519, 5.00233, 4.97331, 4.89564,
                 4.62171, 4.92190, 4.76537, 4.92880, 4.77951),
               c1 = c(-9.26718, -6.78347, -7.29306, -8.19365, -8.19370,
                 -8.39057, -5.32481, -8.30289, -7.61062, -9.32795, -9.31743))
  dimensions <- data.frame(dbh.in = dbh.in,
                           this.order = 1:length(spp),
                           species = spp)
  dimensions <- merge(x=dimensions, y=ht.params, all.x=TRUE, all.y=FALSE)
  dimensions <- dimensions[order(dimensions$this.order, decreasing=FALSE),]
  4.5 + exp(dimensions$c0 + dimensions$c1 / (dimensions$dbh.in + 1))
}

