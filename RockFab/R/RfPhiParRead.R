RfPhiParRead <-
function(rfphi.file){
  my.par <- read.table(file = rfphi.file, header = TRUE, sep = '\t')
  names(my.par) <- c("m.majoraxis", "m.minoraxis", "m.rake.deg")
  my.par$m.eccentricity <- sqrt(1 - (my.par$m.minoraxis^2 / my.par$m.majoraxis^2))
  my.par$m.rake.rad <- my.par$m.rake.deg * (pi / 180)
  my.par$m.ratio <- my.par$m.majoraxis / my.par$m.minoraxis
  my.par$m.area <- pi * (my.par$m.majoraxis / 2) * (my.par$m.minoraxis / 2)
  return(my.par)
}
