vol.fvs.ni.m3 <-
function(spp, dbh.cm, ht.m) 
  vol.fvs.ni.bdft(spp, dbh.cm/2.54, ht.m/0.3048) *
  144 * # Board feet to cubic inches
  2.54^3 / # Cubic inches to cubic centimetres
  100^3  # Cubic centimetres to cubic metres

