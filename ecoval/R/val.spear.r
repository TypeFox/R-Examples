val.spear.create <- function(language="English",dictionaries=NA,col="black")
{
  # ============================================================================
  #
  # References:
  #
  # Langhans, S.D. und Reichert, P. (2011), Einbettung von Verfahren zur Fliess-
  # gewaesserbewertung in ein uebergeordnetes Gewaessermanagementkonzept - 
  # Vorschlaege am Beispiel des Modulstufenkonzepts, 
  # Wasser Energie Luft 103(3), 204-214. 
  #
  # Beketov M.A., Foit K., Schafer R.B., Schriever C.A., Sacchi A., Capri E., 
  # Biggs J., Wells C. & Liess M. (2009). SPEAR indicates pesticide effects in 
  # streams - Comparative use of species- and family-level biomonitoring data. 
  # Environmental Pollution, 157, 1841-1848
  #
  # ============================================================================
  
  # dictionary for node, attribute and attribute level names:
  # =========================================================
  
  dict <- ecoval.dict(language,dictionaries)
  
  # implementation of valuation node:
  # ---------------------------------
  
  spear <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_spear",dict),
      name.attrib = ecoval.translate("A_spear_ind",dict),
      range       = c(0,70),
      x           = c(70, 55, 44, 33, 22, 11, 0),
      u           = c(1.0,1.0,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  return(spear)
}   
