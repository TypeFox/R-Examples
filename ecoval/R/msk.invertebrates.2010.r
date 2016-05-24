msk.invertebrates.2010.create <- function(language     = "English",
                                          dictionaries = NA,
                                          col          = "black",
                                          modify       = FALSE)
{
  # ============================================================================
  #
  # References:
  #
  # Stucki P. (2010), Methoden zur Untersuchung und Beurteilung der Fliess-
  # gewaesser - Makrozoobenthos Stufe F. Bundesamt fuer Umwelt, Bern. 
  # Umwelt-Vollzug Nr. 1026: 61 S.
  # http://www.modul-stufen-konzept.ch
  #
  # Langhans, S.D. und Reichert, P. (2011), Einbettung von Verfahren zur Fliess-
  # gewaesserbewertung in ein uebergeordnetes Gewaessermanagementkonzept - 
  # Vorschlaege am Beispiel des Modulstufenkonzepts, 
  # Wasser Energie Luft 103(3), 204-214. 
  #
  # ============================================================================
  
  # dictionary for node, attribute and attribute level names:
  # =========================================================
  
  dict <- ecoval.dict(language,dictionaries)
  
  # implementation of three (alternative) indices, default is only IBCH:
  # ====================================================================
  
  ibch <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_invertebrates_ibch",dict),
      name.attrib = ecoval.translate("A_invertebrates_ibch_ind",dict),
      range       = c(0,20), 
      x           = c(0,4.5,8.5,12.5,16.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  ibgn <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_invertebrates_ibgn",dict),
      name.attrib = ecoval.translate("A_invertebrates_ibgn_ind",dict),
      range       = c(0,20),
      x           = c(0,4.5,8.5,12.5,16.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1), 
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  makroind <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_invertebrates_makroind",dict),
      name.attrib = ecoval.translate("A_invertebrates_makroind_ind",dict),
      range       = c(1,8),
      x           = c(1,2.5,3.5,4.5,6.5,8),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
 
  invertebrates <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_invertebrates",dict),
      nodes     = list(ibch),
      name.fun  = "utility.aggregate.add",
      par       = c(1),
      required  = FALSE,
      col       = col)
 
  
  invertebrates.mod <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_invertebrates",dict),
      nodes     = list(ibch,ibgn,makroind),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1,1),
      required  = FALSE,
      col       = col)
  
  if ( modify ) invertebrates <- invertebrates.mod
  
  return(invertebrates)
}
