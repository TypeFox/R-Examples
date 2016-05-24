msk.diatoms.2007.create <- function(language     = "English",
                                    dictionaries = NA,
                                    col          = "black")
{
  # ============================================================================
  #
  # References:
  #
  # Huerlimann J., Niederhauser P. (2007), Methoden zur Untersuchung und Beur-
  # teilung der Fliessgewaesser. Kieselalgen Stufe F (flaechendeckend). 
  # Umwelt-Vollzug Nr. 0740. Bundesamt fuer Umwelt, Bern. 130 S
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
  
  # Valuation as a function of the index "diatindex":
  # =================================================
  
  # Huerlimann and Niederhauser (2007), p. 20
  
  diatindex <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_diatoms_diatindex",dict),
      name.attrib = ecoval.translate("A_diatoms_diatindex_ind",dict),
      range       = c(1,8),
      x           = c(1,3.5,4.5,5.5,6.5,8),
      u           = c(1,0.8,0.6,0.4,0.2,0.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  diatoms <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_diatoms",dict),
      nodes     = list(diatindex),
      name.fun  = "utility.aggregate.add",
      par       = c(1),
      required  = FALSE,
      col       = col)
  
  return(diatoms)
}
