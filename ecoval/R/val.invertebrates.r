val.invertebrates.create <- function(language     = "English",
                                     dictionaries = NA,
                                     col          = "black",
                                     modify       = TRUE)
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
  
  # branch with macroinvertebrate indices that indicate organic matter pollution (IBCH,IBGN,Makroindex)
  invertebrates.om <- msk.invertebrates.2010.create(language=language,dictionaries=dictionaries,col=col,modify=modify)
    
  # branch with macroinvertebrate indices that indicate pollution with toxic substances (SPEARpesticides)
  invertebrates.spear <- val.spear.create(language=language,dictionaries=dictionaries,col=col)
  
  invertebrates.tox <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_invertebrates_tox",dict),
      nodes     = list(invertebrates.spear),
      name.fun  = "utility.aggregate.add",
      par       = c(1),
      required  = FALSE,
      col       = col)
  
  invertebrates.integrative <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_invertebrates_int",dict),
      nodes     = list(invertebrates.om,invertebrates.tox),
      name.fun  = "utility.aggregate.mix",
      par       = c(1,1,1,1,0),
      names.par = c("w_invertebrates_om","w_invertebrates_tox",
                    "w_invertebrates_int_add","w_invertebrates_int_min","w_invertebrates_int_geo"),
      required  = FALSE,
      col       = col)
  
  return(invertebrates.integrative)
}
