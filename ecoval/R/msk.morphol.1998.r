# ==========================================
# Create value function for river morphology
# ==========================================

msk.morphol.1998.aggregate <- function(u,par=NA)
{
  if ( length(u) != 2 )
  {
    cat("*** Warning: Morphological aggregation requires 2 values, coverage and uncovered.")
    return(NA)
  }
  
  if ( !is.na(u[1]) & u[1]==0 ) return(0)
  
  return(u[2])
}

msk.morphol.1998.create <- function(language     = "English",
                                    dictionaries = NA,
                                    col          = "black")
{
  # ============================================================================
  #
  # References:
  #
  # Huette, M. and Niederhauser P. (1998), Methoden zur Untersuchung und Beur-
  # teilung der Fliessgewaesser in der Schweiz: Oekomorphologie Stufe F,
  # Mitteilungen zum Gewaesserschutz Nr. 27
  # Bundesamt fuer Umwelt, Wald und Landschaft, BUWAL, Bern
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
  
  # translation of points into values:
  # ==================================
  
  # function to convert assessment points (BUWAL 1998, page 33) into values 
  # (BUWAL 1998, page 34)
  # (value class boundaries (Langhans and Reichert 2011, Fig. 4B): 
  #  0.0-0.3, 0.3-0.6, 0.6-0.8, 0.8-1.0 )
  
  pnt2val <- function(p)
  {
    return(approx(x    = c(0.0,1.5,5.5,9.5,12.0)/4,
                  y    = c(1.0,0.8,0.6,0.3, 0.0),
                  xout = p)$y)
  }
  
  # construction of end nodes:
  # ==========================
  
  # node "Wasserspiegelbreitenvariabilitaet":
  # -----------------------------------------
  
  # Wasserspiegelbreitenvariabilitaet (discrete attribute):
  
  comb <- data.frame(c(ecoval.translate("L_morphol_widthvar_class_high",dict),
                       ecoval.translate("L_morphol_widthvar_class_moderate",dict),
                       ecoval.translate("L_morphol_widthvar_class_none",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_widthvar_class",dict)
  widthvar <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_widthvar",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(2),
                        pnt2val(3)),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  # node "Sohlenverbauung":
  # -----------------------
  
  # Sohlenverbauung definitions:
  
  pnt0          <- 0 
  pnt10         <- 1.5    # midvalues between points
  pnt30_other   <- 2.5    # midvalues between points
  pnt100_riprap <- 2.375  # upper end of the yellow class
  pnt100_other  <- 3.0
  
  # Sohlenverbauung Steinschuettung, Rauhbett (attribute range 0 - 100):
  
  bedmod_cont_riprap <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_morphol_bedmod_cont_riprap",dict),
      name.attrib = ecoval.translate("A_morphol_bedmod_fract_percent",dict),
      range       = c(0,100),
      x           = c(0,10,100),
      u           = c(pnt2val(pnt0),
                      pnt2val(pnt10),
                      pnt2val(pnt100_riprap)),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Sohlenverbauung andere Materialien (attribute range 0 - 100):
  
  bedmod_cont_other <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_morphol_bedmod_cont_other",dict),
      name.attrib = ecoval.translate("A_morphol_bedmod_fract_percent",dict),
      range       = c(0,100),
      x           = c(0,10,30,100),
      u           = c(pnt2val(pnt0),
                      pnt2val(pnt10),
                      pnt2val(pnt30_other),
                      pnt2val(pnt100_other)),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Sohlenverbauung (attribute range 0 - 100):
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bedmod_type_class_riprap1",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap2",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap3",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap4",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap5",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other1",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other2",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other3",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other4",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other5",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bedmod_type_class",dict)
  bedmod_cont <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_bedmod_cont",dict),
      attrib.levels = comb,
      nodes         = list(bedmod_cont_riprap,bedmod_cont_riprap,bedmod_cont_riprap,bedmod_cont_riprap,bedmod_cont_riprap,
                           bedmod_cont_other,bedmod_cont_other,bedmod_cont_other,bedmod_cont_other,bedmod_cont_other),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Sohlenverbauung rip-rap (discrete attribute):
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bedmod_fract_class_0",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_0to10",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_10to30",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_30to60",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_30to100",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_60to100",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_100",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bedmod_fract_class",dict)
  bedmod_discr_riprap <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_bedmod_discr_riprap",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),      #0
                        pnt2val(1),      #0-10
                        pnt2val(2),      #10-30
                        pnt2val(2),      #30-60
                        pnt2val(2),      #30-100
                        pnt2val(2),      #60-100
                        pnt2val(2)),     #100
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Sohlenverbauung other (discrete attribute):
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bedmod_fract_class_0",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_0to10",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_10to30",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_30to60",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_30to100",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_60to100",dict),
                       ecoval.translate("L_morphol_bedmod_fract_class_100",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bedmod_fract_class",dict)
  bedmod_discr_other <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_bedmod_discr_other",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),  #0
                        pnt2val(1),  #0-10
                        pnt2val(2),  #10-30
                        pnt2val(3),  #30-60
                        pnt2val(3),  #30-100
                        pnt2val(3),  #60-100
                        pnt2val(3)), #100
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Sohlenverbauung (discrete attribute):
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bedmod_type_class_riprap1",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap2",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap3",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap4",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_riprap5",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other1",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other2",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other3",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other4",dict),
                       ecoval.translate("L_morphol_bedmod_type_class_other5",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bedmod_type_class",dict)
  bedmod_discr <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_bedmod_discr",dict),
      attrib.levels = comb,
      nodes         = list(bedmod_discr_riprap,bedmod_discr_riprap,bedmod_discr_riprap,bedmod_discr_riprap,bedmod_discr_riprap,
                           bedmod_discr_other,bedmod_discr_other,bedmod_discr_other,bedmod_discr_other,bedmod_discr_other),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Sohlenverbauung:
  
  bedmod <-
    utility.endnode.firstavail.create(
      name.node = ecoval.translate("N_morphol_bedmod",dict),
      nodes     = list(bedmod_cont,bedmod_discr),
      required  = TRUE,
      utility   = FALSE,
      col       = col)
  
  # node "Boeschungsfussverbauung":
  # -------------------------------
  
  # Boeschungsfussverbauung definitions:
  
  pnt0_perm     <- 0 
  pnt10_perm    <- 0.25
  pnt30_perm    <- 1.0
  pnt60_perm    <- 2.0
  pnt100_perm   <- 2.5
  
  pnt0_imperm   <- 0 
  pnt10_imperm  <- 0.5
  pnt30_imperm  <- 1.5
  pnt60_imperm  <- 2.5
  pnt100_imperm <- 3.0
  
  # Boeschungsfussverbauung left permeable (attribute range 0 - 100)
  
  bankmod_left_cont_perm <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_morphol_bankmod_left_cont_perm",dict),
      name.attrib = ecoval.translate("A_morphol_bankmod_fract_left_percent",dict),
      range       = c(0,100),
      x           = c(0,10,30,60,100),
      u           = c(pnt2val(pnt0_perm),
                      pnt2val(pnt10_perm),
                      pnt2val(pnt30_perm),
                      pnt2val(pnt60_perm),
                      pnt2val(pnt100_perm)),
      required    = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung right permeable (attribute range 0 - 100)
  
  bankmod_right_cont_perm <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_morphol_bankmod_right_cont_perm",dict),
      name.attrib = ecoval.translate("A_morphol_bankmod_fract_right_percent",dict),
      range       = c(0,100),
      x           = c(0,10,30,60,100),
      u           = c(pnt2val(pnt0_perm),
                      pnt2val(pnt10_perm),
                      pnt2val(pnt30_perm),
                      pnt2val(pnt60_perm),
                      pnt2val(pnt100_perm)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung left permeable (discrete attribute range)
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_fract_class_0",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_0to10",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_10to30",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_30to60",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_60to100",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_100",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_fract_left_class",dict)
  bankmod_left_discr_perm <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_bankmod_left_discr_perm",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(0),
                        pnt2val(0.5),
                        pnt2val(1.5),
                        pnt2val(2.5),
                        pnt2val(2.5)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung right permeable (discrete attribute range)
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_fract_class_0",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_0to10",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_10to30",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_30to60",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_60to100",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_100",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_fract_right_class",dict)
  bankmod_right_discr_perm <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_bankmod_right_discr_perm",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(0),
                        pnt2val(0.5),
                        pnt2val(1.5),
                        pnt2val(2.5),
                        pnt2val(2.5)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung left impermeable (continuous, attribute range 0 - 100)
  
  bankmod_left_cont_imperm <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_morphol_bankmod_left_cont_imperm",dict),
      name.attrib = ecoval.translate("A_morphol_bankmod_fract_left_percent",dict),
      range       = c(0,100),
      x           = c(0,10,30,60,100),
      u           = c(pnt2val(pnt0_imperm),
                      pnt2val(pnt10_imperm),
                      pnt2val(pnt30_imperm),
                      pnt2val(pnt60_imperm),
                      pnt2val(pnt100_imperm)),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Boeschungsfussverbauung right impermeable (continous, attribute range 0 - 100)
  
  bankmod_right_cont_imperm <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_morphol_bankmod_right_cont_imperm",dict),
      name.attrib = ecoval.translate("A_morphol_bankmod_fract_right_percent",dict),
      range       = c(0,100),
      x           = c(0,10,30,60,100),
      u           = c(pnt2val(pnt0_imperm),
                      pnt2val(pnt10_imperm),
                      pnt2val(pnt30_imperm),
                      pnt2val(pnt60_imperm),
                      pnt2val(pnt100_imperm)),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Boeschungsfussverbauung left impermeable (discrete attribute range)
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_fract_class_0",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_0to10",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_10to30",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_30to60",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_60to100",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_100",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_fract_left_class",dict)
  bankmod_left_discr_imperm <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_bankmod_left_discr_imperm",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(0),
                        pnt2val(1),
                        pnt2val(2),
                        pnt2val(3),
                        pnt2val(3)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung right impermeable (discrete attribute range)
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_fract_class_0",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_0to10",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_10to30",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_30to60",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_60to100",dict),
                       ecoval.translate("L_morphol_bankmod_fract_class_100",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_fract_right_class",dict)
  bankmod_right_discr_imperm <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_bankmod_right_discr_imperm",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(0),
                        pnt2val(1),
                        pnt2val(2),
                        pnt2val(3),
                        pnt2val(3)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung left cont
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_perm_class_perm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm5",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm5",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_perm_left_class",dict)
  bankmod_left_cont <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_bankmod_left_cont",dict),
      attrib.levels = comb,
      nodes         = list(bankmod_left_cont_perm,
                           bankmod_left_cont_perm,
                           bankmod_left_cont_perm,
                           bankmod_left_cont_perm,
                           bankmod_left_cont_perm,
                           bankmod_left_cont_imperm,
                           bankmod_left_cont_imperm,
                           bankmod_left_cont_imperm,
                           bankmod_left_cont_imperm,
                           bankmod_left_cont_imperm),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung right cont
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_perm_class_perm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm5",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm5",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_perm_right_class",dict)
  bankmod_right_cont <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_bankmod_right_cont",dict),
      attrib.levels = comb,
      nodes         = list(bankmod_right_cont_perm,
                           bankmod_right_cont_perm,
                           bankmod_right_cont_perm,
                           bankmod_right_cont_perm,
                           bankmod_right_cont_perm,
                           bankmod_right_cont_imperm,
                           bankmod_right_cont_imperm,
                           bankmod_right_cont_imperm,
                           bankmod_right_cont_imperm,
                           bankmod_right_cont_imperm),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung left discr
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_perm_class_perm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm5",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm5",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_perm_left_class",dict)
  bankmod_left_discr <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_bankmod_left_discr",dict),
      attrib.levels = comb,
      nodes         = list(bankmod_left_discr_perm,
                           bankmod_left_discr_perm,
                           bankmod_left_discr_perm,
                           bankmod_left_discr_perm,
                           bankmod_left_discr_perm,
                           bankmod_left_discr_imperm,
                           bankmod_left_discr_imperm,
                           bankmod_left_discr_imperm,
                           bankmod_left_discr_imperm,
                           bankmod_left_discr_imperm),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung right discr
  
  comb <- data.frame(c(ecoval.translate("L_morphol_bankmod_perm_class_perm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_perm5",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm1",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm2",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm3",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm4",dict),
                       ecoval.translate("L_morphol_bankmod_perm_class_imperm5",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_bankmod_perm_right_class",dict)
  bankmod_right_discr <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_bankmod_right_discr",dict),
      attrib.levels = comb,
      nodes         = list(bankmod_right_discr_perm,
                           bankmod_right_discr_perm,
                           bankmod_right_discr_perm,
                           bankmod_right_discr_perm,
                           bankmod_right_discr_perm,
                           bankmod_right_discr_imperm,
                           bankmod_right_discr_imperm,
                           bankmod_right_discr_imperm,
                           bankmod_right_discr_imperm,
                           bankmod_right_discr_imperm),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Boeschungsfussverbauung left
  
  bankmod_left <-
    utility.endnode.firstavail.create(
      name.node = ecoval.translate("N_morphol_bankmod_left",dict),
      nodes     = list(bankmod_left_cont,bankmod_left_discr),
      required  = FALSE,
      utility   = FALSE,
      col       = col)
  
  # Boeschungsfussverbauung right
  
  bankmod_right <-
    utility.endnode.firstavail.create(
      name.node = ecoval.translate("N_morphol_bankmod_right",dict),
      nodes     = list(bankmod_right_cont,bankmod_right_discr),
      required  = FALSE,
      utility   = FALSE,
      col       = col)
  
  # Boeschungsfussverbauung
  
  bankmod <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_morphol_bankmod",dict),
      nodes     = list(bankmod_left,bankmod_right),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1),
      required  = TRUE,
      col       = col)
  
  # node "Uferbereich":
  # -------------------
  
  # Uferbereichsbreite left width widthvar=high
  
  riparzone_left_width_highvar <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_morphol_riparzone_left_width_highvar",dict),
      name.attrib = c(ecoval.translate("A_morphol_riverbed_width_m",dict),
                      ecoval.translate("A_morphol_riparzone_width_left_m",dict)),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,2,15,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead        = 1,
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Uferbereichsbreite left width widthvar=moderate
  
  riparzone_left_width_moderatevar <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_morphol_riparzone_left_width_moderatevar",dict),
      name.attrib = c(ecoval.translate("A_morphol_riverbed_width_m",dict),
                      ecoval.translate("A_morphol_riparzone_width_left_m",dict)),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,1.5,10,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead        = 1,
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Uferbereichsbreite left width widthvar=none
  
  riparzone_left_width_novar <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_morphol_riparzone_left_width_novar",dict),
      name.attrib = c(ecoval.translate("A_morphol_riverbed_width_m",dict),
                      ecoval.translate("A_morphol_riparzone_width_left_m",dict)),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,1,7.5,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead        = 1,
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Uferbereichsbreite right width widthvar=high
  
  riparzone_right_width_highvar <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_morphol_riparzone_right_width_highvar",dict),
      name.attrib = c(ecoval.translate("A_morphol_riverbed_width_m",dict),
                      ecoval.translate("A_morphol_riparzone_width_right_m",dict)),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,2,15,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead        = 1,
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Uferbereichsbreite right width widthvar=moderate
  
  riparzone_right_width_moderatevar <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_morphol_riparzone_right_width_moderatevar",dict),
      name.attrib = c(ecoval.translate("A_morphol_riverbed_width_m",dict),
                      ecoval.translate("A_morphol_riparzone_width_right_m",dict)),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,1.5,10,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead         = 1,
      utility       = FALSE,
      required    = FALSE,
      col           = col)
  
  # Uferbereichsbreite right width widthvar=none
  
  riparzone_right_width_novar <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_morphol_riparzone_right_width_novar",dict),
      name.attrib = c(ecoval.translate("A_morphol_riverbed_width_m",dict),
                      ecoval.translate("A_morphol_riparzone_width_right_m",dict)),
      ranges      = list(c(0,16),c(0,30)),
      isolines    = list(list(x=c(0,16),y=c(0,0)),
                         list(x=c(0,1,7.5,16),y=c(5,5,15,15)),
                         list(x=c(0,16),y=c(15,15)),
                         list(x=c(0,16),y=c(30,30))),
      u           = c(0.0,0.6,1.0,1.0),
      lead        = 1,
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Uferbereichsbreite left width
  
  comb <- data.frame(c(ecoval.translate("L_morphol_widthvar_class_high",dict),
                       ecoval.translate("L_morphol_widthvar_class_moderate",dict),
                       ecoval.translate("L_morphol_widthvar_class_none",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_widthvar_class",dict)
  riparzone_left_width <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_riparzone_left_width",dict),
      attrib.levels = comb,
      nodes         = list(riparzone_left_width_highvar,
                           riparzone_left_width_moderatevar,
                           riparzone_left_width_novar),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Uferbereichsbreite right width
  
  comb <- data.frame(c(ecoval.translate("L_morphol_widthvar_class_high",dict),
                       ecoval.translate("L_morphol_widthvar_class_moderate",dict),
                       ecoval.translate("L_morphol_widthvar_class_none",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_widthvar_class",dict)
  riparzone_right_width <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_morphol_riparzone_right_width",dict),
      attrib.levels = comb,
      nodes         = list(riparzone_right_width_highvar,
                           riparzone_right_width_moderatevar,
                           riparzone_right_width_novar),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Uferbereichsvegetation left (discrete attribute)
  
  comb <-  data.frame(c(ecoval.translate("L_morphol_riparzone_veg_class_natural",dict),
                        ecoval.translate("L_morphol_riparzone_veg_class_seminatural",dict),
                        ecoval.translate("L_morphol_riparzone_veg_class_artificial",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_riparzone_veg_left_class",dict)
  riparzone_left_veg <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_riparzone_left_veg",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),pnt2val(1.5),pnt2val(3)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Uferbereichsvegetation right (discrete attribute)
  
  comb <-  data.frame(c(ecoval.translate("L_morphol_riparzone_veg_class_natural",dict),
                        ecoval.translate("L_morphol_riparzone_veg_class_seminatural",dict),
                        ecoval.translate("L_morphol_riparzone_veg_class_artificial",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_riparzone_veg_right_class",dict)
  riparzone_right_veg <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_riparzone_right_veg",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),pnt2val(1.5),pnt2val(3)),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Uferbereich left
  
  riparzone_left <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_morphol_riparzone_left",dict),
      nodes     = list(riparzone_left_width,riparzone_left_veg),
      name.fun  = "utility.aggregate.addmin",
      par       = c(0.9,0.1,0.5),
      required  = FALSE,
      col       = col)
  
  # Uferbereich right
  
  riparzone_right <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_morphol_riparzone_right",dict),
      nodes     = list(riparzone_right_width,riparzone_right_veg),
      name.fun  = "utility.aggregate.addmin",
      par       = c(0.9,0.1,0.5),
      required  = FALSE,
      col       = col)
  
  # Uferbereich
  
  riparzone <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_morphol_riparzone",dict),
      nodes     = list(riparzone_left,riparzone_right),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1),
      required  = TRUE,
      col       = col)
  
  # Oekomorphologie
  
  comb <-  data.frame(c(ecoval.translate("L_morphol_coverage_class_covered",dict),
                        ecoval.translate("L_morphol_coverage_class_uncovered",dict)))
  colnames(comb) <- ecoval.translate("A_morphol_coverage_class",dict)
  coverage <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_morphol_coverage",dict),
      attrib.levels = comb,
      u             = c(0,1),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  uncovered <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_morphol_uncovered",dict),
      nodes     = list(widthvar,bedmod,bankmod,riparzone),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1,1,1),
      required  = FALSE,
      col       = col)
  
  morphol <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_morphol",dict),
      nodes     = list(coverage,uncovered),
      name.fun  = "msk.morphol.1998.aggregate",
      par       = NA,
      required  = TRUE,
      col       = col)
  
  return(morphol)
}
