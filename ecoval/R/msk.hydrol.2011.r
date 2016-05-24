msk.hydrol.2011.aggregate <- function(u,par=NA)
{
  # check input:
  
  if ( length(u) != 9 )
  {
    cat("*** Warning: Hydrological aggregation requires 9 attributes")
    return(NA)
  }
  if ( sum(is.na(u)) > 0 ) return(NA)
  
  # calculate aggregated value
  
  u.min <- min(u)
  u.points <- c(0,0.2,0.4,0.6,0.8,1)
  p.points <- c(12,10,6,3,1.5,1)
  #plot(u.points,p.points,type="l") 
  points <- sum(approx(u.points,p.points,u)$y)
  
  u.agg <- u.min   
  if ( u.min < 0.8 )
  {
    if ( u.min >= 0.6 )
    {
      if ( points < 11 ) u.agg <- u.min+0.2
    }
    else
    {
      if ( u.min >= 0.4 )
      {
        if ( points < 13 ) u.agg <- u.min+0.4
        else
        {
          if ( points < 15 ) u.agg <- u.min+0.2
        }
      }
      else
      {
        if ( u.min >= 0.2 )
        {
          if ( points < 17 ) u.agg <- u.min+0.4
          else
          {
            if ( points < 23 ) u.agg <- u.min+0.2
          }
        }
        else
        {
          if ( points < 25 ) u.agg <- u.min+0.4
          else
          {
            if ( points < 31 ) u.agg <- u.min+0.2
          }
        }
      }
    }
  }
  if ( u[7] < u.agg ) u.agg <- u[7]
  if ( sum(u<0.2,na.rm=T) > 1 ) u.agg=0.1
  return(u.agg)
}

msk.hydrol.2011.create <- function(language="English",dictionaries=NA,col="black")
{
  # ============================================================================
  #
  # References:
  #
  # Pfaundler M.,Duebendorfer,C, Zysset, A. (2011): Methoden zur Untersuchung 
  # und Beurteilung der Fliessgewaesser. Hydrologie - Abflussregime Stufe F 
  # (flaechendeckend).
  # Bundesamt fuer Umwelt, Bern. Umwelt-Vollzug Nr. 1107: 113 S.
  #
  # Langhans, S.D. und Reichert, P. (2011), Einbettung von Verfahren zur Fliess-
  # gewaesserbewertung in ein uebergeordnetes Gewaessermanagementkonzept - 
  # Vorschlaege am Beispiel des Modulstufenkonzepts, 
  # Wasser Energie Luft 103(3), 204-214. 
  #
  # ============================================================================
  
  dict <- ecoval.dict(language,dictionaries)
  
  # Attribute Modul Hydrologie (Version Entwurf Oktober 2007)
  # ---------------------------------------------------------
  
  # Mittelwasserabflussverlauf
  
  meandischarge_1 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_1",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_2 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_2",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_3 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_3",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_4 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_4",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_5 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_5",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_6 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_6",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_7 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_7",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(1,31.5,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_8 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_8",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,39.5,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_9 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_9",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,40.4,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_10 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_10",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,42.8,47.7,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_11 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_11",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,44.9,48.2,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_12 <- 
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_12",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,49.1,54.2,60.8,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_13 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_13",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,30,45,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_14 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_14",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,44.5,49.8,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_15 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_15",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,52.8,62,75.1,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  meandischarge_16 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_mean_discharge_16",dict),
      name.attrib = ecoval.translate("A_hydrol_mean_Rrb_percent",dict),
      range       = c(0,100),
      x           = c(0,39.8,45.6,60,85,100),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Mittelwasserabflussverlauf conditional
  
  comb <- data.frame(c(ecoval.translate("L_hydrol_regime_class_1",dict),
                       ecoval.translate("L_hydrol_regime_class_2",dict),
                       ecoval.translate("L_hydrol_regime_class_3",dict),
                       ecoval.translate("L_hydrol_regime_class_4",dict),
                       ecoval.translate("L_hydrol_regime_class_5",dict),
                       ecoval.translate("L_hydrol_regime_class_6",dict),
                       ecoval.translate("L_hydrol_regime_class_7",dict),
                       ecoval.translate("L_hydrol_regime_class_8",dict),
                       ecoval.translate("L_hydrol_regime_class_9",dict),
                       ecoval.translate("L_hydrol_regime_class_10",dict),
                       ecoval.translate("L_hydrol_regime_class_11",dict),
                       ecoval.translate("L_hydrol_regime_class_12",dict),
                       ecoval.translate("L_hydrol_regime_class_13",dict),
                       ecoval.translate("L_hydrol_regime_class_14",dict),
                       ecoval.translate("L_hydrol_regime_class_15",dict),
                       ecoval.translate("L_hydrol_regime_class_16",dict)))
  colnames(comb) <- ecoval.translate("A_hydrol_regime_class",dict)     
  meandischarge <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_hydrol_mean_discharge",dict),
      attrib.levels = comb,
      nodes         = list(meandischarge_1,
                           meandischarge_2,
                           meandischarge_3,
                           meandischarge_4,
                           meandischarge_5,
                           meandischarge_6,
                           meandischarge_7,
                           meandischarge_8,
                           meandischarge_9,
                           meandischarge_10,
                           meandischarge_11,
                           meandischarge_12,
                           meandischarge_13,
                           meandischarge_14,
                           meandischarge_15,
                           meandischarge_16),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Hochwasserh??ufigkeit
  
  floodfrequency_1 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_1",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,2,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_2 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_2",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,2,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_3 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_3",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,2,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_4 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_4",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,2,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_5 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_5",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,2,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_6 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_6",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_7 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_7",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_8 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_8",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,6,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_9 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_9",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,6,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_10 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_10",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_11 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_11",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_12 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_12",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_13 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_13",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,2,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_14 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_14",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_15 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_15",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,3.5,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  floodfrequency_16 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_freq_16",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_freq_pera",dict),
      range       = c(0,20),
      x           = c(0,0.333,0.666,1,6,20),
      u           = c(0,0.2,0.4,0.6,0.8,1.0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  comb <- data.frame(c(ecoval.translate("L_hydrol_regime_class_1",dict),
                       ecoval.translate("L_hydrol_regime_class_2",dict),
                       ecoval.translate("L_hydrol_regime_class_3",dict),
                       ecoval.translate("L_hydrol_regime_class_4",dict),
                       ecoval.translate("L_hydrol_regime_class_5",dict),
                       ecoval.translate("L_hydrol_regime_class_6",dict),
                       ecoval.translate("L_hydrol_regime_class_7",dict),
                       ecoval.translate("L_hydrol_regime_class_8",dict),
                       ecoval.translate("L_hydrol_regime_class_9",dict),
                       ecoval.translate("L_hydrol_regime_class_10",dict),
                       ecoval.translate("L_hydrol_regime_class_11",dict),
                       ecoval.translate("L_hydrol_regime_class_12",dict),
                       ecoval.translate("L_hydrol_regime_class_13",dict),
                       ecoval.translate("L_hydrol_regime_class_14",dict),
                       ecoval.translate("L_hydrol_regime_class_15",dict),
                       ecoval.translate("L_hydrol_regime_class_16",dict)))
  colnames(comb) <- ecoval.translate("A_hydrol_regime_class",dict)     
  floodfrequency <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_hydrol_flood_freq",dict),
      attrib.levels = comb,
      nodes         = list(floodfrequency_1,
                           floodfrequency_2,
                           floodfrequency_3,
                           floodfrequency_4,
                           floodfrequency_5,
                           floodfrequency_6,
                           floodfrequency_7,
                           floodfrequency_8,
                           floodfrequency_9,
                           floodfrequency_10,
                           floodfrequency_11,
                           floodfrequency_12,
                           floodfrequency_13,
                           floodfrequency_14,
                           floodfrequency_15,
                           floodfrequency_16),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Hochwassersaisonalit??t aus Daten
  
  floodseasonality_dat <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_season_dat",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_season_dat_ind",dict),
      range       = c(0,2),
      x           = c(0,0.3,0.6,0.9,1.2,2),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Hochwassersaisonalit??t aus Referenzellipsen       
  
  floodseasonality_ref <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_flood_season_ref",dict),
      name.attrib = ecoval.translate("A_hydrol_flood_season_ref_ind",dict),
      range       = c(0,2),
      x           = c(0,0.25,0.5,0.75,1,2),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Hochwassersaisonalit??t
  
  floodseasonality <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_hydrol_flood_seasonality",dict),
      nodes     = list(floodseasonality_dat,
                       floodseasonality_ref),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Niedrigwasserabfluss
  
  Q347rmax    <- 10000
  CV_Q347r_1  <- 22
  CV_Q347r_2  <- 19
  CV_Q347r_3  <- 19
  CV_Q347r_4  <- 13
  CV_Q347r_5  <- 18
  CV_Q347r_6  <- 22
  CV_Q347r_7  <- 30
  CV_Q347r_8  <- 35
  CV_Q347r_9  <- 38
  CV_Q347r_10 <- 30
  CV_Q347r_11 <- 37
  CV_Q347r_12 <- 38
  CV_Q347r_13 <- 19
  CV_Q347r_14 <- 38
  CV_Q347r_15 <- 34
  CV_Q347r_16 <- 21
  
  # Niedrigwasser, absolute Kriterien an relative Abweichung
  
  lowflow_discharge_reldev <-
    utility.endnode.intpol2d.create(
      name.node     = ecoval.translate("N_hydrol_lowflow_discharge_reldev",dict),
      name.attrib   = c(ecoval.translate("A_hydrol_lowflow_Q347r_lpers",dict),
                        ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict)),
      ranges        = list(c(0,Q347rmax),c(0,100)),
      isolines      = list(list(x=c(0,Q347rmax),y=c(0,0)),
                           list(x=c(0,50,200,500,1000,Q347rmax),y=c(20,20,25,35,45,45)),
                           list(x=c(0,50,200,500,1000,Q347rmax),y=c(40,40,45,55,65,65)),
                           list(x=c(0,50,200,500,1000,Q347rmax),y=c(65,65,70,80,85,85)),
                           list(x=c(0,Q347rmax),y=c(100,100))),
      u             = c(1,0.6,0.4,0.2,0),
      lead          = 1,
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Niedrigwasser, Vergleich von relativer Abweichung mit nat??rlicher Variabilit??t
  
  lowflow_discharge_relvar_1 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_1",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_1),
      x           = c(0,CV_Q347r_1),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_2 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_2",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_2),
      x           = c(0,CV_Q347r_2),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_3 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_3",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_3),
      x           = c(0,CV_Q347r_3),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_4 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_4",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_4),
      x           = c(0,CV_Q347r_4),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_5 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_5",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_5),
      x           = c(0,CV_Q347r_5),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_6 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_6",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_6),
      x           = c(0,CV_Q347r_6),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_7 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_7",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_7),
      x           = c(0,CV_Q347r_7),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_8 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_8",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_8),
      x           = c(0,CV_Q347r_8),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_9 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_9",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_9),
      x           = c(0,CV_Q347r_9),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_10 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_10",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_10),
      x           = c(0,CV_Q347r_10),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_11 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_11",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_11),
      x           = c(0,CV_Q347r_11),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_12 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_12",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_12),
      x           = c(0,CV_Q347r_12),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_13 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_13",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_13),
      x           = c(0,CV_Q347r_13),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_14 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_14",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_14),
      x           = c(0,CV_Q347r_14),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_15 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_15",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_15),
      x           = c(0,CV_Q347r_15),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge_relvar_16 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_discharge_relvar_16",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_deltaQ347rb_percent",dict),
      range       = c(0,CV_Q347r_16),
      x           = c(0,CV_Q347r_16),
      u           = c(1,0.8),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  comb <- data.frame(c(ecoval.translate("L_hydrol_regime_class_1",dict),
                       ecoval.translate("L_hydrol_regime_class_2",dict),
                       ecoval.translate("L_hydrol_regime_class_3",dict),
                       ecoval.translate("L_hydrol_regime_class_4",dict),
                       ecoval.translate("L_hydrol_regime_class_5",dict),
                       ecoval.translate("L_hydrol_regime_class_6",dict),
                       ecoval.translate("L_hydrol_regime_class_7",dict),
                       ecoval.translate("L_hydrol_regime_class_8",dict),
                       ecoval.translate("L_hydrol_regime_class_9",dict),
                       ecoval.translate("L_hydrol_regime_class_10",dict),
                       ecoval.translate("L_hydrol_regime_class_11",dict),
                       ecoval.translate("L_hydrol_regime_class_12",dict),
                       ecoval.translate("L_hydrol_regime_class_13",dict),
                       ecoval.translate("L_hydrol_regime_class_14",dict),
                       ecoval.translate("L_hydrol_regime_class_15",dict),
                       ecoval.translate("L_hydrol_regime_class_16",dict)))
  colnames(comb) <- ecoval.translate("A_hydrol_regime_class",dict)     
  lowflow_discharge_relvar <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_hydrol_lowflow_discharge_relvar",dict),
      attrib.levels = comb,
      nodes         = list(lowflow_discharge_relvar_1,
                           lowflow_discharge_relvar_2,
                           lowflow_discharge_relvar_3,
                           lowflow_discharge_relvar_4,
                           lowflow_discharge_relvar_5,
                           lowflow_discharge_relvar_6,
                           lowflow_discharge_relvar_7,
                           lowflow_discharge_relvar_8,
                           lowflow_discharge_relvar_9,
                           lowflow_discharge_relvar_10,
                           lowflow_discharge_relvar_11,
                           lowflow_discharge_relvar_12,
                           lowflow_discharge_relvar_13,
                           lowflow_discharge_relvar_14,
                           lowflow_discharge_relvar_15,
                           lowflow_discharge_relvar_16),
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  lowflow_discharge <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_hydrol_lowflow_discharge",dict),
      nodes     = list(lowflow_discharge_reldev,lowflow_discharge_relvar),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Niederwassersaisonalit??t aus Daten
  
  lowflow_seasonality_dat <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_season_dat",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_season_dat_ind",dict),
      range       = c(0,2),
      x           = c(0,0.3,0.6,0.9,1.2,2),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Niederwassersaisonalit??t aus Referenzellipsen  
  
  lowflow_seasonality_ref <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_season_ref",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_season_ref_ind",dict),
      range       = c(0,2),
      x           = c(0,0.25,0.5,0.75,1,2),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Niedrigwassersaisonalit??t
  
  lowflow_seasonality <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_hydrol_lowflow_seasonality",dict),
      nodes     = list(lowflow_seasonality_dat,
                       lowflow_seasonality_ref),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Dauer von Niederwasserperioden 
  
  lowflow_duration <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_lowflow_duration",dict),
      name.attrib = ecoval.translate("A_hydrol_lowflow_duration_d",dict),
      range       = c(0,365),
      x           = c(0,20,35,50,65,365),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Schwall-Sunk (2-dim)
  
  hydropeaking <-
    utility.endnode.intpol2d.create(
      name.node     = ecoval.translate("N_hydrol_hydropeaking",dict),
      name.attrib   = c(ecoval.translate("A_hydrol_hydropeaking_stress_ind",dict),
                        ecoval.translate("A_hydrol_hydropeaking_intensity_ind",dict)),
      ranges        = list(c(0,3.25),c(1,7)),
      isolines      = list(list(x=c(0,3.25),y=c(1,1)),
                           list(x=c(0,0.25,1.75,3.25),y=c(2,2,1.36,1.27)),
                           list(x=c(0,0.25,1.75,3.25),y=c(3,3,2.18,1.36)),
                           list(x=c(0,0.25,1.75,3.25),y=c(4.27,4.27,3,1.82)),
                           list(x=c(0,0.25,1.75,3.25),y=c(6,6,4,2.18)),
                           list(x=c(0,3.25),y=c(7,7))),
      u             = c(1,0.8,0.6,0.4,0.2,0),
      lead          = 0,
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Sp??lungen und Entleerungen (2-dim)
  
  flushings <-
    utility.endnode.intpol2d.create(
      name.node     = ecoval.translate("N_hydrol_flushings",dict),
      name.attrib   = c(ecoval.translate("A_hydrol_flushings_stress_ind",dict),
                        ecoval.translate("A_hydrol_flushings_freq_pera",dict)),
      ranges        = list(c(0.1,100),c(0.1,60)),
      isolines      = list(list(x=c(0.1,0.1),y=c(34,0.1)),
                           list(x=c(0.1,0.5,4,9.5,15.8),y=c(34,15,2,0.4,0.1)),
                           list(x=c(0.1,0.5,5,13,23),y=c(57,27,3,0.5,0.1)),
                           list(x=c(0.3,0.5,6,17,33),y=c(60,48,5,0.7,0.1)),
                           list(x=c(0.9,7,24,50),y=c(60,9,0.9,0.1)),
                           list(x=c(100,100),y=c(0.1,100))),
      u             = c(1,0.8,0.6,0.4,0.2,0),
      lead          = 0,
      required      = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Hochwasser durch Regenwassereinleitungen   
  
  rainrunofffloods <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_hydrol_raininput",dict),
      name.attrib = ecoval.translate("A_hydrol_raininput_freq_pera",dict),
      range       = c(0,10),
      x           = c(0,0.2,2,4,8,10),
      u           = c(1,0.8,0.6,0.4,0.2,0),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Modul Hydrologie hierarchische Aggregation
  
  hydrol <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_hydrol",dict),
      nodes     = list(meandischarge,
                       floodfrequency,
                       floodseasonality,
                       lowflow_discharge,
                       lowflow_seasonality,
                       lowflow_duration,
                       hydropeaking,
                       flushings,
                       rainrunofffloods),
      name.fun  = "msk.hydrol.2011.aggregate",
      par       = c(1,1,1,1,1,1,1,1,1),
      required  = FALSE,
      col       = col)
  
  return(hydrol)
}
