msk.nutrients.2010.create <- function(language="English",dictionaries=NA,col="black",modify=F)
{
  # ============================================================================
  #
  # References:
  #
  # Liechti Paul (2010), Methoden zur Untersuchung und Beurteilung der
  # Fliessgewaesser. Chemisch-physikalische Erhebungen, Naehrstoffe.
  # Umwelt-Vollzug Nr. 1005. Bundesamt fuer Umwelt, Bern. 44 S.
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
  
  # construction of nodes for individual chemicals:
  # ===============================================
  
  # valuation parameters according to the pages 18/19 of Liechti (2010):
  # --------------------------------------------------------------------
  
  Z.TP              <-   0.07
  min.TP            <-   0.005
  max.TP            <-   1.5
  
  Z.TP.filt         <-   0.05   # note two exceptions from rule below!
  min.TP.filt       <-   0.005
  max.TP.filt       <-   1.0
  
  Z.PO4             <-   0.04
  min.PO4           <-   0.005
  max.PO4           <-   1.0
  
  Z.TN              <-   7.0    # note exception from rule below
  min.TN            <-   0.5
  max.TN            <-  20
  
  Z.NO3             <-   5.6    # note exception from rule below
  min.NO3           <-   0.25
  max.NO3           <-  15
  
  Z.NO2.Clsmaller10 <-   0.02
  Z.NO2.Cl10to20    <-   0.05
  Z.NO2.Cllarger20  <-   0.10
  min.NO2           <-   0
  max.NO2           <-   0.5
  
  max.ClforNO2      <- 200
  
  Z.NH4.Tsmaller10  <-   0.4    # note exception from rule below
  Z.NH4.Tlarger10   <-   0.2
  min.NH4           <-   0
  max.NH4           <-   2.0
  
  translow.TforNH4  <-   7.5
  transhigh.TforNH4 <-  12.5
  max.TforNH4       <-  25
  
  Z.TOC             <-   5
  min.TOC           <-   0.5
  max.TOC           <-  15
  
  Z.DOC             <-   4
  min.DOC           <-   0.5
  max.DOC           <-  12
  
  Z.BOD5            <-   4
  min.BOD5          <-   1
  max.BOD5          <-  10
  
  # implementation of valuation nodes:
  # ----------------------------------
  
  # Gesamt-P unfiltriert
  
  TP <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_TP",dict),
      name.attrib = ecoval.translate("A_nutrients_TP_mgPperl",dict),
      range       = c(0,max.TP),
      x           = c(0  ,min.TP,4/7*Z.TP,Z.TP,10/7*Z.TP,2*Z.TP,max.TP),         # two exceptions
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Gesamt-P filtriert
  
  TPfilt <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_TPfilt",dict),
      name.attrib = ecoval.translate("A_nutrients_TPfilt_mgPperl",dict),
      range       = c(0,max.TP.filt),
      x           = c(0  ,min.TP.filt,0.5*Z.TP.filt,Z.TP.filt,1.5*Z.TP.filt,2*Z.TP.filt,max.TP.filt),
      u           = c(1.0,    1.0    ,     0.8     ,   0.6   ,     0.4     ,    0.2    ,    0.0    ),
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Ortho-Phosphat 
  
  PO4 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_PO4",dict),
      name.attrib = ecoval.translate("A_nutrients_PO4_mgPperl",dict),
      range       = c(0,max.PO4),
      x           = c(0  ,min.PO4,0.5*Z.PO4,Z.PO4,1.5*Z.PO4,2*Z.PO4,max.PO4),
      u           = c(1.0,  1.0  ,   0.8   , 0.6 ,   0.4   ,  0.2  ,  0.0  ),
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  # Gesamt-N  
  
  TN <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_TN",dict),
      name.attrib = ecoval.translate("A_nutrients_TN_mgNperl",dict),
      range       = c(0,max.TN),
      x           = c(0  ,min.TN,2,Z.TN,1.5*Z.TN,2*Z.TN,max.TN),                   # exception
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Nitrat  
  
  NO3 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_NO3",dict),
      name.attrib = ecoval.translate("A_nutrients_NO3_mgNperl",dict),
      range       = c(0,max.NO3),
      x           = c(0  ,min.NO3,1.5/5.6*Z.NO3,Z.NO3,1.5*Z.NO3,2*Z.NO3,max.NO3),     # exception
      u           = c(1.0,  1.0  ,   0.8   , 0.6 ,   0.4   ,  0.2  ,  0.0  ),
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  # Nitrit
  
  NO2orig <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_nutrients_NO2",dict),
      name.attrib = c(ecoval.translate("A_nutrients_ClforNO2_mgperl",dict),
                      ecoval.translate("A_nutrients_NO2_mgNperl",dict)),
      ranges      = list(c(0,max.ClforNO2),c(0,max.NO2)),
      isolines    = list(list(x=c(0,max.ClforNO2),
                              y=c(0,0)),
                         list(x=c(0,9.95,10.05,19.95,20.05,max.ClforNO2),
                              y=c(0.5*Z.NO2.Clsmaller10,0.5*Z.NO2.Clsmaller10,
                                  2/5*Z.NO2.Cl10to20   ,2/5*Z.NO2.Cl10to20,             # exception
                                  0.5*Z.NO2.Cllarger20 ,0.5*Z.NO2.Cllarger20)),
                         list(x=c(0,9.95,10.05,19.95,20.05,max.ClforNO2),
                              y=c(1.0*Z.NO2.Clsmaller10,1.0*Z.NO2.Clsmaller10,
                                  1.0*Z.NO2.Cl10to20   ,1.0*Z.NO2.Cl10to20,
                                  1.0*Z.NO2.Cllarger20 ,1.0*Z.NO2.Cllarger20)),
                         list(x=c(0,9.95,10.05,19.95,20.05,max.ClforNO2),
                              y=c(1.5*Z.NO2.Clsmaller10,1.5*Z.NO2.Clsmaller10,
                                  1.5*Z.NO2.Cl10to20   ,1.5*Z.NO2.Cl10to20,
                                  1.5*Z.NO2.Cllarger20 ,1.5*Z.NO2.Cllarger20)),
                         list(x=c(0,9.95,10.05,19.95,20.05,max.ClforNO2),
                              y=c(2.0*Z.NO2.Clsmaller10,2.0*Z.NO2.Clsmaller10,
                                  2.0*Z.NO2.Cl10to20   ,2.0*Z.NO2.Cl10to20,
                                  2.0*Z.NO2.Cllarger20 ,2.0*Z.NO2.Cllarger20)),
                         list(x=c(0,max.ClforNO2),
                              y=c(max.NO2,max.NO2))),
      u           = c(1.0,0.8,0.6,0.4,0.2,0.0),
      lead        = 1,
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  NO2modif <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_nutrients_NO2",dict),
      name.attrib = c(ecoval.translate("A_nutrients_ClforNO2_mgperl",dict),
                      ecoval.translate("A_nutrients_NO2_mgNperl",dict)),
      ranges      = list(c(0,max.ClforNO2),c(0,max.NO2)),
      isolines    = list(list(x=c(0,max.ClforNO2),
                              y=c(0,0)),
                         list(x=c(0,10,20,max.ClforNO2),
                              y=c(0.5*Z.NO2.Clsmaller10,0.5*Z.NO2.Clsmaller10,
                                  0.5*Z.NO2.Cllarger20 ,0.5*Z.NO2.Cllarger20)),
                         list(x=c(0,10,20,max.ClforNO2),
                              y=c(1.0*Z.NO2.Clsmaller10,1.0*Z.NO2.Clsmaller10,
                                  1.0*Z.NO2.Cllarger20 ,1.0*Z.NO2.Cllarger20)),
                         list(x=c(0,10,20,max.ClforNO2),
                              y=c(1.5*Z.NO2.Clsmaller10,1.5*Z.NO2.Clsmaller10,
                                  1.5*Z.NO2.Cllarger20 ,1.5*Z.NO2.Cllarger20)),
                         list(x=c(0,10,20,max.ClforNO2),
                              y=c(2.0*Z.NO2.Clsmaller10,2.0*Z.NO2.Clsmaller10,
                                  2.0*Z.NO2.Cllarger20 ,2.0*Z.NO2.Cllarger20)),
                         list(x=c(0,max.ClforNO2),
                              y=c(max.NO2,max.NO2))),
      u           = c(1.0,0.8,0.6,0.4,0.2,0.0),
      lead        = 1,
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  NO2 <- NO2orig
  if ( modify ) NO2 <- NO2modif
    
  # Ammonium
  
  NH4orig <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_nutrients_NH4",dict),
      name.attrib = c(ecoval.translate("A_nutrients_TforNH4_degC",dict),
                      ecoval.translate("A_nutrients_NH4_mgNperl",dict)),
      ranges      = list(c(0,max.TforNH4),c(0,max.NH4)),
      isolines    = list(list(x=c(0,max.TforNH4),
                              y=c(0,0)),
                         list(x=c(0,9.95,10.05,max.TforNH4),
                              y=c(0.2*Z.NH4.Tsmaller10,0.2*Z.NH4.Tsmaller10,   # exception
                                  0.2*Z.NH4.Tlarger10 ,0.2*Z.NH4.Tlarger10)),  # exception
                         list(x=c(0,9.95,10.05,max.TforNH4),
                              y=c(1.0*Z.NH4.Tsmaller10,1.0*Z.NH4.Tsmaller10,
                                  1.0*Z.NH4.Tlarger10 ,1.0*Z.NH4.Tlarger10)),
                         list(x=c(0,9.95,10.05,max.TforNH4),
                              y=c(1.5*Z.NH4.Tsmaller10,1.5*Z.NH4.Tsmaller10,
                                  1.5*Z.NH4.Tlarger10 ,1.5*Z.NH4.Tlarger10)),
                         list(x=c(0,9.95,10.05,max.TforNH4),
                              y=c(2.0*Z.NH4.Tsmaller10,2.0*Z.NH4.Tsmaller10,
                                  2.0*Z.NH4.Tlarger10 ,2.0*Z.NH4.Tlarger10)),
                         list(x=c(0,max.TforNH4),
                              y=c(max.NH4,max.NH4))),
      u           = c(1.0,0.8,0.6,0.4,0.2,0.0),
      lead        = 1,
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  NH4modif <-
    utility.endnode.intpol2d.create(
      name.node   = ecoval.translate("N_nutrients_NH4",dict),
      name.attrib = c(ecoval.translate("A_nutrients_TforNH4_degC",dict),
                      ecoval.translate("A_nutrients_NH4_mgNperl",dict)),
      ranges      = list(c(0,max.TforNH4),c(0,max.NH4)),
      isolines    = list(list(x=c(0,max.TforNH4),
                              y=c(0,0)),
                         list(x=c(0,translow.TforNH4,transhigh.TforNH4,max.TforNH4),
                              y=c(0.2*Z.NH4.Tsmaller10,0.2*Z.NH4.Tsmaller10,           # exception
                                  0.2*Z.NH4.Tlarger10 ,0.2*Z.NH4.Tlarger10)),          # exception
                         list(x=c(0,translow.TforNH4,transhigh.TforNH4,max.TforNH4),
                              y=c(1.0*Z.NH4.Tsmaller10,1.0*Z.NH4.Tsmaller10,
                                  1.0*Z.NH4.Tlarger10 ,1.0*Z.NH4.Tlarger10)),
                         list(x=c(0,translow.TforNH4,transhigh.TforNH4,max.TforNH4),
                              y=c(1.5*Z.NH4.Tsmaller10,1.5*Z.NH4.Tsmaller10,
                                  1.5*Z.NH4.Tlarger10 ,1.5*Z.NH4.Tlarger10)),
                         list(x=c(0,translow.TforNH4,transhigh.TforNH4,max.TforNH4),
                              y=c(2.0*Z.NH4.Tsmaller10,2.0*Z.NH4.Tsmaller10,
                                  2.0*Z.NH4.Tlarger10 ,2.0*Z.NH4.Tlarger10)),
                         list(x=c(0,max.TforNH4),
                              y=c(max.NH4,max.NH4))),
      u           = c(1.0,0.8,0.6,0.4,0.2,0.0),
      lead        = 1,
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  NH4 <- NH4orig
  if ( modify ) NH4 <- NH4modif
  
  # TOC 
  
  TOC <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_TOC",dict),
      name.attrib = ecoval.translate("A_nutrients_TOC_mgCperl",dict),
      range       = c(0,max.TOC),
      x           = c(0  ,min.TOC,0.5*Z.TOC,Z.TOC,1.5*Z.TOC,2*Z.TOC,max.TOC),
      u           = c(1.0,  1.0  ,   0.8   , 0.6 ,   0.4   ,  0.2  ,  0.0  ),
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # DOC
  
  DOC <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_DOC",dict),
      name.attrib = ecoval.translate("A_nutrients_DOC_mgCperl",dict),
      range       = c(0,max.DOC),
      x           = c(0  ,min.DOC,0.5*Z.DOC,Z.DOC,1.5*Z.DOC,2*Z.DOC,max.DOC),
      u           = c(1.0,  1.0  ,   0.8   , 0.6 ,   0.4   ,  0.2  ,  0.0  ),
      utility     = FALSE,
      required    = TRUE,
      col         = col)
  
  # BSB5 
  
  BOD5 <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_nutrients_BOD5",dict),
      name.attrib = ecoval.translate("A_nutrients_BOD5_mgOperl",dict),
      range       = c(0,max.BOD5),
      x           = c(0  ,min.BOD5,0.5*Z.BOD5,Z.BOD5,1.5*Z.BOD5,2*Z.BOD5,max.BOD5),
      u           = c(1.0,   1.0  ,    0.8   ,  0.6 ,    0.4   ,   0.2  ,   0.0  ),
      utility     = FALSE,
      required    = FALSE,
      col         = col)
  
  # Modul Chemie/Naehrstoffe hierarachische Aggregation
  
  nutrients <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_nutrients",dict),
      nodes     = list(TP,
                       TPfilt,
                       PO4,
                       TN,
                       NO3,
                       NO2,
                       NH4,
                       TOC,
                       DOC,
                       BOD5),
      name.fun  = "utility.aggregate.mix",
      par       = c(1,1,1,1,1,1,1,1,1,1,0.5,0.5,0),
      names.par = c("w_nutrients_TP","w_nutrients_TPfilt","w_nutrients_PO4",
                    "w_nutrients_TN","w_nutrients_NO3","w_nutrients_NO2",
                    "w_nutrients_NH4","w_nutrients_TOC","w_nutrients_DOC",
                    "w_nutrients_BOD5",
                    "w_nutrients_add","w_nutrients_min","w_nutrients_geo"),
      required  = TRUE,
      col       = col)
  return(nutrients)
}
