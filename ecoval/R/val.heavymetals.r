val.heavymetals.create <- function(language     = "English",
                                   dictionaries = NA,
                                   col          = "black",
                                   version      = "AWEL")
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
  # AWEL Amt fuer Abfall, Wasser, Energie und Luft, Kanton Zuerich, 
  # Statusbericht 2006: Wasserqualitaet der Seen, Fliessgewaesser und 
  # des Grundwasser im Kanton Zuerich.
  #
  # LAWA (Länderarbeitsgemeinschaft Wasser) 1998: Zielvorgaben zum Schutz oberirdischer 
  # Binnengewässer. Band II: Ableitung und Erprobung von Zielvorgaben zum Schutz  
  # oberiridischer Binnengewässer für die Schwermetalle Blei, Cadmium, Chrom, Kupfer, 
  # Nickel, Quecksilber und Zink. Kulturbuchverlag Berlin GmbH, Berlin.
  #
  # IKSR (Internationale Kommisstion zum Schutz des Rheins) 2009. Bericht Nr. 175, Sedimentmanagementplan Rhein.
  #
  # ============================================================================
  
  # dictionary for node, attribute and attribute level names:
  # =========================================================
  
  dict <- ecoval.dict(language,dictionaries)
  
  # construction of nodes for individual heavy metals:
  # ==================================================
  
  # valuation parameters according to AWEL (2006) based on LAWA or IKSR in mg/kg TS: 
  # --------------------------------------------------------------------------------
  
  Z.Zn              <-   200
  Z.Zn.IKSR         <-   200
  min.Zn            <-   1e-3
  max.Zn            <-   1200
  
  Z.Cu              <-   60
  Z.Cu.IKSR         <-   50
  min.Cu            <-   1e-3
  max.Cu            <-   250
  
  Z.Pb              <-   100
  Z.Pb.IKSR         <-   100  
  min.Pb            <-   1e-3
  max.Pb            <-   800
  
  Z.Cd              <-   1.5    
  Z.Cd.IKSR         <-   1    
  min.Cd            <-   1e-3
  max.Cd            <-   4
  
  Z.Hg              <-   1
  Z.Hg.IKSR         <-   0.5
  min.Hg            <-   1e-3
  max.Hg            <-   3
  
  Z.Ni              <-   50
  Z.Ni.IKSR         <-   50
  min.Ni            <-   1e-3
  max.Ni            <-   150
  
  Z.Cr              <-   100
  # Z.Cr.IKSR inexistent
  min.Cr            <-   1e-3
  max.Cr            <-   300
  
  
  # Version AWEL-LAWA
  # -----------------
  
  # Zinc
  
  Zn <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Zn",dict),
      name.attrib = ecoval.translate("A_heavymetals_Zn_mgkgTS",dict),
      range       = c(0,max.Zn),
      x           = c(0  ,min.Zn,0.5*Z.Zn,Z.Zn,1.5*Z.Zn,2*Z.Zn,max.Zn),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Copper
  
  Cu <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Cu",dict),
      name.attrib = ecoval.translate("A_heavymetals_Cu_mgkgTS",dict),
      range       = c(0,max.Cu),
      x           = c(0  ,min.Cu,0.5*Z.Cu,Z.Cu,1.5*Z.Cu,2*Z.Cu,max.Cu),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Lead
  
  Pb <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Pb",dict),
      name.attrib = ecoval.translate("A_heavymetals_Pb_mgkgTS",dict),
      range       = c(0,max.Pb),
      x           = c(0  ,min.Pb,0.5*Z.Pb,Z.Pb,1.5*Z.Pb,2*Z.Pb,max.Pb),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Cadmium
  
  Cd <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Cd",dict),
      name.attrib = ecoval.translate("A_heavymetals_Cd_mgkgTS",dict),
      range       = c(0,max.Cd),
      x           = c(0  ,min.Cd,0.5*Z.Cd,Z.Cd,1.5*Z.Cd,2*Z.Cd,max.Cd),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Mercury
  
  Hg <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Hg",dict),
      name.attrib = ecoval.translate("A_heavymetals_Hg_mgkgTS",dict),
      range       = c(0,max.Hg),
      x           = c(0  ,min.Hg,0.5*Z.Hg,Z.Hg,1.5*Z.Hg,2*Z.Hg,max.Hg),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)      
  
  # Nickel
  
  Ni <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Ni",dict),
      name.attrib = ecoval.translate("A_heavymetals_Ni_mgkgTS",dict),
      range       = c(0,max.Ni),
      x           = c(0  ,min.Ni,0.5*Z.Ni,Z.Ni,1.5*Z.Ni,2*Z.Ni,max.Ni),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)       
  
  # Chromium
  
  Cr <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Cr",dict),
      name.attrib = ecoval.translate("A_heavymetals_Cr_mgkgTS",dict),
      range       = c(0,max.Cr),
      x           = c(0  ,min.Cr,0.5*Z.Cr,Z.Cr,1.5*Z.Cr,2*Z.Cr,max.Cr),
      u           = c(1.0,  1.0 ,   0.8  , 0.6,   0.4  ,  0.2 ,  0.0 ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)       
  
  # Version IKSR
  # ------------
  
  Zn.IKSR <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Zn",dict),
      name.attrib = ecoval.translate("A_heavymetals_Zn_mgkgTS",dict),
      range       = c(0,12*Z.Zn.IKSR),
      x           = c(0  ,min.Zn,Z.Zn.IKSR,2*Z.Zn.IKSR,4*Z.Zn.IKSR,8*Z.Zn.IKSR,12*Z.Zn.IKSR),
      u           = c(1.0,  1.0 ,   0.8   ,    0.6    ,   0.4     ,  0.2      ,  0.0       ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)       
  
  # Copper
  
  Cu.IKSR <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Cu",dict),
      name.attrib = ecoval.translate("A_heavymetals_Cu_mgkgTS",dict),
      range       = c(0,12*Z.Cu.IKSR),
      x           = c(0  ,min.Cu,Z.Cu.IKSR,2*Z.Cu.IKSR,4*Z.Cu.IKSR,8*Z.Cu.IKSR,12*Z.Cu.IKSR),
      u           = c(1.0,  1.0 ,   0.8   ,    0.6    ,   0.4     ,  0.2      ,  0.0       ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Lead
  
  Pb.IKSR <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Pb",dict),
      name.attrib = ecoval.translate("A_heavymetals_Pb_mgkgTS",dict),
      range       = c(0,12*Z.Pb.IKSR),
      x           = c(0  ,min.Pb,Z.Pb.IKSR,2*Z.Pb.IKSR,4*Z.Pb.IKSR,8*Z.Pb.IKSR,12*Z.Pb.IKSR),
      u           = c(1.0,  1.0 ,   0.8   ,    0.6    ,   0.4     ,  0.2      ,  0.0       ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Cadmium
  
  Cd.IKSR <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Cd",dict),
      name.attrib = ecoval.translate("A_heavymetals_Cd_mgkgTS",dict),
      range       = c(0,12*Z.Cd.IKSR),
      x           = c(0  ,min.Cd,Z.Cd.IKSR,2*Z.Cd.IKSR,4*Z.Cd.IKSR,8*Z.Cd.IKSR,12*Z.Cd.IKSR),
      u           = c(1.0,  1.0 ,   0.8   ,    0.6    ,   0.4     ,  0.2      ,  0.0       ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Mercury
  
  Hg.IKSR <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Hg",dict),
      name.attrib = ecoval.translate("A_heavymetals_Hg_mgkgTS",dict),
      range       = c(0,12*Z.Hg.IKSR),
      x           = c(0  ,min.Hg,Z.Hg.IKSR,2*Z.Hg.IKSR,4*Z.Hg.IKSR,8*Z.Hg.IKSR,12*Z.Hg.IKSR),
      u           = c(1.0,  1.0 ,   0.8   ,    0.6    ,   0.4     ,  0.2      ,  0.0       ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Nickel
  
  Ni.IKSR <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_heavymetals_Ni",dict),
      name.attrib = ecoval.translate("A_heavymetals_Ni_mgkgTS",dict),
      range       = c(0,12*Z.Ni.IKSR),
      x           = c(0  ,min.Ni,Z.Ni.IKSR,2*Z.Ni.IKSR,4*Z.Ni.IKSR,8*Z.Ni.IKSR,12*Z.Ni.IKSR),
      u           = c(1.0,  1.0 ,   0.8   ,    0.6    ,   0.4     ,  0.2      ,  0.0       ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # Chromium inexistent  
  
  
  # hierarchical aggregation
  # ========================
  
  if(version == "AWEL")
  {
    heavymetals <-
      utility.aggregation.create(
        name.node = ecoval.translate("N_heavymetals",dict),
        nodes     = list(Zn,
                         Cu,
                         Pb,
                         Cd,
                         Hg,
                         Ni,
                         Cr),
        name.fun  = "utility.aggregate.mix",
        par       = c(1,1,1,1,1,1,1,1,1,0),
        names.par = c("w_heavymetals_Zn","w_heavymetals_Cu","w_heavymetals_Pb",
                      "w_heavymetals_Cd","w_heavymetals_Hg","w_heavymetals_Ni",
                      "w_heavymetals_Cr",       
                      "w_heavymetals_add","w_heavymetals_min","w_heavymetals_geo"),
        required  = FALSE,
        col         = col)
  }
  
  if(version == "IKSR")
  {
    heavymetals <-
      utility.aggregation.create(
        name.node = ecoval.translate("N_heavymetals",dict),
        nodes     = list(Zn.IKSR,
                         Cu.IKSR,
                         Pb.IKSR,
                         Cd.IKSR,
                         Hg.IKSR,
                         Ni.IKSR),
        name.fun  = "utility.aggregate.mix",
        par       = c(1,1,1,1,1,1,1,1,0),
        names.par = c("w_heavymetals_Zn","w_heavymetals_Cu","w_heavymetals_Pb",
                      "w_heavymetals_Cd","w_heavymetals_Hg","w_heavymetals_Ni",     
                      "w_heavymetals_add","w_heavymetals_min","w_heavymetals_geo"),
        required  = FALSE,
        col       = col)
  }   
  
  if(version!="AWEL" & version !="IKSR") 
  {
    heavymetals=NA
  }  
  return(heavymetals)
}
