msk.physapp.2007.create <- function(language="English",dictionaries=NA,col="black")
{
  # ============================================================================
  #
  # References:
  #
  # Binderheim E., Goeggel W. (2007), Methoden zur Untersuchung und Beurteilung der 
  # Fliessgewaesser. Aeusserer Aspekt. 
  # Umwelt-Vollzug Nr. 0701. Bundesamt fuer Umwelt, Bern. 43 S.
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
  
  # construction of end nodes:
  # ==========================
  
  # Verschlammung Quantifizierung
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_sludge_quant_class_none",dict),
                       ecoval.translate("L_physappearance_sludge_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_sludge_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_sludge_quant_class",dict)
  sludge.quant <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_sludge_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_sludge_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_sludge_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_sludge_origin_class",dict)
  sludge.origin <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_sludge_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  sludge <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_sludge",dict),
      nodes     = list(sludge.quant,sludge.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Truebung
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_turbidity_quant_class_none",dict),
                       ecoval.translate("L_physappearance_turbidity_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_turbidity_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_turbidity_quant_class",dict)
  turbidity.quant <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_turbidity_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_turbidity_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_turbidity_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_turbidity_origin_class",dict)
  turbidity.origin <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_turbidity_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  turbidity <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_turbidity",dict),
      nodes     = list(turbidity.quant,turbidity.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Verfaerbung 
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_discolor_quant_class_none",dict),
                       ecoval.translate("L_physappearance_discolor_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_discolor_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_discolor_quant_class",dict)
  discolor.quant <- 
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_discolor_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_discolor_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_discolor_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_discolor_origin_class",dict)
  discolor.origin <- 
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_discolor_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  discolor <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_discolor",dict),
      nodes     = list(discolor.quant,discolor.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Schaum
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_foam_quant_class_none",dict),
                       ecoval.translate("L_physappearance_foam_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_foam_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_foam_quant_class",dict)
  foam.quant <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_foam_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_foam_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_foam_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_foam_origin_class",dict)
  foam.origin <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_foam_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  foam <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_foam",dict),
      nodes     = list(foam.quant,foam.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Geruch 
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_odor_quant_class_none",dict),
                       ecoval.translate("L_physappearance_odor_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_odor_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_odor_quant_class",dict)
  odor.quant <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_odor_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_odor_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_odor_origin_class_anthropogenic",dict))
  )
  colnames(comb) <- ecoval.translate("A_physappearance_odor_origin_class",dict)
  odor.origin <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_odor_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE, 
      utility       = FALSE,
      col           = col)
  
  odor <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_odor",dict),
      nodes     = list(odor.quant,odor.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Eisensulfid 
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_ironsulfide_quant_class_none",dict),
                       ecoval.translate("L_physappearance_ironsulfide_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_ironsulfide_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_ironsulfide_quant_class",dict)
  ironsulfide.quant <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_ironsulfide_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_ironsulfide_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_ironsulfide_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_ironsulfide_origin_class",dict)
  ironsulfide.origin <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_ironsulfide_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  ironsulfide <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_ironsulfide",dict),
      nodes     = list(ironsulfide.quant,ironsulfide.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Kolmation
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_clogging_quant_class_none",dict),
                       ecoval.translate("L_physappearance_clogging_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_clogging_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_clogging_quant_class",dict)
  clogging.quant <- 
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_clogging_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_clogging_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_clogging_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_clogging_origin_class",dict)
  clogging.origin <- 
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_clogging_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  clogging <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_clogging",dict),
      nodes     = list(clogging.quant,clogging.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Abf??lle 
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_solidwaste_quant_class_none",dict),
                       ecoval.translate("L_physappearance_solidwaste_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_solidwaste_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_solidwaste_quant_class",dict)   
  solidwaste <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_solidwaste",dict),
      attrib.levels = comb,
      u             = c(1.0,0.5,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # Heterotropher Bewuchs 
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_hetcover_quant_class_none",dict),
                       ecoval.translate("L_physappearance_hetcover_quant_class_scattered",dict),
                       ecoval.translate("L_physappearance_hetcover_quant_class_little",dict),
                       ecoval.translate("L_physappearance_hetcover_quant_class_moderate",dict),
                       ecoval.translate("L_physappearance_hetcover_quant_class_high",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_hetcover_quant_class",dict)
  hetcover.quant <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_hetcover_quant",dict),
      attrib.levels = comb,
      u             = c(1.0,0.75,0.5,0.25,0.0),
      required      = TRUE,
      utility       = FALSE,
      col           = col)
  
  comb <- data.frame(c(ecoval.translate("L_physappearance_hetcover_origin_class_natural",dict),
                       ecoval.translate("L_physappearance_hetcover_origin_class_anthropogenic",dict)))
  colnames(comb) <- ecoval.translate("A_physappearance_hetcover_origin_class",dict)
  hetcover.origin <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_physappearance_hetcover_origin",dict),
      attrib.levels = comb,
      u             = c(1.0,0.0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  hetcover <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance_hetcover",dict),
      nodes     = list(hetcover.quant,hetcover.origin),
      name.fun  = "utility.aggregate.max",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # aggregation
  
  physicalappearance <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_physappearance",dict),
      nodes     = list(sludge,
                       turbidity,
                       discolor,
                       foam,
                       odor,
                       ironsulfide,
                       clogging, 
                       solidwaste,
                       hetcover),
      name.fun  = "utility.aggregate.mix",
      par       = c(1,1,1,1,1,1,1,1,1,1,1,0),
      names.par = c("w_physapp_sludge","w_physapp_turbidity","w_physapp_color",
                    "w_physapp_foam","w_physapp_odor","w_physapp_ironsulfide",
                    "w_physapp_clogging","w_physapp_solidwaste","w_physapp_hetcover",
                    "w_physapp_add","w_physapp_min","w_physapp_geo"),
      required  = FALSE,
      col       = col)
  
  return(physicalappearance)
}
