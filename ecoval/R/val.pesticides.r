val.pesticides.create <- function(language="English",dictionaries=NA,col="black")
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
  # Balsiger 2007, Gewaesserbelastung durch Pestizide, Gas Wasser Abwasser 3/2007
  #  
  # Chevre et al. 2006: Pestizide in Schweizer Oberflaechengewaessern. Wirkungs-
  # basierte Qualitaetskriterien. Gas Wasser Abwasser 4/2006. S. 297-307
  #
  # ============================================================================
  
  # dictionary for node, attribute and attribute level names:
  # =========================================================
  
  dict <- ecoval.dict(language,dictionaries)
  
  # implementation of valuation nodes:
  # ----------------------------------
  
  # groups according to mode of action, comparison with acute and chronic environmental quality criteria, 
  # sum of risk quotients for substances with same mode of action, evaluation according to AWEL 2006 based on Chevre et al. 2006
  
  # Mode of actions:
  # AChEI (acetylcholinesterase inhibitors) = organophosphates, e.g. Diazinon, Dimethoat, Pirimicarb
  # PhotosynthInh (photosynthesis inhibitors acting on photosystem II), e.g. Atrazin, Cyanazin, Simazin, Chlortoluron,
  #               Terbuthylazin, Terbutryn, Diuron, Isoproturon, Linuron, Metoxuron
  # VLCFASI (verylong-chain fatty acids (VLCFAs) synthesis inhibitors) = chloroacetanilides, e.g.Alachlor, Dimethenamid, Metazachlor
  # AuxinAct (herbicides that influence the auxin activity of specific plants), e.g. D24, Mcpa, Mecoprop
  
  # photosynthesis inhibitors
  comb <- data.frame(c(ecoval.translate("L_pesticides_photosynthesisinhibitors_class_verygood",dict),
                       ecoval.translate("L_pesticides_photosynthesisinhibitors_class_good",dict),
                       ecoval.translate("L_pesticides_photosynthesisinhibitors_class_moderate",dict),
                       ecoval.translate("L_pesticides_photosynthesisinhibitors_class_poor",dict),
                       ecoval.translate("L_pesticides_photosynthesisinhibitors_class_bad",dict)))
  colnames(comb) <- ecoval.translate("A_pesticides_photosynthesisinhibitors_class",dict)      
  
  photosynthesisinhibitors <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_pesticides_photosynthesisinhibitors",dict),
      attrib.levels = comb,
      u             = c(1,0.75,0.5,0.25,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  #auxinactivity
  comb <- data.frame(c(ecoval.translate("L_pesticides_auxinactivity_class_verygood",dict),
                       ecoval.translate("L_pesticides_auxinactivity_class_good",dict),
                       ecoval.translate("L_pesticides_auxinactivity_class_moderate",dict),
                       ecoval.translate("L_pesticides_auxinactivity_class_poor",dict),
                       ecoval.translate("L_pesticides_auxinactivity_class_bad",dict)))
  colnames(comb) <- ecoval.translate("A_pesticides_auxinactivity_class",dict)      
  
  auxinactivity <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_pesticides_auxinactivity",dict),
      attrib.levels = comb,
      u             = c(1,0.75,0.5,0.25,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  #VLCFASI
  comb <- data.frame(c(ecoval.translate("L_pesticides_VLCFASI_class_verygood",dict),
                       ecoval.translate("L_pesticides_VLCFASI_class_good",dict),
                       ecoval.translate("L_pesticides_VLCFASI_class_moderate",dict),
                       ecoval.translate("L_pesticides_VLCFASI_class_poor",dict),
                       ecoval.translate("L_pesticides_VLCFASI_class_bad",dict)))
  colnames(comb) <- ecoval.translate("A_pesticides_VLCFASI_class",dict)      
  
  VLCFASI <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_pesticides_VLCFASI",dict),
      attrib.levels = comb,
      u             = c(1,0.75,0.5,0.25,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # AChEI
  comb <- data.frame(c(ecoval.translate("L_pesticides_AChEI_class_verygood",dict),
                       ecoval.translate("L_pesticides_AChEI_class_good",dict),
                       ecoval.translate("L_pesticides_AChEI_class_moderate",dict),
                       ecoval.translate("L_pesticides_AChEI_class_poor",dict),
                       ecoval.translate("L_pesticides_AChEI_class_bad",dict)))
  colnames(comb) <- ecoval.translate("A_pesticides_AChEI_class",dict)      
  
  AChEI <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_pesticides_AChEI",dict),
      attrib.levels = comb,
      u             = c(1,0.75,0.5,0.25,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col)
  
  # pesticides hierarachical aggregation
  
  pesticides <- 
    utility.aggregation.create(
      name.node = ecoval.translate("N_pesticides",dict),
      nodes     = list(photosynthesisinhibitors,
                       auxinactivity,
                       VLCFASI,
                       AChEI),
      name.fun  = "utility.aggregate.mix",
      par       = c(1,1,1,1,1,1,0),
      names.par = c("w_pesticides_photosynthesisinhibitors",
                    "w_pesticides_auxinactivity",
                    "w_pesticides_VLCFASI",
                    "w_pesticides_AChEI",
                    "w_pesticides_add",
                    "w_pesticides_min",
                    "w_pesticides_geo"),
      required  = FALSE,
      col       = col)
  
  return(pesticides)
}