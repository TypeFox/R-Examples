msk.fish.2004.create <- function(language="English",dictionaries=NA,col="black")
{
  # ============================================================================
  #
  # References:
  #
  # Schager, E., Peter, A. (2004), Methoden zur Untersuchung und Beurteilung der 
  # Fliessgewaesser Fische Stufe F (flaechendeckend), 
  # Mitteilungenn zum Gewaesserschutz Nr. 44
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
  
  pnt2val <- function(p,s=4)
  {
    return(approx(x    = c(0.0,1.5,5.5,9.5,13.5,16)/(16/s),
                  y    = c(1.0,0.8,0.6,0.4,0.2,0.0),
                  xout = p)$y)
  }
  
  # Artenspektrum (Namen)
  
  comb <- data.frame(spec_div=c(ecoval.translate("L_fish_species_class_natural",dict),
                                ecoval.translate("L_fish_species_class_modified",dict),
                                ecoval.translate("L_fish_species_class_atypical",dict)))
  colnames(comb) <- ecoval.translate("A_fish_species_class",dict)                              
  species <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_species",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0,2),
                        pnt2val(1,2),
                        pnt2val(2,2)),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Dominanzverhaeltnis (Namen)
  
  comb <-  data.frame(dom_spe=c(ecoval.translate("L_fish_dominance_class_typical",dict),
                                ecoval.translate("L_fish_dominance_class_tolerant",dict),
                                ecoval.translate("L_fish_dominance_class_atypical",dict))) 
  colnames(comb) <- ecoval.translate("A_fish_dominance_class",dict)                                 
  dominance <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_dominance",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0,2),
                        pnt2val(1,2),
                        pnt2val(2,2)),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Populationsstruktur Bachforelle: Altersklassen (kontinuierlich) 
  
  popstrucind_popstrucbt_agecoh <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_popstruct_trout_agestructure",dict),
      name.attrib = ecoval.translate("A_fish_popstruct_trout_ratiozeroplustoolder",dict),
      range       = c(0,10),
      x           = c(0,0.4,0.8,1.2,1.59,2.5,10),
      u             = c(pnt2val(4),
                        pnt2val(3.5),
                        pnt2val(2.5),
                        pnt2val(1.5),
                        pnt2val(0.5),
                        pnt2val(0),
                        pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Populationsstruktur Bachforelle: 0+ Fischdichte Mittelland (kontinuierlich) 
  
  popstrucind_popstrucbt_0plusdens_pla <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_popstruct_trout_zeroplusdensity_plateau",dict),
      name.attrib = ecoval.translate("A_fish_popstruct_trout_zeroplusdensity_perha",dict),
      range       = c(0,3500),
      x           = c(0,250,1000,1500,2500,3500),
      u             = c(pnt2val(4),
                        pnt2val(3.5),
                        pnt2val(2.5),
                        pnt2val(1.5),
                        pnt2val(0.5),
                        pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Populationsstruktur Bachforelle: 0+Fischdichte Jura (kontinuierlich)    
  
  popstrucind_popstrucbt_0plusdens_ju <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_popstruct_trout_zeroplusdensity_jura",dict),
      name.attrib = ecoval.translate("A_fish_popstruct_trout_zeroplusdensity_perha",dict),
      range       = c(0,3500),
      x           = c(0,250,1000,1500,2500,3500),
      u           = c(pnt2val(4),
                      pnt2val(3.5),
                      pnt2val(2.5),
                      pnt2val(1.5),
                      pnt2val(0.5),
                      pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Populationsstruktur Bachforelle: 0+Fischdichte Voralpen  (kontinuierlich) 
  
  popstrucind_popstrucbt_0plusdens_preal <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_popstruct_trout_zeroplusdensity_prealps",dict),
      name.attrib = ecoval.translate("A_fish_popstruct_trout_zeroplusdensity_perha",dict),
      range       = c(0,3000),
      x           = c(0,250,500,1000,2000,3000),
      u           = c(pnt2val(4),
                      pnt2val(3.5),
                      pnt2val(2.5),
                      pnt2val(1.5),
                      pnt2val(0.5),
                      pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Populationsstruktur Bachforelle: 0+Fischdichte Alpen (kontinuierlich)  
  
  popstrucind_popstrucbt_0plusdens_al <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_popstruct_trout_zeroplusdensity_alps",dict),
      name.attrib = ecoval.translate("A_fish_popstruct_trout_zeroplusdensity_perha",dict),
      range       = c(0,500),
      x           = c(0,100,200,300,400,500),
      u           = c(pnt2val(4),
                      pnt2val(3.5),
                      pnt2val(2.5),
                      pnt2val(1.5),
                      pnt2val(0.5),
                      pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Populationsstruktur Bachforelle conditional (kontinuierlich und Punktangabe)
  
  comb <- data.frame(c(ecoval.translate("L_region_class_swissplateau",dict),
                       ecoval.translate("L_region_class_swissjura",dict),
                       ecoval.translate("L_region_class_swissprealps",dict),
                       ecoval.translate("L_region_class_swissalps",dict)))
  colnames(comb) <- ecoval.translate("A_region_class",dict)                                                 
  popstrucind_popstrucbt_0plusdens <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_fish_popstruct_trout_zeroplusdensity",dict),
      attrib.levels = comb,
      nodes         = list(popstrucind_popstrucbt_0plusdens_pla,
                           popstrucind_popstrucbt_0plusdens_ju,
                           popstrucind_popstrucbt_0plusdens_preal,
                           popstrucind_popstrucbt_0plusdens_al),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Populationsstruktur Indikatorart 2 (Wanderarten/Aeschen/Kleinfische) (Namen)
  
  comb <- data.frame(c(ecoval.translate("L_fish_popstruct_keyspecies2_class_present",dict),
                       ecoval.translate("L_fish_popstruct_keyspecies2_class_notpresent",dict)))
  colnames(comb) <- ecoval.translate("A_fish_popstruct_keyspecies2_class",dict)   
  popstrucind_popstrucind2 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_popstruct_keyspecies2",dict),
      attrib.levels = comb,
      u             = c(1,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Populationsstruktur Indikatorart 3 (Wanderarten/Aeschen/Kleinfische) (Namen)
  
  comb <- data.frame(c(ecoval.translate("L_fish_popstruct_keyspecies3_class_present",dict),
                       ecoval.translate("L_fish_popstruct_keyspecies3_class_notpresent",dict)))
  colnames(comb) <- ecoval.translate("A_fish_popstruct_keyspecies3_class",dict)   
  popstrucind_popstrucind3 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_popstruct_keyspecies3",dict),
      attrib.levels = comb,
      u             = c(1,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Populationsstruktur Indikatorart 4 (Wanderarten/Aeschen/Kleinfische) (Namen)
  
  comb <- data.frame(c(ecoval.translate("L_fish_popstruct_keyspecies4_class_present",dict),
                       ecoval.translate("L_fish_popstruct_keyspecies4_class_notpresent",dict)))
  colnames(comb) <- ecoval.translate("A_fish_popstruct_keyspecies4_class",dict)   
  popstrucind_popstrucind4 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_popstruct_keyspecies4",dict),
      attrib.levels = comb,
      u             = c(1,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Populationsstruktur Indikatorart 5 (Wanderarten/Aeschen/Kleinfische) (Namen)
  
  comb <- data.frame(c(ecoval.translate("L_fish_popstruct_keyspecies5_class_present",dict),
                       ecoval.translate("L_fish_popstruct_keyspecies5_class_notpresent",dict)))
  colnames(comb) <- ecoval.translate("A_fish_popstruct_keyspecies5_class",dict)   
  popstrucind_popstrucind5 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_popstruct_keyspecies5",dict),
      attrib.levels = comb,
      u             = c(1,0),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Bachforelledichte Mittelland
  
  denskeysp_densbt_pla <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_trout_density_plateau",dict),
      name.attrib = ecoval.translate("A_fish_trout_density_perha",dict),
      range       = c(0,3750),
      x           = c(0,1000,2500,3750),
      u           = c(pnt2val(4),
                      pnt2val(3),
                      pnt2val(1),
                      pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Bachforelledichte Jura      
  
  denskeysp_densbt_ju <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_trout_density_jura",dict),
      name.attrib = ecoval.translate("A_fish_trout_density_perha",dict),
      range       = c(0,5250),
      x           = c(0,1000,3500,5250),
      u             = c(pnt2val(4),
                        pnt2val(3),
                        pnt2val(1),
                        pnt2val(0)),
      required    = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Bachforelledichte Voralpen  
  
  denskeysp_densbt_preal <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_trout_density_prealps",dict),
      name.attrib = ecoval.translate("A_fish_trout_density_perha",dict),
      range       = c(0,3000),
      x           = c(0,500,2000,3000),
      u             = c(pnt2val(4),
                        pnt2val(3),
                        pnt2val(1),
                        pnt2val(0)),
      required    = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Bachforelledichte Alpen  
  
  denskeysp_densbt_al <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_trout_density_alps",dict),
      name.attrib = ecoval.translate("A_fish_trout_density_perha",dict),
      range       = c(0,750),
      x           = c(0,200,500,750),
      u           = c(pnt2val(4),
                      pnt2val(3),
                      pnt2val(1),
                      pnt2val(0)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Bachforellendichte conditional
  
  comb <- data.frame(c(ecoval.translate("L_region_class_swissplateau",dict),
                       ecoval.translate("L_region_class_swissjura",dict),
                       ecoval.translate("L_region_class_swissprealps",dict),
                       ecoval.translate("L_region_class_swissalps",dict)))
  colnames(comb) <- ecoval.translate("A_region_class",dict)   
  denskeysp_densbt <-
    utility.endnode.cond.create(
      name.node     = ecoval.translate("N_fish_trout_density",dict),
      attrib.levels = comb,
      nodes         = list(denskeysp_densbt_pla,
                           denskeysp_densbt_ju,
                           denskeysp_densbt_preal,
                           denskeysp_densbt_al),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Dichte anderer Indikatorarten 2
  
  comb <- data.frame(c(ecoval.translate("L_fish_density_keyspecies2_class_high",dict),
                       ecoval.translate("L_fish_density_keyspecies2_class_moderate",dict),
                       ecoval.translate("L_fish_density_keyspecies2_class_low",dict)))
  colnames(comb) <- ecoval.translate("A_fish_density_keyspecies2_class",dict)      
  denskeysp_denskeysp2 <-  
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_density_keyspecies2",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(2),
                        pnt2val(4)),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Dichte anderer Indikatorarten 3
  
  comb <- data.frame(c(ecoval.translate("L_fish_density_keyspecies3_class_high",dict),
                       ecoval.translate("L_fish_density_keyspecies3_class_moderate",dict),
                       ecoval.translate("L_fish_density_keyspecies3_class_low",dict)))
  colnames(comb) <- ecoval.translate("A_fish_density_keyspecies3_class",dict)      
  denskeysp_denskeysp3 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_density_keyspecies3",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(2),
                        pnt2val(4)),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Dichte anderer Indikatorarten 4
  
  comb <- data.frame(c(ecoval.translate("L_fish_density_keyspecies4_class_high",dict),
                       ecoval.translate("L_fish_density_keyspecies4_class_moderate",dict),
                       ecoval.translate("L_fish_density_keyspecies4_class_low",dict)))
  colnames(comb) <- ecoval.translate("A_fish_density_keyspecies4_class",dict)      
  denskeysp_denskeysp4 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_density_keyspecies4",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(2),
                        pnt2val(4)),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  # Dichte anderer Indikatorarten 5 
  
  comb <- data.frame(c(ecoval.translate("L_fish_density_keyspecies5_class_high",dict),
                       ecoval.translate("L_fish_density_keyspecies5_class_moderate",dict),
                       ecoval.translate("L_fish_density_keyspecies5_class_low",dict)))
  colnames(comb) <- ecoval.translate("A_fish_density_keyspecies5_class",dict)      
  denskeysp_denskeysp5 <-
    utility.endnode.discrete.create(
      name.node     = ecoval.translate("N_fish_density_keyspecies5",dict),
      attrib.levels = comb,
      u             = c(pnt2val(0),
                        pnt2val(2),
                        pnt2val(4)),
      required      = FALSE,
      utility       = FALSE,
      col           = col) 
  
  
  # Deformationen/Anomalien (kontinuierlich)     
  
  deform_anom <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_fish_anomalies",dict),
      name.attrib = ecoval.translate("A_fish_anomalies_fraction",dict),
      range       = c(0,50),
      x           = c(0,1,5,50),
      u           = c(pnt2val(0),
                      pnt2val(1),
                      pnt2val(3),
                      pnt2val(4)),
      required    = FALSE,
      utility     = FALSE,
      col         = col) 
  
  # Aggregationsknoten:
  
  # Artenspektrum und Dominanzverhältnis
  
  specdom <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_fish_communitystructure",dict),
      nodes     = list(species,
                       dominance),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Populationsstruktur der Indikatorarten (Altersklassen, Reproduktion)
  
  popstrucind_popstrucbt <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_fish_popstruct_trout",dict),
      nodes     = list(popstrucind_popstrucbt_agecoh,
                       popstrucind_popstrucbt_0plusdens),
      name.fun  = "utility.aggregate.min",
      par       = c(1,1),
      required  = FALSE,
      col       = col)
  
  # Populationsstruktur der Indikatorarten
  
  populationstructure <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_fish_popstruct",dict),
      nodes     = list(popstrucind_popstrucbt,
                       popstrucind_popstrucind2,
                       popstrucind_popstrucind3,
                       popstrucind_popstrucind4,
                       popstrucind_popstrucind5),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1,1,1,1),
      required  = FALSE,
      col       = col)
  
  # Fischdichte der Indikatorarten
  
  denskeysp <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_fish_density",dict),
      nodes     = list(denskeysp_densbt,
                       denskeysp_denskeysp2,
                       denskeysp_denskeysp3,
                       denskeysp_denskeysp4,
                       denskeysp_denskeysp5),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1,1,1,1),
      required  = FALSE,
      col       = col)
  
  fish <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_fish",dict),
      nodes     = list(specdom,populationstructure,denskeysp,deform_anom),
      name.fun  = "utility.aggregate.add",
      par       = c(1,1,1,1),
      required  = FALSE,
      col       = col)
  
  return(fish)
}
