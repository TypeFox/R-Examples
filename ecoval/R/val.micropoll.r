val.micropoll.create <- function(language="English",dictionaries=NA,col="black")
{
  # Value function for micropollutants with continous inputs (often from point sources) and 
  # micropollutants with episodic inputs (mainly pesticides, often from diffuse sources).
  # The assessment of micropollutants with continous inputs is based on their toxicity regarding 
  # the organism groups fish, invertebrates, and algae as described in Junghans et al 2013.
    
  # References: 
  #
  # Junghans, M. Kunz, P., Werner, I. Toxizitaet von Mischungen, aktuelle praxisorientierte Ansaetze 
  # fuer die Beurteilung von Gewaesserproben. Aqua und Gas 5, 54-61, 2013.
  #
  # Goetz, Ch. Kase, R., Kienle, C., Hollender, J. Mikroverunreinigungen aus kommunalem Abwasser: Kombination von Expositions- 
  # und oekotoxikologischen Effektdaten. Gas Wasser Abwasser 7, 575-585, 2010.
  #
  # Goetz, C.W., R. Kase und J. Hollender. Mikroverunreinigungen - Beurteilungskonzept fuer
  # organische Spurenstoffe aus kommunalem Abwasser". Studie im Auftrag des BAFU. Eawag, Duebendorf, 2010.
  
  # dictionary for node, attribute and attribute level names:
  # =========================================================
  
  dict <- ecoval.dict(language,dictionaries)
  
  # construction of nodes for trophic levels:
  # =========================================
  
  max.TU           <-   20
  u.verygood       <-   0.9 
  
  micropoll_continous_Fish <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_micropollcontinous_fish",dict),
      name.attrib = ecoval.translate("A_micropollcontinous_fish_TU",dict),
      range       = c(0,max.TU),
      x           = c(0  ,0.01,      0.1,1.0,2.0,10 ,max.TU),
      u           = c(1.0,u.verygood,0.8,0.6,0.4,0.2,0.0   ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  micropoll_continous_Invertebrates <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_micropollcontinous_invertebrates",dict),
      name.attrib = ecoval.translate("A_micropollcontinous_invertebrates_TU",dict),
      range       = c(0,max.TU),
      x           = c(0  ,0.01,      0.1,1.0,2.0,10 ,max.TU),
      u           = c(1.0,u.verygood,0.8,0.6,0.4,0.2,0.0   ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  micropoll_continous_Algae <-
    utility.endnode.intpol1d.create(
      name.node   = ecoval.translate("N_micropollcontinous_algae",dict),
      name.attrib = ecoval.translate("A_micropollcontinous_algae_TU",dict),
      range       = c(0,max.TU),
      x           = c(0  ,0.01,      0.1,1.0,2.0,10 ,max.TU),
      u           = c(1.0,u.verygood,0.8,0.6,0.4,0.2,0.0   ),
      required    = FALSE,
      utility     = FALSE,
      col         = col)
  
  # aggregation nodes
  
  micropoll_continous <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_micropollcontinous",dict),
      nodes     = list(micropoll_continous_Fish,
                       micropoll_continous_Invertebrates,
                       micropoll_continous_Algae),
      name.fun  = "utility.aggregate.mix",
      par       = c(1,1,1,0.1,1,0),
      names.par = c("w_micropollcontinous_fish","w_micropollcontinous_invertebrates",
                    "w_micropollcontinous_algae",
                    "w_micropollcontinous_add","w_micropollcontinous_min",
                    "w_micropollcontinous_geo"),
      required  = FALSE, 
      col       = col)   
  
  
  micropoll_episodic <- val.pesticides.create(language=language,dictionaries=dictionaries,col=col)
  
  micropoll <-
    utility.aggregation.create(
      name.node = ecoval.translate("N_micropoll",dict),
      nodes     = list(micropoll_episodic,
                       micropoll_continous),
      name.fun  = "utility.aggregate.mix",
      par       = c(1,1,1,1,0),
      names.par = c("w_micropollepisodic","w_micropollcontinous",
                    "w_micropoll_add","w_micropoll_min","w_micropoll_geo"),
      required  = FALSE, 
      col       = col)
  
  return(micropoll)
}
