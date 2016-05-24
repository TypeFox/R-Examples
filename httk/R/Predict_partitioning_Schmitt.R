# This predicts the coefficient of tissue to FREE plasma fraction via Schmitt's method (2008):  
# fupl: unbound fraction in plasma
# Pow: octonol:water partition coefficient (not log transformed)
# pKa_Donor: compound H dissociation equilibirum constant(s)
# pKa_Accept: compound H association equilibrium constant(s)   
# MA: phospholipid:water distribution coefficient
# KAPPAcell2pu: Ratio of D inside the cell to D in the plasma, as derived from the different pHs and pKas
# Fprotein.plasma: protein fraction in plasma - from Gardner 1980
predict_partitioning_schmitt <- function(chem.name=NULL,
                                         chem.cas=NULL,
                                         species='Human',
                                         default.to.human=F,
                                         parameters=NULL) #Schmitt 2008
{

  if(is.null(parameters)) parameters <- parameterize_schmitt(chem.name=chem.name,chem.cas=chem.cas,species=species,default.to.human=default.to.human)
  parameters$alpha <- 0.001
# For the "rest" tissue containing those tissues in "Physiological Parameter
# Values for PBPK Models" (2004) that are not described by Schmitt (2008)
# we use the average values for the Schmitt (2008) tissues
  tissue.data[tissue.data$Tissue=="rest",2:10]<-apply(tissue.data[1:11,2:10],2,mean)
# Then normalize to one:
  tissue.data[c(1:11,13),2:3] <- tissue.data[c(1:11,13),2:3]/apply(tissue.data[c(1:11,13),2:3],1,sum)

	Ktissue2pu <- list()
	
	# water fraction in plasma:
	FWpl <- 1 - parameters$Fprotein.plasma
	# protein fraction in interstitium:
  FPint <- 0.37 * parameters$Fprotein.plasma
	# water fraction in interstitium:
  FWint <- FWpl
	
	for (this.tissue in tissue.data$Tissue)
	{
		this.row <- tissue.data$Tissue==this.tissue
		
# Tissue-specific cellular/interstial volume fractions:
    # Cellular fraction of total volume:
		Fcell <- as.numeric(tissue.data[this.row,"Fcell"])
		# interstitial fraction of total volume:
		Fint <- as.numeric(tissue.data[this.row,"Fint"])
		if (is.na(Fint)) Fint <- 0
		
# Tissue-specific cellular sub-fractions:
		# water volume fraction:
		FW <- as.numeric(tissue.data[this.row,"FWc"])
		# protein volume fraction:
		FP <-  as.numeric(tissue.data[this.row,"FPc"])

# Tissue-specific cellular lipid sub-sub-fractions:        
		# neutral lipid volume fraction:
		Fn_L <-  as.numeric(tissue.data[this.row,"FLc"]) * as.numeric(tissue.data[this.row,"Fn_Lc"])
		if (is.na(Fn_L)) Fn_L <- 0
		# neutral phospholipid volume fraction:
		Fn_PL <- as.numeric(tissue.data[this.row,"FLc"]) * as.numeric(tissue.data[this.row,"Fn_PLc"])
		if (is.na(Fn_PL)) Fn_PL <- 0
		# acidic phospholipid volume fraction:
		Fa_PL <-  as.numeric(tissue.data[this.row,"FLc"]) * as.numeric(tissue.data[this.row,"Fa_PLc"])
		if (is.na(Fa_PL)) Fa_PL <- 0
		
		# tissue pH
		pH <- as.numeric(tissue.data[this.row,"pH"])
	
 #   # plasma:protein partition coefficient
		KPpl = 1/parameters$Fprotein.plasma*(1/parameters$Funbound.plasma-FWpl)

		# neutral phospholipid:water parition coffficient:
		if (is.null(parameters$MA))
		{
      
			Kn_PL <- 10^(0.999831 - 0.016578*parameters$temperature + 0.881721*log10(parameters$Pow)) # Based on regression to measured MA's in  Schmitt (2008)
		}else if(is.na(parameters$MA)){
 	    Kn_PL <- 10^(0.999831 - 0.016578*parameters$temperature + 0.881721*log10(parameters$Pow)) # Based on regression to measured MA's in  Schmitt (2008)
    }else{
			Kn_PL <- parameters$MA
		}

    # Need to calculate the amount of un-ionized parent:
    ionization <- calc_ionization(pH=pH,pKa_Donor=parameters$pKa_Donor,pKa_Accept=parameters$pKa_Accept)
    fraction_neutral  <- ionization[["fraction_neutral"]]
    fraction_charged <- ionization[["fraction_charged"]]
    fraction_negative <- ionization[["fraction_negative"]]
    fraction_positive <- ionization[["fraction_positive"]]
  
		# Octonol:water distribution coefficient,
    Dow <- calc_dow(parameters$Pow,fraction_neutral=fraction_neutral,alpha=parameters$alpha)

		# neutral lipid:water partition coefficient
		Kn_L <- Dow

		# protein:water partition coefficient:
		KP <- 0.163 + 0.0221*Kn_PL
		
		# acidic phospholipid:water partition coefficient:
		Ka_PL <- Kn_PL * (fraction_neutral + 20*fraction_positive + 0.05*fraction_negative)

		# unbound fraction in interstitium:
		fuint <- 1/(FWint +   FPint/parameters$Fprotein.plasma*(1/parameters$Funbound.plasma - FWpl))
		
		# unbound fraction in cellular space:
		fucell <- 1/(FW + Kn_L*Fn_L+Kn_PL*Fn_PL+Ka_PL*Fa_PL+KP*FP)
		
    KAPPAcell2pu <- calc_dow(parameters$Pow,pH=parameters$plasma.pH,alpha=parameters$alpha,pKa_Donor=parameters$pKa_Donor,pKa_Accept=parameters$pKa_Accept)/Dow
    
    if(this.tissue == 'red blood cells') eval(parse(text=paste("Ktissue2pu[\"Krbc2pu\"] <- ",as.numeric(((Fint/fuint + KAPPAcell2pu*Fcell/fucell))) ,sep='')))
		else eval(parse(text=paste("Ktissue2pu[\"K",this.tissue,"2pu\"] <- ",as.numeric(((Fint/fuint + KAPPAcell2pu*Fcell/fucell))) ,sep='')))
	}
    
 	return(Ktissue2pu)
}