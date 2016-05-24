# Derived traits for Glycan peaks in Plasma for HPLC
# based on paper from 2011.
#
# @references
# Lu et al. 
# "Screening Novel Biomarkers for Metabolic Syndrome by Profiling 
#  Human Plasma N-Glycans in Chinese Han and Croatian Populations"
# \url{http://dx.doi.org/10.1021/pr2004067}
plasma.hplc.derived.traits.2011 <- function(d, print.exp.names=FALSE) {
    if(print.exp.names){
        return("")
    }

	d = d[, c(1:(grep("GP1$", names(d))-1), grep("GP\\d+$|DG\\d+$|sial",names(d)))]
	
	# sialylation of biantennary glycans 
	d$`BAMS`= with(d, (GP7 + GP8)/(DG5 + DG6 + DG7))*100
	d$`BADS`= with(d, (GP9 + GP10 + GP11)/(DG5 + DG6 + DG7))*100
	
	# branching
	d$`BA`= with(d, (DG1 + DG2 + DG3 + DG4 + DG5 + DG6 + DG7))
	d$`TRIA`= with(d, (DG8 + DG9 + DG10))
	d$`TA`= with(d, (DG11 + DG12 + DG13))
	
	# fucosylation
	d$`C-FUC`= with(d, DG6/(DG5 + DG6))*100
	d$`A-FUC`= with(d, DG7/(DG5 + DG7))*100
	
	d$`A2`= with(d, (GP1+DG1)/2)
	
	# galactosylation
	d$`G0`= with(d, (DG1 + DG2))
	d$`G1`= with(d, (DG3 + DG4))
	d$`G2`= with(d, (DG5 + DG6 + DG7))
	d$`G3`= with(d, (GP12 + GP13 + GP14))
	d$`G4`= with(d, (GP15 + GP16))
	
	return(d)
}

