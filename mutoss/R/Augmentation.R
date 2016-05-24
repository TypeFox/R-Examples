
augmentation <- function(adjPValues, newErrorControl, newK, newQ, silent=FALSE) {
	if (newErrorControl == "gFWER") {	
		out = fwer2gfwer(adjPValues, k = newK)
		if(!silent) cat('\n\n\t\tGeneralized Family-Wise Error Rate\n\n')
		return(list(adjPValues=as.numeric(out), rejected=NULL, errorControl = new(Class='ErrorControl',type="gFWER")))
	} else if (newErrorControl == "FDX") {	
		out = fwer2tppfp(adjPValues, q = newQ)
		if(!silent) cat('\n\n\t\tTail Probability of the Proportion of False Positives\n\n')
		return(list(adjPValues=as.numeric(out), rejected=NULL, errorControl = new(Class='ErrorControl',type="FDX")))
	} else if (newErrorControl == "FDR") {	
		out = fwer2fdr(adjPValues, method = "restricted")
		if(!silent) cat('\n\n\t\tFalse Discovery Rate\n\n')
		return(list(adjPValues=as.numeric(out$adjp), rejected=NULL, errorControl = new(Class='ErrorControl',type="FDR")))		
	} else{ 
		if(!silent)cat('\n\n\t\tUnknown newErrorControl method')
	}
}

mutoss.augmentation <- function() { return(new(Class="MutossMethod",
					label="Augmentation MTP adjusted p-values",
					errorControl=c("FWER"),
					callFunction="augmentation",
					output=c("adjPValues", "rejected", "errorControl"),
					info="<h2>Augmentation MTP adjusted p-values</h2>
<p>Wrapper function to the augmentation methods of the multtest package.</p>

<p>The augmentation method turns a vector of p-values which are already adjusted for FWER control
into p-values that are adjusted for gFWER, FDX or FDR. The underlying idea (for gFWER and FDX) 
is that the set of hypotheses rejected at a given level alpha under FWER can be 'augmented'
by rejecting some additional hypotheses while still ensuring (strong) control of the desired weaker type I criterion.
For FDR, it uses the fact that FDX control for q=alpha=1-sqrt(1-beta) entails FDR control at level beta.</p>

<p>Use of these augmentation methods is recommended only in the situation where FWER-controlled p-values are
directly available from the data (using some specific method). When only marginal p-values are available,
it is generally prerefable to use other adjustment methods directly aimed at the intended criterion
(as opposed to first adjust for FWER, then augment)</p>

<p>Note: In the multtest package, two methods ('restricted' and 'conservative') are available for FDR augmentation. 
Here the 'restricted' method is forced for FDR augmentation since it is in fact always valid and better than 
'conservative' (M. van der Laan, personal communication) with respect to power.</p>
<h3>Reference:</h3>
<ul>
<li> S. Dudoit, M.J. van der Laan. \"<i> Multiple Testing Procedures with Applications to Genomics</i>\", Springer, 2008. (chapter 6) </li>\n\
<li></li>
</ul>",
# TODO: <- Add the possibility of filling the rejected slot? (need additional input alpha)
					parameters=list(adjPValues=list(type="numeric"),
								newErrorControl=list(type="character", label="New error control", choices=c("FDR","FDX","gFWER")),
								newK=list(type="numeric", optional=TRUE, label="For gFWER set k"),
                                newQ=list(type="numeric", optional=TRUE, label="For FDX set q")
						
                                    	      ))) }
