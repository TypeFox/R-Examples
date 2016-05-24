
setClass("MutossMethod",
		representation    = representation(
				label			= "character", # this label will be shown in the menus
				errorControl    = "character", # FWER, FWER.weak, FDR, FDX, gFWER, perComparison (?)
				callFunction    = "character", # the function to call
				output          = "character", # this is the character vector of the _possible_ output of the function
				info            = "character", # info text, should contain small description, author, reference etc.
				assumptions     = "character", # assumptions for this method
				parameters		= "list",      # optional description of parameters - see MuToss developer handbook
				misc			= "list"       # a list where you can put all your  miscellaneous stuff
		)
)

bonferroni <- function(pValues, alpha, silent=FALSE) {
	adjPValues=sapply(pValues*length(pValues),function(x){min(x,1)})
	if (missing(alpha)) {
		return(list(adjPValues=adjPValues))
	} else {
	rejected <- (adjPValues<=alpha)
	if (! silent)
	{
		cat("\n\n\t\tBonferroni correction\n\n")
		printRejected(rejected, pValues, adjPValues)
	}
	return(list(adjPValues=adjPValues, rejected=rejected,
					errorControl = new(Class='ErrorControl', type="FWER", alpha=alpha)))
}
}

mutoss.bonferroni <- function() { return(new(Class="MutossMethod",
		label="Bonferroni correction",
		errorControl="FWER",
		callFunction="bonferroni",
		output=c("adjPValues", "rejected", "errorControl"),
		info="<h2>Bonferroni correction</h2>\n\n\
<p>The classical Bonferroni correction outputs adjusted p-values, ensuring strong FWER control under arbitrary
dependence of the input p-values. It simply multiplies each input p-value by the total number of hypotheses
(and ceils at value 1).</p> 
<p>It is recommended to use Holm's step-down instead, which is valid under the exact same assumptions and more powerful.
</p>
<h3>Reference:</h3><ul><li>Bonferroni, C. E. \"<i>Il calcolo delle assicurazioni su gruppi di teste.</i>\" In Studi in Onore del Professore Salvatore Ortu Carboni. Rome: Italy, pp. 13-60, 1935.</li>\n
<li>Bonferroni, C. E. \"<i>Teoria statistica delle classi e calcolo delle probabilita.</i>\" Pubblicazioni del R Istituto Superiore di Scienze Economiche e Commerciali di Firenze 8, 3-62, 1936.</li></ul>",
		parameters=list(pValues=list(type="numeric"), alpha=list(type="numeric", optional=TRUE))
)) }


sidak <- function(pValues, alpha, silent=FALSE) {
	adjPValues <- sapply(1-(1-pValues)^length(pValues), function(x){min(x,1)})
	if (missing(alpha)) {
		return(list(adjPValues=adjPValues))
	} else {
	rejected <- (adjPValues <= alpha)
	if (! silent)
		{
		cat("\n\n\t\tSidak correction\n\n")
		printRejected(rejected, pValues, adjPValues)
	}
	return(list(adjPValues=adjPValues, rejected=rejected,
					errorControl = new(Class='ErrorControl',type="FWER",alpha=alpha)))
			}
}

mutoss.sidak <- function() { return(new(Class="MutossMethod",
					label="Sidak correction",
					errorControl="FWER",
					callFunction="sidak",
					output=c("adjPValues", "rejected", "errorControl"),
					assumptions=c("test independence"),
					info="<h2>Sidak correction</h2>\n\n\
<p>The classical Sidak correction returns adjusted p-values, ensuring strong FWER control under
the assumption of independence of the input p-values. It only uses the fact that the probability of no incorrect
rejection is the product over true nulls of those marginal probabilities (using the assumed independence of p-values).
The procedure is more generally valid for positive orthant dependent test statistics.</p> 
<p>It is recommended to use the step-down version of the Sidak correction instead, 
which is valid under the exact same assumptions and more powerful.</p>
<h3>Reference:</h3>
<ul>
<li>  Sidak, Z. (1967).<i> Rectangular confidence regions for the means of multivariate normal distributions.</i> 
Journal of the American Statistical Association, 62:626-633.</li>
</ul>\n",
					parameters=list(pValues=list(type="numeric"), alpha=list(type="numeric", optional=TRUE))
			)) }


