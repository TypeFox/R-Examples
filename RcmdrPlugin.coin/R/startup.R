# Rcmdr dialogs for the coin package

###### TO DO: 
# to add blocks for all tests cause all have blocks, (2way conting table, max sel stat surv, paired wilcox) - not so usefull, maybe
# maybe optiuni complexe pt teste asympt/approx/exact ...  
# Check dataset whitout factors see error for block? - seems ok

###### Last modified: 
# log 2012.10.29_v1.0-21: broken into smaller files. replaced with onatach. modified namespace. added adjusted p-values for crosstabs for residuals
# log 2011.07.04_v1.0-20: all recalls verified, all help buttons checked. Added Marginal Homogeneity test for variables (not by enter table),replaced pick with select use \n for titles, fixed ansary,fligner to work with ties. added surv test distributions - Works - by Daniel Leucuta 
# log 2011.07.04_v1.0-19: Kruskal Walis added distribution options, Added zero.method for 2 sample Wilcox test - Works - by Daniel Leucuta 
# log 2011.06.18_v1.0-18: modified multiple comparison for friedman test, and Kruskal Walis so that pairwise tests works - Works - by Daniel Leucuta 
# log 2011.06.12_v1.0-17: modified multiple comparison for friedman test, and Kruskal Walis - Works - by Daniel Leucuta 
# log 2011.04.27_v1.0-16: erased old code from RcmdrPlugin.survivalT - Works - by Daniel Leucuta 
# log 2011.04.26_v1.0-15: added surv_test for survival. It workes, including ties! integrated also with survival data definition in RcmdrPlugin.survival - Works - by Daniel Leucuta 
# log 2011.04.26_v1.0-14: added maxstat test for survival. It workes integrated also with survival data definition in RcmdrPlugin.survival - Works - by Daniel Leucuta 
# log 2011.04.03_v1.0-13: added maxstat test. - Works - by Daniel Leucuta 
# log 2011.04.12_v1.0-12: fixed spearman test. - Works - by Daniel Leucuta 
# log 2011.04.03_v1.0-12: added spearman test. - NOT Works - by Daniel Leucuta 
# log 2011.04.03_v1.0-11: added Fligner Killeen test - Works but without ties/CI. Schimbat la toate normal approx cu Monte carlo. Added blocks for median test and normal test. - Works - by Daniel Leucuta 
# log 2011.03.27_v1.0-10: added Ansary Bradley test - Works but without ties/CI. added enter table for Marginal homogeneity test - Works - by Daniel Leucuta 
# log 2011.03.13_v1.0-8: added two-way table ERROR. modified in two/k sample permutation test, added function to allow for it/ Error: exact for k sample doesn t work - by Daniel Leucuta 
# log 2011.02.27_v1.0-7: added Friedman test - merge generic. Verified: OK, except multiple pairwise comparations doesn t work - by Daniel Leucuta
# log 2011.02.17_v1.0-6: added al Walis test - merge generic. minor fixes. Verified: OK, except multiple pairwise comparations doesn t work - by Daniel Leucuta
# log 2011.02.15_v1.0-5: added median and normal test. enabled blocks for Wilcox test. added condition block=group var, and disalbed confint for blocks. added ties for normal test - problem with ties method. added oneway 2 sample permutation test. Verified: OK - by Daniel Leucuta
# log 2011.02.14_v1.0-4: added paired Wilcox test function (exact, approximate, asymptotic, zero.method), added and disabled method for ties for 2 sample Wilcox test. Verified: OK - by Daniel Leucuta
# log 2011.02.14_v1.0-3: added subset box for 2 sample Wilcox test. Verified: OK - by Daniel Leucuta
# log 2011.02.14_v1.0-2: added confidence level for 2 sample Wilcox test. Verified: OK - by Daniel Leucuta

.onAttach <- function(libname, pkgname){
	if (!interactive()) return()
	Rcmdr <- options()$Rcmdr
	plugins <- Rcmdr$plugins
	if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
		Rcmdr$plugins <- c(plugins, pkgname)
		options(Rcmdr=Rcmdr)
		closeCommander(ask=FALSE, ask.save=TRUE)
		Commander()
	}
}
