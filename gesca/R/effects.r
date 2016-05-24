effects <- function (B,C)
{
	#---------------------------------------------------
	# calculate total and indirect effects in measurement and structural models
	# Heungsun Hwang, Sunmee Kim
	# Last revised Aug 27, 2015
	#---------------------------------------------------
	nlv <- nrow(B)
	te_s <- t( solve(diag(1,nlv) - B) - diag(1,nlv))
	ie_s <- te_s - t(B)	
	te_m <- t(t(solve(t(diag(1,nlv) - B),t(C))))
	ie_m <- te_m - t(C)

	output.effects <- list(te_s = te_s, ie_s = ie_s, te_m = te_m, ie_m = ie_m)
	output.effects
}