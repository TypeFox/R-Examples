###
# Example of model selection with models from 'unmarked' package
###
require(MuMIn)
require(unmarked)

opt <- options(width = 110)

# from example(distsamp)
ltUMF <- local({
data(linetran)
dbreaksLine <- c(0, 5, 10, 15, 20)
lengths <- linetran$Length * 1000
with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths, survey = "line", unitsIn = "m")
	})
})

# global model - probably nonsensical:
fmUnmDS <- distsamp( ~ Length + area ~ area + habitat, ltUMF)

# The default null model used for calculating R^2 has a formula ~ 1 ~ 1
msUnmDS <- dredge(fmUnmDS, rank = AIC, extra = "adjR^2")

subset(msUnmDS, delta < 4 | df == min(df))

# Compare with the model selection table from unmarked - the statistics should
# be identical.
# fit the models from the 'top' and a null model.
models <- get.models(msUnmDS, delta < 4 | df == min(df))

modSel(fitList(fits = structure(models, names = model.names(models,
	labels = getAllTerms(fmUnmDS)))), nullmod = "(Null)")
	
options(opt)

########################