cquad <- function(formula, data, index=NULL, model=c("basic","equal","extended","pseudo"), w=rep(1,n), dyn=FALSE){

# INTERFACE FOR CQUAD ACCEPTING A FORMULA IN INPUT
#
# formula : formula with the same sintax as in plm
# data    : data.frame or pdata.frame
# index   : to denote panel structure as in plm
# model   : type of model = "basic", "equal", "extended", "pseudo"
# w       : vector of weights (not used for pseudo version)
# dyn     : for dynamic version (only for basic version)

# preliminaries
	model = match.arg(model)
	if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index)
	if (!inherits(formula, "pFormula")) formula <- pFormula(formula)
	data <- model.frame(formula, data)
	X = model.matrix(formula, data)[,-1]
	yv = pmodel.response(formula, data)
	id = attr(data,"index")[,1] 
	n = length(unique(id))
# call model
	if(model=="basic") out = cquad_basic(id,yv,X,w=w,dyn=dyn)
	if(model=="equal") out = cquad_equ(id,yv,X,w=w)
	if(model=="extended") out = cquad_ext(id,yv,X,w=w)
	if(model=="pseudo") out = cquad_pseudo(id,yv,X)
	out$formula = formula
# return output
	out

}
