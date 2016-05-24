#####################
# extractor functions -- trying to find a way to be consistent with the .lm versions

#	FIXME: Can't find a way to pass m= to mlm.influence, other than providing the infl= argument explicity
#        or calling the .mlm method explicitly
#> hatvalues(Rohwer.mod, m=2)
#Error in UseMethod("hatvalues") : 
#  no applicable method for 'hatvalues' applied to an object of class "c('double', 'numeric')"
# These work:
#> hatvalues.mlm(Rohwer.mod, m=2)
#> hatvalues(Rohwer.mod, infl=mlm.influence(Rohwer.mod,m=2))

#hatvalues <- function (model, ...) 
#	UseMethod("hatvalues")


hatvalues.mlm <- function(model, m=1, infl, ...) 
{
	if (missing(infl)) {
		infl <- mlm.influence(model, m=m, do.coef=FALSE);
	}
	hat <- infl$H
	m <- infl$m
	names(hat) <- if(m==1) infl$subsets else apply(infl$subsets,1, paste, collapse=',')
	hat
}
