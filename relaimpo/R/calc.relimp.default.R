"calc.relimp.default" <-
function (object, x = NULL, ..., type = "lmg", diff = FALSE, rank = TRUE, rela = FALSE, always = NULL,
        groups = NULL, groupnames=NULL, weights=NULL, design=NULL)
{
     ynam = deparse(substitute(object))  ## make labelling of y ok
     do.call("calc.relimp.default.intern", list(object=object, x=x, ..., 
           type=type, diff=diff, rank=rank, rela=rela, always=always,
          groups=groups, groupnames = groupnames, weights=weights, design=design,ynam=ynam))
}