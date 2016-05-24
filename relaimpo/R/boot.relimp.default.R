"boot.relimp.default" <-
function (object, x = NULL, ..., b = 1000, type = "lmg", rank = TRUE, diff = TRUE,
    rela = FALSE, always = NULL, groups = NULL, groupnames = NULL, fixed = FALSE,
    weights = NULL, design = NULL)
{
     ynam = deparse(substitute(object))  ## make labelling of y ok
     do.call("boot.relimp.default.intern", list(object=object, x=x, ...,b=b,
           type=type, diff=diff, rank=rank, rela=rela, always=always,
          groups=groups, groupnames = groupnames, fixed=fixed,
          weights=weights, design=design,ynam=ynam))
}