"cirock" <-
function(nvar,em=.dFvGet()$em,cr=.dFvGet()$cr,iopt=1) {
if (missing(nvar)) messagena("nvar")
vk <- single(1)
f.res <- .Fortran("cirock",
em=to.single(em),
cr=to.single(cr),
nvar=to.integer(nvar),
iopt=to.integer(iopt),
vk=to.single(vk))
list(vk=f.res$vk)
}
