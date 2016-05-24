
 library(aster)
 load("foobar.Rdata")
 aster:::setfam(famlist)
 aster:::mloglhelper(parm.mlh, pred.mlh, fam.mlh, x.mlh, root.mlh,
     modmat.mlh, origin.mlh, deriv.mlh, type.mlh)

