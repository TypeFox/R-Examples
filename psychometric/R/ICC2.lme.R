"ICC2.lme" <-
function (dv, grp, data, weighted=FALSE)
 {
 require(nlme)
 attach(data)
 mod <- lme(dv ~ 1, random=~1|grp, na.action=na.omit)
 detach(data)
 if (!weighted)
  {icc2 <-  mean(GmeanRel(mod)$MeanRel) }
else { icc2 <-  weighted.mean(GmeanRel(mod)$MeanRel, GmeanRel(mod)$GrpSize) }
return(icc2)
 
 }

