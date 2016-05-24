confints.bootpls <- function(bootobject,indices=NULL,typeBCa=TRUE){
nr <- length(bootobject$t0)
if(is.null(indices)){indices <- 1:nr}
ii <- indices[1]
if(typeBCa){
temptemp.ci <- c(boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)],boot.ci(bootobject, conf = 0.95, type = c("bca"), index=ii)$bca[,-c(1,2,3)])
for(ii in indices[-1]){
temptemp.ci <- rbind(temptemp.ci,c(boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)],boot.ci(bootobject, conf = 0.95, type = c("bca"), index=ii)$bca[,-c(1,2,3)]))
}
} else {
temptemp.ci <- c(boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)])
for(ii in indices[-1]){
temptemp.ci <- rbind(temptemp.ci,c(boot.ci(bootobject, conf = 0.95, type = c("norm"), index=ii)$normal[,-1],boot.ci(bootobject, conf = 0.95, type = c("basic"), index=ii)$basic[,-c(1,2,3)],boot.ci(bootobject, conf = 0.95, type = c("perc"), index=ii)$perc[,-c(1,2,3)]))
}
}
attr(temptemp.ci, "typeBCa") <- typeBCa
rownames(temptemp.ci) <- rownames(bootobject$t0)[indices]
return(temptemp.ci)
}
