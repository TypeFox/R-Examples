
a13.eq=function(compteur,ET_Main,biom.mf,fish.m,TopD,FormD,range.TLpred){
  kin=(ET_Main$Kin[compteur] - ET_Main$Fish_mort[compteur]) * (1 + TopD[compteur] * (sum(biom.mf[(compteur + range.TLpred[1]):(compteur + range.TLpred[2])])^(FormD[compteur]) - sum(ET_Main$B[(compteur+ range.TLpred[1]):(compteur + range.TLpred[2])])^(FormD[compteur]))/(sum(ET_Main$B[(compteur + range.TLpred[1]):(compteur + range.TLpred[2])])^(FormD[compteur]))) +  fish.m[compteur]
  return(kin)  
}
a13.eq.ac=function(compteur,ET_Main,biom.mf,fish.m.ac,TopD,FormD,range.TLpred){
  kin.ac=(ET_Main[compteur, "Kin_acc"] - ET_Main[compteur, "Fish_mort_acc"]) * (1 + TopD[compteur] * (sum(biom.mf[(compteur + range.TLpred[1]):(compteur + range.TLpred[2])])^(FormD[compteur]) - sum(ET_Main[(compteur+ range.TLpred[1]):(compteur + range.TLpred[2]), "B"])^(FormD[compteur]))/(sum(ET_Main[(compteur + range.TLpred[1]):(compteur + range.TLpred[2]), "B"])^(FormD[compteur]))) +  fish.m.ac[compteur]
  return(kin.ac)  
}

# function used to compute pB for the higest trophic levels
regPB=function(compteur,pb.mf,TL_out,range.highTL){
    x. <- TL_out[(range.highTL[1]):(range.highTL[2])]
    y <- log(pb.mf[(range.highTL[1]):(range.highTL[2])])
    reg <- coef(lm(y ~ x.))
    reg. <- exp(reg[1] + reg[2] * TL_out[compteur])
    return(reg.)
}
regPB.ac=function(compteur,pb.mf.ac,TL_out,range.highTL){
    x. <- TL_out[(range.highTL[1]):(range.highTL[2])]
    y <- log(pb.mf.ac[(range.highTL[1]):(range.highTL[2])])
    #if(fast){reg <- coef(fastLm(y ~ x.))}else{
    reg <- coef(lm(y ~ x.))
    reg.ac<- exp(reg[1] + reg[2] * TL_out[compteur])
    return(reg.ac)
}