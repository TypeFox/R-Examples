# 4-12-2013 MRC-Epid JHZ

solar_ped <- ped51[,c("id","fid","mid","sex")]
write.table(solar_ped,file="51.ped",col.names=c("id","fa","mo","sex"),
          quote=FALSE, row.names=FALSE,sep=",",na="")
solar_pheno <- ped51[,c("id","qt","bt")]
write.table(solar_pheno,file="51.phen",col.names=c("id","qt","bt"),
          quote=FALSE, row.names=FALSE,sep=",",na="")
