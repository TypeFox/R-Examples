pctAffy <- c(100,70,50,30)
pctIllu <- c(100,79,50,34)
for(tf in 1:4){
    affy <- loadCnRegionData("GSE29172", tumorFraction=pctAffy[tf]/100)
    illu <- loadCnRegionData("GSE11976", tumorFraction=pctIllu[tf]/100)
    nb <- 2*max(table(illu$region))
    tt <- do.call(rbind,
                  lapply(unique(affy$region), function(rr){
                      dat=subset(affy, region==rr)
                      idx=sample(nrow(dat), nb)
                      dat[idx,]
                  })
                  )
    saveObject(tt, sprintf("GenomeWideSNP_6/GSE29172,ASCRMAv2,H1395vsBL1395,%s,cnRegions.xdr",pctAffy[tf]))	
}
