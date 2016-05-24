source("../R/clustering.r")


cat("####################################################################
########################## Test  Clustering ########################
############################ Constructeur ##########################
####################################################################\n")
cleanProg(clustering)
c2a <- clustering(xLongData=ld2,yPartition=p2a)
c2b <- clustering(xLongData=ld2,yPartition=p2b)
c2c <- clustering(xLongData=ld2,yPartition=p2c)
c2d <- clustering(xLongData=ld2,yPartition=p1a)
c2e <- clustering(xLongData=ld2,yPartition=p1b)
c2f <- clustering(xLongData=ld2,yPartition=p1c)

c2an <- clustering(xLongData=ld2n,yPartition=p2a,algorithm=c(imputation="LOCF"))#,trajSizeMin=3)
c2bn <- clustering(xLongData=ld2n,yPartition=p2b,algorithm=c(imputation="LOCB"))#,trajSizeMin=3)
c2cn <- clustering(xLongData=ld2n,yPartition=p2c,algorithm=c(imputation="LI-Global"))#,trajSizeMin=3)
c2dn <- clustering(xLongData=ld2n,yPartition=p1a,algorithm=c(imputation="LI-Local"))#,trajSizeMin=3)
c2en <- clustering(xLongData=ld2n,yPartition=p1b,algorithm=c(imputation="LI-Bissectrice"))#,trajSizeMin=3)
c2fn <- clustering(xLongData=ld2n,yPartition=p1c,algorithm=c(imputation="LI-LOCBF"))#,trajSizeMin=3)

c3a <- clustering(xLongData=ld3,yPartition=p3a)
c3b <- clustering(xLongData=ld3,yPartition=p3b)
c3c <- clustering(xLongData=ld3,yPartition=p3c)
c3d <- clustering(xLongData=ld3,yPartition=p3d)
c3e <- clustering(xLongData=ld3,yPartition=p3e)

c3an <- clustering(xLongData=ld3n,yPartition=p3a,algorithm=c(imputation="LOCF"))
c3bn <- clustering(xLongData=ld3n,yPartition=p3b,algorithm=c(imputation="LOCF"))
c3cn <- clustering(xLongData=ld3n,yPartition=p3c,algorithm=c(imputation="LOCF"))
c3dn <- clustering(xLongData=ld3n,yPartition=p3d,algorithm=c(imputation="LOCF"))
c3en <- clustering(xLongData=ld3n,yPartition=p3e,algorithm=c(imputation="LOCF"))
c3fn <- clustering(xLongData=ld3n,yPartition=p3f,algorithm=c(imputation="LOCF"))
c3gn <- clustering(xLongData=ld3n,yPartition=p3g,algorithm=c(imputation="LOCF"))
c3hn <- clustering(xLongData=ld3n,yPartition=p3h,algorithm=c(imputation="LOCF"))
c3in <- clustering(xLongData=ld3n,yPartition=p3i,algorithm=c(imputation="LOCF"))


c4a <- clustering(xLongData=ld4,yPartition=p4a)
c4b <- clustering(xLongData=ld4,yPartition=p4b)
c4c <- clustering(xLongData=ld4,yPartition=p4c)
c4d <- clustering(xLongData=ld4,yPartition=p4d)
c4e <- clustering(xLongData=ld4,yPartition=p4e)

c6an <- clustering(xLongData=ld6n,yPartition=p6a,algorithm=c(imputation="LI-Bissectrice"))
c6bn <- clustering(xLongData=ld6n,yPartition=p6b,algorithm=c(imputation="LI-Bissectrice"))
c6cn <- clustering(xLongData=ld6n,yPartition=p6c,algorithm=c(imputation="LI-Bissectrice"))
c6dn <- clustering(xLongData=ld6n,yPartition=p6cn,algorithm=c(imputation="LI-Bissectrice"))

c7an <- clustering(xLongData=ld7n,yPartition=p7a)
c7bn <- clustering(xLongData=ld7n,yPartition=p7b)
c7cn <- clustering(xLongData=ld7n,yPartition=p7cn)
c7dn <- clustering(xLongData=ld7n,yPartition=p7d)
c7en <- clustering(xLongData=ld7n,yPartition=p7e)
c7fn <- clustering(xLongData=ld7n,yPartition=p7f)
c7gn <- clustering(xLongData=ld7n,yPartition=p7g)



c8an <- clustering(xLongData=ld8n,yPartition=p8a,algorithm=c(imputation="LI-Bissectrice"))
c8bn <- clustering(xLongData=ld8n,yPartition=p8b,algorithm=c(imputation="LI-Bissectrice"))
c8cn <- clustering(xLongData=ld8n,yPartition=p8c,algorithm=c(imputation="LI-Bissectrice"))
c8dn <- clustering(xLongData=ld8n,yPartition=p8cn,algorithm=c(imputation="LI-Bissectrice"))



