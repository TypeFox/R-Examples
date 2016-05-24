genesFromRandomPathways<-function(seed=1234,GPSrepo,np,ng){
set.seed(seed)
    fr<-GPSrepo$origRepo[[3]]
    p1<-sample(unique(fr[,1]),np)
cat("### randomly selected pathways are: \n")
    cat(paste(as.character(GPSrepo$origRepo[[1]])[p1],"\n"))
g1<-sample(unique(fr[fr[,1]%in%p1,2]),ng)
queryList<-as.character(GPSrepo$origRepo[[2]])[g1]
## cat("### Sigora returns: \n")   
## sigoraRes<-sigora(GPSrepo=GPSrepo,queryList=as.character(queryList),level=4)
## cat("### traditional Overrepresntation Analysis would return: \n")
##oraRes<-ora(queryList,GPSrepo)
## print(oraRes)
res<-list()
res[['genes']]<-queryList
res[['selectedPathways']]<-as.character(GPSrepo$origRepo[[1]])[p1]
##res[['Sigora.Results']]<-sigoraRes
##res[['ORA.Results']]<-oraRes
invisible(res)
}
