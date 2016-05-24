test.idturtle <-
function (data.comp, data.ref, date="date",ID="ID",PL="PL") 
{
  lista<-list("ID","Ident","Years","minPL","PLDif","maxRMSD","minRMSD","NMeasurements","NTurtles")
  compara<-data.frame()
  for (i in 1:length(lista)){
    compara[i]<-vector()
  }
  names(compara)<-lista
  
  idsrepe<-as.factor(data.comp[,c(ID)])
  repes<-data.frame(summary(idsrepe,maxsum=length(idsrepe)))
  
  repes[1]<-row.names(repes)
  repes[2]<-vector()
  for (i in 1:nrow(data.comp)){
    id.comp<-data.comp[i,c(ID)]
    repes[i,2]<-nrow(data.ref[data.ref$ID==id.comp,])
  }
  names(repes)<-list("ID","N")
  
  for(i in 1:nrow(data.comp)){
    print(paste("Reading file ",i," from ",nrow(data.comp),sep=""))
    nombretort<-paste("./turtle_",i,".csv",sep="")
    tort<-read.table(nombretort,header=TRUE,sep=",")
    idcomp<-repes[i,c(ID)]
    datoscomp<-data.ref[data.ref$ID==idcomp,]
    datosref<-data.comp[data.comp$ID==idcomp,]
    fechas<-as.character(datoscomp[,c(date)])
    fechas1<-as.character(datosref[,c(date)])
    fechas<-data.frame(unlist(strsplit(fechas,"-")))
    fechas1<-data.frame(unlist(strsplit(fechas1,"-")))
    names(fechas)<-"Year"
    names(fechas1)<-"Year"
    fechas<-rbind(fechas1,fechas)
    fechas$Year<-as.vector(fechas[,c("Year")])
    fechas$Year<-as.numeric(fechas[,c("Year")])
    fechas<-fechas[fechas$Year>1000,]
    compara[i,1]<-idcomp
    compara[i,3]<-max(fechas)-min(fechas)
    PLs<-rbind(datoscomp,datosref)
    PLs<-PLs[,c(PL)]
    compara[i,4]<-min(PLs)
    compara[i,5]<-max(PLs)-min(PLs)
    compara[i,6]<-max(tort[tort$ID==idcomp,c("RMSD")])
    compara[i,7]<-min(tort[tort$ID!=idcomp,c("RMSD")])
    cont<-0
    for(j in 1:nrow(datoscomp)){
      if (tort[j,1]!=idcomp){
        cont<-cont+1
      }
    }
    if(cont==0){
      compara[i,2]<-"Y"
    }
    else{
      compara[i,2]<-"N"
    }
    grupo<-tort[tort$RMSD<compara[i,6],]
    grupono<-grupo[grupo$ID!=idcomp,]
    compara[i,8]<-nrow(grupono)
    ids<-as.factor(grupono[,c(ID)])
    ids<-data.frame(summary(ids,maxsum=length(ids)))
    names(ids)<-"N"
    ids<-ids[ids$N>=1,]
    compara[i,9]<-(length(ids))
  }
  print("Writing file in hard disc")
  write.table(compara, file="./effectivity.csv", sep=",", row.names=FALSE, quote=FALSE)
  
  print("Done")
  return(compara)
}
