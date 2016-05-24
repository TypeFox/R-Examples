RMSD.idturtle <-
function(data, ID="ID", date="date", PL="PL", GuL="GuL", HumL="HumL", PecL="PecL", AbdL="AbdL",FemL="FemL", AnL="AnL")
{
  datarow<-data[,c(PL,GuL,HumL,PecL,AbdL,FemL,AnL)]
  datarowNA<-na.omit(datarow)
  
  if (nrow(datarow)!=nrow(datarowNA)){
    print("ERROR: You have NA values in important columns of your reference database. Please, correct your data and try again")
    break
  }
  rm(datarowNA)
  datacol<-datarow
  
  matriz<-data.frame()
  for (i in 1:nrow(datarow)){
    matriz[i]<-vector()
  }
  
  dfnames<-data[,c(ID,date)]
  dfnames[,c(date)]<-as.character(dfnames[,c(date)])
  dfnames[,c(ID)]<-as.vector(dfnames[,c(ID)])
  names<-list()
  for(i in 1:nrow(dfnames)){
    names[i]<-paste(dfnames[i,c(ID)],dfnames[i,c(date)],sep="_")
  }
  
  
  datacalc<-data.frame()
  datadif<-vector()
  
  for (i in 1:6){
    datacalc[i]<-vector()
  }
  print("Calculating ratios from raw data")
  
  for(i in 1:nrow(datarow)){
    datacalc[i,1]<-datarow[i,c(GuL)]*10/datarow[i,c(PL)]
    datacalc[i,2]<-datarow[i,c(HumL)]*10/datarow[i,c(PL)]
    datacalc[i,3]<-datarow[i,c(PecL)]*10/datarow[i,c(PL)]
    datacalc[i,4]<-datarow[i,c(AbdL)]*10/datarow[i,c(PL)]
    datacalc[i,5]<-datarow[i,c(FemL)]*10/datarow[i,c(PL)]
    datacalc[i,6]<-datarow[i,c(AnL)]*10/datarow[i,c(PL)]
  }
  
  print("Ratios calculated")
  print("Comparing turtles")
  
  for(i in 1:nrow(datacalc)){
    
    print(paste("Comparing turtle ",i," with the rest",sep=""))
    
    for (j in 1:nrow(datacalc)){
      datadif[1]<-(datacalc[j,1]-datacalc[i,1])^2
      datadif[2]<-(datacalc[j,2]-datacalc[i,2])^2
      datadif[3]<-(datacalc[j,3]-datacalc[i,3])^2
      datadif[4]<-(datacalc[j,4]-datacalc[i,4])^2
      datadif[5]<-(datacalc[j,5]-datacalc[i,5])^2
      datadif[6]<-(datacalc[j,6]-datacalc[i,6])^2
      matriz[j,i]<-sqrt((datadif[1]+datadif[2]+datadif[3]+datadif[4]+datadif[5]+datadif[6])/6)      
    }
  }
  rownames(matriz)<-names
  colnames(matriz)<-names
  print("Done")
  
  return(matriz)
}
