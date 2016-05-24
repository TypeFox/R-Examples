turtle.idturtle <-
function (data.comp, data.ref, date="date", ID="ID", PL="PL", GuL="GuL", HumL="HumL", PecL="PecL", AbdL="AbdL",FemL="FemL", AnL="AnL",lim=10,err=3) 
{
  datapl<-data.ref[,c(ID,date,PL,GuL,HumL,PecL,AbdL,FemL,AnL)]
  dataplpr<-na.omit(datapl)
  
  if (nrow(datapl)!=nrow(dataplpr)){
    print("ERROR: You have NA values in important columns of your reference database. Please, correct your data and try again")
    break
  }
  rm(dataplpr)
  datapl[,c(date)]<-as.character(datapl[,c(date)])
  
  data<-data.comp[,c(ID,date,PL,GuL,HumL,PecL,AbdL,FemL,AnL)]
  datapr<-na.omit(data)
  
  if (nrow(data)!=nrow(datapr)){
    print("ERROR: You have NA values in important columns of your unidentified turtles dataframe. Please, correct your data and try again")
    break
  }
  rm(datapr)
  data[,c(date)]<-as.character(data[,c(date)])
  
  datapl2<-data.frame()
  data2<-data.frame()
  datapl3<-data.frame()
  
  for (i in 1:9){
    datapl2[i]<-vector()
    data2[i]<-vector()
  }
  
  for (i in 1:10){
    datapl3[i]<-vector()
  }
  
  for(i in 1:nrow(datapl)){
    datapl2[i,1]<-datapl[i,c(ID)]
    datapl2[i,2]<-datapl[i,c(date)]
    datapl2[i,3]<-datapl[i,c(PL)]
    datapl2[i,4]<-datapl[i,c(GuL)]*10/datapl[i,c(PL)]
    datapl2[i,5]<-datapl[i,c(HumL)]*10/datapl[i,c(PL)]
    datapl2[i,6]<-datapl[i,c(PecL)]*10/datapl[i,c(PL)]
    datapl2[i,7]<-datapl[i,c(AbdL)]*10/datapl[i,c(PL)]
    datapl2[i,8]<-datapl[i,c(FemL)]*10/datapl[i,c(PL)]
    datapl2[i,9]<-datapl[i,c(AnL)]*10/datapl[i,c(PL)]
  }
  
  print("Reference data prepared")
  
  ident.df<-data.frame()
  for (i in 1:7){
    ident.df[i]<-vector()
  }
  names(ident.df)<-list("TargetID","Rank","CandidateID","Date","PL_New","PL_Dif","RMSE")
  identidad<-ident.df
  
  for(i in 1:nrow(data.comp)){
    data2[i,1]<-data[i,c(ID)]
    data2[i,2]<-data[i,c(date)]
    data2[i,3]<-data[i,c(PL)]
    data2[i,4]<-data[i,c(GuL)]*10/data[i,c(PL)]
    data2[i,5]<-data[i,c(HumL)]*10/data[i,c(PL)]
    data2[i,6]<-data[i,c(PecL)]*10/data[i,c(PL)]
    data2[i,7]<-data[i,c(AbdL)]*10/data[i,c(PL)]
    data2[i,8]<-data[i,c(FemL)]*10/data[i,c(PL)]
    data2[i,9]<-data[i,c(AnL)]*10/data[i,c(PL)]
    
    print(paste("Comparing turtle ",i," with reference data",sep=""))
    
    for (j in 1:nrow(datapl2)){
      datapl3[j,1]<-datapl[j,c(ID)]
      datapl3[j,2]<-datapl[j,c(date)]
      datapl3[j,3]<-data2[i,3]-datapl2[j,3]
      datapl3[j,4]<-(datapl2[j,4]-data2[i,4])^2
      datapl3[j,5]<-(datapl2[j,5]-data2[i,5])^2
      datapl3[j,6]<-(datapl2[j,6]-data2[i,6])^2
      datapl3[j,7]<-(datapl2[j,7]-data2[i,7])^2
      datapl3[j,8]<-(datapl2[j,8]-data2[i,8])^2
      datapl3[j,9]<-(datapl2[j,9]-data2[i,9])^2
      datapl3[j,10]<-sqrt((datapl3[j,4]+datapl3[j,5]+datapl3[j,6]+datapl3[j,7]+datapl3[j,8]+datapl3[j,9])/6)
    }
    names(datapl3)<-list("ID","Date","PLDif","GulDif2","HumDif2","PecDif2","AbdDif2","FemDif2","AnDif2","RMSD")
    
    if(err==9999){}
    else{
      data.ref2<-data.ref[datapl3$PLDif>-err,]
      datapl3<-datapl3[datapl3$PLDif>-err,]
    }
    
    datapl3[,2]<-as.character(datapl3[,2])
    
    #order dataframe by RMSD
    datapl4<-datapl3[order(datapl3[,c("RMSD")]),]
    datapl4[11]<-row.names(datapl4)
    datapl4[,2]<-as.character(datapl4[,2])
    names(datapl4)<-list("ID","Date","PLDif","GulDif^2","HumDif^2","PecDif^2","AbdDif^2","FemDif^2","AnDif^2","RMSD","NRow")
    titulo<-paste("./turtle_",i,".csv",sep="")
    write.table(datapl4, file=titulo, sep=",", row.names=FALSE, quote=FALSE) 
    
    rank<-0
    
    cont<-0
    
    while (rank<lim){
      cont<-cont+1
      if(datapl4[cont,10]!=0){
        if (rank==0){
          rank<-rank+1
          ident.df[1,1]<-data2[i,1]
          ident.df[1,2]<-rank
          ident.df[1,3]<-datapl4[cont,1]
          ident.df[1,4]<-datapl4[cont,2]
          ident.df[1,5]<-data2[i,3]
          ident.df[1,6]<-datapl4[cont,3]
          ident.df[1,7]<-datapl4[cont,10]
        }
        else{
          repetido<-0
          for (j in 1:rank){
            if(datapl4[cont,1]==ident.df[j,3]){
              repetido<-repetido+1
            }
          }
          if (repetido==0){
            rank<-rank+1
            ident.df[rank,1]<-data2[i,1]
            ident.df[rank,2]<-rank
            ident.df[rank,3]<-datapl4[cont,1]
            ident.df[rank,4]<-datapl4[cont,2]
            ident.df[rank,5]<-data2[i,3]
            ident.df[rank,6]<-datapl4[cont,3]
            ident.df[rank,7]<-datapl4[cont,10]
          }
        }
      }
    }
    identidad<-rbind(identidad,ident.df)
  }
  identidad$Date<-as.vector(identidad[,c("Date")])
  print("Done!")
  return(identidad)
}
