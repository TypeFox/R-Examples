lcorrel<-function(DFrame,bw,Coords)
{  
  Dij<- as.matrix(dist(Coords))
  
  Obs<-nrow(DFrame)
  VarNo<-ncol(DFrame)
  VarSQ<-VarNo^2
  
  Ne<-round(Obs*bw,0)
  DFS<-Ne-2
  
  MeltMeNames<-melt(cor(DFrame))
  TNames<-t(paste(MeltMeNames$X1, MeltMeNames$X2, sep="_"))
  
  cor_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  tv_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  sign_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  BF_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  
  wtcor_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  gwsign_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  BF_gw_table_out<-as.data.frame(setNames(replicate(VarSQ,numeric(0), simplify = F), TNames[1:VarSQ]))
  
  DNames<-names(DFrame)
  
  for(m in 1:Obs){
    #Get the data  
    DataSet<-data.frame(DFrame,DNeighbour=Dij[,m])
    
    #Sort by distance
    DataSetSorted<- DataSet[order(DataSet$DNeighbour),]
    
    #Keep Nearest Neighvbours
    SubSet1<-DataSetSorted[1:Ne,]
    
    Kernel_H<-max(SubSet1$DNeighbour)
    Wts<-(1-(SubSet1$DNeighbour/Kernel_H)^2)^2
    
    #Remove Distance to NN columne
    SubSet<-SubSet1[,1:VarNo]
    
    #Calculate Correlation Matrix
    CorMatrix<-cor(SubSet, method="pearson")
    
    WCorMatrix<-wtd.cor(SubSet,weight=Wts)
    
    t_matrix<-matrix(data=NA, nrow=VarNo,ncol=VarNo,dimnames = list(DNames,DNames))
    sig_matrix<-matrix(data=NA,nrow=VarNo,ncol=VarNo,dimnames = list(DNames,DNames))
    
    for(i in 1:VarNo){
      for(j in 1:VarNo){
        if (j!=i){
          t_matrix[i,j]<-(CorMatrix[i,j]*sqrt(DFS))/(sqrt(1-(CorMatrix[i,j]^2)))
          sig_matrix[i,j]<-2*(1-pt(abs(t_matrix[i,j]),df=DFS))
        }
      }
    }
    
    #Store in table
    MeltMeV<-melt(CorMatrix)
    cor_table_out[m,]<-t(MeltMeV$value)
    
    MeltMeTv<-melt(t_matrix)
    tv_table_out[m,]<-t(MeltMeTv$value)
    
    MeltMeSign<-melt(sig_matrix)
    sign_table_out[m,]<-t(MeltMeSign$value)
    
    BF_table_out[m,]<-p.adjust(sign_table_out[m,], method = "bonferroni", n = Obs)
    
    MeltMeWtV<-melt(WCorMatrix$correlation)
    wtcor_table_out[m,]<-t(MeltMeWtV$value)
    
    MeltMeGWSig<-melt(WCorMatrix$p.value)
    gwsign_table_out[m,]<-t(MeltMeGWSig$value)
    
    BF_gw_table_out[m,]<-p.adjust(gwsign_table_out[m,], method = "bonferroni", n = Obs)
    
  }
  lcor.results<-list(LPCC=cor_table_out,LPCC_t=tv_table_out,LPCC_sig=sign_table_out,LPCC_sig_BF=BF_table_out,GWPCC=wtcor_table_out,GWPCC_sig=gwsign_table_out,GWPCC_sig_BF=BF_gw_table_out)
  return(lcor.results)
}