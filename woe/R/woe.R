#' WOE
#'
#' @author
#' Sudarson Mothilal Thoppay
#' @title
#' Weigth of Evidence
#' @description
#' Computes the Weight of Evidence and Information Value between Dependent and Independent variable.
#'
#' @name woe
#' @param Data : Name of Data Set
#' @param Independent : Name of the Independent Variable
#' @param Continuous : True if the variable is continuous, False if variable is Ordinal or Nominal
#' @param Dependent : Name of the Targer Variable
#' @param C_Bin : Count of Bins to be computed
#' @param Good : Which categorical variable do you want to be good
#' @param Bad : Which categorical variable do you want to be bad
#' @return Returns a DataSet with computed WoE and IV values on success or 0 on Failure
#' @note
#' "woe" shows the log-odds ratio between between Goods and Bads.
#' In the Bivalued Dependenet variable, one value represents Goods and others are bads.
#' In Detail with an Example:
#' Let Dependent varaible be ATTRITED (0,1) and Independent variable be TENURE
#' where, 1-Attrited, 0-Non Attrited.
#' If I wish to check WOE and IV of Tenure with ATTRITED to know if Tenure has an effect in getting attrited,
#' Then good would be 1 and bad=0.
#' If I wish to check WOE and IV of Tenure with ATTRITED to know if Tenure has an effect in not getting attrited,
#' Then good would be 0 and bad=1.
#' @examples
#' woe(Data=mtcars,"cyl",FALSE,"am",10,Bad=0,Good=1)
#' woe(Data=mtcars,"mpg",TRUE,"am",10,Bad=0,Good=1)
#'
#' @export


woe<-function(Data,Independent,Continuous,Dependent,C_Bin,Bad,Good)
{
  continuous=Continuous
  C_Bin=C_Bin-1
  success<-0
  success2<-0
  j=0
  for(i in 1:ncol(Data))
  {
    if(colnames(Data[i])==Dependent)
    {
      j=i
      success<-1
      break
    }

  }

  for(Ind in 1:ncol(Data))
  {
    if(colnames(Data[Ind])==Independent)
    {
      j=Ind
      success2<-1
      break
    }

  }
  Ind<-j
  if(success==1 & success2==1)
  {
    CNO_Target=i
    Data[which(Data[,CNO_Target]==Bad),CNO_Target]=0
    Data[which(Data[,CNO_Target]==Good),CNO_Target]=1

    if(Continuous==TRUE)
    {
      BIN<-.sub_woe(Data,Ind,CNO_Target,C_Bin)
      colnames(BIN)<-c("BIN","MIN","MAX","GOOD","BAD","TOTAL","BAD%","GOOD%","TOTAL%","WOE","IV","BAD_SPLIT","GOOD_SPLIT")
      BIN<-BIN[,c(1,2,3,5,4,6,7,8,9,10,11,12,13)]
    }
    else
    {
      BIN<-.sub_woe_ON(Data,Ind,CNO_Target,C_Bin)
      colnames(BIN)<-c("BIN","BAD","GOOD","TOTAL","BAD%","GOOD%","TOTAL%","WOE","IV","BAD_SPLIT","GOOD_SPLIT")
    }
    return(BIN)
  }
  else
  {
    if(success==0)
    {
      print(paste("Variable",Dependent,"is missing in Data set "))
    }
    if(success2==0)
    {
      print(paste("Variable",Independent,"is missing in Data set "))
    }
    if(success==0 & success2==0)
    {
      print(paste("Variables",Independent,",",Dependent,"are missing in Data set "))
    }

    return(0)
  }

}


.sub_woe<-function(Data,CNO_Continuous,CNO_Target,C_Bin){
  DataSet<-data.frame(matrix(ncol=3,nrow=nrow(Data)))
  DataSet$X1<-Data[,CNO_Continuous]
  DataSet$X2<-Data[,CNO_Target]
  #DS<-arrange(DataSet,DataSet$X1)
  DS<-DataSet[order(DataSet$X1),]
  rowno<-nrow(DS)
  BIN<-1
  SHIFT<-round(rowno/C_Bin)
  for(i in 1:rowno){
    if(i<=SHIFT)
    {
      DS[i,3]<-BIN
    }
    else
    {
      SHIFT<-SHIFT+round(rowno/C_Bin)
      BIN<-BIN+1
      DS[i,3]<-BIN
    }
  }
  BINDATA<-data.frame(matrix(ncol=1,nrow=length(table(DS$X3))))
  colnames(BINDATA)[]<-c("Bin")
  for(i in 1:nrow(BINDATA))
  {
    BINDATA$Bin[i]<-i
  }

  BINDATA["mini"]<-NA
  BINDATA["maxi"]<-NA
  BINDATA["total_continuous"]<-NA
  BINDATA["Attrited"]<-NA
  BINDATA["Existed"]<-NA
  BINDATA["Total"]<-NA

  for(i in 1:nrow(BINDATA))
  {
    x<-which(DS$X3==i)
    DS[x,]->temp
    BINDATA$mini[i]<-temp$X1[which.min(temp$X1)]
    BINDATA$maxi[i]<-temp$X1[which.max(temp$X1)]
    BINDATA$total_continuous[i]<-sum(temp$X1)
    BINDATA$Existed[i]<-length(which(temp$X2==0))
    BINDATA$Attrited[i]<-length(which(temp$X2==1))
    BINDATA$Total[i]<-length(x)
  }
  BIN<-BINDATA

  BIN["AVGI"]<-NA
  BIN["P_Existed"]<-NA
  BIN["P_Attrited"]<-NA
  BIN["P_Total"]<-NA
  BIN["woe"]<-NA
  BIN["iv"]<-NA
  data=BIN
  sum_b_existed=sum(BIN$Existed)
  sum_b_attrited=sum(BIN$Attrited)
  sum_b_total=sum(BIN$Total)

  BIN[,8]<-round((BIN[,4]/BIN[,7]),3)
  BIN[,9]<-round((BIN[,6]/sum_b_existed),3)
  BIN[,10]<-round((BIN[,5]/sum_b_attrited),3)
  BIN[,11]<-round((BIN[,7]/sum_b_total),3)
  BIN[,12]<-round(log(BIN[,10]/BIN[,9],base = exp(1)),3)*100
  BIN[,13]<-round(log(BIN[,10]/BIN[,9],base = exp(1))*(BIN[,10]-BIN[,9]),3)

  BIN$EXISTED_PERCENT<-NA
  BIN$ATTR_PERCENT<-NA
  for(i in 1:nrow(BIN))
  {

    BIN$EXISTED_PERCENT[i]<-round(BIN$Existed[i]/BIN$Total[i],3)
    BIN$ATTR_PERCENT[i]<-round(BIN$Attrited[i]/BIN$Total[i],3)

  }

  BIN<-BIN[,c(1,2,3,5,6,7,9,10,11,12,13,14,15)]
  return(BIN)
}



.sub_woe_ON<-function(Data,CNO_Continuous,CNO_Target,C_Bin){

  DataSet<-data.frame(matrix(ncol=2,nrow=nrow(Data)))
  DataSet$X1<-Data[,CNO_Continuous]
  DataSet$X2<-Data[,CNO_Target]
  #DS<-arrange(DataSet,DataSet$X1)
  DS<-DataSet[order(DataSet$X1),]

  rowno<-nrow(DS)
  x<-table(DS[,1],DS[,2])
  DataSet<-data.frame(matrix(ncol=(ncol(x)+1),nrow=nrow(x)))

  for(i in 1:ncol(x))
  {

    for(j in 1:nrow(x))
    {
      DataSet[j,i]<-x[j,i]
    }

  }
  for(i in 1:nrow(x))
  {
    DataSet[i,(ncol(x)+1)]<-rownames(x)[i]
  }
  DataSet$Total<-NA
  DataSet<-DataSet[,c(3,1,2,4)]
  colnames(DataSet)<-c("BUCKET","Existed","Attrited","Total")
  BIN<-DataSet
  BIN[,4]<-BIN[,2]+BIN[,3]

  BIN["P_Existed"]<-NA
  BIN["P_Attrited"]<-NA
  BIN["P_Total"]<-NA
  BIN["woe"]<-NA
  BIN["iv"]<-NA

  sum_b_existed=sum(BIN$Existed)
  sum_b_attrited=sum(BIN$Attrited)
  sum_b_total=sum(BIN$Total)

  for(i in 1:nrow(BIN))
  {
    BIN[i,5]<-round((BIN[i,2]/sum_b_existed),3) #P_Existed<-Existed/sum_E
    BIN[i,6]<-round((BIN[i,3]/sum_b_attrited),3)
    BIN[i,7]<-round((BIN[i,4]/sum_b_total),3)
  }
  for(i in 1:nrow(BIN))
  {
    BIN$woe[i]=round(log(BIN$P_Attrited[i]/BIN$P_Existed[i],base=exp(1)),3)*100
    BIN$iv[i] =round(log(BIN$P_Attrited[i]/BIN$P_Existed[i],base=exp(1))*(BIN$P_Attrited[i]-BIN$P_Existed[i]),3)


  }
  BIN$EXISTED_PERCENT<-NA
  BIN$ATTR_PERCENT<-NA
  for(i in 1:nrow(BIN))
  {

    BIN$EXISTED_PERCENT[i]<-round(BIN$Existed[i]/BIN$Total[i],3)
    BIN$ATTR_PERCENT[i]<-round(BIN$Attrited[i]/BIN$Total[i],3)

  }



  return(BIN)
}


