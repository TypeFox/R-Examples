################################################################################################
################################TMPM Mort Model 1.0.2###########################################
#Libraries##########################Cody Moore##################################################
#library(reshape2)
#marcTable is default injury lexicon unless otherwise specified
#ILex = 9 for ICD-9, 10 = ICD-10, 0 for AIS lexicon

##Function start##
tmpm<-function(Pdat,ILex = 1,ICs = marcTable,Long = FALSE)
{
  xBeta<-NULL #1.0.3, Compatibility for R CMD Check

  if(Long==TRUE) #(c0.0.3 added long format capability)
  {
    Pdat<-dcast(Pdat,Pdat[,1]~Pdat[,2]) #long to wide
    colnames(Pdat)[1]<-"inc_key"
  }
  if (ILex == 9)

  {
    ICodes<-ICs[ICs[,1]=="icdIX",]  #Extract ICD-9 Lexicon

    MortModel<-function(marc1,marc2,marc3,marc4,marc5,S_Region,Interaction) #Model Generation
    {
      xBeta<-(1.406958*marc1)+(1.409992*marc2)+(0.5205343*marc3)+(0.4150946*marc4)+
        (0.8883929*marc5)+(-(0.0890527)*S_Region)+(-(0.7782696*Interaction))-(2.217565)
      return(xBeta)
    }
    cat("Initializing ICD-9 Mortality Model Prediction\n")
  } else if(ILex == 10)
    {
      ICodes<-ICs[ICs[,1]=="icdX",]   #Extract ICD-10 Lexicon

      MortModel<-function(marc1,marc2,marc3,marc4,marc5,S_Region,Interaction) #Model Generation
    {
      xBeta<-(1.406958*marc1)+(1.409992*marc2)+(0.5205343*marc3)+(0.4150946*marc4)+
        (0.8883929*marc5)+(-(0.0890527)*S_Region)+(-(0.7782696*Interaction))-(2.217565)
      return(xBeta)
    }
    cat("Initializing ICD-10 Mortality Model Prediction\n")
  } else
  {
    ICodes<-ICs[ICs[,1]=="ais",]    #Extract AIS Lexicon

    MortModel<-function(marc1,marc2,marc3,marc4,marc5,S_Region,Interaction) #Model Generation
    {
      xBeta<-(1.3138*marc1)+(1.5136*marc2)+(0.4435*marc3)+(0.4240*marc4)+
        (0.6284*marc5)+(-(0.1377)*S_Region)+(-(0.6506*Interaction))-(2.3281)
    }
    return(xBeta)
    cat("Initializing AIS Mortality Model Prediction\n")
  }


  app<-function(x)
  {
    marclist<-ICs[match(x[-1],ICs[,2]),] #Lookup table (v1.0.1)

    marclist<-marclist[order(-marclist[,3]),] #Sort by marc descending order

    TCs<-marclist[1:5,]           #Generation of marc values row1 = marc1..
    TCs[,3][is.na(TCs[,3])]<-0    #Converting marc column NA's to 0's

    RegionCheck<-function(TCs) # same_region
    {
      if(TCs[1,3]!=0 & TCs[2,3]!=0 & TCs[1,4]==TCs[2,4]) #compare top 2 imarcs
      {
        sr<-1 #if !=0 and identical, set to 1
      } else
      {
        sr<-0
      }
      return(sr)
    }

    same_region<-RegionCheck(TCs) # Application of Fx

    #Interaction Term
    Imarc<-TCs[1,3]*TCs[2,3] #Multiply top 2 marcs

    Model_Calc<-MortModel(TCs[1,3],TCs[2,3],TCs[3,3],TCs[4,3],TCs[5,3],same_region,Imarc)  #Calc Model
    probDeath<-pnorm(Model_Calc) #Insert into matrix
    return(probDeath)
  }

  pDeath<-apply(Pdat,1,app) #Apply model calc
  pDeath<-matrix(pDeath,ncol=1)

  tmpm_out<-data.frame(Pdat,pDeath) #Bind to original set
  cat("Mortality Model Prediction Complete\n")

  return(tmpm_out)
}






