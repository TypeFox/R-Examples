# Get the data, either in case order or sorted by the
# sequence variable order. Use only variables, and in the
# order, of those on the MRSP command line.

#Make sure this deals with missing values correctly l_DoScreen4MsngVals

#l_HasSPVar if spvar name was supplied this was true then the data is sorted based on this
#so if a sequencing variable is provided, make sure to sort based on it

#Find default number of permutations for mrsp

sp<-function(data, expon=1,commens=TRUE,
  number.perms,exact=FALSE,save.test,sequence,variables){
   
 
  # deriving some required input from the data matrix
  Call<-match.call()
    if(missing(data)) {if(missing(sequence)) {data <- variables
                         colnames(data)<-Call$variables
                        } else {data <- cbind(variables,sequence)
                             colnames(data)<-cbind(Call$variables,Call$sequence) 
                                }
                         }
    
   if(!is.null(Call$variables)) x<-data[,match(eval(Call$variables),names(data))]
      else x<-data
    x<-as.matrix(x)
   
  ifelse((missing(number.perms)),{do.resample<-FALSE;number.perms<-0},do.resample<-TRUE)
   
       
  if(!is.null(Call$sequence)){
   if(is.character(eval(Call$sequence))) sequen<-order(data[,match(eval(Call$sequence),names(data))])
      else sequen<-order(eval(Call$sequence))
   x<-as.matrix(x[sequen,])
       }
       
  i_NumObs=nrow(x)
  i_NumVars=ncol(x)
   
  #remove incomplete cases and issue a warning
  x<-as.matrix(x[complete.cases(x),])
  comp.cases<-nrow(x)
  if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
    " cases were removed because of missing values."),sep=""))
       i_NumObs<-nrow(x)   
    }
          
      if(missing(save.test)) save.test<-FALSE

      da_STV<-rep(0,times=number.perms)
    
    if(i_NumObs>200 && do.resample==TRUE){
    cat("This may take a while")
    }
   if(i_NumObs<6 & !exact) stop(paste("There are only", i_NumObs, "cases.  There must be at least 6 to do MRSP"))
   if(i_NumObs<2 & exact) stop(paste("There are only", i_NumObs, "cases.  There must be at least 2 to do an EMRSP"))
   if(i_NumObs>2147483647) stop("Number of cases (observations) should not exceed 2147483647 for MRSP")
   if(i_NumVars>65536) stop("Number of variables should not exceed 65536 for MRSP")
    
  # output is initilized here
  d_TestStat=0
  d_ObsDelta=0
  d_ExpDelta=0
  d_VarDelta=0
  d_SkwDelta=0
  d_RhoAgreement=0
  d_PValue=0
  iEr=0

  CommAvgDist=rep(0,times=i_NumVars)
  
   storage.mode(x) <- "double"
   storage.mode(da_STV)<-"double"
   storage.mode(CommAvgDist)="double"
    
    mrspOut=.Fortran("wrapmrsp",
          as.integer(i_NumObs),
          as.integer(i_NumVars),
          as.double(expon),
          as.integer(number.perms),
          as.logical(do.resample),
          as.logical(exact),
          as.logical(commens),
          x,
          as.double(d_TestStat),
          as.double(d_ObsDelta),
          as.double(d_ExpDelta),
          as.double(d_VarDelta),
          as.double(d_SkwDelta),
          as.double(d_RhoAgreement),
          as.double(d_PValue),
          as.integer(iEr),
          as.logical(save.test),
          da_STV,
          CommAvgDist
        )

    mrspOut<-new("MRSPObj",NumObs=mrspOut[[1]],NumVars=mrspOut[[2]],
    DistExp=mrspOut[[3]],NumPerm=mrspOut[[4]],DoResamp=mrspOut[[5]],Exact=mrspOut[[6]],Commens=mrspOut[[7]],inputData=mrspOut[[8]],
    TestStat=mrspOut[[9]],ObsDelta=mrspOut[[10]],ExpectDelta=mrspOut[[11]],DeltaVar=mrspOut[[12]],
    DeltaSkew=mrspOut[[13]],RhoAgreement=mrspOut[[14]],P_value=mrspOut[[15]],
    Call=deparse(Call,width.cutoff=200L),SaveTest=mrspOut[[17]],PermVals=mrspOut[[18]],CommAvgDist=mrspOut[[19]])

mrspOut
}
