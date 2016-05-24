mrbp<-function(variables,group,block,data, expon=1, exact=FALSE,number.perms,commens=TRUE,align=TRUE,save.test=FALSE){
  Call<-match.call()
      if((missing(variables) + missing(group) + missing(block))%in%c(1,2)) stop("either all or none of variables, group, and block must be supplied")      
       if(missing(data)) {data<-cbind(group,variables,block)
                         colnames(data)<-c(Call$group,Call$variables,Call$block)
      }else{ 
        if(missing(variables) & missing(group) & missing(block)){
          data=eval(Call$data)
       }else{
       variables<-as.character(Call$variables)
          if(length(variables)>1) variables<-variables[-c(variables=="c")]
       group<-as.character(Call$group)
       block<-as.character(Call$block)
         if(is.list(variables)) variables<-unlist(variables)
            data<-as.data.frame(data)
            data<-data[match(c(group,block,variables),names(data))]
            } 
       }
# deriving some required input from the data matrix

      if(missing(number.perms)){
        Resample = FALSE
        number.perms = 0
      } else {if(exact) {warning("exact=TRUE cannot be used when number.perms is set, a permutation test will be performed")
                          Resample=TRUE
                          exact=FALSE
            } else Resample = TRUE
          }
  
   x.orig<-data
   data<-data[order(data[,2],data[,1]),]
   group.names<-unique(data[,1])
        #remove incomplete cases
        i_NumObs=nrow(data)
       data<-data[complete.cases(data),]
      comp.cases<-nrow(data)
      if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
        " cases were removed because of missing values."),sep=""))
           i_NumObs<-nrow(data)   
        }
   #x[,1]<-as.numeric(as.factor(x[,1]))
   data<-array(data=(as.matrix(data[3:dim(data)[2]])),dim=c(length(unique(data[,1])),length(unique(data[,2])),(dim(data)[2]-2)))

  i_NumBlocks=dim(data)[2]
  i_NumVars=dim(data)[3]
    if(i_NumVars==1) commens=FALSE
    
  i_NumGrps=dim(data)[1]
  i_NumObs=i_NumGrps*i_NumBlocks
  
      if(i_NumObs>18 && exact==TRUE){
          ans<-1

          while(ans!="y" & ans!="Y" & ans!="yes" & ans!="Yes" & ans!="YES" & ans!="T" & ans!="TRUE"){

          ans<-readline(prompt=paste("Computation times for exact multiresponse randomized block procedure can be high for more than 18 observations\n",
            "Do you wish to continue? (y/n) \n",sep=" "))
            if(ans=="n" | ans=="N" | ans=="no"| ans=="No"| ans=="NO"| ans=="F"| ans=="FALSE") return()

              }
        }
  da_GpVals=vector(mode="numeric",length=i_NumGrps)
  da_BlockV=vector(mode="numeric",length=i_NumObs)
  da_alignVals=matrix(data=0,nrow=i_NumBlocks,ncol=i_NumVars)
  ch_BlkVarNam=""
  CommenAvgDist=rep(0,times=i_NumVars)
  da_STV<-rep(0,times=number.perms)
  
d_ObsDelta=0
d_ExpDelta=0
d_VarDelta=0
d_SkwDelta=0
d_AgreeVal=0
d_StdStat=0
d_PValue=0
iEr=0

  # output is initilized here
  interv=0.0
  d_Delta1=0
  da_YHot=matrix(0.,i_NumVars,i_NumVars)
  d_PValue=0
  iEr=0
  storage.mode(data) <- "double"
  storage.mode(da_YHot) <- "double"
  storage.mode(CommenAvgDist)<-"double"
  storage.mode(da_STV)<-"double"

    if(i_NumGrps < 2) stop("There must be at least 2 grouping variable values in order to do a MRBP")
     if(i_NumBlocks < 2) stop("There must be at least 2 block(s) in the specified blocking variable data")
     if(i_NumBlocks>9 & exact) stop("EMRBP must have a maximum of 9 blocks")

     if(i_NumObs>2147483647) stop("Number of cases (observations) should not exceed 2147483647 for MRBP")
    if(i_NumGrps>2147483647) stop("Number of groups should not exceed 2147483647 for MRBP")
    if(i_NumVars>65536) stop("Number of variables should not exceed 65536 for MRBP")
    if(i_NumBlocks>2147483647) stop("Number of blocks should not exceed 2147483647 for MRBP")
   
mrbpOut=.Fortran("wrapmrbp",
          data,
          as.double(expon),
          as.integer(i_NumVars),
          as.integer(i_NumBlocks),
          as.integer(i_NumGrps),
          as.integer(number.perms),
          as.integer(i_NumObs),
          da_GpVals,
          da_BlockV,
          da_alignVals,
          ch_BlkVarNam,
          as.logical(align),
          as.logical(exact),
          as.logical(Resample),
          as.logical(commens),
          as.double(d_ObsDelta),
          as.double(d_ExpDelta),
          as.double(d_VarDelta),
          as.double(d_SkwDelta),
          as.double(d_AgreeVal),
          as.double(d_StdStat),
          as.double(d_PValue),
          as.integer(iEr),
          CommenAvgDist,
          as.logical(save.test),
          da_STV
        )

mrbpOut<-new("MRBPObj",inputData=as.data.frame(x.orig),DistExp=mrbpOut[[2]],NumVars=mrbpOut[[3]],NumBlocks=mrbpOut[[4]],
          NumGrps=mrbpOut[[5]],NumPerm=mrbpOut[[6]],NumObs=mrbpOut[[7]],
          AlignVals=mrbpOut[[10]],Align=mrbpOut[[12]],Exact=mrbpOut[[13]],Resample=mrbpOut[[14]],Commensurate=mrbpOut[[15]],
          ObsDelta=mrbpOut[[16]], ExpectDelta=mrbpOut[[17]],DeltaVar=mrbpOut[[18]],
          DeltaSkew=mrbpOut[[19]],AgreeVal=mrbpOut[[20]],StdStat=mrbpOut[[21]],P_value=mrbpOut[[22]],
          Call=deparse(Call,width.cutoff=200L),CommenAvgDist=mrbpOut[[24]],
          SaveTest=mrbpOut[[25]],PermVals=mrbpOut[[26]],group.names=as.character(group.names))

mrbpOut
}


