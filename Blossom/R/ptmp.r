ptmp<-function(data,variables,expon=1,exact=FALSE,number.perms,save.test){
# deriving some required input from the data matrix

 Call<-match.call()
   if(missing(data)) {data<-variables
                         colnames(data)<-Call$variables
                         }
   x<-data
   if(!is.null(Call$variables)) data<-data[,match(eval(Call$variables),names(data))]
   #remove incomplete cases and issue a warning
   i_NumObs=nrow(x)
   x<-x[complete.cases(x),]
  comp.cases<-nrow(x)
  if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
    " cases were removed because of missing values."),sep=""))
       i_NumObs<-nrow(x)   
    }
 
 if(dim(x)[2]==3){
     names.used<-names(x)
     group.names<-unique(x[,1])
     block.names<-unique(x[,2])
      if(any(table(x[,2])!=2)) stop("Not all pairs have exactly one match") #check the fortran might actually be able to handle this

     x<-x[order(x[,1],x[,2]),]
     x<-cbind(x[x[,1]==group.names[1],3],x[x[,1]==group.names[2],3])
      } else{
        if(dim(x)[2]==2){
          names.used<-rep("",times=3)
          group.names<-c(1,2)
          block.names<-seq(from=1,to=dim(x)[1])
          } else stop("Input data must have either two or three columns")
          }

      AlignVals=apply(x[,1:2],1,median)

      
 if(missing(number.perms)){
  do.resample<-FALSE;number.perms<-0
  } else {if(exact==TRUE) {warning("Setting number of permutations is not valid for exact ptmp, a permutation test will be performed")
                          do.resample=TRUE
                          exact=FALSE}
          else do.resample<-TRUE
          }

      if(missing(save.test)) save.test<-FALSE

#remove nonzero differences
OrigNumPairs<-nrow(x)
  x<-x[(x[,1]-x[,2])!=0,]
  
NumPairs<-dim(x)[1]
  if(NumPairs>20 && exact==TRUE){
          ans<-1

          while(ans!="y" & ans!="Y" & ans!="yes" & ans!="Yes" & ans!="YES" & ans!="T" & ans!="TRUE"){

          ans<-readline(prompt=paste("Computation times for exact permutation test for matched pairs procedure can be high for more than 20 pairs\n",
            "Do you wish to continue? (y/n) \n",sep=" "))
            if(ans=="n" | ans=="N" | ans=="no"| ans=="No"| ans=="NO"| ans=="F"| ans=="FALSE") return()

              }
        }

da_STV<-rep(0,times=number.perms)
x<-as.matrix(x)
storage.mode(x)<-"double"
storage.mode(da_STV)<-"double"

    if(NumPairs<3) stop(paste("There are only", NumPairs,
                        "nonzero differences in the data set.  There must be at least 3 for PTMP"))
    if(i_NumObs>2147483647) stop("Number of cases (observations) should not exceed 2147483647 for PTMP")

ptmpOut=.Fortran("wrapptmp",
          as.integer(NumPairs),
          as.double(expon),
          as.integer(1),
          as.double(0),
          as.double(x[,1]),
          as.double(x[,2]),
          as.double(0),
          as.double(0),
          as.double(0),
          as.double(0),
          as.double(0),
          as.double(0),
          as.double(0),
          as.logical(do.resample),
          as.logical(exact),
          as.integer(number.perms),
          as.integer(0),
          as.logical(save.test),
          da_STV
          )

ptmpOut<-new("PTMPObj",NumPairs=OrigNumPairs,expon=ptmpOut[[2]],
    Data1=ptmpOut[[5]],Data2=ptmpOut[[6]],ExpectDelta=ptmpOut[[7]],DeltaVar=ptmpOut[[8]],
    DeltaSkew=ptmpOut[[9]],StdStat=ptmpOut[[10]],Rho=ptmpOut[[11]],ObsDelta=ptmpOut[[12]],
    P_value=ptmpOut[[13]],Resample=ptmpOut[[14]],Exact=ptmpOut[[15]],NumPerms=ptmpOut[[16]],
    Call=deparse(Call,width.cutoff=200L),SaveTest=ptmpOut[[18]],PermVals=ptmpOut[[19]],
    NamesUsed=as.character(names.used),GroupNames=as.character(group.names),BlockNames=as.character(block.names),AlignVals=AlignVals)
      
ptmpOut
}


