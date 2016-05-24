mrpp<-function(variables,group,data, expon=1,c.form=1,hotelling=FALSE,commens=TRUE,
      interv=0.0,number.perms,exact=FALSE,has.excess=FALSE,excess.value,max.dist,save.test){
     
      Call<-match.call()
        
      if(missing(data)) {
                         x<-cbind(eval(Call$group),eval(Call$variables))
                         variables<-colnames(eval(Call$variables))
                         group<-colnames(eval(Call$group))
                         if(is.null(group)) group="group"
                           colnames(x)<-c(group,variables)
      }else if(missing(variables) & missing(group)){
       x<-eval(Call$data)
      } else{
         data=eval(Call$data)                   
      variables<-as.character(Call$variables)
        if(length(variables)>1) variables<-variables[-1] 
      group<-as.character(Call$group)
        if((length(variables)==0) | (length(group)==0)) stop("Both of variables and group must be supplied")
        if(is.list(variables)) variables<-unlist(variables)
            x<-data[,match(c(group,variables),names(data))]
      }
  #Switching data to the Fortran format
            x[,1]<-as.character(x[,1])
      if(missing(save.test)) save.test<-FALSE
      
      if(missing(max.dist)) max.dist=0

      if(missing(number.perms)){
        do.resample<-FALSE;number.perms<-0
        } else {if(exact==TRUE) {warning("exact=TRUE cannot be used when number.perms is set, a permutation test will be performed")
                          do.reesample=TRUE
                          exact=FALSE}
          else do.resample<-TRUE
          }
          
      x<-as.data.frame(x)
      x<-x[order(x[,1]),]
      x.input<-x
        i_NumObs=nrow(x)
         x<-x[complete.cases(x),]
        comp.cases<-nrow(x)
        if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
          " cases were removed because of missing values."),sep=""))
        } 
      group.names<-unique(x[,1])
      group.sizes<-as.vector(table(x[,1]))
        
      x[,1]<-as.numeric(as.factor(x[,1]))
      if(has.excess==TRUE | !missing(excess.value)){ #Excess group is assumed to be the largest number
            if(missing(excess.value)) excess.value<-max(x[,1])
            has.excess=TRUE
            excessInd<-x[,1]==excess.value
            x<-rbind(x[!excessInd,],x[excessInd,])
            trueGroups<-x[,1]
            fortranGroupLabs<-rep(seq(1:length(unique(x[,1]))),
                 times=as.vector(table(x[,1]))[match(names(table(x[,1])),unique(x[,1]))])

            x[,1]<-fortranGroupLabs
       } else excess.value<-0
      
      x<-x[,c(2:ncol(x),1)]

      i_NumObs=nrow(x)
      i_NumVars=ncol(x)-1
        if(i_NumVars==1) Commensurate=FALSE
      i_NumGrps=ifelse(has.excess==FALSE,length(unique(x[,ncol(x)])),(length(unique(x[,ncol(x)]))-1))
      ia_GrpSizes=as.vector(c(table(x[,ncol(x)]),0))
      if(has.excess) ia_GrpSizes[length(ia_GrpSizes)]<-ia_GrpSizes[length(ia_GrpSizes)-1]
      da_GpVals<-sort(unique(x[,ncol(x)]))
      ia_GrpValTag<-rep(1,times=ncol(x)) # might want to include the option of excluding some values here
      da_GroupV<-x[,ncol(x)]
      da_STV<-rep(0,times=number.perms)


      # output is initilized here

      da_YHot=matrix(0.0,i_NumVars,i_NumVars)
      da_XI1=rep(0,times=i_NumGrps)
      CommAvgDist=rep(0,times=i_NumVars)

       x<-as.matrix(x)
       storage.mode(x) <- "double"
       storage.mode(da_YHot) <- "double"
       storage.mode(da_XI1) <- "double"
       storage.mode(ia_GrpValTag)<-"integer"
       storage.mode(da_GroupV)<-"double"
       storage.mode(da_STV)="double"
       storage.mode(CommAvgDist)="double"
       storage.mode(da_GpVals)="double"

    if(min(ia_GrpSizes[1:(i_NumGrps+has.excess)]) < 2) stop("Groups must have 2 or more observations to do an MRPP or EMRPP")
    if((i_NumGrps+has.excess) < 2) stop("There must be at least 2 grouping variable values in order to do an MRPP")
    if(i_NumObs<6 & !exact) stop(paste("There are only", i_NumObs, "cases.  There must be at least 6 to do MRPP",
                                          "you could try an EMRPP which only requires 3 cases"))
    if(i_NumObs<3 & exact) stop(paste("There are only", i_NumObs, "cases.  There must be at least 3 to do an EMRPP"))
    if(i_NumObs>2147483647) stop("Number of cases (observations) should not exceed 2147483647 for MRPP or EMRPP")
    if(i_NumGrps>2147483647) stop("Number of groups should not exceed 2147483647 for MRPP or EMRPP")
    if(i_NumVars>65536) stop("Number of variables should not exceed 65536 for MRPP or EMRPP")
    GrpSizeLen<-ifelse(has.excess,i_NumGrps+2,i_NumGrps+1)  
       
    if(!exact){

            mrppOut=.Fortran("wrapmrpp",
              as.integer(i_NumObs),
              as.integer(i_NumVars),
              as.double(expon),
              as.double(max.dist),
              as.integer(c.form),
              as.logical(hotelling),
              as.logical(commens),
              as.integer(i_NumGrps),
              as.double(interv),
              as.integer(ia_GrpSizes),
              as.integer(number.perms),
              as.logical(do.resample),
              x,
              da_XI1,
              as.double(0),
              as.double(0),
              as.double(0),
              as.double(0),
              as.double(0),
              as.double(0),
              da_YHot,
              as.integer(0),
              as.double(excess.value),
              as.logical(has.excess),
              da_GpVals,
              ia_GrpValTag,
              da_GroupV,
              as.logical(save.test),
              da_STV,
              CommAvgDist,
              as.integer(GrpSizeLen)
            )
         
        mrppOut<-new("MRPPObj",NumObs=mrppOut[[1]],NumVars=mrppOut[[2]],
        DistExp=mrppOut[[3]],MaxDist=mrppOut[[4]],CForm=c.form,Hotelling=mrppOut[[6]],Commens=mrppOut[[7]],NumGrps=mrppOut[[8]],Interval=mrppOut[[9]],
        GpSizes=group.sizes,NumPerm=mrppOut[[11]],DoResamp=mrppOut[[12]],inputData=x.input,AvgDist=mrppOut[[14]],
        StandTestStat=mrppOut[[15]],ObsDelta=mrppOut[[16]],ExpectDelta=mrppOut[[17]],DeltaVar=mrppOut[[18]],
        DeltaSkew=mrppOut[[19]],P_value=mrppOut[[20]],YHot=mrppOut[[21]],
        d_ExcessVal=excess.value,l_HasExcess=mrppOut[[24]],da_GpVals=mrppOut[[25]],ia_GrpValTag=mrppOut[[26]],
        da_GroupV=mrppOut[[27]],Call=deparse(Call,width.cutoff=200L),group.names=group.names,SaveTest=mrppOut[[28]],PermVals=mrppOut[[29]],
        CommAvgDist=mrppOut[[30]])


    } else {

      i_ActiveCases=i_NumObs #this might be different for the excess group so I'm leaving this in for now
    # output is initilized here
        
       mrppOut=.Fortran("wrapemrpp",
              x,
              as.integer(ia_GrpSizes),
              as.logical(has.excess),
              as.double(excess.value),
              as.logical(hotelling),
              as.logical(commens),
              as.double(max.dist),
              as.double(interv),
              as.double(expon),
              as.integer(c.form),
              as.integer(i_NumGrps),
              as.integer(i_NumObs),
              as.integer(i_NumVars),
              da_GpVals,
              as.integer(i_ActiveCases),
              da_GroupV,
              as.double(0),
              da_YHot,
              as.double(0),
              as.integer(0),
              CommAvgDist,
              as.integer(GrpSizeLen))
             
              if(is.null(mrppOut[[3]])) mrppOut[[3]]<-FALSE
    mrppOut<-new("EMRPPObj",inputData=x.input,GpSizes=group.sizes,l_HasExcess=mrppOut[[3]],
          d_ExcessVal=excess.value,Hotelling=mrppOut[[5]],Commens=mrppOut[[6]],MaxDist=mrppOut[[7]],Interval=mrppOut[[8]],
          DistExp=mrppOut[[9]],CForm=mrppOut[[10]],NumGrps=mrppOut[[11]],NumObs=mrppOut[[12]],NumVars=mrppOut[[13]],da_GpVals=mrppOut[[14]],
          da_GroupV=mrppOut[[16]],ObsDelta=mrppOut[[17]],YHot=mrppOut[[18]],P_value=mrppOut[[19]],
          Call=deparse(Call,width.cutoff=200L),group.names=group.names,CommAvgDist=mrppOut[[21]])

      }

    mrppOut
  }


