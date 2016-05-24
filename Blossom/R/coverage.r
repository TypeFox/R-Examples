coverage<-function(variable,group,expon=1,interv=0.0,number.perms,exact=FALSE,save.test,data){
      #some default settings for lad regression with no hypothesis testing and such
 Call<-match.call()
          if(missing(save.test)) save.test=FALSE
          if(missing(number.perms)) number.perms=ifelse(exact,0,4000)
       DNAME <- deparse(substitute(variable))
       if(!missing(data)) variable<-data[,match(DNAME,names(data))]
        if(!missing(group)) {GNAME<-deparse(substitute(group))
         if(!missing(data)) group<-data[,match(GNAME,names(data))]
        } else group=rep(1,times=length(variable))
      
        if (any(is.na(group)) | any(is.na(variable)))
            stop("NA's are not allowed")
          
        if (any(diff(c(length(variable), length(group))) !=
            0L))
            stop("variable and groups must have the same length")
        DNAME <- paste(DNAME, ", ", deparse(substitute(group)),sep = "")
       
          #remove incomplete cases
     i_NumObs=length(variable)
       variable<-variable[complete.cases(variable)]
      comp.cases<-length(variable)
      if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
        " cases were removed because of missing values."),sep=""))
           i_NumObs<-nrow(variable)   
        }
        
        group <- factor(group)
        NumGroups <- nlevels(group)
        if(NumGroups==1) {o<-order(variable)
        } else o <- order(group)
        
        variable <- variable[o]
        group <- group[o]
 #probably need something in here for categorical group names to be remembered in Summary, print methods
#check that the formula only has two
#check for missing data and then warn. Offer option to continue
 if(NumGroups==1) GrpSizes<-length(group)
 if(NumGroups>1)   GrpSizes=as.vector(table(group))
 
 if(i_NumObs>2147483647 & NumGroups==1) stop("Number of cases (observations) should not exceed 2147483647 for Kendal-Sherman Goodness of Fit")
    if(i_NumObs>2147483647) stop("Number of cases (observations) should not exceed 2147483647 for G-sample emperical coverage test")
    if(NumGroups>2147483647) stop("Number of groups should not exceed 2147483647 for G-sample emperical coverage test")
    if(min(GrpSizes)<2) stop(paste("Minimum group size is 2.  Your data set contains a group of size", min(GrpSizes)))
    
if(NumGroups==1){
DoArc<-interv!=0
 
    storage.mode(variable)="double"
    d_TObs<-d_Test<-d_VarT<-d_SkwT<-d_Z<-d_PValue<-0
    l_OkZ<-l_OkSkwT<-TRUE
    iEr<-0
    ch_MyMsg=""
    i_NumCases<-length(variable)
    coverageOut=.Fortran("wrapksgf",
              variable,
              as.double(d_TObs),
              as.double(d_Test),
              as.double(d_VarT),
              as.logical(l_OkZ),
              as.double(d_Z),
              as.logical(l_OkSkwT),
              as.double(d_SkwT),
              as.double(d_PValue),
              as.integer(iEr),
              as.double(interv),
              as.logical(DoArc),
              as.integer(i_NumCases))
      CoverageOut<-new("KSGFObj",inputData=coverageOut[[1]],ObsDelta=coverageOut[[2]],ExpectDelta=coverageOut[[3]],
            VarDelta=coverageOut[[4]],StandDelta=coverageOut[[6]],DeltaSkew=coverageOut[[8]],
              P_value=coverageOut[[9]],DoArc=coverageOut[[12]],ArcInterv=coverageOut[[11]],NumCases=coverageOut[[13]],
              Call=deparse(Call,width.cutoff=200L))
               }
               else{
     
      STV<-rep(0,times=(number.perms*save.test))
      GpVals<-as.vector(unique(group))
       
    storage.mode(GrpSizes)<-"integer"
    storage.mode(variable)<-"double"
    storage.mode(GpVals)<-"double"
    storage.mode(STV)="double"
       d_To<-d_Delta<-d_ET<-d_VT<-d_SDVT<-d_Z<-d_SKT<-d_PVal<-d_PZ<-0
       l_OkZ<-l_OkSDVT<-l_OkSKT<-l_OkPVal<-TRUE
       iEr=0
       ch_MyMsg=""
        numObs<-sum(GrpSizes)
        grpLength<-ifelse(exact,NumGroups+1,NumGroups) 
        
 

    coverageOut=.Fortran("wrapCov",
              as.integer(number.perms),
              as.double(expon),
              as.integer(NumGroups),
              GrpSizes,
              variable,
              as.double(d_To),
              as.double(d_ET),
              as.double(d_VT),
              as.logical(l_OkSDVT),
              as.double(d_SDVT),
              as.logical(l_OkZ),
              as.double(d_Z),
              as.logical(l_OkSKT),
              as.double(d_SKT),
              as.logical(l_OkPVal),
              as.double(d_PVal),
              as.double(d_PZ),
              as.integer(iEr),
              as.logical(save.test),
              STV,
              as.integer(NumGroups),
              as.logical(exact),
              GpVals,
              as.integer(numObs),
              as.integer(grpLength),
              as.double(d_Delta))
           
              CoverageOut<-new("CoverageObj",NumPerms=coverageOut[[1]],DistExp=coverageOut[[2]],NumGrps=NumGroups,
              GpSizes=coverageOut[[4]],inputData=variable,ObsDelta=coverageOut[[26]],VarDelta=coverageOut[[8]],ExpectDelta=coverageOut[[7]],
              DeltaSkew=coverageOut[[10]],Z_value=coverageOut[[11]]*coverageOut[[12]],Skt=coverageOut[[14]],
              P_value=coverageOut[[15]]*coverageOut[[16]],PZ=coverageOut[[17]],
              NumObs=length(variable),PermVals=coverageOut[[20]],exact=coverageOut[[22]],group.names=coverageOut[[23]],
              Call=deparse(Call,width.cutoff=200L))
              }

CoverageOut
}