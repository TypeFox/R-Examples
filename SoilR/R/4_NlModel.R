#
# vim:set ff=unix expandtab ts=2 sw=2:
correctnessOfNlModel=function
### The parameters used have a biological meaning, and therefore cannot be arbitrary.
### This functions tests some of the obvious constraints of SOM (Soil Organinc Matter) decomposition models. 
### Up to now these are:
### 1) The inputrates that must be positive for the whole time 
### 
### 2) The compatibility of the time ranges of the supplied functions which is important 

(object)
{   
    times=object@times
    Atm=object@DepComp
    ivList=object@initialValues
    InFluxes=object@inputFluxes
    res=TRUE
     
    tI_min=getTimeRange(InFluxes)["t_min"]
    tI_max=getTimeRange(InFluxes)["t_max"]
    t_min=min(times)
    t_max=max(times)
    if (t_min<tI_min) {
        stop(simpleError("You ordered a timeinterval that starts earlier than the interval your function I(t) (InFluxes) is defined for. \n Have look at the timeMap object of I(t) or the data it is created from")
        )
    }
    if (t_max>tI_max) {
        stop(simpleError("You ordered a timeinterval that ends later than the interval your function I(t) (InFluxes) is defined for. \n Have look at the timeMap object of I(t) or the data it is created from")
        )
    }
    print("tests passed")
    return(res)
}
is.negative=function(number){
   ### the function returns True if the argumente is negative
   return(number<0)
}
### serves as a fence to the interface of SoilR functions. So that later implementations can differ	 
setClass(# NlModel
   Class="NlModel",
   representation=representation(
        times="numeric"
        ,
        DepComp="TransportDecompositionOperator"
        ,
        initialValues="numeric"
        ,
        inputFluxes="BoundInFlux"
        ,
        solverfunc="function"
   )
   , validity=correctnessOfNlModel #set the validating function
)



setMethod(
    f="initialize",
    signature="NlModel",
    definition=function(
        .Object,
        times=c(0,1),
        DepComp=TransportDecompositionOperator.new(
                0,
                1,
                function(t){
                    return(matrix(nrow=1,ncol=1,0))
                },
                function(t){
                    return(matrix(nrow=1,ncol=1,0))
                }
        ) 
        ,
        initialValues=numeric()
        ,
        inputFluxes= BoundInFlux(
            function(t){
                return(matrix(nrow=1,ncol=1,1))
            },
            0,
            1
        )
        ,
        solverfunc=deSolve.lsoda.wrapper
        ,
        pass=FALSE
        ){
       # cat("-initializer at work-\n")
        .Object@times=times
        .Object@DepComp=DepComp
         if (class(inputFluxes)=="TimeMap"){
          warning(
            "The use of object of class TimeMap for the inputFlux argument is deprecated.
            At the moment we cast TimeMap objects to the new class BoundInFlux
            which replaces TimeMap as class of the the inputFlux argument.
            To get rid of this warning adapt your code to use a BoundInFlux instead of a TimeMap.
            Other classes may be implemented in the future." 
            )
            # cast
            inputFluxes<- BoundInFlux(inputFluxes)
         }
        .Object@initialValues=initialValues
        .Object@inputFluxes=inputFluxes
        .Object@solverfunc=solverfunc
        #if (pass==FALSE) validObject(.Object) #call of the ispector if not explicitly disabled
        if (pass==FALSE) correctnessOfNlModel(.Object) #call of the ispector if not explicitly disabled
        return(.Object)
    }
)

#################################################
setMethod(
   f= "getInFluxes",
   signature(object="NlModel"),
   definition=function(object){
       return(object@inputFluxes)
     }
)
#################################################
setMethod(
   f= "plot",
      signature(x="NlModel"),
      definition=function(x){
      ### This function is a stub
      # It only starts the thing ...    
      plot(getTimes(x),getC(x)[,1])
   }
)
#################################################
setMethod(
   f= "print",
      signature(x="NlModel"),
      definition=function(x){
      ### This function is a stub
      # It only starts the thing ...    
      print("Hi there I am the method print for model objects. Change me if you can")
      print(getC(x)[,1])
   }
)
#################################################
setMethod(
   f= "getNumberOfPools",
      signature(object="NlModel"),
      definition=function(object){
      return(length(object@initialValues))
   }
)
#################################################
setMethod(
   f= "summary",
      signature(object="NlModel"),
      definition=function(object){
      ### This function is a stub
      # It only starts the thing ...    
      print("Hi there, I am the method summarize for model objects. 
            I summarize everything:\n
            1.) I have not benn fully implemented yet.\n
            2.) here comes the C stock at least.")
      print(getC(object)[,1])
   }
)
#################################################
setMethod(
   f= "show",
      signature(object="NlModel"),
      definition=function(object){
      ### This function is a stub
      # It only starts the thing ...    
      print("Hi there I am the method show for model objects")
      print(getC(object)[,1])
   }
)
#################################################
setMethod(
   f= "getDecompOp",
      signature= "NlModel",
      definition=function(object){
      ### This method creates a particle simulator 
      return(object@DepComp)
   }
)
#################################################
setMethod(
   f= "getParticleMonteCarloSimulator",
      signature= "NlModel",
      definition=function(object){
      ### This method creates a particle simulator 
      return(new(Class="MCSim",object))
   }
)
#################################################
setMethod(
   f= "getTimes",
      signature= "NlModel",
      definition=function(object){
      ### This functions extracts the times argument from an object of class NlModel
         times=matrix(ncol=1,object@times)
         colnames(times)="times"
      return(times)
   }
)
#################################################
setMethod(
   f= "getInitialValues",
      signature= "NlModel",
      definition=function(object){
      ### This functions extracts the initial values argument from an object of class NlModel
      return(object@initialValues)
   }
)
#################################################
res2fun=function(times,C){
  cn=colnames(C)
  if (length(cn)==0){cn=1:nrow(C)}
  Cs=list()#create an indexed list of functions
  for (i in 1:ncol(C)){
    key=cn[[i]]
    Cs[[key]] <- splinefun(times,C[,i])
  }
  return(Cs)
}
################################################
setMethod(
   f= "getTransferCoefficients",
      signature= "NlModel",
      definition=function # the transfer coefficients for all pool combinations as functions of time
      (object,as.closures=F){
      ### This functions returns the transfer coefficients for all pool combinations. 
      C=getC(object)
      #pp("C",environment())
      t=object@times
      #pp("t",environment())
      DepComp=object@DepComp
      alpha=getTransferCoefficients(DepComp)
      if(length(names(alpha))==1){
        all_tr=matrix(ncol=1,mapply(alpha[[1]],(1:length(t))))
        colnames(all_tr)=names(alpha)
      }else{
        single_tr=function(i){
          l=list()
          for (key in names(alpha)){
            force(key)
            l[[key]]=alpha[[key]](C[i,],t[[i]])
          }
          return(l)
        }
        all_tr=t(mapply(single_tr,(1:length(t))))
      }
      
      if (as.closures){
        return(res2fun(t,all_tr))
      }else{
        return(all_tr)
      }
   }
)
################################################
setMethod(
   f= "getOutputFluxes",
      signature= "NlModel",
      definition=function # complete output of all pools, including the part transfered to other pools.
      (object,as.closures=F){
      ### This functions computes the output flux for all pools. Note that for any given pool not all the output of the pool is released from the system because it migtht as well be fed into other pools. If you are interested what a pool releases from the system use the method \code{\link{getReleaseFlux}}, which internally makes use of this method but afterwards substracts all parts of the outputs  that are fed to other pools.
      C=getC(object)
      t=object@times
      DepComp=object@DepComp
      DotO=getDotOut(DepComp) #this is a function of C and t 
      single_O=function(i){DotO(C[i,],t[[i]])}
      all_O=t(mapply(single_O,(1:length(t))))
      if (as.closures){
        return(res2fun(t,all_O))
      }else{
        return(all_O)
      }
   }
)
#################################################
setMethod(
   f= "computeResults",
      signature= "NlModel",
      definition=function(object){

}
)
#################################################
setMethod(
   f= "getC",
      signature= "NlModel",
      definition=function(object,as.closures=F){
      times=object@times
      DepComp=object@DepComp
      DotO=getDotOut(DepComp) #this is a function of C and t 
      Tr=getTransferMatrix(DepComp) #this is also a function of C and t 
      force(DotO) 
      itm=object@inputFluxes
      inputrates=getFunctionDefinition(itm)
      #print(input)
      DotC=function(C,t){
        return(Tr(C,t)%*%DotO(C,t)+inputrates(t))
        #return(DotO(C,t)+inputrates(t))
      }
      sVmat=matrix(object@initialValues,ncol=1)
      C=solver(times,DotC,sVmat,object@solverfunc) 
      ### A matrix. Every column represents a pool and every row a point in time
      #f=function(i){paste("C",i,sep="")}
      #colnames(C)=sapply((1:ncol(C)),f)
      if (as.closures){
        return(res2fun(times,C))
      }else{
        return(C)
      }
   }
)
#################################################
# overload the [] operator
setMethod("[",signature(x="NlModel",i="character"), #since [] is a already defined generic the names of the arguments are not arbitrary 
        definition=function(x,i){
            getSingleCol=function(slot_name){
                res=""
                #print(paste(sep="",">",slot_name,"<"))
                if(slot_name=="times"){ res=getTimes(x)}
                if(slot_name=="C"){ res=getC(x)}
                #if(slot_name=="ReleaseFlux"){ res=getReleaseFlux(x)}
                #if(slot_name=="AccumulatedRelease"){ res=getAccumulatedRelease(x)}
                #if(res==""){stop(paste("The slot",slot_name,"is not defined"))}
                return(res)
            }
            n=length(i)
            df=getSingleCol(i[1])
            if (n>1){
                for (k in 2:n){
                    df=cbind(df,getSingleCol(i[k]))
                }
            }
            return(df)
        }
)
#################################################
# overload the $ operator
setMethod("$",signature(x="NlModel"), #since $ is a already defined generic the names of the arguments are not arbitrary 
        definition=function(x,name){
            return(x[name])
        }
)

