#############################################################
##### LAD Object
setClass("LADObj",
  representation=list(inputData="matrix",NumObs="numeric",NumVars="numeric",
    theta="numeric", NumPerm="numeric",
    Test="numeric",DoublePermutation="logical",T_o="numeric",
    AsyRankScore="numeric",P_value="numeric",P_valueTN="numeric",
    Betas="vector",RedBetas="vector",SumAbsValRes="numeric", SumAbsValResRed="numeric",
    WtSumAbsDevsFulMod="numeric", WtSumAbsDevsRedMod="numeric",
	  NumIter="numeric",ExitCode="numeric",PermVals="vector",
    HasIntercept="vector",DoAllQuants="numeric",
    DoRankScore="logical",IsOLS="logical",
    NumLaVars="vector",ResRed="vector",Resids="vector",
    Call="call",response="matrix",full.mod.names="character",QuantOut="matrix"),
prototype=list(inputData=matrix(0,0,0),NumObs=0,NumVars=0,
  theta=0, NumPerm=0,
    Test=0,DoublePermutation=FALSE,T_o=0,
    AsyRankScore=0,P_value=0,P_valueTN=0,
    Betas=0,RedBetas=0,SumAbsValRes=0, SumAbsValResRed=0,
	  WtSumAbsDevsFulMod=0, WtSumAbsDevsRedMod=0,
	  NumIter=0,ExitCode=0,PermVals=0,
    HasIntercept=0,DoAllQuants=0,
    DoRankScore=FALSE,IsOLS=FALSE,
    NumLaVars=0,ResRed=0,Resids=0,response=matrix(0,0,0),
    full.mod.names="",QuantOut=matrix(0,0,0))
)

#######################################################
#### LAD Summary

setMethod("summary",
    signature(object = "LADObj"),
    function (object, ...)
    {

    if(object@DoAllQuants==1){
           cat("\tAll Quantile Regression \n")
           cat("\nCall: \n")
               cat(paste(strwrap(as.character(object@Call),width=70),collapse="\n\t"),"\n")
                cat("\n")

          cat("\nSpecification of Analysis:\n")
          cat("\tNumber of Observations:", object@NumObs, "\n")
          cat("\t    Dependent Variable:", colnames(object@response), "\n")
          cat("\t  Number of Parameters:",dim(object@inputData)[2], "\n")
          cat("\t   Number of Solutions:",dim(object@QuantOut)[2], "\n")
          cat("\t   Solution Result Was:", switch(as.character(object@ExitCode),"0"="Successful","2"="Premature End",
                "7"="Too Many Solutions"))
          cat("\n")

           cat("\n\nOutput can be obtained using QuantValues()")
          cat("\n\nColumns are Quantile, Objective Function Solution, Predicted Y at X-Bar,\n",
            "followed by estimates of coefficients for independent variables. \n",
            "Each row describes an estimated conditional quantile function.\n\n")
          return(invisible())
     }
     cv <- c(.001,.01,.05,.1,1)
     name.vect<-colnames(object@inputData)
    if(object@theta>0) cat("\t       Quantile Regression")
    else {if(!object@IsOLS) cat("\t    Least Absolute Deviation Regression (LAD)")
          else cat("\tOrdinary Least Squares Regression")
    }
    if((length(grep("hypothesis.test",object@Call))>0))  {cat("\n\t    Hypothesis Test")
        if(object@DoRankScore==FALSE &  object@IsOLS==FALSE & object@NumLaVars[2]>1) cat(", drop p-q-1 zero residuals")
        if(object@DoRankScore==TRUE) cat(" of Rank Score")
    }
    if(object@DoublePermutation==TRUE){ cat(" with double permutation\n")
                         } else cat("\n")
                         
    cat("\nCall: \n")
       cat(paste(strwrap(deparse(object@Call,width.cutoff=200L),width=70),collapse="\n\t"),"\n")
      cat("\n")

  cat("\nSpecification of Analysis:")

         cat("\n\tNumber of Observations:", object@NumObs, "\n")
           cat("\tResponse Variable     :", colnames(object@response), "\n")
           if(object@theta>0){
           cat("\tFor Quantile          :", object@theta, "\n\n")
            } else cat("\n")

             var.table<-as.data.frame(cbind(colnames(object@inputData),signif(as.data.frame(object@RedBetas),digits=getOption("digits"))))
             colnames(var.table)<-c("  Independent variables","Regression coefficients")
              print(var.table,print.gap=6,quote=FALSE,row.names=FALSE,justify="right")

          if(!object@IsOLS) {cat("\n\tNumber of iterations:",object@NumIter)
            cat("\n\n\tSum of absolute values of the residuals:",object@SumAbsValResRed) }
          if(object@theta>0 & !(length(grep("hypothesis.test",object@Call))>0))
            cat("\n\tWeighted sum of the absolute deviations:",object@WtSumAbsDevsFulMod)
          if(object@theta>0 & (length(grep("hypothesis.test",object@Call))>0))
            cat("\n\tWeighted sum of the absolute deviations:",object@WtSumAbsDevsRedMod)
          if(object@IsOLS)
            cat("\n\n\tSum of squares of the residuals:",object@SumAbsValResRed)

          if(!object@IsOLS) cat("\n\tSolution:", switch(as.character(object@ExitCode),"1"="Successful","2"="Rounding Error",
                  "0"="Non-Unique","7"="Multiple Solutions"))
          cat("\n")
          if((length(grep("hypothesis.test",object@Call))>0) | (object@Call[[1]]=="lad" & object@Test==TRUE)){
          
            cat("\nRegression Evaluation:")
                      if((length(grep("hypothesis.test",object@Call))>0) & object@theta>0)
                      cat(paste("\n\t", object@theta,sep=""),"Quantile Regression Model:")
                      else {if(object@IsOLS==TRUE) cat("\n\tOrdinary Least Squares: ")
                      else cat("\n\tLAD Model: ")}
                      if((length(grep("hypothesis.test",object@Call))>0)){
                          cat("\n\t")
                          cat(colnames(object@response),"~")
                          cat(object@full.mod.names,sep="+")

                      cat("\n\tVersus Hypothesis Model")
                      if((length(grep("hypothesis.test",object@Call))>0) & object@theta>0) 
                      cat(" at Quantile", object@theta)
                      cat(":\n\t")
                          cat(colnames(object@response),"~")
                          cat(name.vect,sep="+")
                        }  else{
                          cat(colnames(object@response),"~")
                          cat(name.vect,sep="+")
                          }
           cat("\n\nTest Summary:")
            cv <- c(.001,.01,.05,.1,1)
                              cat("\n\tNumber of permutations                    : ",object@NumPerm)
                        if((length(grep("hypothesis.test",object@Call))>0)) {
                            if(object@DoRankScore==FALSE) {
                              cat("\n\tObserved Test Statistic                   : ",object@T_o)
                              cat("\n\tP-value of variables in full model but ")
                              cat("\n\tnot reduced model                        : ",object@P_value)
                               switch(min(which(cv>=object@P_value,arr.ind=TRUE)),
                          cat("***\n"),
                          cat("**\n"),
                          cat("*\n"),
                          cat(".\n"),
                          cat("\n")
                        )
                            }
                            if(object@DoRankScore==TRUE) {
                            cat("\n\n\tObserved Rank Score Test Statistic        : ",object@T_o)
                              cat("\n\tP-value of Rank Score Test                : ",object@P_value)
                               switch(min(which(cv>=object@P_value,arr.ind=TRUE)),
                          cat("***\n"),
                          cat("**\n"),
                          cat("*\n"),
                          cat(".\n"),
                          cat("\n")
                        )
                              cat("\n\tAsymptotic Rank Score Statistic           : ",object@AsyRankScore)
                              cat("\n\t(Distributed as Chi-square with degrees of",
                                  "\n\tfreedom equal to difference in number of",
                                  "\n\tparameters between full and reduced models)")
                              cat("\n\tP-Value of Asymptotic RS Stat             : ",object@P_valueTN)
                               switch(min(which(cv>=object@P_valueTN,arr.ind=TRUE)),
                          cat("***\n"),
                          cat("**\n"),
                          cat("*\n"),
                          cat(".\n"),
                          cat("\n")
                        )       
                        Signif <- symnum(cv, corr = FALSE, na = FALSE,
                                          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                          symbols = c("***", "**", "*", ".", " "))
                                cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
                               return(invisible())
                            }
                        }else{
                              cat("\n\tP-value of Full Model                     : ",object@P_value)
                               switch(min(which(cv>=object@P_value,arr.ind=TRUE)),
                            cat("***\n"),
                            cat("**\n"),
                            cat("*\n"),
                            cat(".\n"),
                            cat("\n")
                          )
                        }
                             Signif <- symnum(cv, corr = FALSE, na = FALSE,
                                          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                          symbols = c("***", "**", "*", ".", " "))
                                cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
                                }
    }
)

#################################################################
#### LAD Print


setMethod("print", "LADObj", function(x) {
summary(x)
})

setMethod("show", "LADObj", function(object) {
print(object)
})

# Now setting up just a few generic C style accessor functions as recommended
# in the S4 class documentation

if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

# if (!isGeneric("residuals")) {
#  setGeneric("residuals", function(x){
#  standardGeneric("residuals")
# })
# }
 
# if (!isGeneric("predict")) {
#  setGeneric("predict", function(x){
#  standardGeneric("predict")
# })
# }
 
#if (!isGeneric("coefficients")) {
#  setGeneric("coefficients", function(x){
#  standardGeneric("coefficients")
# })
# }
 

setMethod("pvalue",signature(x = "LADObj"), function(x) {
 x@P_value
 })


setMethod("residuals", signature(object = "LADObj"),
	  function(object, ...)
           object@Resids)


 setMethod("predict", signature(object = "LADObj"), function(object) {
 as.vector(object@response-object@Resids)
 })
 
 setMethod("coefficients",signature(object = "LADObj"), function(object) {
object@Betas
 })


if (!isGeneric("ResampVals")) {
  setGeneric("ResampVals", function(x){
  standardGeneric("ResampVals")
 })
 }

setMethod("ResampVals", signature(x = "LADObj"), function(x) {
 x@PermVals
 })



if (!isGeneric("QuantValues")) {
  setGeneric("QuantValues", function(x){
  standardGeneric("QuantValues")
 })
 }


setMethod("QuantValues",signature(x = "LADObj"), function(x) {
 if(x@DoAllQuants!=1) stop("QuantValues() is not a valid command for the LAD analysis you have performed")
 output<-t(x@QuantOut[c(1,3,2,4:(dim(x@QuantOut)[1]-1+as.numeric(x@HasIntercept[2]))),])
          rownames(output)=rep("",times=dim(output)[1])
          colnames(output)<-c("Quantile", "ObjFuncSol", "PredY_Xbar", paste("b_",colnames(x@inputData),sep=""))
          return(output)
 })


