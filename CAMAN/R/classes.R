setClass("CAMAN.object", representation(dat="matrix", family="character", 
                                        LL="numeric", num.k="numeric", p="numeric", t="numeric",  component.var = "numeric", 
                                        prob="matrix", classification = "numeric", num.obs="numeric", steps = "numeric", 
                                        otherParams = "numeric", BIC = "numeric", VEM_result = "matrix", finalacc = "numeric", cl =  "call", is_metaAnalysis = "numeric"))

setClass("CAMAN.VEM.object", representation(dat="matrix", family="character", 
                                            LL="numeric", num.k="numeric", startk="numeric", p="numeric", t="numeric",  
                                            num.obs="numeric", steps = "numeric", BIC = "numeric", finalacc = "numeric",
                                            otherParams = "numeric", cl =  "call", is_metaAnalysis = "numeric", grid="data.frame", totalgrid="data.frame"))


setClass("CAMAN.glm.object", representation(dat="data.frame", family="character", 
                                            LL="numeric", num.k="numeric", p="numeric", t="numeric",  hetvar = "numeric", prob="matrix", 
                                            classification = "numeric", num.obs="numeric", steps = "numeric", otherParams = "numeric", 
                                            BIC = "numeric", coefMatrix = "data.frame", commonEffect = "numeric", cl =  "call", fittedObs="numeric",
                                            numPara = "numeric", depVar = "character", fixedVar = "character", random = "character",
                                            form="formula", glmModel = "glm", mode="character", residVar= "numeric", idxControl="list", inputData = "data.frame"))



setClass("CAMAN.BIVEM.object", representation(RESULT="matrix", RESULT_uni="matrix",RESULT_meta="matrix",
                                           BIC = "numeric", LL = "numeric"))  


setClass("CAMAN.BIEM.object", representation( RESULT="matrix", Mat="matrix", Z="array",
                                           BIC = "numeric", LL = "numeric",cl="numeric"))   

setClass("CAMAN.BIMIXALG.object", representation(RESULT="matrix", RESULT_uni="matrix",RESULT_meta="matrix", Z="array",
                                           BIC = "numeric", LL = "numeric",cl="numeric",Mat="matrix"
))  
                                           
setClass("CAMAN.BOOT.object", representation( H0="matrix", S1="array",H1="matrix", S2="array",
                                           LL= "numeric", Q95="numeric",Q975="numeric",Q99="numeric"))   

setMethod("show", "CAMAN.object", function(object){
  cat("Computer Assisted Mixture Analysis: \n \n")
  cat("Data consists of", object@num.obs, "observations (rows). \n")
  cat("The Mixture Analysis identified", object@num.k, "component")
  if (length(object@num.k) >0) cat("s")
  cat(" of a", object@family, "distribution: \n \n")
  
  n <- object@num.obs    
  details <- matrix(0, nrow=object@num.k, ncol=2)
  descr_var =""
  if (object@family == "gaussian") {
    tmp <- "mean"
    if (object@is_metaAnalysis == 0) descr_var <- paste("component variance:", object@component.var, "\n")
  }
  else if (object@family == "poisson") tmp <- "lambda"
  else if (object@family == "binomial") tmp <- "prob"
  colnames(details) = c("p", tmp)
  rownames(details) = 1:object@num.k
  details[,1] <- object@p
  details[,2] <- object@t
  cat("DETAILS:\n")
  print(details)
  cat(descr_var,"\n")
  cat("Log-Likelihood:",object@LL,"    ")
  cat("BIC:",object@BIC,"\n")      
})

setMethod("show", "CAMAN.VEM.object", function(object){
  cat("Computer Assisted Mixture Analysis (VEM): \n \n")
  cat("Data consists of", object@num.obs, "observations (rows). \n")
  cat("The VEM-algorithm identified", nrow(object@grid), "grid point")
  if (length(object@num.k) >0) cat("s")
  cat(" with positive support \n \n")
  
  n <- object@num.obs    
  print(object@grid)
  
  cat("Log-Likelihood:",object@LL,"    ")
  cat("BIC:",object@BIC,"\n")      
})


setMethod("show", "CAMAN.glm.object", function(object){
  cat("Computer Assisted Mixture Analysis with covariates: \n \n")
  cat("Data consists of", object@num.obs, "observations (rows). \n")
  cat("The Mixture Analysis identified", object@num.k, "component")
  if (length(object@num.k) >0) cat("s")
  cat(" of a", object@family, "distribution: \n \n")
  
  n <- object@num.obs    
  cat("mixing weights:\n")
  p_tmp <- object@p
  names(p_tmp) = paste("comp.", 1:object@num.k)
  print(p_tmp)
  
  cat("\n Coefficients :\n")
  coefPrint <- round(object@coefMatrix[,1],3)
  names(coefPrint) <- rownames(object@coefMatrix)
  print(coefPrint)
  if (object@family=="gaussian") cat("residual variance:", object@residVar)
  cat("\n")
  
  cat("Log-Likelihood:",object@LL,"    ")
  cat("BIC:",object@BIC,"\n")      
})



setMethod("show","CAMAN.BIVEM.object", function(object){

  cat("Computer Assisted Mixture Analysis (BIVEM): \n \n")

if(length(object@RESULT_meta>1)){

cat("VEM algorithm for diagnostic meta analysis: \n \n ")
colnames(object@RESULT_meta)<- c("Lambda1","Lambda2","Prob")
cat("RESULT_meta: \n \n" ) 
print(object@RESULT_meta[,])}

else if(length(object@RESULT>1)){

cat("Vem for bivariate data: \n \n")  
colnames(object@RESULT)<- c("Lambda_1","Lambda_2","Prob")
cat("RESULT: \n \n" ) 
print(object@RESULT[,])
}

else if(length(object@RESULT_uni>1)){

cat("Vem for univariate data: \n \n")  
colnames(object@RESULT_uni)<- c("Lambda_1","Prob")
cat("RESULT_uni: \n \n" ) 
print(object@RESULT_uni[,])
}

  cat("\n","Log-Likelihood:",object@LL,"\n")
  cat("\n","BIC:",object@BIC,"\n")      
}
)
setMethod("show","CAMAN.BIEM.object", function(object){

cat("Computer Assisted Mixture Analysis (BIEM): \n \n")
if(length(object@Mat>1)){
cat("EM-algorithm for bivariate normally distributed data: \n \n")

colnames(object@RESULT) <- c("Lambda1","Lambda2","Prob","Var1","Var2","Corr")
cat("RESULT: \n \n" )
print(object@RESULT[,])
if (length(object@cl>1))cat("\n","cl:","\n","cl:","Classification of bivariate data with starting values","\n",object@cl)
cat("\n","LL:",object@LL,"\n")
  cat("\n","BIC:",object@BIC,"\n")    
}

else {
cat("EM-algorithm for bivariate diagnostic  meta analysis: \n \n")
colnames(object@RESULT) <- c("Lambda1","Lambda2","Prob")
cat("RESULT: \n \n" )
print(object@RESULT[,]) 
if (length(object@cl>1))cat("\n","cl:","\n","cl:","Classification of meta data with start values","\n",object@cl)

}
  cat("\n","LL:",object@LL,"\n")
  cat("\n","BIC:",object@BIC,"\n")      
})
setMethod("show","CAMAN.BIMIXALG.object", function(object){
cat("Computer Assisted Mixture Analysis (BIMIXALG): \n \n")

if(length(object@RESULT>1)){

cat("Combination of VEM- and EM-algorithm for bivariate normally distributed data: \n \n")
colnames(object@RESULT) <- c("Lambda1","Lambda_2","Prob","Var1","Var2","Corr")
cat("RESULT: \n \n" ) 
print(object@RESULT[,])

if (length(object@cl>1))
cat("\n","cl:","Classification of bivariate data ","\n",object@cl,"\n")
}


else if (length(object@RESULT_meta>1)){
cat("Combination of VEM- and EM-algorithm for bivariate diagnostic  meta analysis: \n \n")
colnames(object@RESULT_meta) <- c("Lambda1","Lambda2","Prob")
cat("RESULT_meta: \n \n" ) 
print(object@RESULT_meta[,]) 


if (length(object@cl>1))cat("\n","cl:","Classification of meta data",object@cl,"\n")}

else if(length(object@RESULT_uni>1)){

cat("Combination of VEM- and EM-algorithm for univariate data : \n \n")
  
colnames(object@RESULT_uni)<- c("Lambda_1","Prob","Var")
cat("RESULT_uni: \n \n" ) 
print(object@RESULT_uni[,])
if (length(object@cl>1))cat("\n","Classification of univariate data ",object@cl,"\n")
}


  cat("\n","LL:",object@LL,"\n")
  cat("\n","BIC:",object@BIC,"\n")      
})



setMethod("show","CAMAN.BOOT.object", function(object){
cat("Computer Assisted Mixture Analysis (BOOT): \n \n")

cat("Parametric Bootstrap: \n \n")
#colnames(object@H0) <- c("Lambda1","Lambda_2","Prob","Var1","Var2","Corr")
print(object@H0[,])

cat("\n","Covariance matrix:","\n")
print(object@S1)

print(object@H1)

 cat("\n","Covariance matrix:","\n")
print(object@S2)

  cat("\n","Log-Likelihood:",object@LL,"\n")
cat("\n","Quantile 0.95 :",object@Q95,"\n")
cat("\n","Quantile 0.975 :",object@Q975,"\n")
cat("\n","Quantile 0.99 :",object@Q99,"\n")
      
}
)
