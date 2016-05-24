Bvs <- function(formula,  fixed.cov=c("Intercept"), data, prior.betas="Robust", prior.models="Constant", n.keep=10, time.test=TRUE, priorprobs=NULL){
  
  #Let's define the result 
  result<- list()
  
  #Get a tempdir as working directory
  wd<- tempdir()
  #remove all the previous documents in the working directory
  unlink(paste(wd,"*",sep="/"))
  

#Set the design matrix if fixed covariates present: 
if(!is.null(fixed.cov)){
 
  #Eval the full model 
  lmfull = lm(formula, data=data, y=TRUE, x=TRUE)
  X.full<- lmfull$x
  namesx<- dimnames(X.full)[[2]]
  
  
  #remove the brackets in "(Intercept)" if present. 
  if(namesx[1]=="(Intercept)"){
    namesx[1]<- "Intercept" #namesx contains the name of variables including the intercept
  }
  
  #Eval of the null model
   response <- strsplit(formula,"~")[[1]][1]
 
  if(length(fixed.cov)==1){
    if("Intercept"%in% fixed.cov){
      formulanull=paste(response,"~1",sep="")
    }
    if(!"Intercept"%in% fixed.cov){
      formulanull=paste(response,"~-1+",fixed.cov,sep="")
    }
  }
  
  
  if(length(fixed.cov)>1){
    if("Intercept"%in% fixed.cov){
      formulanull <- paste(response,"~",paste(fixed.cov[-which(fixed.cov=="Intercept")],collapse="+"),sep="")
    }
    if(!"Intercept"%in% fixed.cov){
      formulanull <- paste(response,"~",paste(fixed.cov,collapse="+"),sep="")
      formulanull=paste(formulanull,"-1",sep="")
    }
  }
  
  lmnull<- lm(formula=formulanull, data=data, y=TRUE, x=TRUE)
  
  #check if null model is contained in the full one:
  namesnull<- dimnames(lmnull$x)[[2]]
 
  
  
  #remove the brackets in "(Intercept)" if present. 
  if(namesnull[1]=="(Intercept)"){
    namesnull[1]<- "Intercept" #namesx contains the name of variables including the intercept
  }
  
  #warning for the use of factors
  if(length(fixed.cov)<length(namesnull)){
    if(!"Intercept"%in%namesnull){
      if("Intercept"%in%namesx){
        stop("When using a factor, Intercept should be included either in both or in non of the models")
      }
    }
      if("Intercept"%in%namesnull){
        if(!"Intercept"%in%namesx){
          stop("When using a factor, Intercept should be included either in both or in non of the models")
        }
    }
    
  }
  
  
  #which variables are wrong
  ausent<-NULL
  j<-0
  for(i in 1:length(namesnull)){
    if(!namesnull[i]%in%namesx){
      ausent<-c(ausent,namesnull[i])
      j<-j+1
    }
  }

  if(j>0){
    stop(paste("object '",ausent,"' not found in the full model\n",sep=""))
  }


 

#There is any variable to select from?
  if(length(namesx)==length(namesnull)){
    stop("The number of fixed covariates is equal to the number of covariates in the full model. No model selection can be done\n")
  }


  #position for fixed variables in the full model 
  fixed.pos<-which(namesx%in%namesnull)

 


  n<- dim(data)[1]

  #the response variable for the C code
  Y<- lmnull$residuals 

  #Design matrix of the null model 
  X0<- lmnull$x
  P0<- X0%*%(solve(t(X0)%*%X0))%*%t(X0)#Intentar mejorar aprovechando lmnull
  knull <- dim(X0)[2]

  #matrix containing the covariates from which we want to select
  X1 <- lmfull$x[,-fixed.pos]
if(dim(X1)[1]<n){
  stop("NA values found for some of the competing variables")
}

  #Design matrix for the C-code
  X<- (diag(n)-P0)%*%X1 #equivalent to X<- (I-P0)X 
  namesx<- dimnames(X)[[2]]
  if(namesx[1]=="(Intercept)"){
    namesx[1]<- "Intercept" #namesx contains the name of variables including the intercept
  }
  
  p<- dim(X)[2]#Number of covariates to select from

  #check if the number of models to save is correct
  if(n.keep>2^(p)){
    warning(paste("The number of models to keep (",n.keep, ") is larger than the total number of models (",2^(p),") and it has been set to ", 2^(p) ,sep=""))
    n.keep <- 2^p
  }

# if(length(fixed.cov)<length(namesnull)){
#   warning("Some of the included covariates are factors. One dummy variable have been included for each level taking the first one (alphabethical order) as the base one.")
# }
}

#If no fixed covariates considered
if(is.null(fixed.cov)){
  #Check that all the fixed covariables are included in the full model 
  lmfull = lm(formula, data, y=TRUE, x=TRUE)
  X.full<- lmfull$x
  namesx<- dimnames(X.full)[[2]]
   #remove the brackets in "(Intercept)" if present. 
  if(namesx[1]=="(Intercept)"){
    namesx[1]<- "Intercept" #namesx contains the name of variables including the intercept
  }
  
  
  X <- lmfull$x
  knull <- 0
  Y <- lmfull$y
  p<- dim(X)[2]
  n<- dim(X)[1]
  #check if the number of models to save is correct
  if(n.keep>2^(p)){
    warning(paste("The number of models to keep (",n.keep, ") is larger than the total number of models (",2^(p),") and it has been set to ", 2^(p) ,sep=""))
    n.keep <- 2^p
  }
}

#write the data files in the working directory
write(Y, ncolumns=1, file=paste(wd,"/Dependent.txt",sep=""))
write(t(X), ncolumns=p, file=paste(wd,"/Design.txt",sep=""))

#prior for betas:
pfb<- substr(tolower(prior.betas),1,1)
#check if the selected option exists
if (pfb!="g" && pfb!="r" && pfb!="z" && pfb!="l" && pfb!="f") stop("I am very sorry: prior for betas no valid\n")

#prior for model space:
pfms<- substr(tolower(prior.models),1,1)
if (pfms!="c" && pfms!="s" && pfms!="u") stop("I am very sorry: prior for model space not valid\n")
if (pfms=="u" && is.null(priorprobs)){stop("A valid vector of prior probabilities must be provided\n")}
if (pfms=="u" && length(priorprobs)!=(p+1)){stop("Vector of prior probabilities with incorrect length\n")}
if (pfms=="u" && sum(priorprobs<0)>0){stop("Prior probabilities must be positive\n")}
if (pfms=="u" && priorprobs[1]==0){stop("Vector of prior probabilities not valid: All the theory here implemented works with the implicit assumption that the null model could be the true model\n")}	
else{
  #The zero here added is for C compatibility
  write(priorprobs, ncolumns=1, file=paste(wd,"/priorprobs.txt",sep=""))
}

#Note: priorprobs.txt is a file that is needed only by the "User" routine. Nevertheless, in order
#to mantain a common unified version the source files of other routines also reads this file
#although they do not use. Because of this we create this file anyway.
if (pfms=="c" | pfms=="s"){
  priorprobs<- rep(0,p+1)
  write(priorprobs, ncolumns=1, file=paste(wd,"/priorprobs.txt",sep=""))	
}

method<- paste(pfb,pfms,sep="")

#Info:
cat("Info. . . .\n")
cat("Most complex model has",p+knull,"covariates\n")
if(!is.null(fixed.cov)){
  if(knull>1){
    cat("From those",knull,"are fixed and we should select from the remaining",p,"\n")
  }
  if(knull==1){
    cat("From those",knull,"is fixed and we should select from the remaining",p,"\n")
  }
  cat(paste(paste(namesx,collapse=", ",sep=""),"\n",sep=""))
}
cat("The problem has a total of", 2^(p), "competing models\n")
cat("Of these, the ", n.keep, "most probable (a posteriori) are kept\n")

#check if the number of covariates is too big. 
if (p>30){stop("Number of covariates too big. . . consider using GibbsBvs\n")}

#The previous test (for time)
estim.time<- 0

if (time.test && p>=18){
  cat("Time test. . . .\n")
aux <- system.time(  result<- switch(method,
                  "gc"=.C("gConst", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "gs"=.C("gSB", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "gu"=.C("gUser", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),	
                  "rc"=.C("RobustConst", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "rs"=.C("RobustSB", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "ru"=.C("RobustUser", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "lc"=.C("LiangConst", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "ls"=.C("LiangSB", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "lu"=.C("LiangUser", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "zc"=.C("ZSConst", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "zs"=.C("ZSSB", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
                  "zu"=.C("ZSUser", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
                          as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
		          "fc"=.C("flsConst", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
		                  as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
		          "fs"=.C("flsSB", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
		                  as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull)),
		          "fu"=.C("flsUser", as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
		                  as.integer(2^(p-1)-1999), as.integer(2^(p-1)+2000), as.character(wd), as.double(estim.time), as.integer(knull))
						  
  ))
  
  estim.time<- result[[8]]*2^(p)/(60*4000) 
  cat("The problem would take ", estim.time, "minutes (approx.) to run\n")
  ANSWER <- readline("Do you want to continue?(y/n) then press enter.\n")
  while (substr(ANSWER, 1, 1) != "n" & substr(ANSWER, 1, 1) !="y"){
    ANSWER <- readline("")
  }
  
  if (substr(ANSWER, 1, 1) == "n")
  {
    return(NULL)
  }
  
}
#if the answer is yes work on the problem
cat("Working on the problem...please wait.\n")

#if the answer is yes work on the problem
result<- switch(method,
                "gc"=.C("gConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
				        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "gs"=.C("gSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "gu"=.C("gUser", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),	
                "rc"=.C("RobustConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "rs"=.C("RobustSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "ru"=.C("RobustUser", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),	
                "lc"=.C("LiangConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "ls"=.C("LiangSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "lu"=.C("LiangUser", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),	
                "zc"=.C("ZSConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "zs"=.C("ZSSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
                "zu"=.C("ZSUser", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
                        as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
		        "fc"=.C("flsConst", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
		                as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
		        "fs"=.C("flsSB", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
		                as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull)),
		        "fu"=.C("flsUser", as.character(""), as.integer(n), as.integer(p), as.integer(n.keep), 
		                as.integer(1), as.integer(2^(p)-1), as.character(wd),  as.double(estim.time), as.integer(knull))
						
)

time <- result[[8]]

#a function to transform the number of the model into a binary number.
integer.base.b_C<-function(x, k){
  #x is the number we want to express in binary
  #k is the number positions we need
  if(x==0)
    return(rep(0,k))
  else{
    ndigits <- (floor(logb(x, base=2))+1)
    res<- rep(0, ndigits)
    for(i in 1:ndigits){#i <- 1
      res[i] <- (x %% 2)
      x <- (x %/% 2)
    }
    return(c(res,rep(0,k-ndigits)))}
}

#read the files given by C
models <- as.vector(t(read.table(paste(wd,"/MostProbModels",sep=""),colClasses="numeric")))
prob <- as.vector(t(read.table(paste(wd,"/PostProb",sep=""),colClasses="numeric")))
incl <- as.vector(t(read.table(paste(wd,"/InclusionProb",sep=""),colClasses="numeric")))
joint <- as.matrix(read.table(paste(wd,"/JointInclusionProb",sep=""),colClasses="numeric"))
dimen <- as.vector(t(read.table(paste(wd,"/ProbDimension",sep=""),colClasses="numeric")))
betahat<- as.vector(t(read.table(paste(wd,"/betahat",sep=""),colClasses="numeric")))


#data.frame with Most probable models
mod.mat <- as.data.frame(cbind(t(rep(0,(p+1)))))

names(mod.mat)<-c(namesx,"prob")
 
N<-n.keep

for(i in 1:N){
   mod.mat[i,1:p]<-integer.base.b_C(models[i],p)
  varnames.aux<-rep("",p)
   varnames.aux[mod.mat[i,1:p]==1]<-"*"
   mod.mat[i,1:p]<-varnames.aux
 }
# 
 mod.mat[,(p+1)]<-prob[]
# 
# 
inclusion <- incl#inclusion probabilities except for the intercept
#names(inclusion) <- namesx 
 
#the final result

result <- list()
# 
result$time <- time #The time it took the programm to finish
result$lmfull <- lmfull # The lm object for the full model
if(!is.null(fixed.cov)){
  result$lmnull <- lmnull # The lm object for the null model
}
 
result$variables <- namesx #The name of the competing variables
result$n <- n #number of observations
result$p <- p #number of competing variables
result$k <- knull
result$HPMbin <- integer.base.b_C(models[1],(p)) #The binary code for the HPM model
names(result$HPMbin) <- namesx
result$modelsprob <- mod.mat #A table with the n.keep most probable models ands its probability
result$inclprob <- inclusion #inclusion probability for each variable 
names(result$inclprob) <- namesx
 
result$jointinclprob <- data.frame(joint[1:p,1:p],row.names=namesx)#data.frame for the joint inclusion probabilities
names(result$jointinclprob) <- namesx 
# 
result$postprobdim <- dimen #vector with the dimension probabilities.
names(result$postprobdim) <- (0:p)+knull #dimension of the true model
# 
#result$betahat <- betahat
#rownames(result$betahat)<-namesx
#names(result$betahat) <- "BetaHat"
result$call <- match.call()
result$method <- "full"
class(result)<- "Bvs"
result 
 
 

}



