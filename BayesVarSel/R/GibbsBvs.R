GibbsBvs <-
function(formula,fixed.cov=c("Intercept"), data,  prior.betas="Robust", prior.models="Constant",  n.iter=1000, init.model="Full", 
        n.burnin=50, n.thin=1, time.test=TRUE, priorprobs=NULL, seed=runif(1,0,16091956)){


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
  }
  
  #write the data files in the working directory
  write(Y, ncolumns=1, file=paste(wd,"/Dependent.txt",sep=""))
  write(t(X), ncolumns=p, file=paste(wd,"/Design.txt",sep=""))
  
  
  
 
#The initial model:
if (is.character(init.model)==TRUE){
	im<- substr(tolower(init.model),1,1)
	if (im!="n" && im!="f" && im!="r"){stop("Initial model not valid\n")}
	if (im=="n"){init.model<- rep(0,p)}
	if (im=="f"){init.model<- rep(1,p)}
	if (im=="r"){init.model<- rbinom(n=p,size=1,prob=.5)}
}
else{
	init.model<- as.numeric(init.model>0)
	if (length(init.model)!=p){stop("Initial model with incorrect length\n")}	
}

write(init.model, ncolumns=1, file=paste(wd,"/initialmodel.txt",sep=""))

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
iter<-n.iter
cat("Of these,", n.burnin+n.iter, "are sampled with replacement\n")

cat("Then,", floor(iter/n.thin), "are kept and used to construct the summaries\n")


#prior for betas:
pfb<- substr(tolower(prior.betas),1,1)
if (pfb!="g" && pfb!="r" && pfb!="z" && pfb!="l" && pfb!="f") stop("I am very sorry: prior for betas not supported\n")
#prior for model space:
pfms<- substr(tolower(prior.models),1,1)
if (pfms!="c" && pfms!="s" && pfms!="u") stop("I am very sorry: prior for model space not valid\n")
	if (pfms=="u" && is.null(priorprobs)){stop("A valid vector of prior probabilities must be provided\n")}
	if (pfms=="u" && length(priorprobs)!=(p+1)){stop("Vector of prior probabilities with incorrect length\n")}
	if (pfms=="u" && sum(priorprobs<0)>0){stop("Prior probabilities must be positive\n")}
	if (pfms=="u" && priorprobs[1]==0){stop("Vector of prior probabilities not valid: All the theory here implemented works with the implicit assumption that the null model could be the true model\n")}
	if (pfms=="u" && priorprobs[sum(init.model)+1]==0){stop("The initial model has zero prior probability\n")}	
	if (pfms=="u"){
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

#The previous test (for time)
estim.time<- 0
if(p<=20){
  warning("The number of variables is small enough to visit every model. Consider Bvs (or pBvs for its parallel version).\n")
}
if (time.test&&p>20){
	cat("Time test. . . .\n")
  

result<- switch(method,
	"gc"=.C("GibbsgConst", as.character(""), as.integer(n), as.integer(p), 
as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"gs"=.C("GibbsgSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"gu"=.C("GibbsgUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"rc"=.C("GibbsRobustConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"rs"=.C("GibbsRobustSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"ru"=.C("GibbsRobustUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),	
	"lc"=.C("GibbsLiangConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"ls"=.C("GibbsLiangSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"lu"=.C("GibbsLiangUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),	
	"zc"=.C("GibbsZSConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"zs"=.C("GibbsZSSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"zu"=.C("GibbsZSUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"fc"=.C("GibbsflsConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"fs"=.C("GibbsflsSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed)),
	"fu"=.C("GibbsflsUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(49),as.character(wd),as.integer(1), as.double(estim.time), as.integer(knull), as.integer(1), as.integer(seed))		
		)
	 
	estim.time<- result[[7]]*(n.burnin+n.iter)/(60*50) 
	cat("The problem would take ", ceiling(estim.time), "minutes (approx.) to run\n")
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

#Call the corresponding function:
result<- switch(method,
	"gc"=.C("GibbsgConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"gs"=.C("GibbsgSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"gu"=.C("GibbsgUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),		
	"rc"=.C("GibbsRobustConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"rs"=.C("GibbsRobustSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"ru"=.C("GibbsRobustUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),		
	"lc"=.C("GibbsLiangConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"ls"=.C("GibbsLiangSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"lu"=.C("GibbsLiangUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),		
	"zc"=.C("GibbsZSConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"zs"=.C("GibbsZSSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"zu"=.C("GibbsZSUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"fc"=.C("GibbsflsConst", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"fs"=.C("GibbsflsSB", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed)),
	"fu"=.C("GibbsflsUser", as.character(""), as.integer(n), as.integer(p), 
		as.integer(floor(n.iter/n.thin)),as.character(wd),as.integer(n.burnin), as.double(estim.time), as.integer(knull), as.integer(n.thin), as.integer(seed))
		)
	   
  time <- result[[7]]
  

#read the files given by C
models <- as.vector(t(read.table(paste(wd,"/MostProbModels",sep=""),colClasses="numeric")))
incl <- as.vector(t(read.table(paste(wd,"/InclusionProb",sep=""),colClasses="numeric")))
joint <- as.matrix(read.table(paste(wd,"/JointInclusionProb",sep=""),colClasses="numeric"))
dimen <- as.vector(t(read.table(paste(wd,"/ProbDimension",sep=""),colClasses="numeric")))
betahat<- as.vector(t(read.table(paste(wd,"/betahat",sep=""),colClasses="numeric")))
allmodels<- as.matrix(read.table(paste(wd,"/AllModels",sep=""),colClasses="numeric"))
allBF<- as.vector(t(read.table(paste(wd,"/AllBF",sep=""),colClasses="numeric")))


#Log(BF) for every model 
modelslBF<- cbind(allmodels, log(allBF))
colnames(modelslBF)<- c(namesx, "logBFi0")


#Highest probability model
mod.mat <- as.data.frame(t(models))


inclusion <- incl
names(inclusion) <- namesx 
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
result$k <- knull#number of fixed covariates
result$HPMbin <- models#The binary code for the HPM model
names(result$HPMbin) <- namesx
#result$modelsprob <- mod.mat
result$modelslogBF <-modelslBF#The binary code for all the visited models (after n.thin is applied) and the correspondent log(BF) 
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
result$method <- "gibbs"
class(result)<- "Bvs"
result 

}