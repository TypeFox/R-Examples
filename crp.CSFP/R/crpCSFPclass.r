setClass(Class="crp.CSFP",                                                                               #define the names and types of attributes of crp.CSFP
	   representation=representation(
							path.in = "character",
							path.out="character",
							port.name="character",
							rating.scale.name="character",
							sec.var.name="character",
							sec.var.est="numeric",
							loss.unit="numeric",
							Niter.max="numeric",
              Niter.max.global="numeric",
							alpha="numeric",
							PLOT.PDF="logical",
							calc.rc="logical",
							rating="numeric",
							rating.PD="numeric",
							rating.SD="numeric",
							NS="numeric",
							NC="numeric",
							sec.var="numeric",
							CP.NR="numeric",
							CP.rating="numeric",
							NEX="numeric",
							LGD="numeric",
							PL="numeric",
							PD="numeric",
							EL="numeric",
							EL.crp="numeric",
							nu="numeric",
							PL.crp="numeric",
							PD.crp="numeric",
							M="numeric",
							mu.k="numeric",
							loss.k="numeric",
              sigma_k="numeric",
							sigma_sqr_div="numeric",
							sigma_sqr_syst="numeric",
							SD="numeric",
							SD.crp="numeric",
							W="matrix",
							alpha.max="numeric",
							a="numeric",
							PDF="numeric",
							CDF="numeric",
							B="matrix",
							loss="numeric",
              PLOT.scale="numeric",
							PLOT.range.x="numeric",
							PLOT.range.y="numeric",
							VaR="numeric",
              EC="numeric",
							ES="numeric",
							VaR.cont="numeric",
							ES.cont="numeric",
							ES.tau.cont="numeric",
							SD.cont="numeric",
							read.OK="logical",
							plausi.OK="logical",
              rc.OK="logical",
              save.memory="logical",
              VaR.pos="numeric",
              alpha.crp="numeric",
              name="character",
              changes.read="logical",
              changes.plausi="logical",
              changes.calc.portfolio.statistics="logical",
              changes.loss="logical",
              changes.measure="logical",
              changes.plot="logical",
              changes.rc.vares="logical",
              changes.rc.sd="logical",
              changes.export="logical",
              file.format="character",
              input="list",
              export.to.file="logical"),
	   prototype=prototype( path.in="",                                                            #define default values
					path.out="",
					port.name="portfolio.csv",
					rating.scale.name="rating_pd.csv",
					sec.var.name="pd_sector_var.csv",
					sec.var.est=5,
					loss.unit=1e6,
					Niter.max=0.9999,
          Niter.max.global=1e5,                
					alpha=c(0.999),
					PLOT.PDF=TRUE,
					calc.rc=FALSE,
					rating=0,
					rating.PD=0,
					rating.SD=0,
					NS=0,
					NC=0,
					sec.var=0,
					CP.NR=0,
					CP.rating=0,
					NEX=0,
					LGD=0,
					PL=0,
					PD=0,
					EL=0,
					EL.crp=0,
					nu=0,
					PL.crp=0,
					PD.crp=0,
					M=0,
					mu.k=0,
					loss.k=0,
					sigma_k=0,
					sigma_sqr_div=0,
					sigma_sqr_syst=0,
					SD=0,
					SD.crp=0,
					W=matrix(),
					alpha.max=0,
					a=0,
					PDF=0,
					CDF=0,
					B=matrix(),
					loss=0,
          PLOT.scale=1e6,                
					PLOT.range.x=c(0,0),
					PLOT.range.y=c(0,0),
					VaR=0,
          EC=0,
					ES=0,
					VaR.cont=0,
					ES.cont=0,
					ES.tau.cont=0,
					SD.cont=0,
          read.OK=FALSE,
          plausi.OK=FALSE,
          rc.OK=FALSE,
          save.memory=FALSE,
          VaR.pos=0,
          alpha.crp=0,
          name="MyModel",
          changes.read=FALSE,
          changes.plausi=FALSE,
          changes.calc.portfolio.statistics=FALSE,
          changes.loss=FALSE,
          changes.measure=FALSE,
          changes.plot=FALSE,
          changes.rc.vares=FALSE,
          changes.rc.sd=FALSE,
          changes.export=FALSE,
          file.format="csv",
          input=list(),
          export.to.file=FALSE)
)

setGeneric("read",function(this) standardGeneric("read"))
setMethod("read",c("crp.CSFP"),function(this) {
#
#       <read>      Import portfolio and risk data
#
#       Last Modified:  24/06/2013
#       Author:         Matthias Fischer & Kevin Jakob        
#

   ERROR=0
# --------------------------------------------------------------------------------------
#               Importing rating, probability of default and standard deviation
# ---------------------------------------------------------------------------------------
  
  cat("Importing risk information (Rating, PD, SD)....\n")
  if(nrow(this@input$rating.scale)>0)
    RISK.MATRIX=this@input$rating.scale
  else if(file.exists(paste(this@path.in,this@rating.scale.name,sep=""))==TRUE){                                  # check if the file exists
    if(this@file.format=="csv")
      RISK.MATRIX=read.csv(paste(this@path.in,this@rating.scale.name,sep=""),header=TRUE)
    else if(this@file.format=="csv2")
      RISK.MATRIX=read.csv2(paste(this@path.in,this@rating.scale.name,sep=""),header=TRUE)
  }
  else{
    cat(paste("ERROR: The ",this@rating.scale.name," file is missing!!!\n",sep=""))
    ERROR=ERROR+1
  }
  if(ERROR==0){
    if(is.null(RISK.MATRIX$RATING)){
      cat("The RATING-column in ",this@rating.scale.name," is missing.")
      ERROR=ERROR+1
    }
    else
      this@rating=RISK.MATRIX$RATING                                                                 # assigne data to attributes
    if(is.null(RISK.MATRIX$PD)){
      cat("The PD-column in ",this@rating.scale.name," is missing.\n")
      ERROR=ERROR+1
    }
    else
      this@rating.PD=RISK.MATRIX$PD
    if(this@sec.var.est<5){
      if(is.null(RISK.MATRIX$SD)){
        cat("The SD-column in ",this@rating.scale.name," is missing.\n")
        ERROR=ERROR+1
      }
      else
        this@rating.SD=RISK.MATRIX$SD
    }
    remove(RISK.MATRIX)
    this@input$rating.scale=data.frame()
  }
    
   


# --------------------------------------------------------------------------------------
#               Importing portfolio data
# ---------------------------------------------------------------------------------------

  cat("Importing portfolio data....\n")
  if(nrow(this@input$portfolio)>0)
    PORTFOLIO=this@input$portfolio
  else if(file.exists(paste(this@path.in,this@port.name,sep=""))==TRUE){                                  # check if file exists
    if(this@file.format=="csv")
   	  PORTFOLIO=read.csv(paste(this@path.in,this@port.name,sep=""),header=TRUE)
    else if(this@file.format=="csv2")
      PORTFOLIO=read.csv2(paste(this@path.in,this@port.name,sep=""),header=TRUE)
  }
  else{
   cat(paste("ERROR: The ",this@port.name," file is missing!!!\n"),sep="")
   ERROR=ERROR+1
  }
  if(ERROR==0){
    this@NS=length(PORTFOLIO)-6                                                                      # assigne data to attributes

    if(is.null(PORTFOLIO$CPnumber)){
      cat("The CPnumber-column in ",this@port.name," is missing.\n")
      ERROR=ERROR+1
    }
    else{
      this@CP.NR=PORTFOLIO$CPnumber
      this@NC=length(PORTFOLIO$CPnumber)
    }
    if(is.null(PORTFOLIO$rating)){
      cat("The rating-column in ",this@port.name," is missing.\n")
      ERROR=ERROR+1
    }
    else
   	  this@CP.rating=PORTFOLIO$rating
   	if(is.null(PORTFOLIO$exposure)){
      cat("The exposure-column in ",this@port.name," is missing.\n")
      ERROR=ERROR+1
    }
    else
   	  this@NEX=PORTFOLIO$exposure
    if(is.null(PORTFOLIO$lgd)){
      cat("The lgd-column in ",this@port.name," is missing.\n")
      ERROR=ERROR+1
    }
    else
   	  this@LGD=PORTFOLIO$lgd
    TEMP=(do.call(cbind,PORTFOLIO))[,-(1:6)]                                                         # extract weight matrix
    if(is.null(TEMP)){
      cat("The sector-weights in ",this@port.name," are missing.\n")
      ERROR=ERROR+1
    }
    else{
      remove(PORTFOLIO)
      this@input$portfolio=data.frame()
      if(this@NS==1)
        TEMP=matrix(ncol=1,TEMP)                                                                       # convert to matrix
      this@W=cbind(1-apply(TEMP,1,sum),TEMP)                                                           # attache idiosyncratic weights as first column
      remove(TEMP)
    }   
  }
   
   
# --------------------------------------------------------------------------------------
#               Importing sector information
# ---------------------------------------------------------------------------------------

   if(this@sec.var.est==5){                                                                          # sector variances are just needed if sec.var.est=5, otherwise they are calculated out of rating.SD
   	 cat("Importing PD sector variances....\n")   
     if(nrow(this@input$sec.var)>0)
       sec.var=this@input$sec.var
     else if(file.exists(paste(this@path.in,this@sec.var.name,sep=""))==TRUE){                             # check if file exists
        if(this@file.format=="csv")
   		    sec.var=read.csv(paste(this@path.in,this@sec.var.name,sep=""),header=TRUE)
        else if(this@file.format=="csv2")
          sec.var=read.csv2(paste(this@path.in,this@sec.var.name,sep=""),header=TRUE)
     }
   	 else{
   	   cat(paste("ERROR: The ",this@sec.var.name," file is missing!!!\n",sep=""))
   	   ERROR=ERROR+1
   	 }
    if(ERROR==0){
      if(is.null(sec.var$Var)){
        cat("The Var-column in ",this@sec.var.name," is missing.\n")
        ERROR=ERROR+1
      }
      else if(length(sec.var$Var)<this@NS){
        cat("ERROR: Number of sector variances specified in ",this@sec.var.name," (",length(sec.var$Var),") is smaller then the number of sectors (",this@NS,").\n",sep="")
        ERROR=ERROR+1
      }
      else{
   		  this@sec.var=sec.var$Var
   		  remove(sec.var)
        this@input$sec.var=data.frame()
      }
    }
   }
   if(this@sec.var.est<5)
	   this@sec.var=0

   if(ERROR>0)                                                                                       # set state read.ok
     this@read.OK=FALSE
   else
     this@read.OK=TRUE

   this@changes.read=FALSE
   gc()
   return(this)
}
)

setGeneric("plausi",function(this) standardGeneric("plausi"))
setMethod(f="plausi",signature=c("crp.CSFP"),definition=function(this){
#
#       <plausi>    Verify that the input data makes sence
#
#       Last Modified:  24/06/2013
#       Author:         Dr. Matthias Fischer, Stefan Kolb & Kevin Jakob
#

  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")

  ERROR=0
  if(length(this@rating.PD)>1){                                                                     # check that higher rating has higher probability of default
    for(i in 1:(length(this@rating.PD)-1)){
      if(this@rating.PD[i+1]-this@rating.PD[i]<0){
        cat("ERROR: Higher Rating cannot have smaller Probability of Default!!!\n")
        ERROR=ERROR+1
      }
    }
  }

  if(sum(this@NEX>=0)!=length(this@NEX)){                                                            # check if exposure is not negative
    cat("ERROR: Exposure has to be not negative. \n")
    ERROR=ERROR+1
  }
  if(any(this@LGD<0) || any(this@LGD>1))
    stop("LGDs are supposed to be between 0 and 1!\n")
  
  if(this@NS>1){                                                                                     # check if sector weights are plausible
    if(sum(rowSums(this@W[,-1])<=1)!=nrow(this@W[,-1])|sum(this@W[,-1]>=0)!=nrow(this@W[,-1])*ncol(this@W[,-1])){
      cat("ERROR: Sum of weights and have to be between 0 and 1 for each counterparty!\n")
      ERROR=ERROR+1
     }
  }
  else if(this@NS==1){
	  if(any(this@W>=0 & this@W<=1)==FALSE){
	    cat("ERROR: Sum of weights have to be between 0 and 1 for each counterparty!\n")
	    ERROR=ERROR+1
	  }
  }

  if(this@sec.var.est<5){
    if(ERROR>0)
      this@plausi.OK=FALSE
    else
      this@plausi.OK=TRUE
    gc()
    this@changes.plausi=FALSE
	  return(this)                                                                                     # return if the variance is calculated and not from an input file
  }                                                       

  if(sum(this@sec.var>=0)!=length(this@sec.var)){                                                    # check if sector variances are postive
    cat("ERROR: Sector Variance has to be positive\n")
    ERROR=ERROR+1
  }
  if(ERROR>0)
    this@plausi.OK=FALSE
  else
    this@plausi.OK=TRUE
  gc()
  this@changes.plausi=FALSE
  return(this)                                                                        
}
)

setGeneric("calc.portfolio.statistics",function(this) standardGeneric("calc.portfolio.statistics"))
setMethod(f="calc.portfolio.statistics",signature=c("crp.CSFP"),definition=function(this){
#
#       <calc.portfolio.statistics>              calculate key numbers for portfolio
#
#       Last Modified:  24/06/2013
#       Author:         Matthias Fischer & Kevin Jakob        
#

  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(this@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
 
  i=which(this@rating==0)
  if(length(i)>0){
    this@rating=this@rating[-i]
    this@rating.PD=this@rating.PD[-i]
    this@rating.SD=this@rating.SD[-i]
  }
  
  emptyCP=0
  l=0
  for(i in 1:this@NC){
    if(this@CP.rating[i]==0){
      l=l+1
      emptyCP[l]=i
    }
    else if(this@NEX[i]==0 || this@rating.PD[this@CP.rating[i]]==0 || this@LGD[i]==0){
      l=l+1
      emptyCP[l]=i
    }
  }
  if(l>0){
    this@CP.NR=this@CP.NR[-emptyCP]
    this@CP.rating=this@CP.rating[-emptyCP]
    this@NEX=this@NEX[-emptyCP]
    this@LGD=this@LGD[-emptyCP]
    this@W=this@W[-emptyCP,]
    this@NC=this@NC-l
  }
  
  cat(this@NS," sectors ...",sep="")
  cat(this@NC," counterparties (",l," removed) ... OK\n",sep="")

	this@PL=this@NEX*this@LGD                                                                         # potential loss
  this@PD=this@rating.PD[this@CP.rating]                                                            # probability of default
  this@nu=crp.round(this@PL/this@loss.unit)                                                         # integer muliples of loss.unit representing poZential loss
	this@nu=ifelse(this@nu==0,1,this@nu)                                                               # should not be 0
	this@PL.crp=this@nu*this@loss.unit                                                                     # discretized potential loss
  this@PD.crp=this@PD*(this@PL/this@PL.crp)                                                                # transformed PD
  
 	cat("Calculating sector information...\n")

  for(k in 1:this@NS){
    TEMP=this@W[,1+k]
    SELECTION=which(TEMP>0)
    if(length(SELECTION)>0 && this@sec.var.est!=5){
      this@mu.k[k]=sum(TEMP[SELECTION]*this@PD.crp[SELECTION])                                             # formula (39), CSFP 1997, p. 43
      this@loss.k[k]=sum(TEMP[SELECTION]*this@PD.crp[SELECTION]*this@PL.crp[SELECTION])                           # required for risk contributions
      if(this@sec.var.est==1)
        this@sigma_k[k]=sum(TEMP[SELECTION]*this@rating.SD[this@CP.rating[SELECTION]])                 # (2.21) p. 17 in Gundlach/Lehrbass (2003)
      else if(this@sec.var.est==2)
        this@sigma_k[k]=sum(TEMP[SELECTION]*this@rating.SD[this@CP.rating[SELECTION]])/this@mu.k[k]    # (2.21) p. 17 in Gundlach/Lehrbass (2003)
      else if(this@sec.var.est==3)
        this@sigma_k[k]=sum(sqrt(TEMP[SELECTION])*this@rating.SD[this@CP.rating[SELECTION]])           # (++++) p. 18 in Gundlach/Lehrbass (2003)
      else if(this@sec.var.est==4)
        this@sigma_k[k]=sum(sqrt(TEMP[SELECTION])*this@rating.SD[this@CP.rating[SELECTION]])/this@mu.k[k]     # divided by mu.k
    }
    else if(length(SELECTION)==0 && this@sec.var.est!=5){
      this@loss.k[k]=NA
      this@mu.k[k]=NA
      this@sigma_k=NA
    }
    if(this@sec.var.est==5){
      if(length(SELECTION)>0){
        this@mu.k[k]=sum(TEMP[SELECTION]*this@PD.crp[SELECTION])                                             # formula (39), CSFP 1997, p. 43
        this@loss.k[k]=sum(TEMP[SELECTION]*this@PD.crp[SELECTION]*this@PL.crp[SELECTION])
      }  
      this@sigma_k[k]=sqrt(as.numeric(this@sec.var[k]))                                              # from an external file
    }
    else
        cat("Wrong specification for the sector variance!\n")

  }
  
  remove(TEMP)
  remove(SELECTION)

  this@sigma_sqr_syst=sum(this@sigma_k^2*this@loss.k*this@loss.k)  			                                   # portfolio standard diviation
  
  cat("Sector information completed...\n")
  this@sigma_sqr_div=sum(this@PL.crp^2*this@PD.crp)                                                             # diversifible risk
  cat("Diversifible risk: ",fo(this@sigma_sqr_div),"  Systematic risk: ",fo(this@sigma_sqr_syst),"\n")
    
  this@EL=sum(this@PL*this@PD)                                                                     # expected loss
  this@SD=sqrt(this@sigma_sqr_div+this@sigma_sqr_syst)                                                     # standard deviation
  
  cat("Calculating portfolio statistics....\n")
  cat("Portfolio net exposure:",fo(sum(as.numeric(this@NEX))),"\n")
  cat("Portfolio potential loss:",fo(sum(as.numeric(this@PL))),"\n")
  cat("Portfolio expected loss:",fo(this@EL),"\n")
  cat("Portfolio standad deviation:",fo(this@SD),"\n")
  cat("Max. exposure band per CP: ",max(this@nu),"\n",sep="")
  cat("Discretization completed...\n")

  if(this@save.memory){                                                                              # PD0 and PL0 are not needed anymore
    this@PL=0
    this@PD=0
  }
  this@changes.calc.portfolio.statistics=FALSE
  gc()
	return(this)
}
)



setGeneric("loss.dist",function(this) standardGeneric("loss.dist"))
setMethod(f="loss.dist",signature=c("crp.CSFP"),definition=function(this){
#
#       <loss> Calculation of the portfolio loss distribution 
#                               modified Nested evaluation according to Haaf, Rei?, Schoenmakers (2003)
#                               Chapter 5 in Gundlach/Lehrbass (2003), p.74-75
#       Author:         Dr. Matthias Fischer & Kevin Jakob
#       Last Modified:  24/06/2013
#
  
  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(this@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
  if(this@changes.calc.portfolio.statistics)
    cat("WARNING: Inputparameter effecting 'calc.portfolio.statistics()' have been changed, without running 'calc.portfolio.statistics()' and the following routines afterwards.\n")  
  
  if(deparse(substitute(this))!="this")                                                              # synchronize model name to the calling model 
    this@name=deparse(substitute(this))                                                              # if it is not called from crp.CSFP() 
   
  if(this@Niter.max==0 && this@alpha.max>0 && this@alpha.max<1){			# calculate CDF until alpha.max level is reached
	  this@M=min(this@Niter.max.global,sum(this@nu))						       # loss distribution should be calculated till certain confidenc level.                                                                  			                                                               
    cat("Calculate the loss distribution till ",this@alpha.max,"-confidence level is reached.",sep="","\n")  # It has no effect if Niter.max > 1 is set manuelly.
  }
	else{
	  this@M=this@Niter.max
		this@alpha.max=1
    cat("Number of probabilities to be calculated (Niter.max): ",this@Niter.max,sep="","\n")
	} 
  cat("Loss unit: ",fo(this@loss.unit),"\n",sep="")  
   
      
  A=matrix(0,ncol=(this@M+1),nrow=this@NS)                                                           # Make all needed variables local to this 
  B=A                                                                                                # function instead of going up the environment 
  a=rep(0,this@M+1)                                                                                  # every calling of 'this' . This will increase 
  PDF=rep(0,this@M+1)                                                                                # speed but needs more memory.
  CDF=rep(0,this@M+1)
  sigma_k=this@sigma_k
  mu.k=this@mu.k
  nu=this@nu
  W=this@W
  PD=this@PD.crp

##############################################################################################################
#       1. Calculation of PDF
##############################################################################################################

        
  if(sum(W[,2:(this@NS+1)])==0){                                                                     # just idiosyncratic
    a[1]=-sum(W[,1]*PD)
    PDF[1]=exp(a[1])
    CDF[1]=PDF[1]

    for(j in 1:this@M){
      if(CDF[j]<1){
        SELECT=which(nu==j)
        sum1=0
        if(length(SELECT)>0)
  	      sum1=sum(W[SELECT,1]*PD[SELECT])
        a[j+1]=sum1
        PDF[j+1]=sum((1:j)/j*a[2:(j+1)]*PDF[j:1])
        CDF[j+1]=CDF[j]+PDF[j+1]
      }
    }
  }

  else{                                                                                              # standard case
    for(k in 1:this@NS){
	    A[k,1]=1+sigma_k[k]^2*mu.k[k]
	    B[k,1]=-log(A[k,1])
	  }
	  a[1]=-sum(W[,1]*PD)+sum(B[,1]/(sigma_k[1:this@NS]^2))
	  PDF[1]=exp(a[1])
	  CDF[1]=PDF[1]
	
	  j=0
	  while(j<this@M && CDF[j+1]<this@alpha.max){
		  j=j+1
		  if(j%%1e3==0)
  	    cat(" it=",j," CDF=",CDF[j],";")
      SELECT=which(nu==j)
		  for(k in 1:this@NS){
		    if(length(SELECT)>0)
			    A[k,j+1]=sigma_k[k]^2*sum(W[SELECT,k+1]*PD[SELECT])
		    if(j+1==2)
			    B[k,2]=A[k,2]/A[k,1]
		    else
			    B[k,j+1]=1/A[k,1]*(A[k,j+1]+(1/j*sum((1:(j-1))*B[k,2:j]*A[k,j:2])))
		  }
		
		  sum1=0
      if(length(SELECT)>0)
        sum1=sum(W[SELECT,1]*PD[SELECT])
      a[j+1]=sum1+sum(B[,j+1]/(sigma_k[1:this@NS]^2))
      PDF[j+1]=sum((1:j)/j*a[2:(j+1)]*PDF[j:1])
		  CDF[j+1]=CDF[j]+PDF[j+1]
	  }
    cat("\n")

    remove(list=c("A","SELECT","W","sigma_k","mu.k","PD","nu"))	
	  this@PDF=PDF
    remove(PDF)
	  this@CDF=CDF
    remove(CDF)
	  this@a=a
    remove(a)
	  this@B=B
	  remove(B)      
      
  }
      
  length(this@CDF)=j+1                                                                               # cut off the variables to the needed length
  length(this@a)=j+1
	length(this@PDF)=j+1
	this@B=matrix(this@B[,-((j+1):(this@M+1))],nrow=this@NS)

##############################################################################################################
#       3. Summarize
##############################################################################################################

  this@loss=(0:j)*this@loss.unit                                                                     # create losses
  this@alpha.max=this@CDF[length(this@CDF)]
  this@M=j
  cat("Reached level of confidence: ",this@alpha.max," ( iterations actually done: ",this@M," )\n",sep="")
  cat("Calculation completed...\n")
  
  if(this@save.memory && !this@calc.rc){                                                      # delete data a and B if not needed anymore
    this@B=matrix()
    this@a=0
  }
  if(this@save.memory){                                                                              # don't store loss and CDF permanently in 
    this@loss=0                                                                                      # save.memory mode
    this@CDF=0
  }

  this@changes.loss=FALSE
  gc()
  return(this)        
}
)



setGeneric("measure",function(this) standardGeneric("measure"))
setMethod(f="measure",signature=c("crp.CSFP"),definition=function(this){
#
#       <measure>   Calculate portfolio measures out of the calculated distribution.
#
#       Last Modified:  24/06/2013
#       Author:         Dr. Matthias Fischer, Stefan Kolb & Kevin Jakob
#
#

  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(this@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
  if(this@changes.calc.portfolio.statistics)
    cat("WARNING: Inputparameter effecting 'calc.portfolio.statistics()' have been changed, without running 'calc.portfolio.statistics()' and the following routines afterwards.\n")
  if(this@changes.loss)
    cat("WARNING: Inputparameter effecting 'loss()' have been changed, without running 'loss()' and the following routines afterwards.\n") 
  
  if(length(this@loss)==1 && this@save.memory)                                                       # recreate losses if not available
    this@loss=this@loss.unit*(0:this@M)
  if(length(this@CDF)==1 && this@save.memory)                                                        # recalculate CDF if not available
    this@CDF=cumsum(this@PDF)
            
  cat("Calculating portfolio statistics....\n")
        
  this@EL.crp=sum(as.numeric(this@loss*this@PDF))                                                    # crp. expected loss
  cat("CR+ portfolio expected loss:",fo(this@EL.crp),"\n")
  cat("Expected loss difference:",fo(this@EL.crp-this@EL),"\n")
  cat("CR+ Exceedance Probability of the expected loss:",1-this@CDF[floor(this@EL/this@loss.unit)],"\n")

  this@SD.crp=sqrt(sum(as.numeric(this@PDF*(this@loss-this@EL.crp)^2)))                              # crp. standard deviation
  cat("CR+ portfolio standard deviation:",fo(this@SD.crp),"\n")  
  
  L.alpha=length(this@alpha)
  this@VaR=rep(0,L.alpha)
  this@EC=rep(0,L.alpha)
  this@VaR.pos=rep(0,L.alpha)
  
  for(i in 1:L.alpha){
	  if(this@alpha.max>this@alpha[i]){					                                                       # check if alpha-quantile exists
      this@VaR.pos[i]=min(which(this@CDF>this@alpha[i]))                                                # save position of alpha-quantile and VaR
      this@VaR[i]=round(this@loss[this@VaR.pos[i]],0)
      cat("CR+ portfolio Value-at-risk(",this@alpha[i],"): ",fo(this@VaR[i]),sep="","\n")
		}
		else{
		  cat("WARNING: Value-at-risk(",this@alpha[i],") is not available, increase loss.unit or Niter.max",sep="","\n")
			this@VaR[i]=0
		}
	}
  this@EC=this@VaR-this@EL.crp
  for(i in 1:L.alpha){
    if(this@EC[i]>0)
      cat("CR+ portfolio economic capital(",this@alpha[i],"): ",fo(this@EC[i]),sep="","\n")
  }
  for(i in 1:L.alpha){
	  if(this@alpha.max>this@alpha[i]){                                                                # check if alpha-quantiles exists
      if(this@VaR.pos[i]==1)
        this@ES[i]=this@EL                                                                           # calculate ES
      else
        this@ES[i]=(this@EL-sum(this@loss[1:(this@VaR.pos[i]-1)]*this@PDF[1:(this@VaR.pos[i]-1)]))/(1-this@CDF[this@VaR.pos[i]-1])   # calculate ES 
      cat("CR+ portfolio Expected Shortfall(",this@alpha[i],"): ",fo(this@ES[i]),sep="","\n")        # "backwards" because of EL/EL.crp difference
		}
	  else{
		  cat("WARNING: Expected Shortfall (",this@alpha[i],") is not available, increase loss.unit or Niter.max",sep="","\n")
		  this@VaR.pos[i]=0
		  this@ES[i]=0
		}
	}
  this@alpha.crp=this@CDF[this@VaR.pos]                                                                 # save crp. alphas for further calculations
  ES.EL=(this@EL-sum(this@loss[1:floor(this@EL/this@loss.unit)]*this@PDF[1:floor(this@EL/this@loss.unit)])) /(1-this@CDF[floor(this@EL/this@loss.unit)])
  cat("CR+ portfolio mean expected loss exceedance: ",fo(ES.EL),sep="","\n")

  if(this@save.memory){
    this@loss=0
    this@CDF=0
  }
  this@changes.measure=FALSE 
  gc()
	return(this)
}
)

setMethod("plot","crp.CSFP",function(x,y=""){
#
#       <plot>      Plotting the calculated loss distribution together with key numbers from measure
#       Author:         Dr. Matthias Fischer & Kevin Jakob
#       Last Modified:  24/06/2013
#

  if(x@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(x@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
  if(x@changes.calc.portfolio.statistics)
    cat("WARNING: Inputparameter effecting 'calc.portfolio.statistics()' have been changed, without running 'calc.portfolio.statistics()' and the following routines afterwards.\n")
  if(x@changes.loss)
    cat("WARNING: Inputparameter effecting 'loss()' have been changed, without running 'loss()' and the following routines afterwards.\n")
  if(x@changes.measure)
    cat("WARNING: Inputparameter effecting 'measure()' have been changed, without running 'measure()' and the following routines afterwards.\n")
  
  if(deparse(substitute(x))!="x")                                                              # synchronize model name to the calling model
    x@name=deparse(substitute(x))                                                              # if it is not called from crp.CSFP()

  if(length(x@loss)==1 && x@save.memory)                                                       # recreate losses if not available
    x@loss=x@loss.unit*(0:x@M)
  if(length(x@CDF)==1 && x@save.memory)                                                        # recalculate CDF if not available
    x@CDF=cumsum(x@PDF)

  if(x@PLOT.range.x[1]<1 && x@PLOT.range.x[1]>0)                                          	   # convert plot intervall of confidence levels
	  x@PLOT.range.x[1]=x@loss[min(which(x@CDF>x@PLOT.range.x[1]))]	                       # into an intervall of absolute losses
	if(x@PLOT.range.x[2]<1 && x@PLOT.range.x[2]>0)
	  x@PLOT.range.x[2]=x@loss[min(which(x@CDF>x@PLOT.range.x[2]))]
	 
  unit=fo(x@PLOT.scale)
  PLOT.TITLE.XLAB=paste("Loss in ",unit,sep="")
  PLOT.TITLE.YLAB="Probability"
  PLOT.TITLE=x@name
  title=("CSFP-model")
	  	  
  if(all((x@PLOT.range.x==c(0,0) & x@PLOT.range.y==c(0,0))==TRUE))                             # distinguish if x-/y ranges are given or choosen by R automatically
    plot(x@loss/x@PLOT.scale,x@PDF,type="l",main=paste("Portfolio Credit Loss (",PLOT.TITLE,")\n",title,sep=""),xlab=PLOT.TITLE.XLAB,ylab=PLOT.TITLE.YLAB,cex.axis=1.3,cex.main=1.3,lwd=3,cex.lab=1.3)
	else if(any((x@PLOT.range.x!=c(0,0))==TRUE))
  	plot(x@loss/x@PLOT.scale,x@PDF,type="l",main=paste("Portfolio Credit Loss (",PLOT.TITLE,")\n",title,sep=""),xlab=PLOT.TITLE.XLAB,ylab=PLOT.TITLE.YLAB,cex.axis=1.3,cex.main=1.3,lwd=3,cex.lab=1.3,xlim=x@PLOT.range.x/x@PLOT.scale)
	else if(any((x@PLOT.range.y!=c(0,0))==TRUE))
  	plot(x@loss/x@PLOT.scale,x@PDF,type="l",main=paste("Portfolio Credit Loss (",PLOT.TITLE,")\n",title,sep=""),xlab=PLOT.TITLE.XLAB,ylab=PLOT.TITLE.YLAB,cex.axis=1.3,cex.main=1.3,lwd=3,cex.lab=1.3,ylim=x@PLOT.range.y)
	else 
  	plot(x@loss/x@PLOT.scale,x@PDF,type="l",main=paste("Portfolio Credit Loss (",PLOT.TITLE,")\n",title,sep=""),xlab=PLOT.TITLE.XLAB,ylab=PLOT.TITLE.YLAB,cex.axis=1.3,cex.main=1.3,lwd=3,cex.lab=1.3,ylim=x@PLOT.range.y,xlim=x@PLOT.range.x/x@PLOT.scale)

  nVaR=0
  nES=0
  if(x@EL.crp>0)                                                                                  # check if EL, VaR and ES exists
    lines=c(x@EL.crp)
  for(i in 1:length(x@VaR)){
    if(x@VaR[i]>0){
      lines=c(lines,x@VaR[i])
      nVaR=i
    }
  }
  for(i in 1:length(x@ES)){
    if(x@ES[i]>0){
      lines=c(lines,x@ES[i])
      nES=i
    }
  }
  
  abline(v=c(x@EL.crp,x@VaR,x@ES)/x@PLOT.scale,col=c("green",rep("blue",nVaR),rep("red",nES)))   # add lines for EL, VaR and ES
  legend("topright",c("EL.crp","VaR","ES"),col=c("green","blue","red"),lwd=1)                        # add legend
  
  if(x@save.memory){                                                                              # dont store loss and CDF permanently in
    x@loss=0                                                                                      # save.memory mode
    x@CDF=0
  }
  x@changes.plot=FALSE
  gc()
  return(x)
}
)


setGeneric("rc.vares",function(this) standardGeneric("rc.vares"))
setMethod(f="rc.vares",signature=c("crp.CSFP"),definition=function(this){
#
#       <rc.vares>    Calculation of contributions to VaR and Expected Shortfall
#
#       Last Modified:  24/06/2013
#       Author:         Dr. Matthias Fischer, Stefan Kolb & Kevin Jakob
#

  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(this@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
  if(this@changes.calc.portfolio.statistics)
    cat("WARNING: Inputparameter effecting 'calc.portfolio.statistics()' have been changed, without running 'calc.portfolio.statistics()' and the following routines afterwards.\n")
  if(this@changes.loss)
    cat("WARNING: Inputparameter effecting 'loss()' have been changed, without running 'loss()' and the following routines afterwards.\n")
  if(this@changes.measure)
    cat("WARNING: Inputparameter effecting 'measure()' have been changed, without running 'measure()' and the following routines afterwards.\n") 
  
  l=length(this@VaR.pos)
	
  if(!(this@VaR.pos[l]>0 && this@ES[l]>0) || length(this@CDF)<=this@VaR.pos){                                                             # check if VaR and ES exists
	  cat("VaR and ES contributions are not available. \n")
		this@rc.OK=FALSE
		return(this)
	}
	if(length(this@a)==0 || all(dim(this@B)==1)){                                                      # check if required data exists
    cat("Necessary data for calculation of risk contributions (VaR/ES) are not availabel, probably you`re running in save.memory mode. \n Change save.memory to FALSE or calc.rc to TRUE, run loss and start calculation of risk contributions (VaR/ES) again.\n")
    return(this)
  }
  
  if(l>1)
	  alpha=this@alpha[l]                                                                              # calculate RC just for last alpha
  cat("Calculating VaR and ES contributions....\n")
         
  
  Z=length(this@PDF)                                                                                 # initialization
  VaR.K=numeric(this@NC)
  VaR.K=this@VaR.pos[l]-this@PL.crp/this@loss.unit
  alpha.M=matrix(0,nrow=this@NS+1,ncol=this@VaR.pos[l])
  PLH=matrix(0,nrow=this@NS+1,ncol=Z)
  PLH[1,]=this@PDF
  PDVAR=this@PDF[this@VaR.pos[l]]
  VaR.cont=numeric(this@NC)                                                                          # Make all needed variables local to this
  ES.cont=numeric(this@NC)                                                                           # function instead of going up the environment
  a=this@a                                                                                           # every calling of 'this' . This will increase
  B=this@B                                                                                           # speed.
  VaR.pos=this@VaR.pos
  PL=this@PL.crp
  PD=this@PD.crp
  NS=this@NS
  W=this@W
  NC=this@NC                                                                                         # probability of the VaR
      
  if (NS==1)
    B=matrix(nrow=1,B)
  a1=a[1]+B[,1]                                                                                      # for details please see SAS documentation or p.87-90 in "Credit Risk in the Banking Industry".
  PLH[2:(NS+1),1]=exp(a1)                                                                            # calculation of the "distorted" PDF for each sector 
  for(h in 2:(NS+1)){
    alpha.M[h,]=a[1:VaR.pos[l]]+B[h-1,1:VaR.pos[l]]
    for(j in 2:VaR.pos[l])
      PLH[h,j]=sum((1:(j-1))/(j-1)*alpha.M[h,2:j]*PLH[h,(j-1):1])
  }
  CLH=PLH
  for(x in 1:(NS+1))                                                                                 # calculation of the CDF
    CLH[x,]=cumsum(PLH[x,])
  alpha=this@alpha.crp[l]-this@PDF[this@VaR.pos[l]]
  for(i in 1:NC){                                                                                    # Verify if the exposure is equal or bigger than the VaR
    if(VaR.K[i]<=0){
      VaR.cont[i]=0
      ES.cont[i]=(PL[i]*PD[i])/(1-alpha)
    }
    else if(VaR.K[i]==1){
      VaR.cont[i]=(PL[i]*PD[i]*sum(W[i,]*PLH[1:(NS+1),VaR.K[i]]))/PDVAR
      ES.cont[i]=(PL[i]*PD[i])/(1-alpha)
    }
    else{
      VaR.cont[i]=(PL[i]*PD[i]*sum(W[i,]*PLH[1:(NS+1),VaR.K[i]]))/PDVAR
      ES.cont[i]=(PL[i]*PD[i]*sum(W[i,]*(1-CLH[1:(NS+1),(VaR.K[i]-1)])))/(1-alpha)
    }
  }
      
  ES.alt=this@ES[l]
  VaR=this@VaR[l]
  ES=this@ES[l]
  VaR.pos=this@VaR.pos[l]
  alpha=this@alpha.crp[l]-this@PDF[VaR.pos]
  while(VaR<ES){                                                                                     # calculation of expected shortfall TAU
    VaR.pos=VaR.pos-1
    alpha=alpha-this@PDF[VaR.pos]
    ES=(this@EL-sum((0:(VaR.pos-2))*this@loss.unit*this@PDF[1:(VaR.pos-1)]))/(1-alpha)
  }
        
  cat("CR+ Expected Shortfall TAU(",alpha+this@PDF[VaR.pos],"):",fo(ES),"\n")
        
  VaR.K=VaR.pos-this@PL.crp/this@loss.unit                                                                  # calculation of the contributions for ES TAU 
  ES.tau.cont=ES.cont                                                                                # (analog to ES)
  for(i in 1:NC){
    if(VaR.K[i]<=1)
       ES.tau.cont[i]=(PL[i]*PD[i])/(1)
    else           
    	 ES.tau.cont[i]=(PL[i]*PD[i]*sum(W[i,]*(1-CLH[1:(NS+1),(VaR.K[i]-1)])))/(1-alpha)
  }
        
  cat("Scale Factor for TAU",sum(VaR.cont)/sum(ES.tau.cont),"\n")
  ES.tau.cont=ES.tau.cont*sum(VaR.cont)/sum(ES.tau.cont)  		                                       # upscaling
      
  if(round(sum(VaR.cont),0)==round(VaR,0))
	   cat("Sum Check VaR: OK\n")
  else
	   cat("Sum Check VaR: NOT OK, Difference:",fo(sum(VaR.cont)-VaR),"\n")
  if(round(sum(ES.cont),0)==round(ES.alt,0))
	   cat("Sum Check ES: OK\n")
  else
	   cat("Sum Check ES: NOT OK, Difference:",fo(sum(this@ES.cont)-ES.alt),"\n")
  
  remove(list=c("a","B","VaR.pos","PL","PD","W")) 
  
  this@VaR.cont=VaR.cont
  remove(VaR.cont)
  this@ES.cont=ES.cont
  remove(ES.cont)
  this@ES.tau.cont=ES.tau.cont
  
  remove(list=c("ES.tau.cont","alpha.M","PLH","CLH","VaR.K"))
  this@rc.OK=TRUE
	cat("Done....\n")
  this@changes.rc.vares=FALSE
  gc()
	return(this)
}
)

setGeneric("rc.sd",function(this) standardGeneric("rc.sd"))
setMethod(f="rc.sd",signature="crp.CSFP",definition=function(this){
#
#       <rc.sd>     Calculation of SD Contributions
#
#       Last Modified:  24/06/2013
#       Author:         Dr. Matthias Fischer, Stefan Kolb & Kevin Jakob
#

  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(this@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
  if(this@changes.calc.portfolio.statistics)
    cat("WARNING: Inputparameter effecting 'calc.portfolio.statistics()' have been changed, without running 'calc.portfolio.statistics()' and the following routines afterwards.\n")
  if(this@changes.loss)
    cat("WARNING: Inputparameter effecting 'loss()' have been changed, without running 'loss()' and the following routines afterwards.\n")
  if(this@changes.measure)
    cat("WARNING: Inputparameter effecting 'measure()' have been changed, without running 'measure()' and the following routines afterwards.\n")   
  
  cat("Calculating SD contributions....\n")

      
  NC=this@NC                                                                                         # Make all needed variables local to this
  NS=this@NS                                                                                         # function instead of going up the environment
  PL=this@PL.crp                                                                                         # every calling of 'this' . This will increase
  PD=this@PD.crp                                                                                         # speed.
  W=this@W
  SD.cont=as.numeric(this@NC)
  SD=this@SD
  sigma_k=this@sigma_k
  loss.k=this@loss.k
    
  if(sum(W[,2:(NS+1)])==0){                                                                          # just ideosyncratic factor exists
    for(i in 1:NC)
      SD.cont[i]=(PL[i]*PD[i]/SD)*PL[i]
  }
  else{
    for(i in 1:NC){                                                                                  # calculation based on formula (2.27) in 
      k=sum(sigma_k^2*W[i,-1]*loss.k)                                                                   # Gundlach/Lehrbass (2003)
      SD.cont[i]=(PL[i]*PD[i]/SD)*(PL[i]+k)
    }
  }
  if(round(sum(SD.cont),0)==round(SD,0))
	  cat("Sum Check SD: OK\n")
  else
	  cat("Sum Check SD: NOT OK, Difference:",fo(sum(SD.cont)-SD),"\n")

  remove(list=c("PL","PD","W","sigma_k","loss.k"))
  this@SD.cont=SD.cont
  remove(SD.cont)    
    
  cat("Done....\n")
  this@changes.rc.sd=FALSE
  gc()
  return(this)
}
)


setGeneric("export",function(this) standardGeneric("export"))
setMethod(f="export",signature=c("crp.CSFP"),definition=function(this){
#
#       <export>    Export the calculated Risk Contributions to a CSV file              
#
#       Last Modified:  24/06/2013
#       Author:         Dr. Matthias Fischer, Stefan Kolb & Kevin Jakob
#
#
  if(nchar(this@path.out)==0)
    stop("path.out is unspecified!")

  if(this@changes.read)
    cat("WARNING: Inputparameter effecting 'read()' have been changed, without running 'read()' and the following routines afterwards.\n")
  if(this@changes.plausi)
    cat("WARNING: Inputparameter effecting 'plausi()' have been changed, without running 'plausi()' and the following routines afterwards.\n")
  if(this@changes.calc.portfolio.statistics)
    cat("WARNING: Inputparameter effecting 'calc.portfolio.statistics()' have been changed, without running 'calc.portfolio.statistics()' and the following routines afterwards.\n")
  if(this@changes.loss)
    cat("WARNING: Inputparameter effecting 'loss()' have been changed, without running 'loss()' and the following routines afterwards.\n")
  if(this@changes.measure)
    cat("WARNING: Inputparameter effecting 'measure()' have been changed, without running 'measure()' and the following routines afterwards.\n") 
  if(this@changes.rc.vares)
    cat("WARNING: Inputparameter effecting 'rc.vares()' have been changed, without running 'rc.vares()' and the following routines afterwards.\n") 
  if(this@changes.rc.sd)
    cat("WARNING: Inputparameter effecting 'rc.sd()' have been changed, without running 'rc.sd()' and the following routines afterwards.\n")   
  
  if(deparse(substitute(this))!="this")                                                              # synchronize model name to the calling model
    this@name=deparse(substitute(this))                                                              # if it is not called from crp.CSFP()

  if(length(this@loss)==1 && this@save.memory)                                                       # recreate losses if not available
    this@loss=this@loss.unit*(0:this@M)
  if(length(this@CDF)==1 && this@save.memory)                                                        # recalculate CDF if not available
    this@CDF=cumsum(this@PDF)
  
  cat("Exporting loss distribution...\n")
  FILE.TITLE="lossdist"
  FILE.name=paste(this@path.out,this@name,"/",FILE.TITLE,".csv",sep="")
  dir.create(paste(this@path.out,this@name,"/",sep=""),showWarnings=FALSE)
  if(this@file.format=="csv")
    write.csv(matrix(cbind(this@loss,this@PDF,this@CDF),nrow=this@M+1,ncol=3,dimnames=list(1:length(this@PDF),c("Loss","PDF","CDF"))), file=FILE.name,row.names=FALSE)
  else if(this@file.format=="csv2")
    write.csv2(matrix(cbind(this@loss,this@PDF,this@CDF),nrow=this@M+1,ncol=3,dimnames=list(1:length(this@PDF),c("Loss","PDF","CDF"))), file=FILE.name,row.names=FALSE)
  
  cat("Writing summary file...\n")
  this=write.summary(this)                                                                       # write input/output data to file
  
  if(this@calc.rc && this@rc.OK){
    cat("Exporting risk contributions...\n")
    if(this@NS==1)
      COL.nameS=c("CP.NR","PL","CP.rating","Sektor:1","RC_EL","SD","SD/SD_Gesamt","VaR","VaR/VaR_Gesamt","TAU","TAU/(SUMME TAU)","ES","ES/ES_Gesamt")
    else
      COL.nameS=c("CP.NR","PL","CP.rating","Sektor:1",as.character(c(2:this@NS)),"RC_EL","SD","SD/SD_Gesamt", "VaR","VaR/VaR_Gesamt","TAU","TAU/(SUMME TAU)","ES","ES/ES_Gesamt")
    RC=matrix(,nrow=this@NC,ncol=12+this@NS,dimnames=list(1:this@NC,COL.nameS))                                                           # collect risk contribution in one matrix
    RC[,1]=this@CP.NR
    RC[,2]=this@PL.crp
    RC[,3]=this@CP.rating
    RC[,4:(this@NS+3)]=this@W[,2:(this@NS+1)]
    RC[,this@NS+4]=this@PL.crp*this@PD.crp
    RC[,this@NS+5]=this@SD.cont
    RC[,this@NS+6]=this@SD.cont/sum(this@SD.cont)
    RC[,this@NS+7]=this@VaR.cont
    RC[,this@NS+8]=this@VaR.cont/sum(this@VaR.cont)
    RC[,this@NS+9]=this@ES.tau.cont
    RC[,this@NS+10]=this@ES.tau.cont/sum(this@ES.tau.cont)
    RC[,this@NS+11]=this@ES.cont
    RC[,this@NS+12]=this@ES.cont/sum(this@ES.cont)
    
    
    FILE.TITLE="RC"# File name of loss distribution
    dir.create(paste(this@path.out,this@name,"/",sep=""),showWarnings=FALSE)
    if(this@file.format=="csv")
      write.csv(RC, file=paste(this@path.out,this@name,"/",FILE.TITLE,".csv",sep=""),row.names=FALSE)
    else if(this@file.format=="csv2")
      write.csv2(RC, file=paste(this@path.out,this@name,"/",FILE.TITLE,".csv",sep=""),row.names=FALSE)
    cat("Done....\n")
    remove(RC)
  }
  
  if(this@save.memory){                                                                              # dont store risk contributions after export 
    this@ES.cont=0                                                                                   # in save.memory mode
    this@VaR.cont=0
    this@ES.tau.cont=0
    this@SD.cont=0
    this@loss=0
    this@CDF=0
  }
	this@changes.export=FALSE
  gc()
	return(this)
}
)


setGeneric("crp.CSFP",function(this,skip.read) standardGeneric("crp.CSFP"))
setMethod(f="crp.CSFP",signature=c("crp.CSFP","logical"),definition=function(this,skip.read){
  #
  #       <crp.CSFP>    Main routine of portfolio model              
  #
  #       Last Modified:  24/06/2013
  #       Author:         Kevin Jakob
  #
  
  if(deparse(substitute(this))!="this")                                                              # synchronize model name to the calling model    
    this@name=deparse(substitute(this))                                                                # synchronize model name to the calling model

  org.save.memory=this@save.memory                                                                   # save original state of save.memory and 
  this@save.memory=FALSE                                                                             # change to FALSE during calculations in order 
                                                                                                     # to increase performance
  if(!skip.read){
    this=read(this)
    if(!this@read.OK)
      return(this)
  }
  this=plausi(this)
  if(!this@plausi.OK)
    return(this)
  this=calc.portfolio.statistics(this)
  this=loss.dist(this)
  this=measure(this)
  if(this@PLOT.PDF)
    this=plot(this)
  else
    cat("You dont want to plot the PDF, otherwise set attribute PLOT.PDF=TRUE (default case)!\n")
  if(this@calc.rc){
    this=rc.vares(this)
    this=rc.sd(this)
  }
  else
    cat("You dont want to calculate the Risk Contributions otherwise set attribute calc.rc=TRUE!\n")
  
  if(this@export.to.file)
    this=export(this)
  this@save.memory=org.save.memory                                                                   # restore original save.memory state
    
  if(this@save.memory && !this@calc.rc){
    this@B=matrix()
    this@a=0
    this@ES.cont=0
    this@VaR.cont=0
    this@ES.tau.cont=0
    this@SD.cont=0
  }
  if(this@save.memory){
    this@loss=0
    this@CDF=0
    this@PD=0
    this@PL=0
  }
  gc()
  return(this)
}
)

setMethod(f="crp.CSFP",signature=c("crp.CSFP","missing"),definition=function(this){
  if(deparse(substitute(this))!="this")                                                              # synchronize model name to the calling model    
    this@name=deparse(substitute(this)) 
  crp.CSFP(this,skip.read=FALSE)})

setMethod(f="summary",signature=c("crp.CSFP"),definition=function(object){
    S=list(name=object@name,sec.var.est=object@sec.var.est,
           loss.unit=object@loss.unit,Niter.max=object@Niter.max,NC=object@NC,
           NS=object@NS,EL=object@EL,EL.crp=object@EL.crp,SD=object@SD,
           SD.crp=object@SD.crp,SIGMA.DIV=object@sigma_sqr_div,
           SIGMA.SYST=object@sigma_sqr_syst,alpha.max=object@alpha.max,
           PL=sum(object@PL),alpha=object@alpha,VaR=object@VaR,EC=object@EC,
           ES=object@ES)
    return(S)
}
)

setGeneric("write.summary",function(this) standardGeneric("write.summary"))
setMethod(f="write.summary",signature=c("crp.CSFP"),definition=function(this){
    
  if(nchar(this@path.out)==0)
    stop("path.out is unspecified!")
  if(deparse(substitute(this))!="this")
        this@name=deparse(substitute(this))
  dir.create(paste(this@path.out,this@name,"/",sep=""),showWarnings=FALSE)
  if(this@file.format=="csv")
    write.csv(summary(this),paste(this@path.out,this@name,"/summary.csv",sep=""),row.names=FALSE)
  else if(this@file.format=="csv2")
    write.csv2(summary(this),paste(this@path.out,this@name,"/summary.csv",sep=""),row.names=FALSE)

  return(this)
}
)

setGeneric("integrity.check",function(method,this) standardGeneric("integrity.check"))
setMethod(f="integrity.check",signature=c("character","crp.CSFP"),definition=function(method,this){

  if(method=="read"){
    if(this@changes.read)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="plausi"){
    if(this@changes.read || this@changes.plausi)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="calc.portfolio.statistics"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="loss"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics || this@changes.loss)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="measure"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics || this@changes.loss || this@changes.measure)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="plot"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics || this@changes.loss || this@changes.measure || this@changes.plot)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="rc.vares"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics || this@changes.loss || this@changes.measure || this@changes.rc.vares)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="rc.sd"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics || this@changes.loss || this@changes.measure || this@changes.rc.sd)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  else if(method=="export"){
    if(this@changes.read || this@changes.plausi || this@changes.calc.portfolio.statistics || this@changes.loss || this@changes.measure || this@changes.rc.vares || this@changes.rc.sd || this@changes.export)
      cat("WARNING: Inputparameter possibly effecting this variable have been changed without running all necessary routines afterwards.\n ")
  }
  return()
}
)

setGeneric("set.changes",function(method,this) standardGeneric("set.changes"))
setMethod(f="set.changes",signature=c("character","crp.CSFP"),definition=function(method,this){
  
  if(method=="read"){
    this@changes.read=TRUE
    this@changes.plausi=TRUE
    this@changes.calc.portfolio.statistics=TRUE
    this@changes.loss=TRUE
    this@changes.measure=TRUE
    this@changes.plot=TRUE
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.rc.vares=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.rc.sd=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0) || !(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="plausi"){
    this@changes.plausi=TRUE
    this@changes.calc.portfolio.statistics=TRUE
    this@changes.loss=TRUE
    this@changes.measure=TRUE
    this@changes.plot=TRUE
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.rc.vares=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.rc.sd=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0) || !(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="calc.portfolio.statistics"){
    this@changes.calc.portfolio.statistics=TRUE
    this@changes.loss=TRUE
    this@changes.measure=TRUE
    this@changes.plot=TRUE
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.rc.vares=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.rc.sd=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0) || !(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="loss"){
    this@changes.loss=TRUE
    this@changes.measure=TRUE
    this@changes.plot=TRUE
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.rc.vares=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.rc.sd=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0) || !(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="measure"){
    this@changes.measure=TRUE
    this@changes.plot=TRUE
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.rc.vares=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.rc.sd=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0) || !(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="plot"){
    this@changes.plot=TRUE
  }
  if(method=="rc.vares"){
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.rc.vares=TRUE
    if(!(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="rc.sd"){
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.rc.sd=TRUE
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0))
      this@changes.export=TRUE
  }
  if(method=="export"){
    if(!(length(this@SD.cont)==1 && this@SD.cont[1]==0) || !(length(this@VaR.cont)==1 && this@VaR.cont[1]==0))
      this@changes.export=TRUE
  }
  return(this)
})


setMethod("show", signature = "crp.CSFP", definition = function(object){
  
  if(any(object@changes.read,object@changes.plausi,object@changes.calc.portfolio.statistics,object@changes.loss,object@changes.measure,object@changes.plot,object@changes.rc.vares,object@changes.rc.sd,object@changes.export))
    cat("WARNING: Input parameters have been changed without recalculation of the model. Displayed results may be incorrect.\n")
  S=summary(object)
  cat("Model name:",S$name,"\n")
  cat("Number of counterparties:",S$NC,"\n")
  cat("Number of sectors:",S$NS,"\n")
  cat("Expected loss:",fo(S$EL),"(analytical);",fo(S$EL.crp),"(discretized)\n")
  cat("      SD loss:",fo(S$SD),"(analytical);",fo(S$SD.crp),"(discretized)\n")
  cat("Systematic SD:",fo(sqrt(S$SIGMA.SYST)),"; Diversifiable SD:",fo(sqrt(S$SIGMA.DIV)),"\n")
  temp=data.frame(S$alpha,fo(S$EC),fo(S$VaR),fo(S$ES))
  colnames(temp)=c("alpha     ","EC        ","VaR       ","ES   ")
  write.matrix(temp)
  cat("Loss unit:",fo(S$loss.unit),"\n")
  if(S$Niter.max<1)
    cat("Maximum level of CDF:",S$Niter.max,"\n")
  else
    cat("Maximum number of iterations for loss distribution:",fo(S$Niter.max),"\n")
  cat("Calculation of risk contributions:",object@calc.rc,"\n")
  
})




  

setGeneric("path.in",function(this) standardGeneric("path.in"))
setMethod(f="path.in",signature=c("crp.CSFP"),definition=function(this) this@path.in)
setGeneric("path.in<-",function(this,value) standardGeneric("path.in<-"))
setReplaceMethod(f="path.in",signature=c("crp.CSFP","character"),definition=function(this,value){this@path.in=value;this=set.changes("read",this);this})

setGeneric("path.out",function(this) standardGeneric("path.out"))
setMethod(f="path.out",signature=c("crp.CSFP"),definition=function(this) this@path.out)
setGeneric("path.out<-",function(this,value) standardGeneric("path.out<-"))
setReplaceMethod(f="path.out",signature=c("crp.CSFP","character"),definition=function(this,value){this@path.out=value;this=set.changes("export",this);this})

setGeneric("port.name",function(this) standardGeneric("port.name"))
setMethod(f="port.name",signature=c("crp.CSFP"),definition=function(this) this@port.name)
setGeneric("port.name<-",function(this,value) standardGeneric("port.name<-"))
setReplaceMethod(f="port.name",signature=c("crp.CSFP","character"),definition=function(this,value){this@port.name=value;this=set.changes("read",this);this})

setGeneric("name",function(this) standardGeneric("name"))
setMethod(f="name",signature=c("crp.CSFP"),definition=function(this) this@name)
setGeneric("name<-",function(this,value) standardGeneric("name<-"))
setReplaceMethod(f="name",signature=c("crp.CSFP","character"),definition=function(this,value){this@name=value;this=set.changes("read",this);this})

setGeneric("rating.scale.name",function(this) standardGeneric("rating.scale.name"))
setMethod(f="rating.scale.name",signature=c("crp.CSFP"),definition=function(this) this@rating.scale.name)
setGeneric("rating.scale.name<-",function(this,value) standardGeneric("rating.scale.name<-"))
setReplaceMethod(f="rating.scale.name",signature=c("crp.CSFP","character"),definition=function(this,value){this@rating.scale.name=value;this=set.changes("read",this);this})

setGeneric("sec.var.name",function(this) standardGeneric("sec.var.name"))
setMethod(f="sec.var.name",signature=c("crp.CSFP"),definition=function(this) this@sec.var.name)
setGeneric("sec.var.name<-",function(this,value) standardGeneric("sec.var.name<-"))
setReplaceMethod(f="sec.var.name",signature=c("crp.CSFP","character"),definition=function(this,value){this@sec.var.name=value;this=set.changes("read",this);this})

setGeneric("sec.var.est",function(this) standardGeneric("sec.var.est"))
setMethod(f="sec.var.est",signature=c("crp.CSFP"),definition=function(this) this@sec.var.est)
setGeneric("sec.var.est<-",function(this,value) standardGeneric("sec.var.est<-"))
setReplaceMethod(f="sec.var.est",signature=c("crp.CSFP","numeric"),definition=function(this,value){    if(value==5 || this@sec.var.est==5)
                                                                                                         this=set.changes("read",this)
                                                                                                       else
                                                                                                         this=set.changes("calc.portfolio.statistics",this)
                                                                                                       this@sec.var.est=value
                                                                                                       this})

setGeneric("loss.unit",function(this) standardGeneric("loss.unit"))
setMethod(f="loss.unit",signature=c("crp.CSFP"),definition=function(this) this@loss.unit)
setGeneric("loss.unit<-",function(this,value) standardGeneric("loss.unit<-"))
setReplaceMethod(f="loss.unit",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(value<=0)
    stop("Loss unit has to be positive!\n")
  this@loss.unit=value
  this=set.changes("calc.portfolio.statistics",this)
  this})

setGeneric("Niter.max",function(this) standardGeneric("Niter.max"))
setMethod(f="Niter.max",signature=c("crp.CSFP"),definition=function(this) this@Niter.max)
setGeneric("Niter.max<-",function(this,value) standardGeneric("Niter.max<-"))
setReplaceMethod(f="Niter.max",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(!value>0)
    stop("Niter.max has to be positive!\n")
  this@Niter.max=value
  this=set.changes("loss",this)
  this
  })

setGeneric("alpha",function(this) standardGeneric("alpha"))
setMethod(f="alpha",signature=c("crp.CSFP"),definition=function(this) this@alpha)
setGeneric("alpha<-",function(this,value) standardGeneric("alpha<-"))
setReplaceMethod(f="alpha",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(any(value<=0) || any(value>=1))
    stop("alpha has to be between 0 and 1!\n")
  this@alpha=value
  this=set.changes("measure",this)
  this})

setGeneric("PLOT.PDF",function(x,...) standardGeneric("PLOT.PDF"))
setMethod(f="PLOT.PDF",signature="crp.CSFP",definition=function(x,...){return(x@PLOT.PDF)})
setGeneric("PLOT.PDF<-",function(x,value) standardGeneric("PLOT.PDF<-"))
setReplaceMethod(f="PLOT.PDF",signature=c("crp.CSFP","logical"),definition=function(x,value){x@PLOT.PDF=value;x=set.changes("plot",x);x})

setGeneric("export.to.file",function(this) standardGeneric("export.to.file"))
setMethod(f="export.to.file",signature=c("crp.CSFP"),definition=function(this) this@export.to.file)
setGeneric("export.to.file<-",function(this,value) standardGeneric("export.to.file<-"))
setReplaceMethod(f="export.to.file",signature=c("crp.CSFP","logical"),definition=function(this,value){this@export.to.file=value;this=set.changes("export",this);this})

setGeneric("calc.rc",function(this) standardGeneric("calc.rc"))
setMethod(f="calc.rc",signature=c("crp.CSFP"),definition=function(this) this@calc.rc)
setGeneric("calc.rc<-",function(this,value) standardGeneric("calc.rc<-"))
setReplaceMethod(f="calc.rc",signature=c("crp.CSFP","logical"),definition=function(this,value){this@calc.rc=value;this=set.changes("rc.vares",this);this})

setGeneric("save.memory",function(this) standardGeneric("save.memory"))
setMethod(f="save.memory",signature=c("crp.CSFP"),definition=function(this) this@save.memory)
setGeneric("save.memory<-",function(this,value) standardGeneric("save.memory<-"))
setReplaceMethod(f="save.memory",signature=c("crp.CSFP","logical"),definition=function(this,value){this@save.memory=value;this=set.changes("loss",this);this})

setGeneric("PLOT.range.x",function(x,...) standardGeneric("PLOT.range.x"))
setMethod(f="PLOT.range.x",signature=c("crp.CSFP"),definition=function(x,...) x@PLOT.range.x)
setGeneric("PLOT.range.x<-",function(x,value) standardGeneric("PLOT.range.x<-"))
setReplaceMethod(f="PLOT.range.x",signature=c("crp.CSFP","numeric"),definition=function(x,value){x@PLOT.range.x=value;x=set.changes("plot",x);x})

setGeneric("PLOT.range.y",function(x,...) standardGeneric("PLOT.range.y"))
setMethod(f="PLOT.range.y",signature=c("crp.CSFP"),definition=function(x,...) x@PLOT.range.y)
setGeneric("PLOT.range.y<-",function(x,value) standardGeneric("PLOT.range.y<-"))
setReplaceMethod(f="PLOT.range.y",signature=c("crp.CSFP","numeric"),definition=function(x,value){x@PLOT.range.y=value;x=set.changes("plot",x);x})

setGeneric("PLOT.scale",function(x,...) standardGeneric("PLOT.scale"))
setMethod(f="PLOT.scale",signature=c("crp.CSFP"),definition=function(x,...) x@PLOT.scale)
setGeneric("PLOT.scale<-",function(x,value) standardGeneric("PLOT.scale<-"))
setReplaceMethod(f="PLOT.scale",signature=c("crp.CSFP","numeric"),definition=function(x,value){x@PLOT.scale=value;x=set.changes("plot",x);x})

setGeneric("file.format",function(this) standardGeneric("file.format"))
setMethod(f="file.format",signature=c("crp.CSFP"),definition=function(this) this@file.format)
setGeneric("file.format<-",function(this,value) standardGeneric("file.format<-"))
setReplaceMethod(f="file.format",signature=c("crp.CSFP","character"),definition=function(this,value){
  if(!(value=="csv" || value=="csv2")){
   cat("Wrong specification of file.format. Choose between csv and csv2.\n")
   return(this)
  }
  this@file.format=value;this=set.changes("read",this);this})


setGeneric("NS",function(this) standardGeneric("NS"))
setMethod(f="NS",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("read",this)
  this@NS})

setGeneric("NC",function(this) standardGeneric("NC"))
setMethod(f="NC",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("read",this)
  this@NC})

setGeneric("sec.var",function(this) standardGeneric("sec.var"))
setMethod(f="sec.var",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("plausi",this)
  this@sec.var})
setGeneric("sec.var<-",function(this,value) standardGeneric("sec.var<-"))
setReplaceMethod(f="sec.var",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(length(value)!=this@NS)
    stop("Length of new value (",length(value),") is unequal number of sectors in the model (",this@NS,")!\n")
  if(any(value<=0))
    stop("Sector variances have to be positive")
  this@sec.var=value
  this=set.changes("plausi",this)
  this})
  

setGeneric("NEX",function(this) standardGeneric("NEX"))
setMethod(f="NEX",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("plausi",this)
  this@NEX})
setGeneric("NEX<-",function(this,value) standardGeneric("NEX<-"))
setReplaceMethod(f="NEX",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(any(value<0))
    stop("NEX has to be non negative!\n")
  if(length(value)!=this@NC)
    stop("Length of new value (",length(value),") and number of counterparties (",this@NC,") are unequal!\n")
  this@NEX=value
  this=set.changes("plausi",this)
  this})


setGeneric("LGD",function(this) standardGeneric("LGD"))
setMethod(f="LGD",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@LGD})
setGeneric("LGD<-",function(this,value) standardGeneric("LGD<-"))
setReplaceMethod(f="LGD",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(any(value<0) || any(value>1))
    stop("LGDs are supposed to be between 0 and 1!\n")
  if(length(value)!=this@NC)
    stop("Length of new value (",length(value),") and number of counterparties (",this@NC,") are unequal!\n")
  this@LGD=value
  this=set.changes("calc.portfolio.statistics",this)
  this})


setGeneric("PL",function(this) standardGeneric("PL"))
setMethod(f="PL",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@PL})

setGeneric("PD",function(this) standardGeneric("PD"))
setMethod(f="PD",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@PD})

setGeneric("EL",function(this) standardGeneric("EL"))
setMethod(f="EL",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@EL})

setGeneric("nu",function(this) standardGeneric("nu"))
setMethod(f="nu",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@nu})

setGeneric("PL.crp",function(this) standardGeneric("PL.crp"))
setMethod(f="PL.crp",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@PL.crp})

setGeneric("PD.crp",function(this) standardGeneric("PD.crp"))
setMethod(f="PD.crp",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@PD.crp})

setGeneric("CP.rating",function(this) standardGeneric("CP.rating"))
setMethod(f="CP.rating",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@CP.rating})
setGeneric("CP.rating<-",function(this,value) standardGeneric("CP.rating<-"))
setReplaceMethod(f="CP.rating",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(length(value)!=this@NC)
    stop("Length of new value (",length(value),") and number of counterparties (",this@NC,") are unequal!\n")
  for(i in 1:this@NC){
    if(!is.element(value[i],this@rating))
      stop("Counterparties ratings does not correspond to rating possible rating grades!\n")
  }
  this@CP.rating=value
  this=set.changes("calc.portfolio.statistics",this)
  this})


setGeneric("CP.NR",function(this) standardGeneric("CP.NR"))
setMethod(f="CP.NR",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@CP.NR})
setGeneric("CP.NR<-",function(this,value) standardGeneric("CP.NR<-"))
setReplaceMethod(f="CP.NR",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(length(value)!=this@NC)
    stop("Length of new value (",length(value),") and number of counterparties (",this@NC,") are unequal!\n")
  this@CP.NR=value
  this=set.changes("calc.portfolio.statistics",this)
  this})


setGeneric("W",function(this) standardGeneric("W"))
setMethod(f="W",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("plausi",this)
  this@W})
setGeneric("W<-",function(this,value) standardGeneric("W<-"))
setReplaceMethod(f="W",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(nrow(value)!=this@NC)
    stop("Number of rows is unequal to number of counterparties in the model!\n")
  if(ncol(value)!=this@NS)
    stop("number of columns is unequal to number of sectors in the model!\n ")
  if(this@NS>1){                                                                                     # check if sector weights are plausible
    if(sum(rowSums(value[,-1])<=1)!=nrow(value[,-1])|sum(value[,-1]>=0)!=nrow(value[,-1])*ncol(value[,-1]))
      stop("Sum of weights have to be between 0 and 1 for each counterparty!\n")
  }
  else if(this@NS==1){
    if(any(value>=0 & value<=1)==FALSE)
      stop("Sum of weights have to be between 0 and 1 for each counterparty!\n")
  }  
  this@W=value
  this=set.changes("plausi",this)
  this})


setGeneric("rating.PD",function(this) standardGeneric("rating.PD"))
setMethod(f="rating.PD",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("plausi",this)
  this@rating.PD})
setGeneric("rating.PD<-",function(this,value) standardGeneric("rating.PD<-"))
setReplaceMethod(f="rating.PD",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(length(value)!=length(this@rating))
    cat("Warning: Slots rating and rating.PD are of different length! Please adjust rating \n")
  if(length(value)!=length(rating.SD) && this@sec.var.est!=5)
    cat("Warning: Slots rating and rating.SD are of different length! Please adjust rating.SD \n")
  if(any(diff(value)<0))
    stop("PDs have to be in ascending order!\n")
  this@rating.PD=value
  this=set.changes("plausi",this)
  this})

setGeneric("rating",function(this) standardGeneric("rating"))
setMethod(f="rating",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@rating})
setGeneric("rating<-",function(this,value) standardGeneric("rating<-"))
setReplaceMethod(f="rating",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(length(value)!=length(this@rating.PD))
    cat("Warning: Slots rating and rating.PD are of different length! Please adjust rating.PD \n")
  if(length(value)!=length(rating.SD) && this@sec.var.est!=5)
    cat("Warning: Slots rating and rating.SD are of different length! Please adjust rating.SD \n")
  if(any(diff(value)<0))
    stop("Rating classes have to be in ascending order!\n")
  this@rating=value
  this=set.changes("calc.portfolio.statistics",this)
  this})


setGeneric("rating.SD",function(this) standardGeneric("rating.SD"))
setMethod(f="rating.SD",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@rating.SD})
setGeneric("rating.SD<-",function(this,value) standardGeneric("rating.SD<-"))
setReplaceMethod(f="rating.SD",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(length(value)!=length(this@rating))
    cat("Warning: Slots rating and rating.PD are of different length! Please adjust rating \n")
  if(length(value)!=length(rating.PD) && this@sec.var.est!=5)
    cat("Warning: Slots rating and rating.SD are of different length! Please adjust rating.PD \n")  
  this@rating.SD=value
  this=set.changes("calc.portfolio.statistics",this)
  this})


setGeneric("M",function(this) standardGeneric("M"))
setMethod(f="M",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@M})

setGeneric("mu.k",function(this) standardGeneric("mu.k"))
setMethod(f="mu.k",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@mu.k})

setGeneric("loss.k",function(this) standardGeneric("loss.k"))
setMethod(f="loss.k",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@loss.k})

setGeneric("sigma_k",function(object,...) standardGeneric("sigma_k"))
setMethod(f="sigma_k",signature=c("crp.CSFP"),definition=function(object,...){
  integrity.check("calc.portfolio.statistics",object)
  object@sigma_k})

setGeneric("sigma_sqr_div",function(object,...) standardGeneric("sigma_sqr_div"))
setMethod(f="sigma_sqr_div",signature=c("crp.CSFP"),definition=function(object,...){
  integrity.check("calc.portfolio.statistics",object)
  object@sigma_sqr_div})

setGeneric("sigma_sqr_syst",function(object,...) standardGeneric("sigma_sqr_syst"))
setMethod(f="sigma_sqr_syst",signature=c("crp.CSFP"),definition=function(object,...){
  integrity.check("calc.portfolio.statistics",object)
  object@sigma_sqr_syst})

setGeneric("SD",function(this) standardGeneric("SD"))
setMethod(f="SD",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("calc.portfolio.statistics",this)
  this@SD})

setGeneric("alpha.max",function(this) standardGeneric("alpha.max"))
setMethod(f="alpha.max",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@alpha.max})
setGeneric("alpha.max<-",function(this,value) standardGeneric("alpha.max<-"))
setReplaceMethod(f="alpha.max",signature=c("crp.CSFP","numeric"),definition=function(this,value){
  if(value<=0 || value>1)
    stop("alpha.max has to be between 0 and 1!\n")  
  this@alpha.max=value
  this=set.changes("loss",this)
  this})


setGeneric("a",function(this) standardGeneric("a"))
setMethod(f="a",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@a})

setGeneric("PDF",function(this) standardGeneric("PDF"))
setMethod(f="PDF",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@PDF})

setGeneric("CDF",function(this) standardGeneric("CDF"))
setMethod(f="CDF",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@CDF})

setGeneric("B",function(this) standardGeneric("B"))
setMethod(f="B",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@B})

setGeneric("loss",function(this) standardGeneric("loss"))
setMethod(f="loss",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("loss",this)
  this@loss})

setGeneric("VaR",function(this) standardGeneric("VaR"))
setMethod(f="VaR",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@VaR})

setGeneric("EC",function(this) standardGeneric("EC"))
setMethod(f="EC",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@EC})

setGeneric("ES",function(this) standardGeneric("ES"))
setMethod(f="ES",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@ES})

setGeneric("VaR.cont",function(this) standardGeneric("VaR.cont"))
setMethod(f="VaR.cont",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("rc.vares",this)
  this@VaR.cont})

setGeneric("ES.cont",function(this) standardGeneric("ES.cont"))
setMethod(f="ES.cont",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("rc.vares",this)
  this@ES.cont})

setGeneric("ES.tau.cont",function(this) standardGeneric("ES.tau.cont"))
setMethod(f="ES.tau.cont",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("rc.vares",this)
  this@ES.tau.cont})

setGeneric("EL.crp",function(this) standardGeneric("EL.crp"))
setMethod(f="EL.crp",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@EL.crp})

setGeneric("SD.crp",function(this) standardGeneric("SD.crp"))
setMethod(f="SD.crp",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@SD.crp})

setGeneric("VaR.pos",function(this) standardGeneric("VaR.pos"))
setMethod(f="VaR.pos",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@VaR.pos})

setGeneric("alpha.crp",function(this) standardGeneric("alpha.crp"))
setMethod(f="alpha.crp",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("measure",this)
  this@alpha.crp})

setGeneric("SD.cont",function(this) standardGeneric("SD.cont"))
setMethod(f="SD.cont",signature=c("crp.CSFP"),definition=function(this){
  integrity.check("rc.sd",this)
  this@SD.cont})

setGeneric("changes.read",function(this) standardGeneric("changes.read"))
setMethod(f="changes.read",signature=c("crp.CSFP"),definition=function(this){
  this@changes.read
  })

setGeneric("changes.plausi",function(this) standardGeneric("changes.plausi"))
setMethod(f="changes.plausi",signature=c("crp.CSFP"),definition=function(this){
  this@changes.plausi
  })

setGeneric("changes.calc.portfolio.statistics",function(this) standardGeneric("changes.calc.portfolio.statistics"))
setMethod(f="changes.calc.portfolio.statistics",signature=c("crp.CSFP"),definition=function(this){
  this@changes.calc.portfolio.statistics
  })

setGeneric("changes.loss",function(this) standardGeneric("changes.loss"))
setMethod(f="changes.loss",signature=c("crp.CSFP"),definition=function(this){
  this@changes.loss
  })

setGeneric("changes.measure",function(this) standardGeneric("changes.measure"))
setMethod(f="changes.measure",signature=c("crp.CSFP"),definition=function(this){
  this@changes.measure
  })

setGeneric("changes.plot",function(this) standardGeneric("changes.plot"))
setMethod(f="changes.plot",signature=c("crp.CSFP"),definition=function(this){
  this@changes.plot
  })

setGeneric("changes.rc.vares",function(this) standardGeneric("changes.rc.vares"))
setMethod(f="changes.rc.vares",signature=c("crp.CSFP"),definition=function(this){
  this@changes.rc.vares
  })

setGeneric("changes.rc.sd",function(this) standardGeneric("changes.rc.sd"))
setMethod(f="changes.rc.sd",signature=c("crp.CSFP"),definition=function(this){
  this@changes.rc.sd
  })

setGeneric("changes.export",function(this) standardGeneric("changes.export"))
setMethod(f="changes.export",signature=c("crp.CSFP"),definition=function(this){
  this@changes.export
  })

