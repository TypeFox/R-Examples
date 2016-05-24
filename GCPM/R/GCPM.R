setClass(Class="GCPM",                                                                               #define the names and types of attributes of GCPM
	   representation=representation(
							model.type="character",
              default="character",
              link.function="character",
              loss.unit="numeric",
							NS="numeric",
							NC="numeric",
              name="character",
							NR="numeric",
							EAD="numeric",
							LGD="numeric",
							PL="numeric",
							PD="numeric",
              business="factor",
              country="factor",
							EL.analyt="numeric",
							EL="numeric",
							nu="numeric",
							PL.disc="numeric",
							PD.disc="numeric",
              sec.var="numeric",
              sector.names="factor",
							SD.div="numeric",
							SD.syst="numeric",
							SD.analyt="numeric",
							SD="numeric",
							W="matrix",
              idiosyncr="numeric",
							alpha.max="numeric",
							a="numeric",
							PDF="numeric",
							CDF="numeric",
							B="matrix",
							loss="numeric",
              random.numbers="matrix",
              LHR="numeric",
              N="numeric",
              scenarios="numeric",
              seed="numeric",
              loss.thr="numeric",
              sim.losses="numeric",
              CP.sim.losses="matrix",
              max.entries="numeric"),
	   prototype=prototype( 
          model.type="CRP",
          default="Poisson",
          link.function="CRP",
          loss.unit=1e6,               
					NS=0,
					NC=0,
          name="",
					NR=0,
					EAD=0,
					LGD=0,
					PL=0,
					PD=0,
          business=factor(),
          country=factor(),
					EL.analyt=0,
					EL=0,
					nu=0,
					PL.disc=0,
					PD.disc=0,
          sec.var=0,
          sector.names=factor(),
					SD.div=0,
					SD.syst=0,
					SD.analyt=0,
					SD=0,
					W=matrix(),
          idiosyncr=0,
					alpha.max=0,
					a=0,
					PDF=0,
					CDF=0,
					B=matrix(),
					loss=0,
          random.numbers=matrix(),
          LHR=0,
          N=1e4,
          scenarios=1,
          seed=0,
          loss.thr=0,
          sim.losses=0,
          CP.sim.losses=matrix(),
          max.entries=1e6)
)

setGeneric("read",function(this,portfolio) standardGeneric("read"))
setMethod("read",c("GCPM","data.frame"),function(this,portfolio) {

#portfolio data
#####################################
message("Importing portfolio data....")
this@NC=nrow(portfolio)
  if(is.null(portfolio$Number)){
    warning("The Number-column is missing in the portfolio\n")
    this@NR=1:this@NC
  }
  else
    this@NR=as.numeric(portfolio$Number)

  if(is.null(portfolio$Name)){
    warning("The Name-column is missing in the portfolio!")
    this@name=paste("CP",this@NR)
  }
  else
    this@name=as.character(portfolio$Name)

  if(is.null(portfolio$Default)){
    stop("The Default-column is missing in the portfolio.\n")
  }
  else{
    if(any(!is.element(portfolio$Default,c("Bernoulli","Poisson"))))
      stop("Entries of Default column in the portfolio have to be 'Bernoulli' or 'Poisson'")
    this@default=as.character(portfolio$Default)
  }

  if(is.null(portfolio$PD)){
    stop("The PD-column is missing in the portfolio.\n")
  }
  else
 	  this@PD=as.numeric(portfolio$PD)
 	
  if(is.null(portfolio$EAD)){
    stop("The EAD-column is missing in the portfolio.\n")
  }
  else
 	  this@EAD=as.numeric(portfolio$EAD)
  
  if(is.null(portfolio$LGD)){
    stop("The LGD-column is missing in the portfolio.\n")
  }
  else
 	  this@LGD=as.numeric(portfolio$LGD)

  if(is.null(portfolio$Country)){
    warning("The Country-column is missing in the portfolio.")
    this@country=as.factor(rep(1,this@NC))
  }
  else
    this@country=as.factor(portfolio$Country)
  
  if(is.null(portfolio$Business)){
    warning("The Business-column is missing in the portfolio.")
    this@business=as.factor(rep(1,this@NC))
  }
  else
    this@business=as.factor(portfolio$Business)
  temp=(do.call(cbind,portfolio))[,-(1:8)]
  temp=as.matrix(temp)
  if(ncol(temp)==1)
    colnames(temp)=colnames(portfolio)[9]
  this@NS=ncol(temp)
  this@sector.names=as.factor(colnames(temp))
  #this@sector.names=this@sector.names[match(levels(this@sector.names),this@sector.names)]
  temp=temp[,match(this@sector.names,colnames(temp))]
  temp=as.matrix(temp)
  if(ncol(temp)==1)
    colnames(temp)=colnames(portfolio)[9]
  if(is.null(temp)){
    stop("The sector-weights are missing in the portfolio.\n")
  }
  else{
    temp=matrix(as.numeric(temp),nrow=this@NC,ncol=this@NS,dimnames = list(1:nrow(temp),colnames(temp)))
    
    remove(portfolio)
    if(this@NS==1)
      temp=matrix(ncol=1,temp)                                                                       # convert to matrix
    this@W=temp
    remove(temp)
  }
  
  #exclude counterparties without EAD, pd==0 or lgd==0
  noEAD=which(this@EAD==0)
  noLGD=which(this@LGD==0)
  noPD=which(this@PD<=0)
  if(this@link.function=="CRP")
    PD1=which(this@PD>=1 & this@default=="Bernoulli")
  else if(this@link.function=="CM")
    PD1=which(this@PD>=1)
  else
    stop("Wrong specification of link.function. Choose between CRP and CM.")
  emptyCP=union(noEAD,union(noLGD,union(noPD,PD1)))
  if(length(emptyCP)>0){
    this@NR=this@NR[-emptyCP]
    this@name=this@name[-emptyCP]
    this@business=this@business[-emptyCP]
    this@country=this@country[-emptyCP]
    this@PD=this@PD[-emptyCP]
    this@EAD=this@EAD[-emptyCP]
    this@LGD=this@LGD[-emptyCP]
    this@W=this@W[-emptyCP,]
    this@NC=this@NC-length(emptyCP)
  }
  
  message(this@NS," sectors ...")
  message(this@NC," counterparties (",length(emptyCP)," removed due to EAD=0 (",length(noEAD),"), lgd=0 (",length(noLGD),"), pd<=0 (",length(noPD),") pd>=1 (",length(PD1),"))\n")
  gc()
  return(this)
}
)

setGeneric("plausi",function(this) standardGeneric("plausi"))
setMethod(f="plausi",signature=c("GCPM"),definition=function(this){ 
  
  # sector variances in case of CRP
  ###################################  
  if(this@model.type=="CRP"){   
    if(length(this@sec.var)<this@NS){
      stop("Number of sector variances specified (",length(this@sec.var),") is smaller then the number of sectors within the portfolio (",this@NS,").\n",sep="")
    }
    else{
      temp=which(!is.na(match(this@sector.names,gsub(" ",".",names(this@sec.var)))))
      if(length(temp)<this@NS)
        stop("Not all sectors names in the portfolio can be found as sector names")
      else{
        this@sec.var=this@sec.var[temp]
        this@sec.var=this@sec.var[match(this@sector.names,names(this@sec.var))]
      }
    }
  }  
  
  # Random numbers (sector realisations) in case of simulative model
  #########################################
   if(this@model.type=="simulative"){
    if(ncol(this@random.numbers)<this@NS){
      stop("Number of columns of  random numbers (",ncol(this@random.numbers),") is smaller then the number of sectors within the portfolio (",this@NS,").\n",sep="")
    }
    else{
      temp=which(!is.na(match(this@sector.names,levels(as.factor(colnames(this@random.numbers))))))
      if(length(temp)<this@NS)
        stop("Not all sectors names in the portfolio can be found as colnames for random numbers.")
      else if(this@NS>1){
        this@random.numbers=as.matrix(this@random.numbers[,temp])
        this@random.numbers=as.matrix(this@random.numbers[,match(this@sector.names,colnames(this@random.numbers))])
      }
      else{
        this@random.numbers=as.matrix(this@random.numbers[,temp])
        colnames(this@random.numbers)=this@sector.names
      }
    }
  }
  
  if(any(this@EAD<0))                                                                               # check if EAD is not negative 
    stop("EAD can not be negative. \n")

  if(any(this@LGD<0) || any(this@LGD>1))
    stop("LGDs are supposed to be between 0 and 1!\n")
  
  if(this@model.type=="CRP"){
    if(this@NS>1){                                                                                     # check if sector weights are plausible
      if(any(rowSums(this@W[,-1])>1) | any(this@W[,-1]<0)){
        stop("Sum of sector weights have to be between 0 and 1 for each counterparty!\n")
      }
    }
    else if(this@NS==1){
  	  if(any(this@W<0 | this@W>1)){
  	    stop("Sum of sector weights have to be between 0 and 1 for each counterparty!\n")
  	  }
    }                                                    
    if(any(this@sec.var<=0)){                                                                         # check if secor variances are positive    stop("Sector variances have to be positive\n")
      stop("Sector variances have to be positive!")
    }
  }
  gc()
  return(this)                                                                        
}
)
setGeneric("sector.info",function(this) standardGeneric("sector.info"))
setMethod(f="sector.info",signature=c("GCPM"),definition=function(this){
  if(this@model.type=="CRP"){
    mu.k=loss.k=numeric(this@NS)
    for(k in 1:this@NS){
      select=which(this@W[,k]>0)
      if(length(select)>0){
        mu.k[k]=sum(this@W[select,k]*this@PD.disc[select])                                             # formula (39), CSFP 1997, p. 43
        loss.k[k]=sum(this@W[select,k]*this@PD.disc[select]*this@PL.disc[select])
      }
      else{
        mu.k[k]=loss.k[k]=0
      }
    }
    remove(select)    
  }
  else{
    stop("No sector info available for simulative model")
  }
  return(list(mu=mu.k,loss=loss.k))
})

setGeneric("calc.portfolio.statistics",function(this) standardGeneric("calc.portfolio.statistics"))
setMethod(f="calc.portfolio.statistics",signature=c("GCPM"),definition=function(this){

	this@PL=this@EAD*this@LGD                                                                         # potential loss
  this@nu=gcpm.round(this@PL/this@loss.unit)                                                         # integer muliples of loss.unit representing poZential loss
	this@nu=ifelse(this@nu==0,1,this@nu)                                                               # should not be 0
	this@PL.disc=this@nu*this@loss.unit                                                                     # discretized potential loss
  this@PD.disc=this@PD*(this@PL/this@PL.disc)                                                                # transformed PD
  this@EL.analyt=sum(this@PL*this@PD)                                                                     # expected loss
  
  message("Portfolio statistics....")
  message("Loss unit: ",fo(this@loss.unit))  
  message("Portfolio EAD:",fo(sum(as.numeric(this@EAD))))
  message("Portfolio potential loss:",fo(sum(as.numeric(this@PL))))
  message("Portfolio expected loss:",fo(this@EL.analyt),"(analytical)")
  
  if(this@model.type=="CRP"){
    SI=sector.info(this)
    this@SD.syst=sqrt(sum(this@sec.var*SI$loss*SI$loss))  			                                   # portfolio standard diviation
    this@SD.div=sqrt(sum(this@PL.disc^2*this@PD.disc))                                              # diversifible risk
    
    message("Diversifible risk: ",fo(this@SD.div),"  Systematic risk: ",fo(this@SD.syst))
    this@SD.analyt=sqrt((this@SD.div)^2+(this@SD.syst)^2)
    message("Portfolio standad deviation:",fo(this@SD.analyt),"(analytical)")
  }
  gc()
	return(this)
}
)

setGeneric("loss.dist.sim",function(this,Ncores) standardGeneric("loss.dist.sim"))
setMethod(f="loss.dist.sim",signature=c("GCPM","numeric"),definition=function(this,Ncores){

LHR=(this@LHR[this@scenarios][1:this@N]/sum(this@LHR[this@scenarios][1:this@N]))[1:this@N]                                     #normalizing LHR to sum=1
sim.losses=rep(0,this@N)

#C++ simulation
if(this@NS==1)
  Sigma=matrix(1,1,1)
else
  Sigma=cor(this@random.numbers[this@scenarios,])
default.distr.int=numeric(this@NC)
for(i in 1:this@NC)
  default.distr.int[i]=which(this@default[i]==c("Bernoulli","Poisson"))

link.function.int=which(this@link.function==c("CRP","CM"))
if(this@link.function=="CRP")
  W=cbind(this@idiosyncr,this@W)
else if(this@link.function=="CM")
  W=this@W
W=as.matrix(W)
calc.rc=ifelse(is.finite(this@loss.thr),1,0)
message("Starting simulation (",this@N,"simulations )")

Ncores=max(round(Ncores),1)
if(Ncores>1){
  if(Ncores>=detectCores()){
    Ncores=detectCores()-1
    warning("Ncores reduced to ", Ncores," due to system restrictions.")
  }
  message("Parallel computing on ",Ncores," cores (no progress bar)")
  
  do.tasks<-function(X,default.distr.int,link.function.int,random.numbers,Sigma,W,
                     PD.disc,PL.disc,calc.rc,loss.thr,max.entries,Ncores,seed){
    library(GCPM)
    set.seed(seed = ifelse(is.na(seed),as.integer(Sys.time()),seed))
    N=nrow(random.numbers)
    start=end=rep(0,Ncores)
    start[1]=1
    end[1]=ceiling(N/Ncores)
    for(i in 2:Ncores){
      start[i]=end[i-1]+1
      end[i]=start[i]+ceiling(N/Ncores)
    }
    end=pmin(end,N)
    start=pmin(start,N)
    ret=cpploss(default.distr.int,link.function.int,as.matrix(random.numbers[start[X]:end[X],]),
                as.matrix(Sigma),as.matrix(W),PD.disc,PL.disc,calc.rc,loss.thr,max.entries)
    return(ret)    
  }
  
  cluster <- makeCluster(Ncores)
  retPar=parLapplyLB(cluster,X = 1:Ncores,fun = do.tasks,default.distr.int,link.function.int,as.matrix(this@random.numbers[this@scenarios,]),as.matrix(Sigma),as.matrix(W),this@PD.disc,this@PL.disc,calc.rc,this@loss.thr,this@max.entries,Ncores,this@seed)
  stopCluster(cluster)
  
  #merge simulated losses
  sim.losses=retPar[[1]]$simlosses
  if(calc.rc){
    CP.sim.losses=retPar[[1]]$CPsimlosses
    loss.szenarios=retPar[[1]]$lossszenarios
  }
  for(i in 2:Ncores){ 
    if(calc.rc){
      CP.sim.losses=cbind(CP.sim.losses,retPar[[i]]$CPsimlosses)
      loss.szenarios=c(loss.szenarios,retPar[[i]]$lossszenarios+length(sim.losses))
    }
    sim.losses=c(sim.losses,retPar[[i]]$simlosses)
  }  
}
else{ #single core
  set.seed(seed = ifelse(is.na(this@seed),as.integer(Sys.time()),this@seed))
  ret=cpploss(default.distr.int,link.function.int,as.matrix(this@random.numbers[this@scenarios,]),as.matrix(Sigma), as.matrix(W),this@PD.disc,this@PL.disc,calc.rc,this@loss.thr,this@max.entries)
  sim.losses=ret$simlosses
  
  if(sim.losses[length(sim.losses)]==-1)
    stop("C++ simulation was terminated by the user\n")
  if(calc.rc){
    CP.sim.losses=ret$CPsimlosses
    loss.szenarios=ret$lossszenarios
  }
}

if(calc.rc){
  if(length(loss.szenarios)>0){
    if(any(loss.szenarios==-1)){
      warning("Number of loss scenarios exceeded memory limit. Increase loss.thr or max.entries.\nRisk contributions are not available.\n")
      calc.rc=FALSE
    }
  }
  else{
    warning("No loss scenarios stored because of too high value for loss.thr.\nRisk contributions are not available.\n")
    calc.rc=FALSE
  }
}
this@sim.losses=sim.losses
if(calc.rc)
  this@CP.sim.losses=matrix(CP.sim.losses,ncol=ncol(CP.sim.losses),nrow=this@NC,dimnames=list(this@NR,loss.szenarios))
message("Simulation finished\n")
gc()
##########################################
message("Calculating loss distribution...")
if(abs(this@loss.unit)>=1e-8){
  sim.losses=round(sim.losses/this@loss.unit)*this@loss.unit
  if(calc.rc)
    this@CP.sim.losses=round(this@CP.sim.losses/this@loss.unit)*this@loss.unit
}
tab=table(sim.losses)
loss.levels=sort(as.numeric(names(tab)))
this@loss=loss.levels
this@sim.losses=sim.losses
if(all(LHR==LHR[1]))
  PDF=as.numeric(tab)*LHR[1]
else{
  PDF=rep(0,length(this@loss))
  for(i in 1:length(loss.levels)){
    select=which(abs(sim.losses-loss.levels[i])<1e-8)
    PDF[i]=sum(LHR[select])
    sim.losses=sim.losses[-select]
    LHR=LHR[-select]
  }
}
this@PDF=PDF
this@CDF=cumsum(this@PDF)
this@alpha.max=this@CDF[length(this@CDF)]
gc()
return(this)
})

setGeneric("loss.dist.crp",function(this) standardGeneric("loss.dist.crp"))
setMethod(f="loss.dist.crp",signature=c("GCPM"),definition=function(this){

  if(!all(this@default=="Poisson")){
    warning("If model.type==CRP, only Poisson default distribution will be used")
    this@default=rep("Poisson",this@NC)
  }
  message("Calculate the loss distribution till ",this@alpha.max,"-confidence level is reached.",sep="")  
  M=1e5 #max number of exposure bands
  A=matrix(0,ncol=M,nrow=this@NS)                                                           # Make all needed variables local to this 
  B=A                                                                                                # function instead of going up the environment 
  a=rep(0,M)                                                                                  # every calling of 'this' . This will increase 
  PDF=rep(0,M)                                                                                # speed but needs more memory.
  CDF=rep(0,M)
  sigma.k=sqrt(this@sec.var)
  SI=sector.info(this)
  loss.k=SI$loss
  mu.k=SI$mu
  
  nu=this@nu
  W=cbind(this@idiosyncr,this@W)
  PD=this@PD.disc

##############################################################################################################
#       1. Calculation of PDF
##############################################################################################################

  if(sum(W[,2:(this@NS+1)])==0){                                                                     # just idiosyncratic
    a[1]=-sum(W[,1]*PD)
    PDF[1]=exp(a[1])
    CDF[1]=PDF[1]

    while(j<M && CDF[j]<this@alpha.max){
      select=which(nu==j)
      sum1=0
      if(length(select)>0)
	      sum1=sum(W[select,1]*PD[select])
      a[j+1]=sum1
      PDF[j+1]=sum((1:j)/j*a[2:(j+1)]*PDF[j:1])
      CDF[j+1]=CDF[j]+PDF[j+1]
      j=j+1
    }
  }
  else{                                                                          # standard case
    for(k in 1:this@NS){
	    A[k,1]=1+sigma.k[k]^2*mu.k[k]
	    B[k,1]=-log(A[k,1])
	  }
	  a[1]=-sum(W[,1]*PD)+sum(B[,1]/(sigma.k[1:this@NS]^2))
	  PDF[1]=exp(a[1])
	  CDF[1]=PDF[1]
	
	  j=0
    pb=txtProgressBar(style = 3,min = 0.1)
	  while(j<M && CDF[j+1]<this@alpha.max){
		  j=j+1
      select=which(nu==j)
		  for(k in 1:this@NS){
		    if(length(select)>0)
			    A[k,j+1]=sigma.k[k]^2*sum(W[select,k+1]*PD[select])
		    if(j+1==2)
			    B[k,2]=A[k,2]/A[k,1]
		    else
			    B[k,j+1]=1/A[k,1]*(A[k,j+1]+(1/j*sum((1:(j-1))*B[k,2:j]*A[k,j:2])))
		  }		
      sum1=ifelse(length(select>0),sum(W[select,1]*PD[select]),0)
		  a[j+1]=sum1+sum(B[,j+1]/(sigma.k[1:this@NS]^2))
      PDF[j+1]=sum((1:j)/j*a[2:(j+1)]*PDF[j:1])
		  CDF[j+1]=CDF[j]+PDF[j+1]
      if(j%%1e2==0){
        setTxtProgressBar(pb,value = 10^(CDF[j+1]-this@alpha.max))
      }
	  }
    message("")
	  this@PDF=PDF[1:(j+1)]
	  this@CDF=CDF[1:(j+1)]
	  this@a=a[1:(j+1)]
	  this@B=B[,-((j+1):M)]
    remove(list=c("A","select","W","sigma.k","mu.k","PD","PDF","CDF","a","B"))
  }
##############################################################################################################
#       3. Summarize
##############################################################################################################

  this@loss=(0:j)*this@loss.unit                                                                     # create losses
  this@alpha.max=this@CDF[length(this@CDF)]
  message("Calculation completed...")
  message("Reached level of confidence: ",this@alpha.max," ( iterations actually done: ",j," )\n")
  gc()
  return(this)        
}
)

setGeneric("measure",function(this,alpha) standardGeneric("measure"))
setMethod(f="measure",signature=c("GCPM","missing"),definition=function(this,alpha){
  return(measure(this,numeric()))
  })

setMethod(f="measure",signature=c("GCPM","numeric"),definition=function(this,alpha){

  message("Calculating risk measures from loss distribution....")        
  this@EL=sum(as.numeric(this@loss*this@PDF))                                                    # expected loss
  message("Expected loss from loss distribution: ",fo(this@EL)," (deviation from EL calculated from portfolio data: ",round(this@EL/this@EL.analyt-1,4)*100,"%)")  
  alpha.EL=this@CDF[min(which(num.geq(this@loss,this@EL.analyt)))]
  message("Exceedance Probability of the expected loss:",1-alpha.EL)
  message("Portfolio mean expected loss exceedance: ",fo(ES(this,alpha.EL)))
  this@SD=sqrt(sum(as.numeric(this@PDF*(this@loss-this@EL)^2)))                              # standard deviation
  message("Portfolio loss standard deviation:",fo(this@SD),"\n")  
  
  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  
  if(length(alpha)==0)
    return(this)
    
  L.alpha=length(alpha)
  VaR=VaR.pos=EC=ES=numeric(L.alpha)
  for(i in 1:L.alpha){
	  if(num.geq(this@alpha.max,alpha[i])){					                                                       # check if alpha-quantile exists
      VaR.pos[i]=index(this,alpha[i])
      VaR[i]=VaR(this = this,alpha[i])
		}
		else
		  warning(paste("Value-at-risk(",alpha[i],") is not available, increase alpha.max",sep=""))
    EC[i]=VaR[i]-this@EL.analyt
    
	  if(num.geq(this@alpha.max,alpha[i])){                                                                # check if alpha-quantiles exists
      if(VaR.pos[i]==1)
        ES[i]=this@EL.analyt                                                                           # calculate ES
      else
        ES[i]=ES(this,alpha[i]) 
		}
	  else
		  warning(paste("Expected Shortfall (",alpha[i],") is not available, increase alpha.max",sep=""))
    
    #output
		message("Portfolio Economic Capital(",alpha[i],"): ",fo(EC[i]))
		message("Portfolio Value-at-risk(",alpha[i],"): ",fo(VaR[i]))
		message("Portfolio Expected Shortfall(",alpha[i],"): ",fo(ES[i]),"\n")        
	}  
  gc()
	return(this)
}
)

setMethod("plot","GCPM",function(x,y,alpha=numeric(),nbins=100){

  if(missing(y))
    y=1
  xlab=paste("Loss in ",fo(y),sep="")
  ylab="Probability"
  main="Portfolio Loss Distribution"
  if(x@model.type=="CRP")
    plot(x@loss/y,x@PDF,type="l",main=main,xlab=xlab,ylab="Probability",cex.axis=1.3,cex.main=1.3,lwd=3,cex.lab=1.3)  
  else{
    breaks=seq(0,max(x@loss)+1,length.out = nbins+1)
    pdf=loss=numeric(nbins)
    for(i in 2:length(loss)){
      select=which(x@loss>=breaks[i-1] & x@loss<breaks[i])
      pdf[i-1]=sum(x@PDF[select])
      loss[i-1]=(breaks[i-1]+breaks[i])/2
    }
    plot(loss/y,pdf,type="l",main=main,xlab=xlab,ylab="Probability",cex.axis=1.3,cex.main=1.3,lwd=3,cex.lab=1.3) 
  }
  
  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  if(length(alpha)==0){
    abline(v=x@EL/y,col="green")   # add lines for EL
    legend("topright",c("EL"),col=c("green"),lwd=1)                        # add legend
  }
  else{
    VaR=ES=alpha.VaR=alpha.ES=c()
    for(i in 1:length(alpha)){
      if(!is.na(VaR(x,alpha[i]))){
        VaR=c(VaR,VaR(x,alpha[i]))
        alpha.VaR=c(alpha.VaR,alpha[i])
      }
      if(!is.na(ES(x,alpha[i]))){
        ES=c(ES,ES(x,alpha[i]))
        alpha.ES=c(alpha.ES,alpha[i])
      }
    }
    abline(v=x@EL/y,col="green")   # add lines for EL
    
    abline(v=c(VaR,ES)/y,col=c(rep("blue",length(VaR)),rep("red",length(ES))))   # add lines for VaR and ES
    legend("topright",c("EL","VaR","ES"),col=c("green","blue","red"),lwd=1)                        # add legend
    text(x = VaR/y,y=0.2*par("yaxp")[2],labels = alpha.VaR,cex = 0.7,font=2)
    text(x = ES/y,y=0.5*par("yaxp")[2],labels = alpha.ES,cex = 0.7,font=2)
    
    #add labels with numeric values
    text(x = x@EL/y,y=0.15*par("yaxp")[2],labels = fo(x@EL),cex = 0.65)
    text(x = VaR/y,y=0.15*par("yaxp")[2],labels = fo(VaR),cex = 0.65)
    text(x = ES/y,y=0.45*par("yaxp")[2],labels = fo(ES),cex = 0.65)
  }
  gc()
  invisible(x)
}
)

setGeneric("EC.VaR.ES.cont",function(this,alpha,type) standardGeneric("EC.VaR.ES.cont"))
setMethod(f="EC.VaR.ES.cont",signature=c("GCPM","numeric","missing"),definition=function(this,alpha,type){
  return(EC.VaR.ES.cont(this,alpha,type="all"))
})

setMethod(f="EC.VaR.ES.cont",signature=c("GCPM","numeric","character"),definition=function(this,alpha,type){

  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  if(length(alpha)==0)
    return(matrix(NA,nrow=this@NC,ncol=3))
  
  L.alpha=length(alpha)
  RC=matrix(NA,nrow=this@NC,ncol=3*L.alpha)
  for(i in 1:L.alpha){
    if(this@model.type=="CRP"){
      temp=EC.VaR.ES.cont.single.crp(this,alpha[i],type)
      RC[,i]=temp[,1]
      RC[,i+L.alpha]=temp[,2]
      RC[,i+2*L.alpha]=temp[,3]
    }
    else if(this@model.type=="simulative"){
      temp=EC.VaR.ES.cont.single.sim(this,alpha[i],type)
      RC[,i]=temp[,1]
      RC[,i+L.alpha]=temp[,2]
      RC[,i+2*L.alpha]=temp[,3]
    }
    else
      stop("Wrong specification of model.type. Choose between CRP and simulative!")
  }
  colnames(RC)=c(paste("EC", alpha),paste("VaR", alpha),paste("ES", alpha))
  rownames(RC)=this@NR
  return(RC)
}
)

setGeneric("EC.VaR.ES.cont.single.sim",function(this,alpha,type) standardGeneric("EC.VaR.ES.cont.single.sim"))
setMethod("EC.VaR.ES.cont.single.sim",signature=c("GCPM","numeric","missing"),definition=function(this,alpha,type){
  return(EC.VaR.ES.cont.single.sim(this,alpha,"all"))
})
setMethod("EC.VaR.ES.cont.single.sim",signature=c("GCPM","numeric","character"),definition=function(this,alpha,type){
  
  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  if(length(alpha)==0)
    return(matrix(NA,nrow=this@NC,ncol=3))
  
  EC.cont=ES.cont=VaR.cont=rep(NA,this@NC)
  loss.szenarios=as.numeric(colnames(this@CP.sim.losses))
  min.loss=min(apply(this@CP.sim.losses,MARGIN=2,sum))
  LHR=(this@LHR[this@scenarios][1:this@N]/sum(this@LHR[this@scenarios][1:this@N]))[1:this@N]

  if(type=="ES" || type=="all"){
    loss=VaR(this,alpha)
    if(loss<min.loss){
      warning("Not enough szenarios available. Increase level alpha (currently ",alpha,") for ES contributions or lower loss.thr!")
    }
    else{
      select=which(num.geq(apply(this@CP.sim.losses,2,sum),loss))
      ES.cont=(this@CP.sim.losses[,select]%*%LHR[loss.szenarios[select]])/(sum(LHR[loss.szenarios[select]]))
    }
  }
  
  if(type=="VaR" || type=="all"){
    tau=search.tau(this,VaR(this,alpha))
    loss=VaR(this,tau)
    if(loss<min.loss){
      warning("Not enough szenarios available. Increase level alpha (currently ",alpha,") for VaR contributions or lower loss.thr!")
    }
    else{
      select=which(num.geq(apply(this@CP.sim.losses,2,sum),loss))
      VaR.cont=(this@CP.sim.losses[,select]%*%LHR[loss.szenarios[select]])/(sum(LHR[loss.szenarios[select]]))
    }
  }
  
  if(type=="EC" || type=="all"){
    tau=search.tau(this,VaR(this,alpha)-this@EL.analyt)
    loss=VaR(this,tau)
    if(loss<min.loss){
      warning("Not enough szenarios available. Increase level alpha (currently ",alpha,") for VaR contributions or lower loss.thr!")
    }
    else{     
      select=which(num.geq(apply(this@CP.sim.losses,2,sum),loss))
      EC.cont=(this@CP.sim.losses[,select]%*%LHR[loss.szenarios[select]])/(sum(LHR[loss.szenarios[select]]))
    }
  }
  
  if(any(type==c("all","VaR","EC"))){
    dev.EC=0
    dev.VaR=0
    if(type=="all" || type=="EC"){
      if(!any(is.na(EC.cont)))
        dev.EC=round((sum(EC.cont)-VaR(this,alpha)+this@EL.analyt)/(VaR(this,alpha)-this@EL.analyt)*100,2)
    }
    if(type=="all" || type=="VaR"){
      if(!any(is.na(VaR.cont)))
        dev.VaR=round((sum(VaR.cont)-VaR(this,alpha))/VaR(this,alpha)*100,2)
    }
    
    if(abs(dev.EC)>0.01 || abs(dev.VaR)>0.01){
      message("Deviation between VaR, EC contributions and VaR, EC caused by discontinuity of loss distribution:")
      if((type=="all" || type=="EC") && abs(dev.EC)>0.01)
        message("Sum-check (alpha = ",alpha," ) EC-cont: ",round((sum(EC.cont)-VaR(this,alpha)+this@EL.analyt)/(VaR(this,alpha)-this@EL.analyt)*100,2),"% deviation")                                                                                                                                                                  
      if((type=="all" || type=="VaR") && abs(dev.VaR)>0.01)
        message("Sum-check (alpha = ",alpha," ) VaR-cont: ",round((sum(VaR.cont)-VaR(this,alpha))/VaR(this,alpha)*100,2),"% deviation")
    }
  }
  temp=cbind(EC.cont,VaR.cont,ES.cont)
  return(temp)
})


setGeneric("EC.VaR.ES.cont.single.crp",function(this,alpha,type) standardGeneric("EC.VaR.ES.cont.single.crp"))
setMethod(f="EC.VaR.ES.cont.single.crp",signature=c("GCPM","numeric","missing"),definition=function(this,alpha,type){  
  return(EC.VaR.ES.cont.single.crp(this,alpha,"all"))
})

setMethod(f="EC.VaR.ES.cont.single.crp",signature=c("GCPM","numeric","character"),definition=function(this,alpha,type){  

  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  if(length(alpha)==0)
    return(matrix(NA,nrow=this@NC,ncol=3))
  
  VaR=VaR(this,alpha)
  ES=ES(this,alpha)
  EC.start=index(this,search.tau(this,VaR-this@EL.analyt))
  VaR.start=index(this,search.tau(this,VaR))
  ES.start=index(this,alpha)
  
  EC.K=EC.start-this@nu
  VaR.K=VaR.start-this@nu
  ES.K=ES.start-this@nu
  
  if(!(EC.start>0 && ES>0) || length(this@CDF)<=ES.start){                                                             # check if VaR and ES exists
	  warning("VaR and ES contributions are not available. \n")
    return(matrix(NA,nrow = this@NC,ncol=3))
	}
	if(length(this@a)==0 || all(dim(this@B)==1)){                                                      # check if required data exists
    warning("Necessary data for calculation of risk contributions (VaR/ES) are not availabel, probably you`re running in save.memory mode. \n Change save.memory to FALSE or calc.rc to TRUE, run loss and start calculation of risk contributions (VaR/ES) again.\n")
    return(matrix(NA,nrow = this@NC,ncol=3))
  }
  
  EC.cont=VaR.cont=ES.cont=numeric(this@NC)                                                                          # Make all needed variables local to this
  a=this@a                                                                                           # every calling of 'this' . This will increase
  B=this@B                                                                                           # speed.
  PL=this@PL.disc
  PD=this@PD.disc
  NS=this@NS
  W=cbind(this@idiosyncr,this@W)
  NC=this@NC    
  
  alpha.M=matrix(0,nrow=this@NS+1,ncol=ES.start)
  PLH=matrix(0,nrow=this@NS+1,ncol=length(this@CDF))
  PLH[1,]=this@PDF
  
  if(NS==1)
    B=matrix(nrow=1,B)
  a1=a[1]+B[,1]                                                                                      # for details please see p.87-90 in "Credit Risk in the Banking Industry".
  PLH[2:(NS+1),1]=exp(a1)                                                                            # calculation of the "distorted" PDF for each sector 
  for(h in 2:(NS+1)){
    alpha.M[h,]=a[1:ES.start]+B[h-1,1:ES.start]
    for(j in 2:ES.start)
      PLH[h,j]=sum((1:(j-1))/(j-1)*alpha.M[h,2:j]*PLH[h,(j-1):1])
  }
  CLH=PLH
  for(x in 1:(NS+1))                                                                                 # calculation of the CDF
    CLH[x,]=cumsum(PLH[x,])
  
  alpha.EC=ifelse(EC.start>1,this@CDF[EC.start-1],0)
  alpha.VaR=ifelse(VaR.start>1,this@CDF[VaR.start-1],0)
  alpha.ES=ifelse(ES.start>1,this@CDF[ES.start-1],0)
  
  if(type=="all" || type=="EC"){
    if(any(EC.K<=1))
      warning("Loss level alpha=",alpha," is to low for analytical risk-contributions.\n Counterparty contributions to EC may be overestimated.")
    for(i in 1:NC){                                                                                    
      if(EC.K[i]<=1)
        EC.cont[i]=(PL[i]*PD[i])/(1-alpha.EC)
      else
        EC.cont[i]=(PL[i]*PD[i]*sum(W[i,]*(1-CLH[1:(NS+1),(EC.K[i]-1)])))/(1-alpha.EC)
    }
  }
  
  if(type=="all" || type=="VaR"){
    if(any(VaR.K<=1))
      warning("Loss level alpha=",alpha," is to low for analytical risk-contributions.\n Counterparty contributions to VaR may be overestimated.")
    for(i in 1:NC){
      if(VaR.K[i]<=1)
        VaR.cont[i]=(PL[i]*PD[i])/(1-alpha.VaR)
      else
        VaR.cont[i]=(PL[i]*PD[i]*sum(W[i,]*(1-CLH[1:(NS+1),(VaR.K[i]-1)])))/(1-alpha.VaR)
    }
  }
  
  if(type=="all" || type=="ES"){
    if(any(ES.K<=1))
      warning("Loss level alpha=",alpha," is to low for analytical risk-contributions.\n Counterparty contributions to ES may be overestimated.")
    for(i in 1:NC){
      if(ES.K[i]<=1)
        ES.cont[i]=(PL[i]*PD[i])/(1-alpha.ES)
      else
        ES.cont[i]=(PL[i]*PD[i]*sum(W[i,]*(1-CLH[1:(NS+1),(ES.K[i]-1)])))/(1-alpha.ES)
    }
  }  

  temp=cbind(EC.cont,VaR.cont,ES.cont)
  gc()
	return(temp)
}
)

setGeneric("SD.cont",function(this) standardGeneric("SD.cont"))
setMethod(f="SD.cont",signature=c("GCPM"),definition=function(this){

  if(!this@model.type=="CRP")
    stop("Contributions to portolio standard deviation are only available for analytical CRP models.")
  
  NC=this@NC                                                                                         # Make all needed variables local to this
  NS=this@NS                                                                                         # function instead of going up the environment
  PL=this@PL.disc                                                                                         # every calling of 'this' . This will increase
  PD=this@PD.disc                                                                                         # speed.
  W=cbind(this@idiosyncr,this@W)
  SD.cont=as.numeric(this@NC)
  SD.analyt=this@SD.analyt
  sigma.k=sqrt(this@sec.var)
  loss.k=sector.info(this)$loss
    
  if(sum(W[,2:(NS+1)])==0){                                                                          # just ideosyncratic factor exists
    for(i in 1:NC)
      SD.cont[i]=(PL[i]*PD[i]/SD.analyt)*PL[i]
  }
  else{
    for(i in 1:NC){                                                                                  # calculation based on formula (2.27) in 
      k=sum(sigma.k^2*W[i,-1]*loss.k)                                                                   # Gundlach/Lehrbass (2003)
      SD.cont[i]=(PL[i]*PD[i]/SD.analyt)*(PL[i]+k)
    }
  }

  remove(list=c("PL","PD","W","sigma.k"))
  gc()
  names(SD.cont)=this@NR
  return(SD.cont)
}
)

setGeneric("export",function(this,path.out,file.format,alpha) standardGeneric("export"))
setMethod(f="export",signature=c("GCPM","missing","character","missing"),definition=function(this,path.out,file.format,alpha){
  export(this,"",file.format,numeric())
})
setMethod(f="export",signature=c("GCPM","character","missing","missing"),definition=function(this,path.out,file.format,alpha){
  export(this,path.out,"",numeric())
})
setMethod(f="export",signature=c("GCPM","missing","missing","missing"),definition=function(this,path.out,file.format,alpha){
  export(this,"","",numeric())
})
setMethod(f="export",signature=c("GCPM","character","character","missing"),definition=function(this,path.out,file.format,alpha){
  export(this,path.out,file.format,numeric())
})
setMethod(f="export",signature=c("GCPM","missing","character","numeric"),definition=function(this,path.out,file.format,alpha){
  export(this,"",file.format,alpha)
})
setMethod(f="export",signature=c("GCPM","character","missing","numeric"),definition=function(this,path.out,file.format,alpha){
  export(this,path.out,"",alpha)
})
setMethod(f="export",signature=c("GCPM","missing","missing","numeric"),definition=function(this,path.out,file.format,alpha){
  export(this,"","",alpha)
})
setMethod(f="export",signature=c("GCPM","character","character","numeric"),definition=function(this,path.out,file.format,alpha){

  return.files=(nchar(path.out)==0)
  if(nchar(file.format)==0 && !return.files)
    stop("file format is not specified!")
  
  if(!return.files){
    if(substr(path.out,nchar(path.out),nchar(path.out))!="/" && substr(path.out,nchar(path.out),nchar(path.out))!="\\")
      path.out=paste(path.out,"/",sep="")
    if(!any(file.format==c("csv","csv2")))
      stop("Wrong specification of file format.Choose between csv and csv2.")
  }
    
  ret=list()
  temp=matrix(cbind(this@loss,this@PDF,this@CDF),nrow=length(this@CDF),ncol=3,dimnames=list(1:length(this@PDF),c("Loss","PDF","CDF")))
  if(!return.files){
    message("Exporting loss distribution...")
    file.name=paste(path.out,"lossdist.csv",sep="")
    dir.create(paste(path.out,sep=""),showWarnings=FALSE)
    if(file.format=="csv")
      write.csv(temp, file=file.name,row.names=FALSE)
    else if(file.format=="csv2")
      write.csv2(temp, file=file.name,row.names=FALSE)
    
  }
  else
    ret$loss.dist=temp
  
  if(!return.files){
    message("Writing summary file...")
    this=write.summary(this,path.out,file.format)                                                                       # write input/output data to file
  }
  else
    ret$summary=summary(this)
  
  
  if(this@model.type=="CRP"){
    RC=data.frame(NR=this@NR,name=this@name,business=this@business,country=this@country,
                default=this@default,PL.disc=this@PL.disc,PD.disc=this@PD.disc,EL.disc=this@PD.disc*this@PL.disc,SD=SD.cont(this))
  }
  else{
    RC=data.frame(NR=this@NR,name=this@name,business=this@business,country=this@country,
                  default=this@default,PL.disc=this@PL.disc,PD.disc=this@PD.disc,EL.disc=this@PD.disc*this@PL.disc,SD=rep(NA,this@NC))
  } 
  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  if(length(alpha)>0){
    col.names=c(paste("EC",alpha),paste("VaR",alpha),paste("ES",alpha),this@sector.names)
    L.alpha=length(alpha)
    temp=matrix(NA,nrow=this@NC,ncol=3*L.alpha+this@NS)                                                           # collect risk contribution in one matrix
    temp[,1:(3*L.alpha)]=EC.VaR.ES.cont(this,alpha)
    temp[,(1+3*L.alpha):(3*L.alpha+this@NS)]=this@W
    dimnames(temp)=list(1:this@NC,col.names)
  }
  else{
    temp=this@W
    colnames(temp)=this@sector.names
  }
  RC=cbind(RC,temp)
  
  if(!return.files){
    dir.create(paste(path.out,sep=""),showWarnings=FALSE)
    if(file.format=="csv")
      write.csv(RC, file=paste(path.out,"RC.csv",sep=""),row.names=FALSE)
    else if(file.format=="csv2")
      write.csv2(RC, file=paste(path.out,"RC.csv",sep=""),row.names=FALSE)
    remove(RC)
  }
  else
    ret$RC=RC
  gc()
  if(!return.files)
	  return(this)
  else
    return(ret)
}
)


setGeneric("index",function(this,alpha) standardGeneric("index"))
setMethod("index",signature=c("GCPM","numeric"),definition=function(this,alpha){
  ret=numeric(length(alpha))
  for(i in 1:length(alpha)){
    if(num.geq(max(this@CDF),alpha[i]))
      ret[i]=min(which(num.geq(this@CDF,alpha[i])))
    else{
      warning("Index is not available, returning Inf. Decrease alpha!")
      ret[i]=Inf
    }
  }
  return(ret)
})

setGeneric("EC",function(this,alpha) standardGeneric("EC"))
setMethod("EC",signature=c("GCPM","missing"),definition=function(this,alpha){
  stop("Level alpha required.")
})
setMethod("EC",signature=c("GCPM","numeric"),definition=function(this,alpha){
  return(VaR(this,alpha)-this@EL)
})

setGeneric("VaR",function(this,alpha) standardGeneric("VaR"))
setMethod("VaR",signature=c("GCPM","missing"),definition=function(this,alpha){
  stop("Level alpha required.")
})
  
setMethod("VaR",signature=c("GCPM","numeric"),definition=function(this,alpha){
  if(any(alpha<=0 | alpha>=1))
    stop("Level alpha has to be in (0,1)")
  index=index(this,alpha)
  ret=numeric(length(alpha))
  for(i in 1:length(alpha)){
    if(index[i]<Inf)
      ret[i]=this@loss[index[i]]
    else{
      warning("VaR is not available, returning Inf. Decrease alpha!")
      ret[i]=Inf
    }
  }
  return(ret)
})

setGeneric("ES",function(this,alpha) standardGeneric("ES"))
setMethod("ES",signature=c("GCPM","missing"),definition=function(this,alpha){
  stop("Level alpha required.")
})

setMethod("ES",signature=c("GCPM","numeric"),definition=function(this,alpha){
  if(any(alpha<=0 | alpha>=1))
    stop("Level alpha has to be in (0,1)")
  index=index(this,alpha)
  ret=numeric(length(alpha))
  for(i in 1:length(alpha)){
    if(index[i]<Inf){
      if(this@model.type=="CRP"){
        if(index[i]>1){
          select=1:(index[i]-1)
          ret[i]=(this@EL.analyt-sum(this@PDF[select]*this@loss[select]))/(1-this@CDF[index[i]-1])
        }
        else
          ret[i]=this@EL.analyt
      }
      else if(this@model.type=="simulative"){
        if(index[i]>1){
          select=index[i]:length(this@CDF)
          ret[i]=sum(this@PDF[select]*this@loss[select])/(1-this@CDF[index[i]-1])
        }
        else
          ret[i]=this@EL.analyt
      }
      else
        stop("Wrong specification of model.type! Choose between CRP and simulative.")
    }
    else{
      warning("ES is not available, returning Inf. Decrease alpha!")
      ret[i]=Inf
    }
  }
  return(ret)
})

setGeneric("search.tau",function(this,loss) standardGeneric("search.tau"))
setMethod("search.tau",signature=c("GCPM","numeric"),definition=function(this,loss){
  if(loss>tail(this@loss,1)){ #KJA spaeter wieder einschalten
    warning("Desired loss level exceeds maximum value of realized losses, returning Inf!")
    return(Inf)
  }
  ret=numeric(length(loss))
  
  for(i in 1:length(loss)){
    if(num.geq(max(this@loss),loss))
      max.ind=min(which(num.geq(this@loss,loss),arr.ind = TRUE))
    else
      max.ind=length(this@loss)
    min.ind=1
    obj=Inf
    obj.old=Inf
    
    while(abs(max.ind-min.ind)>1){
      x=floor((min.ind+max.ind)/2)
      obj.old=obj
      obj=ES(this,this@CDF[x])
      if(obj>loss)
        max.ind=x
      else if(obj<loss)
        min.ind=x
      else if(obj==loss)
        ret[i]=this@CDF[x]
    }
    if(abs(ES(this,this@CDF[min.ind])-loss)>abs(ES(this,this@CDF[max.ind])-loss))
      ret[i]=this@CDF[max.ind]
    else
      ret[i]=this@CDF[min.ind]
  }
  return(ret)
})

setGeneric("analyze",function(this,portfolio,alpha,Ncores) standardGeneric("analyze"))
setMethod(f="analyze",signature=c("GCPM","data.frame","missing","missing"),definition=function(this,portfolio,alpha,Ncores){
  return(analyze(this,portfolio,numeric(),1))
})
setMethod(f="analyze",signature=c("GCPM","data.frame","numeric","missing"),definition=function(this,portfolio,alpha,Ncores){
  return(analyze(this,portfolio,alpha,1))
})
setMethod(f="analyze",signature=c("GCPM","data.frame","missing","numeric"),definition=function(this,portfolio,alpha,Ncores){
  return(analyze(this,portfolio,numeric(),Ncores))
})

setMethod(f="analyze",signature=c("GCPM","data.frame","numeric","numeric"),definition=function(this,portfolio,alpha,Ncores){

  if(any(alpha<=0 | alpha>=1))
    warning("Level alpha has to be in (0,1)")
  select=which(alpha>0 & alpha<1)
  if(length(select)>0)
    alpha=alpha[select]
  else
    alpha=numeric()
  
  this=read(this,portfolio)
  this=plausi(this)
  this=calc.portfolio.statistics(this)
  if(this@model.type=="CRP")
    this=loss.dist.crp(this)
  else if(this@model.type=="simulative")
    this=loss.dist.sim(this,Ncores)
  else
    stop("Wrong speficification of model.type. Choose beween CRp and simulative.")
  this=measure(this,alpha)
  gc()
  return(this)
}
)

  
setMethod(f="summary",signature=c("GCPM"),definition=function(object){
  S=list(model.type=object@model.type,link.function=object@link.function,loss.unit=object@loss.unit,
         NC=object@NC,NS=object@NS,EL.analyt=object@EL.analyt,EL=object@EL,SD.analyt=object@SD.analyt,
         SD=object@SD,SD.div=object@SD.div,SD.syst=object@SD.syst,
         alpha.max=object@alpha.max,EAD=sum(object@EAD),PL=sum(object@PL.disc))
  return(S)
})

setGeneric("write.summary",function(this,path.out,file.format) standardGeneric("write.summary"))
setMethod(f="write.summary",signature=c("GCPM","character","character"),definition=function(this,path.out,file.format){
  if(nchar(path.out)==0)
    stop("path.out is unspecified!")
  if(nchar(file.format)==0)
    stop("file format is unspecified!")
  
  if(substr(path.out,nchar(path.out),nchar(path.out))!="/" && substr(path.out,nchar(path.out),nchar(path.out))!="\\")
    path.out=paste(path.out,"/",sep="")
  if(!any(file.format==c("csv","csv2")))
    stop("Wrong specification of file format.Choose between csv and csv2.")
  
  dir.create(paste(path.out,sep=""),showWarnings=FALSE)
  if(file.format=="csv")
    write.csv(summary(this),paste(path.out,"summary.csv",sep=""),row.names=FALSE)
  else if(file.format=="csv2")
    write.csv2(summary(this),paste(path.out,"summary.csv",sep=""),row.names=FALSE)
  return(this)
}
)


setMethod("show", signature = "GCPM", definition = function(object){

  S=summary(object)
  message("Model type: ",S$model.type)
  message("Link function: ",S$link.function)  
  message("Number of counterparties: ",S$NC)
  message("Number of sectors: ",S$NS)
  message("Sum EAD: ",fo(S$EAD))
  message("Sum PL: ",fo(S$PL))  
  message("Expected loss: ",fo(S$EL)," (loss distribution); ",fo(S$EL.analyt)," (analytical)")
  if(S$model.type=="CRP"){
    message("      SD loss: ",fo(S$SD)," (loss distribution); ",fo(S$SD.analyt)," (analytical)")
    message("Systematic SD: ",fo(S$SD.syst),"; Diversifiable SD:",fo(S$SD.div))
  }  
  else if(S$model.type=="simulative")
    message("      SD loss: ",fo(S$SD)," (loss distribution);")
  message("Loss unit: ",fo(S$loss.unit))
  message("Maximum level of CDF: ",S$alpha.max)  
})

setGeneric("model.type",function(this) standardGeneric("model.type"))
setMethod(f="model.type",signature=c("GCPM"),definition=function(this) this@model.type)

setGeneric("default",function(this) standardGeneric("default"))
setMethod(f="default",signature=c("GCPM"),definition=function(this) this@default)

setGeneric("link.function",function(this) standardGeneric("link.function"))
setMethod(f="link.function",signature=c("GCPM"),definition=function(this) this@link.function)

setGeneric("name",function(this) standardGeneric("name"))
setMethod(f="name",signature=c("GCPM"),definition=function(this) this@name)

setGeneric("sec.var",function(this) standardGeneric("sec.var"))
setMethod(f="sec.var",signature=c("GCPM"),definition=function(this) this@sec.var)

setGeneric("business",function(this) standardGeneric("business"))
setMethod(f="business",signature=c("GCPM"),definition=function(this) this@business)

setGeneric("country",function(this) standardGeneric("country"))
setMethod(f="country",signature=c("GCPM"),definition=function(this) this@country)

setGeneric("sector.names",function(this) standardGeneric("sector.names"))
setMethod(f="sector.names",signature=c("GCPM"),definition=function(this) this@sector.names)

setGeneric("idiosyncr",function(this) standardGeneric("idiosyncr"))
setMethod(f="idiosyncr",signature=c("GCPM"),definition=function(this) this@idiosyncr)

setGeneric("random.numbers",function(this) standardGeneric("random.numbers"))
setMethod(f="random.numbers",signature=c("GCPM"),definition=function(this) this@random.numbers)

setGeneric("LHR",function(this) standardGeneric("LHR"))
setMethod(f="LHR",signature=c("GCPM"),definition=function(this) this@LHR)

setGeneric("N",function(this) standardGeneric("N"))
setMethod(f="N",signature=c("GCPM"),definition=function(this) this@N)

setGeneric("loss.thr",function(this) standardGeneric("loss.thr"))
setMethod(f="loss.thr",signature=c("GCPM"),definition=function(this) this@loss.thr)

setGeneric("loss.unit",function(this) standardGeneric("loss.unit"))
setMethod(f="loss.unit",signature=c("GCPM"),definition=function(this) this@loss.unit)

setGeneric("seed",function(this) standardGeneric("seed"))
setMethod(f="seed",signature=c("GCPM"),definition=function(this) this@seed)

setGeneric("NS",function(this) standardGeneric("NS"))
setMethod(f="NS",signature=c("GCPM"),definition=function(this){
  this@NS})

setGeneric("NC",function(this) standardGeneric("NC"))
setMethod(f="NC",signature=c("GCPM"),definition=function(this){
  this@NC})

setGeneric("EAD",function(this) standardGeneric("EAD"))
setMethod(f="EAD",signature=c("GCPM"),definition=function(this){
  this@EAD})

setGeneric("LGD",function(this) standardGeneric("LGD"))
setMethod(f="LGD",signature=c("GCPM"),definition=function(this){
  this@LGD})

setGeneric("EL.analyt",function(this) standardGeneric("EL.analyt"))
setMethod(f="EL.analyt",signature=c("GCPM"),definition=function(this){
  this@EL.analyt})

setGeneric("PL",function(this) standardGeneric("PL"))
setMethod(f="PL",signature=c("GCPM"),definition=function(this){
  this@PL.disc})

setGeneric("PD",function(this) standardGeneric("PD"))
setMethod(f="PD",signature=c("GCPM"),definition=function(this){
  this@PD.disc})


setGeneric("NR",function(this) standardGeneric("NR"))
setMethod(f="NR",signature=c("GCPM"),definition=function(this){
  this@NR})

setGeneric("W",function(this) standardGeneric("W"))
setMethod(f="W",signature=c("GCPM"),definition=function(this){
  this@W})

setGeneric("SD.div",function(this) standardGeneric("SD.div"))
setMethod(f="SD.div",signature=c("GCPM"),definition=function(this){
  this@SD.div})

setGeneric("SD.syst",function(this) standardGeneric("SD.syst"))
setMethod(f="SD.syst",signature=c("GCPM"),definition=function(this){
  this@SD.syst})

setGeneric("SD.analyt",function(this) standardGeneric("SD.analyt"))
setMethod(f="SD.analyt",signature=c("GCPM"),definition=function(this){
  this@SD.analyt})

setGeneric("alpha.max",function(this) standardGeneric("alpha.max"))
setMethod(f="alpha.max",signature=c("GCPM"),definition=function(this){
  this@alpha.max})

setGeneric("PDF",function(this) standardGeneric("PDF"))
setMethod(f="PDF",signature=c("GCPM"),definition=function(this){
  this@PDF})

setGeneric("CDF",function(this) standardGeneric("CDF"))
setMethod(f="CDF",signature=c("GCPM"),definition=function(this){
  this@CDF})

setGeneric("loss",function(this) standardGeneric("loss"))
setMethod(f="loss",signature=c("GCPM"),definition=function(this){
  this@loss})

setGeneric("EC.cont",function(this,alpha) standardGeneric("EC.cont"))
setMethod(f="EC.cont",signature=c("GCPM","numeric"),definition=function(this,alpha){
  temp=matrix(NA,ncol=length(alpha),nrow=this@NC)
  colnames(temp)=paste("EC",alpha)
  for(i in 1:length(alpha))
    temp[,i]=EC.VaR.ES.cont(this,alpha[i],type="EC")[,1]
  return(temp)
})

setGeneric("VaR.cont",function(this,alpha) standardGeneric("VaR.cont"))
setMethod(f="VaR.cont",signature=c("GCPM","numeric"),definition=function(this,alpha){
  temp=matrix(NA,ncol=length(alpha),nrow=this@NC)
  colnames(temp)=paste("VaR",alpha)
  for(i in 1:length(alpha))
    temp[,i]=EC.VaR.ES.cont(this,alpha[i],type="VaR")[,2]
return(temp)})

setGeneric("ES.cont",function(this,alpha) standardGeneric("ES.cont"))
setMethod(f="ES.cont",signature=c("GCPM","numeric"),definition=function(this,alpha){
  temp=matrix(NA,ncol=length(alpha),nrow=this@NC)
  colnames(temp)=paste("ES",alpha)
  for(i in 1:length(alpha))
    temp[,i]=EC.VaR.ES.cont(this,alpha[i],type="ES")[,3]
return(temp)})


setGeneric("EL",function(this) standardGeneric("EL"))
setMethod(f="EL",signature=c("GCPM"),definition=function(this){
  this@EL})

setGeneric("SD",function(this) standardGeneric("SD"))
setMethod(f="SD",signature=c("GCPM"),definition=function(this){
  this@SD})
