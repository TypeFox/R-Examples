## ============== setD ====================================
setD <- function(nmax,SL.p,immunity.p=0,risk.type=2,target.EC.p=10,
                 plot=FALSE,alpha.p=5,beta.p=20,
                 print.result="02.sample size.txt")
{
  SL <- SL.p/100
  target.EC <- target.EC.p/100
  immunity  <- immunity.p/100 
  alpha     <- alpha.p/100
  beta      <- beta.p/100

  if (risk.type == 1) {pATarget <- target.EC}
  if (risk.type == 2) {pATarget <- SL + target.EC}
  if (risk.type == 3) {pATarget <- SL + (1-immunity-SL)*target.EC}

  nmin  <- ceiling(1/pATarget)  
  if (nmax < nmin) stop(paste("nmax muss >",nmin,"sein"))

  res           <- matrix(NA,ncol=1,nrow=8)
  rownames(res) <- c("number.organisms","spontaneous.lethality",
                     "immunity",
                     "delta.to.zero",                    
                     "risk.type","target.EC",
                     "exact.alpha","exact.beta")

  pA <- pATarget
  deltan <- ceiling((nmax-nmin)/10)
  n1 <- nmin
  n2 <- 0
  ok <- FALSE
  while((n1 < nmax) & !ok)
  { 
    n2 <- n1 + deltan - 1
    nseq    <- seq(n1,n2, by=1)
    param <- matrix(c(rep(pA,times=length(nseq)),
                    nseq),
                    byrow=TRUE,nrow=2) 

    plan    <- apply(param,MARGIN=2,FUN=OptNBinoA,        
                     p0=SL,alpha=alpha,beta=beta,nmin=nmin)
    plan <- as.matrix(plan[ ,plan["beta ok?", ] == 1])       
    ok   <- length(plan) > 0 #
    n1   <- n2 + 1
  }

  if (length(plan) > 0)
  {
     res["number.organisms",1]      <- plan["n",1]
     res["spontaneous.lethality",1] <- plan["p0",1]*100 
     res["immunity",1]              <- immunity.p
     res["delta.to.zero",1]         <- plan["pA",1]*100
     res["risk.type",1]             <- risk.type
     res["target.EC",1]             <- target.EC.p
     res["exact.alpha",1]           <- round(plan["alpha eff.",1]*100,digits=2)
     res["exact.beta",1]            <- round(plan["beta eff.",1]*100,digits=2)
   }else
   {
     stop("\n","*****  the chosen nmax is too small for the EC value (target.EC) *****\n",
          "***************  please raise nmax or target.EC *********************\n") 
   }


# --- results ----------
 if (print.result!=FALSE)
 {
    sink(print.result,append = FALSE)

    cat("[Part 2: Planning sample size]","\n","\n")
    cat("  [2]",format(Sys.time(),"date: %a., %d.%b.%Y / time: %H:%M:%S"),"\n","\n")
    cat("   ********************************************************","\n",
        "  ** sample size in each treatment per substance and run *","\n",
        "  ********************************************************","\n","\n")
    cat("   number of organisms (n)       = ",res["number.organisms",1],"\n",
        "  risk type                     = ",risk.type,"\n",
        "  verifiable difference (delta) = ",target.EC.p,"%","\n",
        "  immunity                      = ",immunity,"%","\n",
        "  spontaneous lethality         = ",res["spontaneous.lethality",1],"%","\n",
        "  effective alpha               = ",res["exact.alpha",1],"%","\n",
        "  effective beta                = ",res["exact.beta",1],"%","\n","\n",
        "  ********************************************************","\n")
    sink()
    cat(" ---- results are written to '",print.result,"' ---- ","\n")   
 }


# --- plots -----
 if(plot=="single" | plot=="all")
 {  maxk <- res["number.organisms",1]
                       
    dev.new()
    par(las=1,xpd=FALSE)
    n.target  <- maxk    
    nseq.plot <- 1:n.target

    cdf.H0 <- pbinom(nseq.plot,n.target,SL)
    cdf.HA <- pbinom(nseq.plot,n.target,pA)

    plot(nseq.plot ,cdf.H0,
         type="s",col="white",lwd=2,
         xlab=paste("number of dead individuals (of",n.target,"individuals)"),
         ylab="probability of lethality",
         ylim=c(0,1), main="spontaneous vs. treated group")

    abline(h=0)
    abline(h=1)

    axis(2,at=seq(0.1,0.9,by=0.2),labels=FALSE)
    rect(0,(1-alpha),maxk,1,col="mistyrose",density=30,
         border="transparent")
    rect(0,(beta),maxk,0,col="mistyrose",density=30,
         border="transparent")
    lines(nseq.plot ,cdf.H0,type="s",col="gray60",lwd=2)
    lines(nseq.plot ,cdf.HA,type="s",col="gray10",lwd=2)
    lines(c(0,maxk),c(1-alpha,1-alpha),col="lightpink2",lty=5)
    lines(c(0,maxk),c(beta,beta),col="lightpink2",lty=5)
    lines(c(plan["k krit",1],plan["k krit",1]),c(0,1),col="green3",lty=1,lwd=1)

    text(maxk-(maxk/4.5),1-alpha-0.025,col="lightpink2",
         labels=paste("limit of alpha  = ",100*alpha,"%"),
         pos=4,offset=0,cex=0.75)
    text(maxk-(maxk/5),beta+0.025,col="lightpink2",labels=
         paste("limit of beta = ",100*beta,"%"),pos=4,offset=0,cex=0.75)
    text(maxk-(maxk/3.5),0.6,col="gray60",
         labels=paste("-- control group\n","   (spontaneous lethality)"),
         pos=4,offset=0,cex=0.75)
    text(maxk-(maxk/3.5),0.55,col="gray10",labels=paste("-- treated group"),
         pos=4,offset=0,cex=0.75)

    text(plan["k krit",1]+(maxk/15),0.3,col="green4",labels=paste("alpha = ",
         res["exact.alpha",1],"%"),pos=4,cex=0.75)
    text(plan["k krit",1]+(maxk/15),0.25,col="green4",labels=paste("beta = ",
         res["exact.beta",1],"%"),pos=4,cex=0.75)

 }

  
 if(plot=="all")
 {
    roundSL.all <- round(SL, digits=2)
    pA.all      <- seq(roundSL.all,1, by=0.01)
    nseq.all    <- seq(nmin,nmax, by=1)
    param.all <- matrix(c(rep(pA.all,
                    each=length(nseq.all)),
                    rep(nseq.all,times=length(pA.all))),
                    byrow=TRUE,nrow=2) 
    plan.all  <- apply(param.all,MARGIN=2,FUN=OptNBinoA,       
                       p0=SL,alpha=alpha,beta=beta,nmin=nmin)
    plan.all <- plan.all[ ,plan.all["beta ok?", ] == 1]
    plan.all <- unique(plan.all, MARGIN=2)
    xplotmax <- max(plan.all[2,])+5  
    lborder  <- -xplotmax/10
    plan.b <- cbind(c(rev(plan.all["n",]),plan.all["n",1]+5),
                    c(rev(plan.all["pA",]),plan.all["pA",1]-0.01))
    colnames(plan.b) <- c("n","pA")

    dev.new()
    par(las=1,bg="white")
    plot(c(lborder,-1),c(SL,SL),type="n", 
         xlim=c(lborder,xplotmax),ylim=c(0,1),xlab="number of individuals",
         ylab="verifiable response [%]",yaxt="n",main="sample size")
    rg <- seq(1,(nrow(plan.b)-1),by=1)
    for (i in rg)
    { rect(plan.b[i,"n"],plan.b[i,"pA"],plan.b[i+1,"n"],1,
           col="limegreen",border="limegreen")
      rect((plan.b[i,"n"]),0,(plan.b[i+1,"n"]),
           plan.b[i,"pA"],col="red",border="red")
    }
    lines(c(lborder,xplotmax),c(pATarget,pATarget),lty=2,col="gray20")
    lines(c(res[1,1],res[1,1]),c(0,(pATarget+0.1)),lty=2,col="gray20")
    rect(lborder,0,xplotmax,SL, col="gray50", 
         density=30, border="gray60")
    text(lborder,(SL+0.02),labels="spontaneous lethality",col="gray40",
         cex=1.01,pos=4)
    text(lborder,(pATarget+0.02),labels="target.EC",col="gray20",
         cex=1.01,pos=4)
    axis(1,seq(0,xplotmax,by=1),labels=FALSE,tck=-0.01,col.ticks="gray70")
    axis(1,seq(0,xplotmax,by=10),labels=FALSE,col.ticks="black")
    axis(2,axTicks(2),labels=(axTicks(2)*100))
    axis(2,seq(0,1,by=0.01),labels=FALSE,tck=-0.01,col.ticks="gray70")
    axis(2,seq(0,1,by=0.1),labels=FALSE,col.ticks="black")
    legend((xplotmax/2),0.99,c("no identifiable differences",
           "differences are identifiable"),pch=22,col=c("black","black"),
           pt.bg=c("red","limegreen"), bty="n")
    res <- plan.all
 }

return(res)
                     
}

