##  ============== doseD ==========================
doseD <- function(DP,immunity.p=0,SL.p=0,target.EC.p=10,nconc=8,
                  text=TRUE,risk.type=2,print.result="03.dosestrategy.txt")
{

  logit <- function(p,pNull=1.e-8)
  {  
    if (pNull <= p &  p <= (1-pNull) ) logit <- log(p/(1-p))
    if (p < pNull)                     logit <- log(pNull/(1-pNull)) 
    if (p > (1-pNull))                 logit <- log((1-pNull)/pNull)
    return(logit)
  }
 

 dNull <- 0.1            
 pNull <- 1.e-12         
 n        <- DP[ ,2]
 r        <- DP[ ,3]
 d        <- DP[ ,4]
 unit     <- DP[1,5]
 immunity <- immunity.p/100
 SL       <- SL.p/100

 if(SL==0)
 { SL.ini <- sum(r[d==0])/sum(n[d==0])}else { SL.ini <- SL}

  I <- nconc	
 delta.ini    <- c(0.10,0.50,0.90) 
 label.delta  <- expression("EC"[10],"EC"[50],"EC"[90])
 target.ini   <- sort(target.EC.p/100)

 riskname <- c("total","added risk","extra risk")
 name.risk <- riskname[risk.type]

 if(risk.type==1)
 { label.delta <- label.delta[which(1-immunity > delta.ini &
                              delta.ini > SL.ini)]
   delta.ini   <- delta.ini[delta.ini > SL.ini & 
                            delta.ini < 1-immunity]
   delta       <- delta.ini 
   ok <- ((0 < target.ini) & (target.ini <= (1 - immunity)))
   if(any(!ok) )
   { cat("\n","###  invalid target.EC [%] found  ###\n")
     print(100*target.ini[!ok])
     cat("\n")
     target.ini  <- target.ini[ok]
   }
   target.tot  <- target.ini               
   target      <- target.tot               
   error       <- (length(target) == 0)
 }
 
 if(risk.type==2)
 { label.delta <- label.delta[which((delta.ini + SL.ini) < (1-immunity))]
   delta.ini   <- delta.ini[(delta.ini + SL.ini) < (1-immunity)]
   delta       <- SL.ini+delta.ini 
   ok <- ((0 < target.ini) & ((target.ini+SL.ini) <= (1-immunity)))
   if(any(!ok) )
   { cat("\n","###  invalid target.EC [%] found  ###\n")
     print(100*target.ini[!ok])
     cat("\n")
     target.ini <- target.ini[ok]
   }
   target.tot  <- SL.ini+target.ini        
   target      <- target.tot               
   error       <- (length(target.tot) == 0)
 }

 if(risk.type==3)
 { delta    <- SL.ini+(1-SL.ini-immunity)*delta.ini 
   ok <- ((0 < target.ini) & (target.ini < 1))
   if(any(!ok) )
   { cat("\n","###  invalid target.EC [%] found  ###\n")
     print(100*target.ini[!ok])
     cat("\n")
     target.ini <- target.ini[ok]
   }
   target.tot  <- SL.ini+(1-SL.ini-immunity)*target.ini 
   target      <- target.tot               
   error       <- (length(target.tot) == 0)
 }

 if (error) {cat("\n###  NO valid target values entered  ###\n\n") }

 tau  <- target
 lt   <- length(target)
 ld   <- length(delta)
 l.delta      <- matrix(delta,ncol=1)
 l.target     <- matrix(target,ncol=1)
 EC.target    <- rep(NA,times=lt)
 target.index <- rep(NA,times=lt)
 TCIO         <- rep(NA,times=lt)
 TCIU         <- rep(NA,times=lt)

 # ---  technical settings -----------

 dplot       <- d
 dplot[d==0] <- dNull
 logx        <- log(dplot)
 col.delta <- c("aquamarine2","cadetblue2","cyan3","blue3","blue4")

 logxgrid <- seq(log(dNull),max(c(logx,log(1.1*max(DP[ ,4])))),length=100)
 logxpred <- data.frame(cbind(rep(1,times=length(logxgrid)),logxgrid))
 colnames(logxpred) <- c("(Intercept)","logx")
 dpred <- data.frame(cbind(rep(1,times=length(logxpred[ ,2])),
                     exp(logxpred[ ,2])))
 colnames(dpred) <- c("(Intercept)","d")

 p.emp    <- r/n
 p.ini    <- (p.emp - SL.ini)/(1-SL.ini-immunity)
 okp      <- (0 < p.ini) & (p.ini < 1)
 if(length(p.ini[okp])> 2)
 { p.ini    <- p.ini[okp]
   beta.ini <- lm(-log((1-p.ini)/(p.ini)) ~ logx[okp])
 } else
 { logxpoint <- logx
   p.up   <- (p.ini>1)| (p.ini==1)
   logxpoint1 <- min(logx[p.up])
   p.down <- (p.ini<0)| (p.ini==0)
   logxpoint2 <- max(logx[p.down])
   if(logxpoint2 ==0) { logxpoint2 <- mean(logx[p.down])}
   if(!is.character(logxpoint1) & !is.character(logxpoint1))
   {  stop("----------------------------------------------------------\n",
           "   not enough observations available for further estimations\n",
           "       please repeat the pretest\n",
           "  ----------------------------------------------------------\n\n")
   }else
   { lmma <- data.frame(cbind(rbind(p1=0.999, p2=0.001),
                        rbind(logx1=logxpoint1, logx2=logxpoint2)))
     beta.ini <- lm(-log((1-lmma[,1])/(lmma[,1])) ~ lmma[,2])
   }
 }

 theta.ini <- c(beta.ini$coefficients)
 p.ini <- DW2s(theta.ini,dpred[ ,"d"],SL.ini,immunity,pNull=pNull)
 
 M2Sn <- EstM2Sn(theta.ini,n,r,d,dNull,dpred,logxpred,tau=tau,
                 SL=SL.ini,immunity=immunity,pNull=pNull)

 dgrid    <- M2Sn$dgrid
 logdgrid <- M2Sn$logdgrid

 dev.new()
 par(las=1)
 logxplotmin <- dNull
 logxplotmax <- max(dgrid[ ,"d"])
 yplotmin <- 0
 yplotmax <- 1
 xlabel <- paste("concentration  [",unit,"]") 

 plot(dplot,p.emp,pch=16,type="p",log="x",col="black",
      xlim=c(logxplotmin,logxplotmax),ylim=c(yplotmin,yplotmax),
      xlab=xlabel,ylab="effect [%]",xaxt="n",yaxt="n"
      )
 old <- axTicks(1)
 new <- log10(old)
 axis(1,at=old,labels=10^new)
 axis(2,at=axTicks(2),labels=100*axTicks(2))
 lines(dgrid$d,M2Sn$predp.hut, col="blue",lwd=2)
 lines(dgrid[ ,"d"],M2Sn$predp.ciu,type="l",col="blue",lty=2)
 lines(dgrid[ ,"d"],M2Sn$predp.cio,type="l",col="blue",lty=2)

 if(risk.type==1)
 { if (text)
   { text((logxplotmax-logxplotmax/1.5),0.4,
          labels=paste("target effects:"),col= "green3")
     text((logxplotmax-logxplotmax/1.5),0.36,
          labels=paste("(total risk)"), col= "green3",cex=0.9)
   }
   for (i in 1:lt)
   { lines(exp(c(M2Sn$logxi.hut[i],M2Sn$logxi.hut[i])),
          c(0,(target[i])+0.2),col="green3",lty=1,lwd=1)
     if(text)
     { EC.text <-  paste(100*target.ini[i],"%\n")
       text((logxplotmax-logxplotmax/1.35),(0.34-i/25),
            labels=EC.text, col= "green3",cex=0.9)
     }
   }
 } 

 if(risk.type==2)
 { if (text)
   {text((logxplotmax-logxplotmax/1.5),0.4,
          labels=paste("target effects:"),col= "mediumorchid")
     text((logxplotmax-logxplotmax/1.5),0.36,
          labels=paste("(added risk)"), col= "mediumorchid",cex=0.9)
   }
   for (i in 1:lt)
   { lines(exp(c(M2Sn$logxi.hut[i],M2Sn$logxi.hut[i])),
          c(0,(target[i])+0.2),col="mediumorchid",lty=1,lwd=1)
     if(text)
     { EC.text <-  paste(100*target.ini[i],"%\n")
       text((logxplotmax-logxplotmax/1.35),(0.34-i/25),
            labels=EC.text, col= "mediumorchid",cex=0.9)
     }
   }
 } 

 if(risk.type==3)
 { if (text)
   { text((logxplotmax-logxplotmax/1.5),0.4,
          labels=paste("target effects:"),col= "red2")
     text((logxplotmax-logxplotmax/1.5),0.36,
          labels=paste("(extra risk)"), col= "red2",cex=0.9)
   }
   for (i in 1:lt)
   {lines(exp(c(M2Sn$logxi.hut[i],M2Sn$logxi.hut[i])),
          c(0,(target[i])+0.2),col="red2",lty=1,lwd=1)
    if(text)
    {EC.text <-  paste(100*target.ini[i],"%\n")
     text((logxplotmax-logxplotmax/1.35),(0.34-i/25),
          labels=EC.text, col= "red2",cex=0.9)
    }
   }
 }

# ---- dose positions ------------------------
 a.hut <- M2Sn$theta.hut[1]
 b.hut <- M2Sn$theta.hut[2]

 abline(h=0,lty=1,col="black") 
 abline(h=delta,lty=3,col="gray60")
 ECdelta   <- exp(-(log((1-immunity-delta)/(delta-SL.ini))+ a.hut)/b.hut)
 text((logxplotmax),(delta-0.015),labels=label.delta, col= "gray50",cex=0.65)
 abline(h=SL,lty=3,col="green4")
 text((logxplotmax-(logxplotmax/log(4))),(SL-0.015),
       labels="spontaneous lethality",col= "green4",cex=0.65)
 legend("topleft", c("observed points","estimated curve","confidence interval"),
         lwd=c(0,2,2), lty=c(-1,1,5),pch=c(16,-1,-1),
         col=c("black","blue","blue"),inset=0.01,bg="white")
 
 ECdelta.ciu <- apply(l.delta,MARGIN=1,FUN=InterpolxA,
                      x=dpred[ ,"d"],y=M2Sn$predp.ciu)
 ECdelta.cio <- apply(l.delta,MARGIN=1,FUN=InterpolxA,
                      x=dpred[ ,"d"],y=M2Sn$predp.cio)

 Design <- data.frame(effect=rep(NA,times=I),concentration=rep(NA,times=I),
                      CIO=rep(NA,times=I),CIU=rep(NA,times=I))
 Design[1,1] <- immunity.p 
 Design[1,1] <- paste(SL.p,"% SL (fixed)")
 Design[1,2] <- 0
 Design[2:(ld+1),1] <- paste(delta.ini*100,"% (fixed)")
 Design[2:(ld+1),2] <- ECdelta
 Design[2:(ld+1),4] <- ECdelta.ciu[1:ld]
 Design[2:(ld+1),3] <- ECdelta.cio[1:ld]
 Design.ini <- Design
 
 EC.target[1:lt] <- exp(-(log( (1-immunity-target[1:lt])/
                         (target[1:lt]-SL.ini))+
                          a.hut)/b.hut)

 M2Sn99 <- EstM2Sn(theta.ini,n,r,d,dNull,dpred,logxpred,tau=tau,
                   SL=SL.ini,immunity=immunity,pNull=pNull,alpha=0.01)

 TCIO95 <- apply(l.target,MARGIN=1,FUN=InterpolxA,x=dpred[ ,"d"],y=M2Sn$predp.cio)
 TCIU95 <- apply(l.target,MARGIN=1,FUN=InterpolxA,x=dpred[ ,"d"],y=M2Sn$predp.ciu)
 TCIO99 <- apply(l.target,MARGIN=1,FUN=InterpolxA,x=dpred[ ,"d"],y=M2Sn99$predp.cio)
 TCIU99 <- apply(l.target,MARGIN=1,FUN=InterpolxA,x=dpred[ ,"d"],y=M2Sn99$predp.ciu)
 difference <- data.frame(EC.target = EC.target[1:lt],
                          FCIO95 = TCIO95[1:lt]-EC.target[1:lt], 
                          FCIU95 = TCIU95[1:lt]-EC.target[1:lt],
                          FCIO99 = TCIO99[1:lt]-EC.target[1:lt], 
                          FCIU99 = TCIU99[1:lt]-EC.target[1:lt]
                          )

 nt     <- nrow(difference)
 nfrei  <- nconc-ld-1     

 q <- nfrei / nt

 if (q-trunc(q) < 1.e-6)
 { ntarget <- rep(q,times=nt)
   Design <- plan(difference,ntarget)
 }  else
 { if (q < 1) 
   {  ntarget          <- rep(0,times=nt)
      ntarget[1:nfrei] <- rep(1,times=nfrei)  
      Design <- plan(difference,ntarget)
   }  else
   {  ntarget          <- rep(floor(q),times=nt)
      nrest <- nfrei - sum(ntarget)
      ntarget <- ntarget + c(rep(1,times=nrest),rep(0,times=nt-nrest))  
      Design <- plan(difference,ntarget)
   }
 }

 ndiv <- length(Design.ini[!is.na(Design.ini[ ,"concentration"]),
                          "concentration"])
 ndi <- nrow(Design.ini)
 Jt    <- ncol(Design) 
 i     <- ndiv   
 for (it in 1:lt)
 {  for (jt in 1:Jt)
    { if (!is.na(Design[it,jt]))
      { i <- i + 1 
        Design.ini[i,"effect"]        <- paste(target.EC.p[it],"% (target)")
        Design.ini[i,"concentration"] <- Design[it,jt]
      }
    }
 }

 new.Design <- Design.ini[order(Design.ini[,2]),1:2]
 new.unit   <- rep(unit,times=nrow(new.Design))
 new.Design <- cbind(new.Design,new.unit)
 new.Design[,2] <- round(new.Design[,2], digits=2)
 dose.scheme <- data.frame("concentration" = new.Design[ ,"concentration"],
                           "unit"          = new.Design[ ,"new.unit"],
                           "effect"        = new.Design[ ,"effect"])

 if (print.result!=FALSE)
 {
    sink(print.result,append = FALSE)
    cat("[Part 3: dose strategy]","\n","\n")
    cat("  [3]",format(Sys.time(),"date: %a., %d.%b.%Y / time: %H:%M:%S"),"\n","\n")
    cat("=========== dose strategy  =====================================",   
        "\n    spontaneous lethality =",round((100*SL.ini),digits=2),"%",
        "\n    immunity              =",immunity.p,"%",
        "\n    risk type             =",name.risk,
        "\n    target effect(s)      =",target*100,"% effect \n","\n"
        )

    if(I>3) 
      { cat(" ***********************************************************\n")
                  print(dose.scheme)
        cat(" ***********************************************************\n")
      }else
      { cat(" ***********************************************************\n",
            "  less points for testing than adequate","\n",
            " ***********************************************************\n")
      }
    sink()

    cat(" ---- results are written to '",print.result,"' ---- ","\n")   
 } 

return(dose.scheme)
}
