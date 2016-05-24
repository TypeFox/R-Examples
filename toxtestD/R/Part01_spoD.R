##  ============== spoD ==========================
spoD <- function(n=500,SL.p=5,SLmin=NA,SLmax=NA,bio.sd.p=2.008,
                 maxCI = 2.5,analysis=FALSE,SLdataset=NA,
                 print.result="01_spontaneous lethality.txt")
{

##  ============== spoD_A ========================================
spoD_A <- function(SLdataset)
{
   TV    <- nrow(SLdataset)
   SL    <- SLdataset   
   sum.n <- sum(SLdataset[ ,"n"])
   sum.r <- sum(SLdataset[ ,"bearer"])

   p.i      <- SLdataset[ ,"bearer"]/
               SLdataset[ ,"n"]
   mue.p    <- sum(SLdataset[ ,"n"]* p.i)/
               sum.n                      
   mSL      <- round(mue.p*100,digits=2) 
   sig.p    <- sqrt(TV/(TV-1) * 
                   (SLdataset[ ,"n"]/sum.n) %*% ((p.i-mue.p)^2) )
   sdSL     <- round(sig.p*100,digits=2) 

# --- CIs -------------------------
  CI <- ConfA(c(sum.n,sum.r))   
  ciu <- round(CI["CIlow",], digits=2)
  cio <- round(CI["CIup",], digits=2)

if (TV > 1) 
{
  Chi.test <- chisq.test(cbind( SL[ ,"n"]-SL[ ,"bearer"],SL[ ,"bearer"]),
                         simulate.p.value=TRUE) 
  if(Chi.test["p.value"]<= 0.05)
  {  
    res.txt <- paste(" the standard deviation is composed of the normal random variation",
                     " and a time depending biological variation","\n")
  } else
  {  res.txt <- paste(" the standard deviation results only from the normal random variation\n")
  }
}

# --- results ----------------------
  if (print.result!=FALSE)
  {
    sink(file=print.result,append=TRUE)

    cat("[Part 1b: analysis]","\n","\n")
    cat("  [1]",format(Sys.time(),"date: %a., %d.%b.%Y / time: %H:%M:%S"),"\n","\n")
    cat("   --------------------------------------------------------------","\n",
        "  -------          spontaneous lethality               ----------","\n",
        "  --------------------------------------------------------------","\n")
    cat("      mean rate of the spontaneous lethality: ",mSL,"%\n")
    cat("      95%  confidence interval:              ",ciu,"-",cio,"%\n",
        "      standard deviation:                    ",sdSL,"%\n",
        "  --------------------------------------------------------------","\n","\n")
    cat(" ",res.txt,"\n",
        "  --------------------------------------------------------------","\n","\n")
    sink()

    cat(" ---- results are written to '",print.result,"' ---- ","\n")    
  }
  result       <- c(mue.p*100,CI["CIlow",],CI["CIup",],sig.p*100)
 names(result) <- c("SL","CIlo","CIup","sdSL")
 return(result)
}


##  ============== spoD_P ====================================================
spoD_P <- function(n=n,SL.p=NA,SLmin=NA,SLmax=NA,maxCI=2.5,
                   bio.sd.p=bio.sd.p,plot=FALSE)
{ 
 SL.txt    <- paste("01_spontaneous lethality.txt",sep="")
 SL.name   <- c("LfdNr","total.number","death.number","sponLeth",
                "CIlow","CIup","CIwidth","CIso","CIsu")  
 maxn       <- n
 ok <-  (!is.na(SL.p) &  is.na(SLmin) &  is.na(SLmax) )  |  
        ( is.na(SL.p) & !is.na(SLmin) & !is.na(SLmax) )  |
        (!is.na(SL.p) & !is.na(SLmin) & !is.na(SLmax) )   

 if (!ok) 
 { result <- NA
   stop("Warning: incorrect submission of rates (SL.p, SLmin or SLmax) - proceeding is terminated!")
 }
 if(!is.na(SL.p) & is.na(SLmin) & is.na(SLmax) )  
 { targetSL <- SL.p
   targetmin  <- floor(SL.p)
   targetmax  <- ceiling(SL.p)
   if(targetmin == targetmax)
   {  targetmin  <- SL.p-0.5
      targetmax  <- SL.p+0.5
   }  
 } else 
 { if(is.na(SL.p) & !is.na(SLmin) & !is.na(SLmax) )
   { targetSL <- mean(c(SLmin,SLmax))
     targetmin <- SLmin
     targetmax <- SLmax
   } else
   { if(!is.na(SL.p) & !is.na(SLmin) & !is.na(SLmax) )   
     { targetSL <- SL.p
       targetmin  <- SLmin
       targetmax  <- SLmax     
     } else
     { cat("Warning: incorrect submission of rates (SL.p, SLmin or SLmax): ",
            SL.p,SLmin,SLmax,"\n")
     }
   }
 }

##  -------------- pretest ------------------------------------------
 if(( maxn*(targetmax/100) ) < 1 )
 { result <- NA 
   cat("the combination of spontaneous lethality rate and the maximum total number (n)\n", 
        "does not allow any observation (minimum is one deathd organism)\n")
   stop("Warning: proceeding is terminated!")
 }

##  -------------- CIs and F-distribution --------------
 max.rea   <- ceiling(maxn*targetmax/100)  
 min.nmax  <- floor(1/(targetmax/100))
 span.rea  <- maxn-min.nmax-1
 Go <- (targetmax/100)+0.0001
 Gu <- (targetmin/100)-0.0001

 confid <- matrix(c(rep(min.nmax:maxn,each=(max.rea)),rep(1:max.rea,times=span.rea)),
                   byrow=TRUE,nrow=2)
 p.confid <- confid[2, ]/confid[1, ]
 confid <- confid[ ,p.confid <= Go & Gu < p.confid]
   
 Spontan <- apply(confid,MARGIN=2,FUN=ConfA)
 rownames(Spontan) <-SL.name[-1]

 Spontan <-  Spontan[,Spontan["CIso",] < maxCI]
 Spontan <-  Spontan[,Spontan["CIsu",] < maxCI]
 Spontan["sponLeth",] <- round(Spontan["sponLeth",], digits = 2)
 u.SL <- unique(Spontan["sponLeth",])
 l.SL <- length(u.SL)
 if(l.SL == 0)
 {
   result <- NA
   stop("Warning: proceeding is terminated!\n" 
        ,"the chosen nmax is too small\n"
        ,"there are no results within maxCI = ",maxCI,"\n")
 }

 SP   <- matrix(NA,ncol=7,nrow=l.SL)
 colnames(SP) <- SL.name[1:7]
 SP[ ,"LfdNr"] <- seq(1,l.SL,by=1) 

 i <- 1
 for (AS in u.SL)
 { sponLeth.eq.AS <- Spontan["sponLeth",]== AS
   ok.spv      <- min(Spontan["total.number",sponLeth.eq.AS])
   ok.sp.index <- which(Spontan["total.number",] == ok.spv & 
                        sponLeth.eq.AS)
   SP[i,"total.number"]  <- Spontan["total.number",ok.sp.index] 
   SP[i,"death.number"]     <- Spontan["death.number",ok.sp.index]
   SP[i,"sponLeth"]     <- Spontan["sponLeth",ok.sp.index]
   SP[i,"CIlow"]           <- round(Spontan["CIlow",ok.sp.index], digits=3)
   SP[i,"CIup"]            <- round(Spontan["CIup",ok.sp.index], digits=3)
   SP[i,"CIwidth"]         <- round(Spontan["CIwidth",ok.sp.index], digits=3)
   i <- i+1
 }

 delta.rate <- 0.05
 SP.1 <- SP[which( abs(SP[ ,"sponLeth"]-targetSL) < delta.rate ), ]
 SP.2 <- SP.1[which(     abs(SP.1[ ,"sponLeth"]-targetSL) == 
                     min(abs(SP.1[ ,"sponLeth"]-targetSL)) ),"LfdNr"]
 ntarget  <- SP[SP.2,"total.number"]
 SP.3        <- SP[    SP[ ,"total.number"] == 
                   max(SP[ ,"total.number"]),  ]
 SP.3        <- rbind(SP.3,SP.3)
 nSLmax      <- unique( SP.3[ ,"total.number"] )
 maxSL       <- max(SP.3[ ,"sponLeth"])
 bio.sd      <- bio.sd.p/100


##  --------  partition ---------
  TVr.opt <- TV.model(targetSL=targetSL,ntarget=ntarget,
                      sigma.bio.wahr=bio.sd,plot=FALSE)
  TVopt   <- TVr.opt["TV"]
  nTVopt  <- TVr.opt["nTV"]
  TVmax.opt   <- TV.model(targetSL=maxSL,ntarget=nSLmax,
                      sigma.bio.wahr=bio.sd,plot=FALSE)
  optmaxTV    <- TVmax.opt["TV"]
  nTVoptmaxTV <- TVmax.opt["nTV"]


# -----  results -----
  result        <- c(targetSL,ntarget,TVopt,nTVopt,maxSL,
                     targetmin,targetmax,nSLmax,optmaxTV,nTVoptmaxTV)
  names(result) <- c("targetSL","ntarget","optnum","nopt","maxSL",
                     "Intmin","Intmax","maxn","optmax","noptmax")

  if (print.result!=FALSE)
  {
  sink(file=print.result,append=TRUE)

  cat("[Part 1a: planning spontaneous lethality]","\n","\n")
  cat("  [1]",format(Sys.time(),"date: %a., %d.%b.%Y / time: %H:%M:%S"),"\n","\n")
  cat("   --------------------------------------------------------------\n",
      "  ------  Determination of the spontaneous lethality    ---------\n",
      "  --------------------------------------------------------------\n\n")
  cat("      >>>>>>  spontaneous lethality of:        ",targetSL,"%    <<<<<<\n",
      "    total number of test organisms (n): ",ntarget,"\n",
      "    proposed partitioning: ",TVopt, "separated tests","\n",
      "    (with ",nTVopt,"organisms respectively)\n\n\n")

  cat("        >>>>>>  spontaneous lethality between: ",targetmin,"-",targetmax,"%  <<<<<<\n",
      "    total number of test organisms (n) raises to:  ",nSLmax,"\n",
      "    proposed partitioning: ",optmaxTV, "separated tests","\n",
      "    (with ",nTVoptmaxTV,"organisms respectively)\n\n",
      "  --------------------------------------------------------------\n","\n")
  sink()
  cat(" ---- results are written to '",print.result,"' ---- ","\n")    
  }
  return(result)
}

#==============================================================
#==============================================================

 if(analysis==FALSE)
 { result <- spoD_P(n=n,SL.p=SL.p,SLmin=SLmin,SLmax=SLmax,
                    bio.sd.p=bio.sd.p,maxCI=maxCI)
 } 
 if(analysis==TRUE)
 { result <- spoD_A(SLdataset=SLdataset)
 }

 if(any(is.na(result))) 
 { stop("[spontaneous lethality] please choose a higher number of 
        organisms")
 }
 return(result)
}

############################################################################

