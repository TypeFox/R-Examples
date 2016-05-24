`HWTernaryPlot` <- function(X, n=NA, addmarkers=TRUE, newframe=TRUE, hwcurve=TRUE, vbounds=TRUE, mafbounds=FALSE, mafvalue=0.05, axis=0, region=1, vertexlab=colnames(X), alpha = 0.05, vertex.cex = 1, pch = 19, cc = 0.5, markercol = "black", markerbgcol= "black", cex=0.75, axislab ="", verbose=FALSE, markerlab=NULL, markerpos=NULL, mcex=1, connect = FALSE, curvecols=rep("black",5) , signifcolour=TRUE, curtyp = "solid", ssf = "max", pvaluetype = "dost", ...)
  {
# plot a ternary diagram that represents all rows of X as points. Acceptance regions for various tests for HWE
#    can be added.
    if(is.vector(X)) {
      if(length(X)!=3) {
        stop("X must have three elements")
      }
      else {
        X <- matrix(X,ncol=3,dimnames=list(c("1"),names(X)))
      }
    }
    X <- as.matrix(X)
    nr <- nrow(X)
    nc <- ncol(X)
    if(any(X<0)) stop("X must be non-negative")
    if(nc != 3) stop("X must have three columns")
    if(is.na(n)) {
      if((sum(apply(X,1,sum))==nr))  # data are compositions
        stop("argument n (the sample size) should be supplied")
      else { # raw counts
        ssf <- match.fun(ssf)
        vsums <- as.matrix(apply(X,1,sum),ncol=1)
        n <- apply(vsums,2,ssf)
        Xr <- X
        if (nrow(X) == 1) 
            Xcom <- X/sum(X)
        else {
            Xcom <- HWClo(X)
        }
      }
    }
    else {
      if((sum(apply(X,1,sum))==nr)) {  # data are compositions
        Xr <- round(n*X)
        Xcom <- X
      }
      else { # raw counts
        Xr <- X
        if (nrow(X) == 1) 
            Xcom <- X/sum(X)
        else {
            Xcom <- HWClo(X)
        }
      }
    }
    chiquant <- qchisq(1-alpha,1)
    r <- sqrt(chiquant/n)
    k <- cc/n
    M <- matrix(c(-1/sqrt(3),0,0,1,1/sqrt(3),0),ncol=2,byrow=T)
    nsignif <- NA
    
    markerq <- (Xcom[,2]+2*Xcom[,3])/2

      if(newframe) {
       opar <- par(pty="m",xpd=TRUE)
       on.exit(par(opar))

       plot(M[,1],M[,2], type="n", axes=FALSE, xlab="", ylab="", pch=19, asp=1, cex.main=2, ... )
       polygon(M)
       
       eps <- 0.04 * vertex.cex
       Mlab <- M + matrix(c(-eps,0,0,eps,eps,0),ncol=2,byrow=T)
       text(Mlab[,1],Mlab[,2], vertexlab, cex=vertex.cex)

       text(0,-0.1,axislab,cex=vertex.cex)
     }

       if (axis==1) {
          AXA <- rbind(c(0,0.5,0.5),c(1,0,0))
          AXA <- AXA%*%M
          lines(AXA[,1],AXA[,2],...)      
       }

       if (axis==2) {
          AXAB <- rbind(c(0.5,0,0.5),c(0,1,0))
          AXAB <- AXAB%*%M
          lines(AXAB[,1],AXAB[,2],...)
       }

       if (axis==3) {
          AXB <- rbind(c(0.5,0.5,0),c(0,0,1))
          AXB <- AXB%*%M
          lines(AXB[,1],AXB[,2],...)
       }


    if(hwcurve) {
         p <- seq(0,1,by=0.005)
         HW <- cbind(p^2,2*p*(1-p),(1-p)^2)
         HWc <- HW%*%M
         points(HWc[,1],HWc[,2],type="l",col=curvecols[1])
       }

       minp  <- sqrt(5/n)
       minpt <- 2*(minp-0.5)/sqrt(3)
       maxp  <- 1-sqrt(5/n)
       maxpt <- 2*(maxp-0.5)/sqrt(3)

       inrange <- sum((markerq >= minp & markerq <= maxp))
       percinrange <- round(100*inrange/nr,digits=2)

       ind1 <- markerq==1
       ind0 <- markerq==0
       nfixed <- sum(ind1) + sum(ind0)
    
       D <- 0.5*(Xcom[,2] - 2*(1-markerq)*markerq)
       Dpos <- sum(Xcom[,2] > 2*(1-markerq)*markerq)
       Dneg <- sum(Xcom[,2] < 2*(1-markerq)*markerq)
       Dzer <- sum(Xcom[,2] == 2*(1-markerq)*markerq)
       Dtot <- Dpos+Dneg


#       cat("D>0:",Dpos,"D<0:",Dneg,"D=0:",Dzer,"nfix:",nfixed,"Tot:",Dtot,"mD:",mean(D),"medD:",median(D),"\n")
#       cat("D>0:",round(100*Dpos/Dtot,digits=5),"D<0:",round(100*Dneg/Dtot,digits=5),"D=0:",
#           round(100*Dzer/Dtot,digits=5),"\n")
       
       if(vbounds) {
         if (n >= 20) {
            lines(c(minpt,minpt),c(0,2*minp),lty="dashed")
            lines(c(maxpt,maxpt),c(0,2-2*maxp),lty="dashed")
          }
       }

       if(mafbounds) {
          minaf <- mafvalue
          minaft <- 2*(minaf-0.5)/sqrt(3)
          maxaf <- 0.95
          maxaft <- 2*(maxaf-0.5)/sqrt(3)
          lines(c(minaft,minaft),c(0,2*minaf),lty="dashed")
          lines(c(maxaft,maxaft),c(0,2-2*maxaf),lty="dashed")
       }    

       if(region==1) { # simple hw ci
          HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
          HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
       }

       if(region==2) { # all curves for hw with cc

         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

         # upper curve for D<0

         DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

         # lower curve for D>0

         DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
      
         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

         if(verbose) {
             cat("D<0 LL",round(DnegLL,digits=2),"\n")
             cat("D<0 UL",round(DnegUL,digits=2),"\n")
             cat("D>0 LL",round(DposLL,digits=2),"\n")
             cat("D>0 UL",round(DposUL,digits=2),"\n")
         }
       }
       if(region==3) { # only limits for D>0, chisq with cc
         # lower curve for D>0

         DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
      
         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

       }
       if(region==4) { # only limits for D<0, chisq with cc
         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

         # upper curve for D<0

         DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
       }

       if(region==5) { # all limits

         HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
         HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])

         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

         # upper curve for D<0

         DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

         # lower curve for D>0

         DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
      
         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
         
       }

       if(region==6) { # lower for D<0, upper for D>0

#         HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
#         HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])

         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
         
         
       }

       if(region==7) { # For Haldane's Exact test

          Crit <- CritSam(n,alphalimit=alpha,pvaluetype=pvaluetype)$Xn
          Critcar <- Crit%*%M        # cartesian coordinates
          points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l",lty=curtyp)

          Crit <- CritSam(n,Dpos=FALSE,alphalimit=alpha,pvaluetype=pvaluetype)$Xn
          Critcar <- Crit%*%M        # cartesian coordinates
          points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l",lty=curtyp)

       }

       if(addmarkers) {
          Xc <- Xcom%*%M # cartesian coordinates

          if (signifcolour==TRUE) { 
  
             if (region == 1) {
                chi.stats <- numeric(nr)           
                chisq.crit <- qchisq(1-alpha,1)
                chisq.stats <- HW.chi.mat(Xr) 
                markerbgcol <- rep("green",nr)             
                markerbgcol[chisq.stats > chisq.crit] <- "red" # look if chisquare is too large.
                markercol <- rep("green",nr)             
                markercol[chisq.stats > chisq.crit] <- "red"
                nsignif <- sum(chisq.stats > chisq.crit)

             }
              
             if (region == 2) {
                pvals <- numeric(nr)
                for (i in 1:nr)
                   pvals[i] <- HWChisq(Xr[i,],cc=0.5,verbose=FALSE)$pval
                markerbgcol <- rep("green",nr)             
                markerbgcol[pvals<alpha] <- "red"
                markercol <- rep("green",nr)             
                markercol[pvals<alpha] <- "red"
                nsignif <- sum(pvals<alpha)

             }

             if (region == 7) {
                pvals <- numeric(nr)
                for (i in 1:nr) {
                    
                   x <- Xr[i,]
                   pvals[i] <- HWExact(Xr[i,],alternative="two.sided",verbose=FALSE,pvaluetype)$pval
                   markerbgcol <- rep("green",nr)             
                   markerbgcol[pvals<alpha] <- "red"
                   markercol <- rep("green",nr)             
                   markercol[pvals<alpha] <- "red"
                   nsignif <- sum(pvals<alpha)
             }

             }


          }

          if (connect)
              points(Xc[,1],Xc[,2],pch=pch,col=curvecols[1],cex=cex,type="l",lty=curtyp)
          else
              points(Xc[,1],Xc[,2],pch=pch,bg=markerbgcol,col=markercol,cex=cex)
          text(Xc[,1],Xc[,2],markerlab,cex=mcex,pos=markerpos)
       }

       return(list(minp=minp,maxp=maxp,inrange=inrange,percinrange=percinrange,nsignif=nsignif))
  }

