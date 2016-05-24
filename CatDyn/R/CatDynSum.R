CatDynSum <-
function(x, season, meths)
   {
    if(class(x) != "list")
      {stop("'x' should be a list with each component an output of function CatDynFit(), class 'catdyn'")}
    if(sum(sapply(x,class) == "catdyn") != length(x))
      {stop("'x' should be a list with each component an output of function CatDynFit(), class 'catdyn'")}
    if(length(meths) != length(x))
      {stop("One numerical method for each component of 'x' should be provided")}
    if(unique(sapply(1:length(x), function(u) length(x[[u]]$Data$Properties$Fleets$Fleet)) == 1)) #1 Fleet
      {
       fleet.name  <- unique(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet));
       effort.unit <- sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Units);
       maxP        <- abs(max(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Type)))
       res         <- data.frame(matrix(0,nrow=length(x),ncol=12+3*maxP+7))
       res$X1      <- fleet.name
       res$X2      <- effort.unit
       res$X3      <- rep(season,length(x))
       res$X4      <- paste(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Type),"P",sep="")
       res$X5      <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Distr)
       res$X6      <- meths
       res$X7      <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$AIC)
       res$X8      <- sapply(1:length(x), function(u) max(abs(x[[u]]$Model[[meths[u]]]$num.grads)))
       res$X9      <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par["M"]))
       res$X10     <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev["M"])
       res$X11     <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par["N0"])
       res$X12     <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev["N0"])
       if(maxP == 0)
         {
          res$X13 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("k.",fleet.name,sep="")]))
          res$X14 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("k.",fleet.name,sep="")])
          res$X15 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("alpha.",fleet.name,sep="")]))
          res$X16 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("alpha.",fleet.name,sep="")])
          res$X17 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("beta.",fleet.name,sep="")]))
          res$X18 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("beta.",fleet.name,sep="")])
          res$X19 <- 0
          names(res) <- c("Fleet","Effort","Season","Model","Distribution","Method","AIC","Max.Abs.Grads.","M","SE.M","N0","SE.N0",
                          paste(c("k.","SE.k.","alpha.","SE.alpha.","beta.","SE.beta."),fleet.name,sep=""),"Sel.Model")
         }
       else
         {
          for(i in 1:length(x))
            {
             if(abs(x[[i]]$Model[[meths[i]]]$Type!=0))
               {
                for(j in 1:abs(x[[i]]$Model[[meths[i]]]$Type))
                  {
                   res[i,13+3*(j-1)] <- x[[i]]$Model[[meths[i]]]$bt.par[paste("P",as.character(j),".",fleet.name,sep="")]
                   res[i,14+3*(j-1)] <- x[[i]]$Model[[meths[i]]]$bt.stdev[paste("P",as.character(j),".",fleet.name,sep="")]
                   res[i,15+3*(j-1)] <- x[[i]]$Model[[meths[i]]]$Dates[1+j]
                  }
               }
            }
          res[,3*maxP+13] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("k.",fleet.name,sep="")]))
          res[,3*maxP+14] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("k.",fleet.name,sep="")])
          res[,3*maxP+15] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("alpha.",fleet.name,sep="")]))
          res[,3*maxP+16] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("alpha.",fleet.name,sep="")])
          res[,3*maxP+17] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("beta.",fleet.name,sep="")]))
          res[,3*maxP+18] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("beta.",fleet.name,sep="")])
          res[,3*maxP+19] <- 0
          names(res) <- c("Fleet","Effort","Season","Model","Distribution","Method","AIC","Max.Abs.Grads.","M","SE.M","N0","SE.N0",
                          paste(paste(rep(c("P","SE.P","ts.P"),maxP),sort(rep(1:maxP,3)),rep(".",maxP),sep=""),fleet.name,sep=""),
                          paste(c("k.","SE.k.","alpha.","SE.alpha.","beta.","SE.beta."),fleet.name,sep=""),
                          "Sel.Model")
         } 
      }   
    else if(anyDuplicated(t(sapply(1:length(x), function(u) unique(x[[u]]$Data$Properties$Fleets$Fleet)))) == 2) #2 Fleets
      {
       maxP        <- max(sapply(1:length(x), function(u) sum(x[[u]]$Model[[meths[u]]]$Type)))
       maxPF1      <- max(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Type)[1,])
       maxPF2      <- max(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Type)[2,])
       res         <- data.frame(matrix(0,nrow=length(x),ncol=12+3*maxP+13))
       for(i in 1:length(x))
         {
          res$X1[i] <- paste(t(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet))[i,1],
                             t(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet))[i,2],sep=",")
          res$X2[i] <- paste(t(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Units))[i,1],
                             t(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Units))[i,2],sep=",")
         }
       fleet.name  <- c(unique(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet)[1,]),
                        unique(sapply(1:length(x), function(u) x[[u]]$Data$Properties$Fleets$Fleet)[2,]))
       res$X3      <- rep(season,length(x))
       for(i in 1:length(x))
         {
          res$X4[i] <- paste(x[[i]]$Model[[meths[i]]]$Type[1],"P", x[[i]]$Model[[meths[i]]]$Type[2],"P",sep="")
          res$X5[i] <- paste(t(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Distr))[i,1],
                             t(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$Distr))[i,2],sep=",")
         }
       res$X6     <- meths
       res$X7     <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$AIC)
       res$X8     <- sapply(1:length(x), function(u) max(abs(x[[u]]$Model[[meths[u]]]$num.grads)))
       res$X9     <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par["M"]))
       res$X10    <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev["M"])
       res$X11    <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par["N0"])
       res$X12    <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev["N0"])
       if(maxP == 0)
         {
          res$X13 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("k.",fleet.name[1],sep="")]))
          res$X14 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("k.",fleet.name[1],sep="")])
          res$X15 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("alpha.",fleet.name[1],sep="")]))
          res$X16 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("alpha.",fleet.name[1],sep="")])
          res$X17 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("beta.",fleet.name[1],sep="")]))
          res$X18 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("beta.",fleet.name[1],sep="")])
          res$X19 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("k.",fleet.name[2],sep="")]))
          res$X20 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("k.",fleet.name[2],sep="")])
          res$X21 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("alpha.",fleet.name[2],sep="")]))
          res$X22 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("alpha.",fleet.name[2],sep="")])
          res$X23 <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("beta.",fleet.name[2],sep="")]))
          res$X24 <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("beta.",fleet.name[2],sep="")])
          res$X25 <- 0
          names(res) <- c("Fleets","Effort","Season","Model","Distribution","Method","AIC","Max.Abs.Grads.","M","SE.M","N0","SE.N0",
                          paste(c("k.","SE.k.","alpha.","SE.alpha.","beta.","SE.beta."),fleet.name[1],sep=""),
                          paste(c("k.","SE.k.","alpha.","SE.alpha.","beta.","SE.beta."),fleet.name[2],sep=""),
                          "Sel.Model")
         }
       else
         {
          for(i in 1:length(x))
            {
             if(x[[i]]$Model[[meths[i]]]$Type[1]!=0) 
               {
                for(j in 1:x[[i]]$Model[[meths[i]]]$Type[1])
                  {
                   res[i,13+3*(j-1)] <- x[[i]]$Model[[meths[i]]]$bt.par[paste("P",as.character(j),".",fleet.name[1],sep="")]
                   res[i,14+3*(j-1)] <- x[[i]]$Model[[meths[i]]]$bt.stdev[paste("P",as.character(j),".",fleet.name[1],sep="")]
                   res[i,15+3*(j-1)] <- x[[i]]$Model[[meths[i]]]$Dates[1+j]
                  }
                } 
             if(x[[i]]$Model[[meths[i]]]$Type[2]!=0) 
               {
                for(k in 1:x[[i]]$Model[[meths[i]]]$Type[2])
                  {
                   res[i,13+3*(maxPF1)+3*(k-1)] <- x[[i]]$Model[[meths[i]]]$bt.par[paste("P",as.character(k),".",fleet.name[2],sep="")]
                   res[i,14+3*(maxPF1)+3*(k-1)] <- x[[i]]$Model[[meths[i]]]$bt.stdev[paste("P",as.character(k),".",fleet.name[2],sep="")]
                   res[i,15+3*(maxPF1)+3*(k-1)] <- x[[i]]$Model[[meths[i]]]$Dates[1+x[[i]]$Model[[meths[i]]]$Type[1]+k]
                  }
               }
            }
          res[,3*maxP+13] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("k.",fleet.name[1],sep="")]))
          res[,3*maxP+14] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("k.",fleet.name[1],sep="")])
          res[,3*maxP+15] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("alpha.",fleet.name[1],sep="")]))
          res[,3*maxP+16] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("alpha.",fleet.name[1],sep="")])
          res[,3*maxP+17] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("beta.",fleet.name[1],sep="")]))
          res[,3*maxP+18] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("beta.",fleet.name[1],sep="")])
          res[,3*maxP+19] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("k.",fleet.name[2],sep="")]))
          res[,3*maxP+20] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("k.",fleet.name[2],sep="")])
          res[,3*maxP+21] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("alpha.",fleet.name[2],sep="")]))
          res[,3*maxP+22] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("alpha.",fleet.name[2],sep="")])
          res[,3*maxP+23] <- unlist(sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.par[paste("beta.",fleet.name[2],sep="")]))
          res[,3*maxP+24] <- sapply(1:length(x), function(u) x[[u]]$Model[[meths[u]]]$bt.stdev[paste("beta.",fleet.name[2],sep="")])
          res[,3*maxP+25] <- 0
          names(res) <- c("Fleets","Effort","Season","Model","Distribution","Method","AIC","Max.Abs.Grads.","M","SE.M","N0","SE.N0",
                          paste(paste(rep(c("P","SE.P","ts.P"),maxPF1),sort(rep(1:maxPF1,3)),rep(".",maxPF1),sep=""),fleet.name[1],sep=""),
                          paste(paste(rep(c("P","SE.P","ts.P"),maxPF2),sort(rep(1:maxPF2,3)),rep(".",maxPF2),sep=""),fleet.name[2],sep=""),
                          paste(c("k.","SE.k.","alpha.","SE.alpha.","beta.","SE.beta."),fleet.name[1],sep=""),
                          paste(c("k.","SE.k.","alpha.","SE.alpha.","beta.","SE.beta."),fleet.name[2],sep=""),
                          "Sel.Model")
         } 
      } 
    else
      {
       stop("All components of 'x' should have been built with data from the same fleets with the number of fleets > 0 and < 3")
      }
    return(res);
   }
