MirSEA <- function(
input.ds, 
input.cls, 
p_value,
p2miR,
reshuffling.type = "miR.labels", 
nperm = 1000, 
weighted.score.type = 1, 
ms.size.threshold.min = 10, 
ms.size.threshold.max = 500) {


  # Read input data matrix
      content <- input.ds
      content <- content[-(1:2)]
      col.names <- noquote(unlist(strsplit(content[1], "\t")))
      col.names <- col.names[c(-1, -2)]
      num.cols <- length(col.names)
      content <- content[-1]
      num.lines <- length(content)


      row.nam <- vector(length=num.lines, mode="character")
      row.des <- vector(length=num.lines, mode="character")
      m <- matrix(0, nrow=num.lines, ncol=num.cols)

      for (i in 1:num.lines) {
         line.list <- noquote(unlist(strsplit(content[i], "\t")))
         row.nam[i] <- noquote(line.list[1])
         row.des[i] <- noquote(line.list[2])
         line.list <- line.list[c(-1, -2)]
         for (j in 1:length(line.list)) {
            m[i, j] <- as.numeric(line.list[j])
         }
      }
      dataset <- data.frame(m)
      names(dataset) <- col.names
      row.names(dataset) <- row.nam
     
  
  miR.labels <- row.names(dataset)
  sample.names <- names(dataset)
  A <- data.matrix(dataset)
  dim(A) 
  cols <- length(A[1,])
  rows <- length(A[,1])

##################################################
# Read input class.labels
  cls.cont <- input.cls
  num.lines <- length(cls.cont)
      class.list <- unlist(strsplit(cls.cont[[3]], " "))
      s <- length(class.list)
      t <- table(class.list)
      l <- length(t)
      class.phen <- vector(length=l, mode="character")
      phen.label <- vector(length=l, mode="numeric")
      class.labels <- vector(length=s, mode="numeric")
      for (i in 1:l) {
         class.phen[i] <- noquote(names(t)[i])
         phen.label[i] <- i - 1
      }
      for (i in 1:s) {
         for (j in 1:l) {
             if (class.list[i] == class.phen[j]) {
                class.labels[i] <- phen.label[j]
             }
         }
      }
  

 col.index <- order(class.labels, decreasing=FALSE)
 class.labels <- class.labels[col.index]
 sample.names <- sample.names[col.index]
 for (j in 1:rows) {
    A[j, ] <- A[j, col.index]
 }
 names(A) <- sample.names
###################################################
#creat p_value weighting matrix and p2miR profile
 w=1-p_value
  miR.names <- colnames(p_value)
  pathwaynames <- rownames(p_value)
      temp <- p2miR
     
      max.Nm <-nrow(temp)
      temp.size.M <- vector(length = max.Nm, mode = "numeric") 
      for (i in 1:max.Nm) {
          temp.size.M[i] <- length(temp[i,which(temp[i,]!="")]) - 2
      }

      max.size.M <- max(temp.size.M)      
      ms <- matrix(rep("null", max.Nm*max.size.M), nrow=max.Nm, ncol= max.size.M)
      temp.names <- vector(length = max.Nm, mode = "character")
      temp.desc <- vector(length = max.Nm, mode = "character")
      ms.count <- 1
      for (i in 1:max.Nm) {
          miR.set.size <-  temp.size.M[i] 
          ms.line <- temp[i,which(temp[i,]!="")]         
          miR.set.name <- ms.line[1] 
          miR.set.desc <- ms.line[2] 
          miR.set.tags <- ms.line[-(1:2)]
          existing.set <- is.element(miR.set.tags, miR.labels)
          set.size <- length(existing.set[existing.set == TRUE])
          if ((set.size < ms.size.threshold.min) || (set.size > ms.size.threshold.max)) next
          temp.size.M[ms.count] <- set.size
          ms[ms.count,] <- c(miR.set.tags[existing.set], rep("null", max.size.M - temp.size.M[ms.count])) 
          temp.names[ms.count] <- miR.set.name
          temp.desc[ms.count] <- miR.set.desc
          ms.count <- ms.count + 1
      } 
      Nm <- ms.count - 1
      ms.names <- vector(length = Nm, mode = "character")
      ms.desc <- vector(length = Nm, mode = "character")
      size.M <- vector(length = Nm, mode = "numeric") 
      ms.names <- temp.names[1:Nm]
      ms.desc <- temp.desc[1:Nm] 
      size.M <- temp.size.M[1:Nm]

  N <- length(A[,1])


 

  all.ms.descs <- vector(length = Nm, mode ="character") 
 
  
  for (i in 1:Nm) {
        all.ms.descs[i] <- ms.desc[i]
     }
  
  Obs.indicator <- matrix(nrow= Nm, ncol=N)
  
  Obs.ES <- vector(length = Nm, mode = "numeric")
  Obs.arg.ES <- vector(length = Nm, mode = "numeric")
  Obs.ES.norm <- vector(length = Nm, mode = "numeric")
  obs.s2n <- vector(length=N, mode="numeric")
  correl.s2n <- vector(length=N, mode="numeric")
  obs.correl.s2n <- vector(length=N, mode="numeric")
  obs.s2n.matr <- matrix(nrow = Nm, ncol = N)
  signal.strength <- vector(length=Nm, mode="numeric")
  tag.frac <- vector(length=Nm, mode="numeric")
  miR.frac <- vector(length=Nm, mode="numeric")
  
  correl.matrix <- matrix(nrow = N, ncol = nperm)
  obs.correl.matrix <- matrix(nrow = N, ncol = nperm)

   nperm.per.call <- 100
   n.groups <- nperm %/% nperm.per.call
   n.rem <- nperm %% nperm.per.call
   n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
   n.ends <- cumsum(n.perms)
   n.starts <- n.ends - n.perms + 1

   if (n.rem == 0) {
     n.tot <- n.groups
   } else {
     n.tot <- n.groups + 1
   }

 for (nk in 1:n.tot) {
   call.nperm <- n.perms[nk]
   O <- S2N(A, class.labels, miR.labels, call.nperm)
   gc()

   
   correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
   obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
    rm(O)
 }


 for (i in 1:Nm) {
       obs.s2n <- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
	   miR.set <- ms[i,ms[i,] != "null"]
	   miR.set2 <- vector(length=length(miR.set), mode = "numeric")
       miR.set2 <- match(miR.set, miR.labels)
	   msnam <- match(ms.names[i],pathwaynames)
	   for(s in 1:length(miR.set)){
         ww<-w[msnam,match(miR.set[s],miR.names)]
		 obs.s2n[miR.set2[s]]<-obs.s2n[miR.set2[s]]*(1+ww)
		 }
	   obs.index <- order(obs.s2n, decreasing=TRUE) 
       miR.list2 <- obs.index
	   obs.s2n <- sort(obs.s2n, decreasing=TRUE)
		obs.s2n.matr[i,]<-obs.s2n
	
	   MirSEA.results <- EnrichmentScore(miR.list=miR.list2, miR.set=miR.set2, weighted.score.type=weighted.score.type, correl.vector = obs.s2n.matr[i,])
      
       Obs.ES[i] <- MirSEA.results$ES
       Obs.arg.ES[i] <- MirSEA.results$arg.ES
   
       Obs.indicator[i,] <- MirSEA.results$indicator
       if (Obs.ES[i] >= 0) {  # compute signal strength
           tag.frac[i] <- sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.M[i]
           miR.frac[i] <- Obs.arg.ES[i]/N
       } else {
           tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.M[i]
           miR.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
       }
       signal.strength[i] <- tag.frac[i] * (1 - miR.frac[i]) * (N / (N - size.M[i]))
   }
     
    phi <- matrix(nrow = Nm, ncol = nperm)
   phi.norm <- matrix(nrow = Nm, ncol = nperm)
   

    if (reshuffling.type == "sample.labels") {
for (i in 1:Nm) {
        miR.set <- ms[i,ms[i,] != "null"]
        miR.set2 <- vector(length=length(miR.set), mode = "numeric")
        miR.set2 <- match(miR.set, miR.labels)
		msnam <- match(ms.names[i],pathwaynames)
		for (r in 1:nperm) { 
             correl.s2n <- correl.matrix[,r]                     
		     for(s in 1:length(miR.set)){  
                 ww<-w[msnam,match(miR.set[s],miR.names)]
				 correl.s2n[miR.set2[s]]<-correl.s2n[miR.set2[s]]*(1+ww)
			}
		    miR.list2 <- order(correl.s2n,decreasing=TRUE)
			correl.s2n<-sort(correl.s2n,decreasing=TRUE)	
			
			MirSEA.results <- EnrichmentScore2(miR.list=miR.list2, miR.set=miR.set2, weighted.score.type=weighted.score.type, correl.vector=correl.s2n)   
           
            phi[i, r] <- MirSEA.results$ES
        }
       
        gc()
     }
	}else if (reshuffling.type == "miR.labels") {
     for (i in 1:Nm) {
        miR.set <- ms[i,ms[i,] != "null"]
        miR.set2 <- vector(length=length(miR.set), mode = "numeric")
        miR.set2 <- match(miR.set, miR.labels)
		msnam <- match(ms.names[i],pathwaynames)
	for (r in 1:nperm) {
		    obs.s2n <- obs.correl.matrix[,1]
            reshuffled.miR.labels <- sample(1:rows)
			 for(s in 1:length(miR.set)){
			  ww<-w[msnam,match(miR.set[s],miR.names)]
			    xx<-match(miR.set2[s],reshuffled.miR.labels)
				obs.s2n[xx]<-obs.s2n[xx]*(1+ww)
		    }
			obs.s2n<-sort(obs.s2n,decreasing=TRUE)
            
            MirSEA.results <- EnrichmentScore2(miR.list=reshuffled.miR.labels, miR.set=miR.set2, weighted.score.type=weighted.score.type, correl.vector=obs.s2n)   
            
            phi[i, r] <- MirSEA.results$ES
			
        }
         
     }
 }
    
    p.vals <- matrix(0, nrow = Nm, ncol = 2)
    for (i in 1:Nm) {
      if (Obs.ES[i] >= 0) {
         p.vals[i, 1] <-  sum(phi[i,] >= Obs.ES[i])/length(phi[i,])
         p.vals[i, 1] <-  signif(p.vals[i, 1], digits=5)
      } else {
         p.vals[i, 1] <-  sum(phi[i,] <= Obs.ES[i])/length(phi[i,])
         p.vals[i, 1] <-  signif(p.vals[i, 1], digits=5)
      }
   }
   
   
  p.vals[,2]<-p.adjust(p.vals[,1],"fdr")
 
 
 for (i in 1:Nm) {
         pos.phi <- NULL
         neg.phi <- NULL
         for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
               pos.phi <- c(pos.phi, phi[i, j]) 
            } else {
               neg.phi <- c(neg.phi, phi[i, j]) 
            }
         }
         pos.m <- mean(pos.phi)
         neg.m <- mean(abs(neg.phi))
		 pos.phi <- pos.phi/pos.m
         neg.phi <- neg.phi/neg.m
         for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                phi.norm[i, j] <- phi[i, j]/pos.m
            } else {
                phi.norm[i, j] <- phi[i, j]/neg.m
            }
          }
      
          if (Obs.ES[i] >= 0) {
             Obs.ES.norm[i] <- Obs.ES[i]/pos.m
          } else {
             Obs.ES.norm[i] <- Obs.ES[i]/neg.m
          }
   }
   
 
 # Produce results report

       Obs.ES <- signif(Obs.ES, digits=5)
       Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
       p.vals <- signif(p.vals, digits=4)
       signal.strength <- signif(signal.strength, digits=3)
       tag.frac <- signif(tag.frac, digits=3)
       miR.frac <- signif(miR.frac, digits=3)
       

       report <- data.frame(cbind(ms.names, size.M, all.ms.descs, Obs.ES, Obs.ES.norm, p.vals[,1], p.vals[,2], tag.frac, miR.frac, signal.strength))
       names(report) <- c("Pathway", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", "FDR q-val",  "Tag \\%", "Mir \\%", "Signal")

       report2 <- report
       report.index2 <- order(Obs.ES.norm, decreasing=TRUE)
       for (i in 1:Nm) {
           report2[i,] <- report[report.index2[i],]
       }   
       report3 <- report
       report.index3 <- order(Obs.ES.norm, decreasing=FALSE)
       for (i in 1:Nm) {
           report3[i,] <- report[report.index3[i],]
       }   
       phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
       phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
       report.phen1 <- report2[1:phen1.rows,]
       report.phen2 <- report3[1:phen2.rows,]
	
	result<-list(report.phen1 = report.phen1,report.phen2 = report.phen2)
	return(result)
}