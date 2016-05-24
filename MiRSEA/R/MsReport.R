
###################################################################################      
MsReport<- function(
     MsNAME = "",
     input.ds, 
     input.cls, 
     p_value,
	 p2miR,
     weighted.score.type = 1 
    ) {


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
  cols <- length(A[1,])
  rows <- length(A[,1])

##################################################
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
     
     phen1 <- class.phen[1]
     phen2 <- class.phen[2]

 col.index <- order(class.labels, decreasing=FALSE)
 class.labels <- class.labels[col.index]
 sample.names <- sample.names[col.index]
 for (j in 1:length(A[,1])) {
    A[j, ] <- A[j, col.index]
 }
 names(A) <- sample.names
   N <- length(A[,1])

 #    
 ###################################################
#get p_value weighting matrix and p2miR profile
  w=1-p_value
  miR.names <- colnames(p_value)
  pathwaynames <- rownames(p_value)
      temp <- p2miR
     
	  i<-match(MsNAME,pathwaynames)
	  if(!is.na(i)){
      temp.size.M<-length(temp[i,which(temp[i,]!="")]) - 2
      ms.line <- temp[i,which(temp[i,]!="")] 
      miR.set.name <- ms.line[1]
      miR.set.tags <- ms.line[-(1:2)]
	  existing.set <- is.element(miR.set.tags, miR.labels)
	  size.M <- length(existing.set[existing.set == TRUE])
	  ms <- miR.set.tags[existing.set]



  Obs.indicator <- vector(length = N, mode = "numeric")
  Obs.RES <- vector(length = N, mode = "numeric")
  obs.s2n <- vector(length=N, mode="numeric")

    O <- S2N(A, class.labels, miR.labels,1)
    obs.s2n <- O$obs.s2n.matrix[,1]
    rm(O)
 
	   miR.set <- vector(length=length(ms), mode = "numeric")
       miR.set <- match(ms, miR.labels)
	 
	   for(s in 1:length(ms)){
         ww<-w[i,match(ms[s],miR.names)]
		 obs.s2n[miR.set[s]]<-obs.s2n[miR.set[s]]*(1+ww)
		 }
	   obs.index <- order(obs.s2n, decreasing=TRUE) 
       obs.s2n <- sort(obs.s2n, decreasing=TRUE)
	   obs.miR.labels<-  miR.labels[obs.index]  
	
       PathwayResult <- EnrichmentScore(miR.list=obs.index, miR.set=miR.set, weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
      
       Obs.ES <-  PathwayResult$ES
       Obs.arg.ES <-  PathwayResult$arg.ES
       Obs.RES <-  PathwayResult$RES
       Obs.indicator <-  PathwayResult$indicator

            kk <- 1
            miR.number <- vector(length = size.M, mode = "character")
            miR.names <- vector(length = size.M, mode = "character")
      
            miR.list.loc <- vector(length = size.M, mode = "numeric")
            core.enrichment <- vector(length = size.M, mode = "character")
            miR.s2n <- vector(length = size.M, mode = "numeric")
            miR.RES <- vector(length = size.M, mode = "numeric")
            rank.list <- seq(1, N)
           
            if (Obs.ES >= 0) {
              set.k <- seq(1, N, 1)
              
            } else {
              set.k <- seq(N, 1, -1)
            
            }

            for (k in set.k) {
               if (Obs.indicator[k] == 1) {
                  miR.number[kk] <- kk
                  miR.names[kk] <- obs.miR.labels[k]
             
                  miR.list.loc[kk] <- k
                  miR.s2n[kk] <- signif(obs.s2n[k], digits=3)
                  miR.RES[kk] <- signif(Obs.RES[k], digits = 3)
                  if (Obs.ES >= 0) {
                     core.enrichment[kk] <- ifelse(miR.list.loc[kk] <= Obs.arg.ES, "YES", "NO")
                  } else {
                     core.enrichment[kk] <- ifelse(miR.list.loc[kk] > Obs.arg.ES, "YES", "NO")
                  }
                  kk <- kk + 1
               }
            }

       PathwayReport <- data.frame(cbind(miR.number, miR.names, miR.list.loc, miR.s2n, miR.RES, core.enrichment))
       names( PathwayReport) <- c("#", "MiR","LIST LOC", "TW-SCORE", "RES", "CORE_ENRICHMENT")
     result1<-PathwayReport
	 result2<-list(N,t(Obs.RES),t(obs.s2n),Obs.ES,size.M,Obs.arg.ES,t(Obs.indicator),phen1,phen2,t(obs.index),t(obs.miR.labels),MsNAME)
	 result<-list("Msreport"=result1,"miRList"=result2)
	 return(result)
	}
	
	}