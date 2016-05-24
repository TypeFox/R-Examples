################################################################################################
########################################## P A G I #############################################
################################################################################################
# Main PAGI Analysis Function that implements the entire methodology
#library(igraph)

initializePAGI<-function(){
   utils::data("PAGIData",package="PAGI")
  
}





PAGI.Main <- function(
dataset, 
class.labels, 
nperm = 100, 
p.val.threshold=-1,
FDR.threshold = 0.01, 
gs.size.threshold.min = 25, 
gs.size.threshold.max = 500 ) 
{

  if(!exists("PAGIData")) initializePAGI()

 pathway.db<-get("pathway.db",envir=PAGIData)
 netWorkdata<-get("netWorkdata",envir=PAGIData)
  print("Running PAGI Analysis...")

# Start of PAGI methodology 
  gene.labels <- as.character(dataset[,1])
  sample.names <- names(dataset)[2:length(dataset[1,])]
  A <- data.matrix(dataset[2:length(dataset[1,])])
  row.names(A)<-gene.labels
  cols <- length(A[1,])
  rows <- length(A[,1])

# sort samples according to phenotype
  class.labels[which(class.labels==unique(class.labels)[1])]<-rep(0,length(which(class.labels==unique(class.labels)[1])))
  class.labels[which(class.labels==unique(class.labels)[2])]<-rep(1,length(which(class.labels==unique(class.labels)[2])))
  col.index <- order(class.labels, decreasing=F)
  class.labels <- as.character(class.labels[col.index])
  sample.names <- sample.names[col.index]
  for (j in 1:rows) {
    A[j, ] <- A[j, col.index]
  }
    names(A) <- sample.names
 
# Read input pathway database
      temp<-pathway.db
      max.Ng <- length(temp)
      temp.size.G <- vector(length = max.Ng, mode = "numeric") 
      for (i in 1:max.Ng) {
          temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
      }

      max.size.G <- max(temp.size.G)      
      gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
      temp.names <- vector(length = max.Ng, mode = "character")
      temp.desc <- vector(length = max.Ng, mode = "character")
      gs.count <- 1
      for (i in 1:max.Ng) {
          gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
          gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
          gene.set.name <- gs.line[1] 
          gene.set.desc <- gs.line[2] 
          gene.set.tags <- vector(length = gene.set.size, mode = "character")
          for (j in 1:gene.set.size) {
              gene.set.tags[j] <- gs.line[j + 2]
          } 
          existing.set <- is.element(gene.set.tags, gene.labels)
          set.size <- length(existing.set[existing.set == T])
          if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
          temp.size.G[gs.count] <- set.size
          gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
          temp.names[gs.count] <- gene.set.name
          temp.desc[gs.count] <- gene.set.desc
          gs.count <- gs.count + 1
      } 
      Ng <- gs.count - 1
      gs.names <- vector(length = Ng, mode = "character")
      gs.desc <- vector(length = Ng, mode = "character")
      size.G <- vector(length = Ng, mode = "numeric") 
      gs.names <- temp.names[1:Ng]
      gs.desc <- temp.desc[1:Ng] 
      size.G <- temp.size.G[1:Ng]
      N <- length(A[,1])
      Ns <- length(A[1,])
      all.gene.symbols <- vector(length = N, mode ="character") 
      all.gs.descs <- vector(length = Ng, mode ="character")  
      for (i in 1:N) {       
        all.gene.symbols[i] <- gene.labels[i]
      } 
      for (i in 1:Ng) {
        all.gs.descs[i] <- gs.desc[i]
      }  
  Obs.indicator <- matrix(nrow= Ng, ncol=N)
  Obs.RES <- matrix(nrow= Ng, ncol=N)
  Obs.ES <- vector(length = Ng, mode = "numeric")
  Obs.arg.ES <- vector(length = Ng, mode = "numeric")

# PAGI methodology

  obs.Test<- vector(length=N, mode="numeric")
  obs.Weight<- vector(length=N, mode="numeric")
  Weight<- vector(length=N, mode="numeric")
  obs.s2n <- vector(length=N, mode="numeric")
  signal.strength <- vector(length=Ng, mode="numeric")
  tag.frac <- vector(length=Ng, mode="numeric")
  gene.frac <- vector(length=Ng, mode="numeric")
  obs.Test<-abs(calTest(A, class.labels)$T)
  obs.Weight<-CalGIF2(obs.Test, netWorkdata)
  Weight<-as.vector(obs.Weight)+1
  rm(A)
  obs.s2n <- obs.Test^Weight
  obs.index <- order(obs.s2n, decreasing=T)            
  obs.s2n   <- sort(obs.s2n, decreasing=T)            
  obs.gene.labels <- gene.labels[obs.index]        
  obs.gene.symbols <- all.gene.symbols[obs.index]       
  gc()  
  gene.list2 <- obs.index
  for (i in 1:Ng) {
       gene.set <- gs[i,gs[i,] != "null"]
       gene.set2 <- vector(length=length(gene.set), mode = "numeric")
       gene.set2 <- match(gene.set, gene.labels)
       PAGI.results <- PAGI.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2,  correl.vector = obs.s2n)
       Obs.ES[i] <- PAGI.results$ES
       Obs.arg.ES[i] <- PAGI.results$arg.ES
       Obs.RES[i,] <- PAGI.results$RES
       Obs.indicator[i,] <- PAGI.results$indicator
       tag.frac[i] <- sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.G[i]
       gene.frac[i] <- Obs.arg.ES[i]/N
       signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * (N / (N - size.G[i]))
   }
   phi <- matrix(nrow = Ng, ncol = nperm)  
        for (i in 1:Ng) {
        gene.set <- gs[i,gs[i,] != "null"]
        gene.set2 <- vector(length=length(gene.set), mode = "numeric")
        gene.set2 <- match(gene.set, gene.labels)
          for (r in 1:nperm) {
            reshuffled.gene.labels <- sample(1:rows)
            PAGI.results <- PAGI.EnrichmentScore2(gene.list=reshuffled.gene.labels, gene.set=gene.set2, correl.vector=obs.s2n)   
            phi[i, r] <- PAGI.results$ES
          }
        gc()
       }
# Find nominal p-values       
p.vals <- matrix(0, nrow = Ng, ncol = 2)
   for (i in 1:Ng) {
      pos.phi <- NULL      
      for (j in 1:nperm) {        
            pos.phi <- c(pos.phi, phi[i, j])         
      }
      ES.value <- Obs.ES[i]
      p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=7)    
   }
# Compute FDRs 
      p.vals[, 2]<-p.adjust(p.vals[, 1],method="fdr")
# Produce results report
       result<-list()
       Obs.ES <- signif(Obs.ES, digits=5)
       signal.strength <- signif(signal.strength, digits=3)
       tag.frac <- signif(tag.frac, digits=3)
       gene.frac <- signif(gene.frac, digits=3)
       report <- data.frame(cbind(gs.names, size.G, all.gs.descs, Obs.ES,  p.vals[,1], p.vals[,2],  tag.frac, gene.frac, signal.strength))
       names(report) <- c("Pathway Name", "SIZE", "PathwayID", "Pathway Score",  "NOM p-val", "FDR q-val",  "Tag \\%", "Gene \\%", "Signal")
       report1 <- report
       report.index2 <- order(p.vals[,2], decreasing=F)
       for (i in 1:Ng) {
           report1[i,] <- report[report.index2[i],]
       }   
       result[[1]]<-report1  
	   names(result)[1]<-"SUMMARY.RESULTS"
	   
	   n<-1
	   report2<-list()
       for (i in 1:Ng) {
          if ((p.vals[i,1] <= p.val.threshold) ||(p.vals[i, 2] <= FDR.threshold)) {
# produce report per pathway
            kk <- 1
            gene.number <- vector(length = size.G[i], mode = "character")
            gene.names <- vector(length = size.G[i], mode = "character")
			temp.gene.obs.test<-vector(length = size.G[i], mode = "numeric")
			gene.obs.pvalue<-vector(length = size.G[i], mode = "numeric")
            gene.obs.test<-vector(length = size.G[i], mode = "numeric")
			gene.obs.weight<-vector(length = size.G[i], mode = "numeric")
            gene.list.loc <- vector(length = size.G[i], mode = "numeric")
            core.enrichment <- vector(length = size.G[i], mode = "character")
            rank.list <- seq(1, N)
            set.k <- seq(1, N, 1)  
            loc <- match(i, report.index2)
            for (k in set.k) {
               if (Obs.indicator[i, k] == 1) {
                  gene.number[kk] <- kk
                  gene.names[kk] <- obs.gene.labels[k]
                  gene.list.loc[kk] <- k
                  temp.gene.obs.test[kk]<-signif(obs.Test[obs.gene.labels[k]],digits=3)
				  gene.obs.pvalue[kk]<-signif(pt(temp.gene.obs.test[kk],length(class.labels)-2,lower.tail = FALSE),digits=3)
				  gene.obs.test[kk]<-paste(temp.gene.obs.test[kk],"(",gene.obs.pvalue[kk],")")
				  gene.obs.weight[kk]<-signif(obs.Weight[obs.gene.labels[k]],digits=3)                 
                  core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= Obs.arg.ES[i], "YES", "NO")
                  kk <- kk + 1
               }
            }

       gene.report <- data.frame(cbind(gene.number, gene.names,  gene.list.loc, gene.obs.test,gene.obs.weight, core.enrichment))
       names(gene.report) <- c("#", "GENE SYMBOL", "LIST LOC", "Tscore(p-value)","GIF", "CORE_ENRICHMENT")
       report2[[n]]<-gene.report
       names(report2)[n]<-gs.names[i]
	   n<-n+1
     } 
   
  }

  result[[2]]<-report2
  names(result)[2]<-"pathway.report"
  return(result)

}  

###########################################################	
##Calculate GIF by RandomWalk
CalGIF<-function(dataset,class.labels){

     if(!exists("PAGIData")) initializePAGI()

 pathway.db<-get("pathway.db",envir=PAGIData)
 netWorkdata<-get("netWorkdata",envir=PAGIData)
     
     gene.labels <- as.character(dataset[,1])
     sample.names <- names(dataset)[2:length(dataset[1,])]
     A <- data.matrix(dataset[2:length(dataset[1,])])
     row.names(A)<-gene.labels
     cols <- length(A[1,])
     rows <- length(A[,1])
     class.labels[which(class.labels==unique(class.labels)[1])]<-rep(0,length(which(class.labels==unique(class.labels)[1])))
     class.labels[which(class.labels==unique(class.labels)[2])]<-rep(1,length(which(class.labels==unique(class.labels)[2])))
     col.index <- order(class.labels, decreasing=F)
     class.labels <- as.character(class.labels[col.index])
     sample.names <- sample.names[col.index]
     for (j in 1:rows) {
      A[j, ] <- A[j, col.index]
     }
     names(A) <- sample.names
	 obs.Test<-abs(calTest(A, class.labels)$T)
     igraphM<-graph.data.frame(netWorkdata, directed=FALSE, vertices=NULL)
     ##generation of VertexWeight vector
     Ve<-V(igraphM)##取得igraphM这个图对象的顶点序列
     VertexWeight<-numeric(length=length(Ve))##定义顶点的权重向量（初始值为0）
     names(VertexWeight)<-get.vertex.attribute(igraphM, "name", index=V(igraphM))##将igraphM这个图的顶点名字赋给VertexWeight向量中元素的名字
     SeedGenes<-intersect(names(obs.Test),names(VertexWeight))##差异表达基因映射到网络中作为seednodes
     VertexWeight[SeedGenes] <- obs.Test[SeedGenes]##修改初始顶点权值向量中种子节点的权值得igraphM中每个节点的权值向量VertexWeight
	 gc()
	 geneNetworkVertexW<-RandomWalk2igraph(igraphM,VertexWeight,EdgeWeight=FALSE)
     names(geneNetworkVertexW)<-names(VertexWeight)
     W<-((geneNetworkVertexW-min(geneNetworkVertexW))/(max(geneNetworkVertexW)-min(geneNetworkVertexW)))
	 Weight<-rep(0,length(obs.Test))
	 names(Weight)<-names(obs.Test)
	 Weight[intersect(names(obs.Test),names(W))]<-W[intersect(names(obs.Test),names(W))]
	 rm(igraphM)
	 gc()
	 return(Weight)
}
###############################################################################################
###############################################################################################
###############################################################################################


# start of auxiliary functions
##RandomWalk function
RandomWalk2igraph<-function(igraphM,VertexWeight,EdgeWeight=TRUE,gamma=0.7){
      if(EdgeWeight==TRUE){
        adjM<-get.adjacency(igraphM,attr="weight") # convert igraph object to a weight matrix
      }
      if(EdgeWeight==FALSE){
        adjM<-get.adjacency(igraphM) # convert igraph object to a conventional matrix
      }
      res<-rw(adjM,VertexWeight,gamma)
      return(drop(res))
}

rw<-function(W,p0,gamma) {
      p0<-t(p0)
      p0 <- p0/sum(p0)
      PT <- p0   
      k <- 0
      delta <- 1
      Ng <- dim(W)[2]
      for (i in 1:Ng) {
        sumr<-sum(W[i,])
        if(sumr==0){
        W[i,] <-numeric(length=length(W[i,]))
        }
        if(sumr>0){
        W[i,] <- W[i,]/sum(W[i,])
        }
     }
     W<-as.matrix(W)
     W <- t(W)
   
     while(delta>1e-10) {
      PT1 <- (1-gamma)*W
      PT2 <- PT1 %*% t(PT)
      PT3 <- (gamma*p0)
      PT4 <- t(PT2) + PT3
      delta <- sum(abs(PT4 - PT))
      PT <- PT4
      k <- k + 1
    }
    PT<-t(PT)
    rownames(PT)<-NULL
    return(PT)
}
#################################################
##Calculate T-test
calTest <- function(A, class.labels) {
      m <- length(ind1 <- which(class.labels == unique(class.labels)[1]))
      n <- length(ind2 <- which(class.labels == unique(class.labels)[2]))
      inData1 <- A[, ind1]
      inData2 <- A[, ind2]
      rmean1 <- rowMeans(inData1)
      rmean2 <- rowMeans(inData2)
      ss1 <- rowSums((inData1 - rmean1)^2)
      ss2 <- rowSums((inData2 - rmean2)^2)
      tt <- (m + n - 2)^0.5 * (rmean2 - rmean1)/((1/m + 1/n) * (ss1 + ss2))^0.5
      return(list(T = tt, df = m + n - 2))
    }
###########################################################	
##Calculate GIF by RandomWalk
CalGIF2<-function(obs.Test, netWorkdata){

     if(!exists("PAGIData")) initializePAGI()

     pathway.db<-get("pathway.db",envir=PAGIData)
     netWorkdata<-get("netWorkdata",envir=PAGIData)
     
     igraphM<-graph.data.frame(netWorkdata, directed=FALSE, vertices=NULL)
     ##generation of VertexWeight vector
     Ve<-V(igraphM)##取得igraphM这个图对象的顶点序列
     VertexWeight<-numeric(length=length(Ve))##定义顶点的权重向量（初始值为0）
     names(VertexWeight)<-get.vertex.attribute(igraphM, "name", index=V(igraphM))##将igraphM这个图的顶点名字赋给VertexWeight向量中元素的名字
     SeedGenes<-intersect(names(obs.Test),names(VertexWeight))##差异表达基因映射到网络中作为seednodes
     VertexWeight[SeedGenes] <- obs.Test[SeedGenes]##修改初始顶点权值向量中种子节点的权值得igraphM中每个节点的权值向量VertexWeight
	 gc()
	 geneNetworkVertexW<-RandomWalk2igraph(igraphM,VertexWeight,EdgeWeight=FALSE)
     names(geneNetworkVertexW)<-names(VertexWeight)
     W<-((geneNetworkVertexW-min(geneNetworkVertexW))/(max(geneNetworkVertexW)-min(geneNetworkVertexW)))
	 Weight<-rep(0,length(obs.Test))
	 names(Weight)<-names(obs.Test)
	 Weight[intersect(names(obs.Test),names(W))]<-W[intersect(names(obs.Test),names(W))]
	 rm(igraphM)
	 gc()
	 return(Weight)
}
###########################################################
PAGI.EnrichmentScore <- function(gene.list, gene.set, correl.vector = NULL) {  
   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   correl.vector <- abs(correl.vector)
   sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
   norm.tag    <- 1.0/sum.correl.tag
   norm.no.tag <- 1.0/Nm
   RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES) 
   ES <- signif(max.ES, digits = 5)
   arg.ES <- which.max(RES)
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}
##################################################################################
PAGI.EnrichmentScore2 <- function(gene.list, gene.set, correl.vector = NULL) {  
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   loc.vector <- vector(length=N, mode="numeric")
   peak.res.vector <- vector(length=Nh, mode="numeric")
   valley.res.vector <- vector(length=Nh, mode="numeric")
   tag.correl.vector <- vector(length=Nh, mode="numeric")
   tag.diff.vector <- vector(length=Nh, mode="numeric")
   tag.loc.vector <- vector(length=Nh, mode="numeric")
   loc.vector[gene.list] <- seq(1, N)
   tag.loc.vector <- loc.vector[gene.set]
   tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
   tag.correl.vector <- correl.vector[tag.loc.vector]
   tag.correl.vector <- abs(tag.correl.vector)
   norm.tag <- 1.0/sum(tag.correl.vector)
   tag.correl.vector <- tag.correl.vector * norm.tag
   norm.no.tag <- 1.0/Nm
   tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
   tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
   tag.diff.vector <- tag.diff.vector * norm.no.tag
   peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)   
   max.ES <- max(peak.res.vector)
   ES <- signif(max.ES, digits=5)
   return(list(ES = ES))
}
# end of auxiliary functions
