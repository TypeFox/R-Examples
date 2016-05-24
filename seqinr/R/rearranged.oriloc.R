#############################################################################
#                        Rearranged oriloc                                  #  
#                                                                           #
# Detection of replication-associated effects on base composition asymmetry #
#                                                                           #
#############################################################################

rearranged.oriloc <- function(
    seq.fasta = "ftp://pbil.univ-lyon1.fr/pub/seqinr/data/ct.fasta" ,
    g2.coord = "ftp://pbil.univ-lyon1.fr/pub/seqinr/data/ct.coord" )
  {
    
    seq.fasta <- read.fasta(seq.fasta)[[1]]
    g2 <- readLines(g2.coord)
    g2 <- lapply(g2, function(x) unlist(strsplit(x, split=" ")))
    g2 <- lapply(g2, function(x) x[which(x!="")])
        
    start <-as.numeric(unlist(lapply(g2, function(x) x[2]))) 
    end <- as.numeric(unlist(lapply(g2, function(x) x[3])))
    strand <- rep("forward",length(start))
    strand[which(end<start)] <- "reverse"


    meancoord <- (start+end)/2

    ncds <- length(start)

    gcskew <- numeric(ncds)
    atskew <- numeric(ncds)
    
    for(i in seq_len(ncds)){

      cds=seq.fasta[start[i]:end[i]]
      
      g3=sum(cds[seq(from=3,to=length(cds),by=3)]%in%c("g","G"))
      c3=sum(cds[seq(from=3,to=length(cds),by=3)]%in%c("c","C"))
      t3=sum(cds[seq(from=3,to=length(cds),by=3)]%in%c("t","T"))
      a3=sum(cds[seq(from=3,to=length(cds),by=3)]%in%c("a","A"))

      gcskew[i]=(g3-c3)/(g3+c3)
      atskew[i]=(a3-t3)/(a3+t3)
      if(is.na(gcskew[i])){
        gcskew[i]=0
      }
      if(is.na(atskew[i])){
        atskew[i]=0
      }
      
    }

    neworder=c(which(strand=="forward"),which(strand=="reverse"))

    meancoord.rear <- seq_len(ncds)
    gcskew.rear=gcskew[neworder]
    atskew.rear=atskew[neworder]
    strand.rear=strand[neworder]
    meancoord.real=meancoord[neworder]
    data=data.frame(meancoord.rear,gcskew.rear,atskew.rear,strand.rear,neworder,meancoord.real)

    colnames(data)=c("meancoord.rearr","gcskew.rearr","atskew.rearr","strand.rearr","order","meancoord.real")

    return(data)
    
    
  }
