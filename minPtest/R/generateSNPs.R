generateSNPs <-
function(n,gene.no,block.no,block.size,p.same,p.different=NULL,p.minor,n.sample,SNPtoBETA){
  call <- match.call()
  if(length(p.minor)!=block.no) { stop("length of minor allel frequency vector p.minor should equal the number of blocks block.no") }
  snp.no <- gene.no*block.no*block.size
  if(max(SNPtoBETA[,1])>snp.no){
    stop("SNP items for the offset beta should not exceed the number of actual SNPs")
  }
  if(length(p.same)>1){
    if(!is.null(p.different)){
      warning("p.differnent is not used if p.same is a vector")
    }
    if(length(p.same)!=(block.no*block.size)){
      stop("length of the vector p.same should equal the number of SNPs per gene")
    }
    p.vec <- p.same
  }else{
    if(is.null(p.different)){
      stop("p.different is missing if p.same is a numeric value")
    }
    p.vec <- rep(c(p.different,rep(p.same,block.size - 1)),block.no)
  }
  minor <- lapply(seq_along(p.minor),function(i){
    x <- rep(p.minor[i],block.size)
    x
  })
  p.minor.vec <- as.vector(do.call(cbind, minor))
  sim.block <- function(N,minor.allel,same){
    same[1] <- 0
    not.same <- rbinom(N*length(same),1,rep(1-same,N))
    matrix(rbinom(sum(not.same),1,rep(minor.allel,N)[not.same==1])[cumsum(not.same)],N,length(same),byrow=TRUE)
  }
  cov.continuous <- sample(c(20:80),size=n,replace=TRUE)
  sample.cov.continuous <- sample(cov.continuous,size=n/2)
  beta.cov.continuous <- 0.01
  beta.cov.binary <- 0.1
  beta.snp <- rep(0,length=snp.no)
  beta.snp[SNPtoBETA[,1]] <- SNPtoBETA[,2]
  sim.temp <- lapply(1:(n/2),function(i){
    sim.snp <- matrix(unlist(lapply(1:gene.no,function(j){
      x1 <- sim.block(N=n.sample,minor.allel=p.minor.vec,same=p.vec)
      x2 <- sim.block(N=n.sample,minor.allel=p.minor.vec,same=p.vec)
      x <- x1+x2
    })),ncol=snp.no,nrow=n.sample)
    cov.binary <- rbinom(1,1,0.6)
    sim.x <- cbind(sim.snp,rep(sample.cov.continuous[i],length=n.sample),rep(cov.binary,length=n.sample))
    sim.beta <- c(beta.snp,beta.cov.continuous,beta.cov.binary)
    nue <- exp(sim.x%*%sim.beta)/(1+exp(sim.x%*%sim.beta))
    sim.y <- rbinom(n.sample,1,nue)
    y_case <- which(sim.y==1)
    case <- sample(y_case,1)
    y_control <- which(sim.y==0)
    control <- sample(y_control,1)
    matchset <- i
    n_case <- c(sim.y[case],sim.x[case,],matchset)
    n_control <- c(sim.y[control],sim.x[control,],matchset)
    dat <- c(n_case,n_control)
  })
  dim.sim.temp <- snp.no+4
  data.temp <- lapply(sim.temp, matrix, ncol=dim.sim.temp, byrow=TRUE)
  sim.data <- do.call(rbind, data.temp)
  SNPs <- paste("SNP",1:snp.no,sep="")
  colnames(sim.data) <- c("response",SNPs,"cov.continuous","cov.binary","matchset")
  rownames(sim.data) <- c(1:n)
  GENEs.temp <- paste("G",1:gene.no,sep="")
  g <- lapply(seq_along(GENEs.temp),function(i){
    g.temp <- rep(GENEs.temp[i],length=snp.no/gene.no)
  })
  Genes <- as.vector(do.call(cbind, g))
  y <- as.vector(sim.data[,"response"])
  x <- sim.data[,SNPs]	
  SNPtoGene <- cbind(SNPs,Genes)
  cov <- cbind(sim.data[,"cov.continuous"],sim.data[,"cov.binary"])
  colnames(cov) <- c("cov.continuous","cov.binary")
  matchset <- as.vector(sim.data[,"matchset"])
  fit <- list(call=call,snp.no=snp.no,sim.data=sim.data,
              y=y,x=x,SNPtoGene=SNPtoGene,cov=cov,matchset=matchset)
  class(fit) <- "generateSNPs"
  fit
}

print.generateSNPs <- function(x, ...){
  cat("Call: ",deparse(x$call),"\n")
  cat("\nSimulated data set containing", x$snp.no, "SNPs, two matching covariates and a matchset column (containing the matchset numbers).", "\n")
  cat("Output y, x, SNPtoGene, cov and matchset can directly be used for the minPtest function.", "\n\n")
  invisible()
}

