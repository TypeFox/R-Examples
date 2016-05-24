design.score <- function(Mi,model,ploidy,min.MAF,max.geno.freq){
  n <- length(Mi)
  freq <- mean(Mi,na.rm=T)/ploidy
  if (min(freq,1-freq) < min.MAF){
    return(NULL)
  } else {
    if (model=="additive") {
      geno.freq <- table(Mi)/n
      if (max(geno.freq) <= max.geno.freq) {
        return(matrix(Mi))
      } else {
        return(NULL)
      }
    } 
    if (model=="diplo-additive") {
      Mi[which((Mi>0)&(Mi<ploidy))] <- ploidy/2
      Mi <- Mi/(ploidy/2)
      geno.freq <- table(Mi)/n
      if (max(geno.freq) <= max.geno.freq) {
        return(matrix(Mi))
      } else {
        return(NULL)
      }
    }
    if (model=="diplo-general") {
      Mi[(Mi>0)&(Mi<ploidy)] <- ploidy/2
      Mi <- Mi/(ploidy/2)
      geno.freq <- table(Mi)/n
      if (max(geno.freq)<=max.geno.freq) {
        tmp <- model.matrix(~x,data.frame(x=factor(Mi)))[,-1]
        if (is.null(dim(tmp))) {
          return(matrix(tmp))
        } else {
          return(tmp)
        }
      } else {
        return(NULL)
      }
    }
    if (length(grep("dom",model,fixed=T))>0) {
      tmp <- strsplit(model,split="-",fixed=T)[[1]]
      dom.order <- as.integer(tmp[1])
      if (tmp[3]=="alt") {
        Mi <- ifelse(Mi>=dom.order,1,0)
      } else {
        Mi <- ifelse(Mi<=ploidy-dom.order,0,1)
      }
      geno.freq <- table(Mi)/n
      if (max(geno.freq) <= max.geno.freq) {
        return(matrix(Mi))
      } else {
        return(NULL)
      }
    }
    if (model=="general") {
      geno.freq <- table(Mi)/n
      if (max(geno.freq)<=max.geno.freq) {
        tmp <- model.matrix(~x,data.frame(x=factor(Mi)))[,-1]
        if (is.null(dim(tmp))) {
          return(matrix(tmp))
        } else {
          return(tmp)
        }
      } else {
        return(NULL)
      }
    }
  }
}
