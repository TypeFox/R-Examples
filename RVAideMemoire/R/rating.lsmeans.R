rating.lsmeans <- function(lsm,type=c("prob","cumprob","class1","class2"),level=0.9) {
  name <- deparse(substitute(lsm))
  if (is.list(lsm)) {lsm <- lsm$lsmeans}
  if (!"cut" %in% names(lsm@levels)) {
    stop(paste0("no 'cut' in ",name,
	", use  ~...|cut  as formula in lsmeans()"))
  }
  type <- match.arg(type)
  summ <- as.data.frame(summary(lsm,type="response"))
  levs.cut <- levels(summ$cut)
  ncut <- nlevels(summ$cut)
  levs.cut.rep <- rle(as.character(summ$cut))$lengths
  levs.cut.all <- unique(unlist(strsplit(levs.cut,split="\\|")))
  levs.cut2 <- levs.cut.all[-length(levs.cut.all)]
  col.cut <- grep("cut",colnames(summ))
  levs <- apply(as.data.frame(summ[,1:(col.cut-1)]),1,function(x) paste0(x,collapse=":"))
  levs.unique <- levs[1:levs.cut.rep[1]]
  summ2 <- data.frame(fac=levs,cut=rep(levs.cut2,levs.cut.rep),
    cumprob=summ$cumprob)
  summ2$fac <- factor(summ2$fac,levels=levs.unique)
  summ2$cut <- factor(summ2$cut,levels=levs.cut.all)
  prob <- tapply(summ2$cumprob,list(summ2$fac,summ2$cut),identity)
  prob[,ncol(prob)] <- rep(1,nrow(prob))
  for (i in 2:ncol(prob)) {
    prob[,i] <- prob[,i]-rowSums(as.matrix(prob[,1:(i-1)],nrow=nrow(prob)))
  }
  if (type=="prob") {
    result <- summ[,1:(col.cut-1)]
    if (col.cut==2) {
	result <- as.data.frame(result)
	to.bind <- as.data.frame(result[1:levs.cut.rep[1],1])
	colnames(result) <- colnames(to.bind) <- colnames(summ)[1]
    } else {
	to.bind <- result[1:levs.cut.rep[1],]
    }
    result <- rbind(result,to.bind)
    result$Rating <- factor(rep(levs.cut.all,each=levs.cut.rep[1]),levels=levs.cut.all)
    result$Prob <- as.vector(prob)
  } else if (type=="cumprob") {
     result <- summ[,1:(col.cut-1)]
    if (col.cut==2) {
	result <- as.data.frame(result)
	to.bind <- as.data.frame(result[1:levs.cut.rep[1],1])
	colnames(result) <- colnames(to.bind) <- colnames(summ)[1]
    } else {
	to.bind <- result[1:levs.cut.rep[1],]
    }
    result <- rbind(result,to.bind)
    result$Rating <- factor(rep(levs.cut.all,each=levs.cut.rep[1]),levels=levs.cut.all)
    result$Cumprob <- c(summ2$cumprob,rep(1,levs.cut.rep[1]))
  } else if (type=="class1") {
    result <- summ[1:levs.cut.rep[1],1:(col.cut-1)]
    if (col.cut==2) {
	result <- as.data.frame(result)
	colnames(result) <- colnames(summ)[1]
    }
    result$Rating <- apply(prob,1,function(x) colnames(prob)[which.max(x)])
  } else if (type=="class2") {
    cumprob <- t(apply(prob,1,cumsum))
    result <- summ[1:levs.cut.rep[1],1:(col.cut-1)]
    if (col.cut==2) {
	result <- as.data.frame(result)
	colnames(result) <- colnames(summ)[1]
    }
    result$Rating <- apply(cumprob,1,function(x) colnames(cumprob)[which(x>=level)[1]])
  }
  return(result)
}
