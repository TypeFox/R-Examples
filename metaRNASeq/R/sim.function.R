sim.function <-
function(param, dispFuncs, nrep=4, classes=NULL, inter.sd = 0.30) {
  
	if(is.null(classes)) {
		classes=list(c(rep(1,nrep),rep(2,nrep)), c(rep(1,nrep),rep(2,nrep))) }
	sim.means=param[,c("mucond1","mucond2")] 
	
	vecstudies <- unlist(lapply(1:length(classes), FUN=function(x) rep(x,length(classes[[x]]))))
	vecstudies12 <- vecstudies
	vecstudies12[which(vecstudies %in% seq(1,19,by=2))]=1        
	vecstudies12[which(vecstudies %in% seq(2,20,by=2))]=2

  varinter <- vector("list", length(classes))
	varinter <- lapply(varinter, function(x) matrix(rnorm(dim(param)[1]*2, mean = 0, sd = inter.sd),
		ncol = 2))

	for(i in 1:length(varinter)) {
		varinter[[i]][which(param$DE == FALSE),] <- 
			cbind(varinter[[i]][which(param$DE == FALSE),1],
			varinter[[i]][which(param$DE == FALSE),1])
	}
	
	calc.size1 <- function(mean) 1 / (dispFuncs[[1]][1] + dispFuncs[[1]][2] / mean)
	calc.size2 <- function(mean) 1 / (dispFuncs[[2]][1] + dispFuncs[[2]][2] / mean) 
	
	matsim <- mapply(function(condition, study) {
			mean <- sim.means[,condition] * exp(varinter[[study]][,condition])

				size <- get(paste("calc.size", study, sep = ""))(rowMeans(sim.means));
				size <- ifelse(mean > 10, size, 10e10)
			
			rnbinom(dim(param)[1], mu = mean, size = size)}, 
		unlist(classes), vecstudies12)

	colnames(matsim)=paste("study",vecstudies,"cond",unlist(classes),sep="")
	matsim
}

extractfromsim <- function(matsim,studyname)
{
  study <- matsim[,grep(studyname,colnames(matsim))]
  study.conds <- gsub(studyname,"",colnames(study))
  study.conds <- as.factor(study.conds)
  levels(study.conds) <- c("untreated","treated")
  colnames(study) <- paste("rep",1:dim(study)[2],sep="")
  pheno <- data.frame(study=rep(studyname,dim(study)[2]),condition=study.conds,row.names=colnames(study))
  simstudy <- list("study"=study,"pheno"=pheno)
  simstudy
}