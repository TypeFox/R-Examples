OUshifts <- function(y, phy, method=c("mbic","aic","bic","saic","sbic"), nmax, check.pruningwise = TRUE) {
	if (!inherits(phy, "phylo")) 
      	stop("object \"phy\" is not of class \"phylo\".")
	if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")
  
  origphy = phy
	if (check.pruningwise) phy = reorder(phy, "pruningwise")
  
	method = match.arg(method)
	options(warn = -1)
	n <- length(phy$tip.label)
	N <- dim(phy$edge)[1]
	ROOT <- n + 1L
	anc <- phy$edge[, 1]
	des <- phy$edge[, 2]
	el <- phy$edge.length
	v <- matrix(0,n + phy$Nnode,n)
	
	if (is.null(names(y))) 
	  warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
	else{
	  ordr = match(phy$tip.label,names(y))
	  if (sum(is.na(ordr))>0)
	    stop("the row names in the data do not match the tip labels in the tree.\n")
	  y = y[ordr]
	}
  
  for (i in 1:N) {
		if (des[i] <= n) v[des[i],des[i]] = 1
		v[anc[i],] = v[anc[i],]+v[des[i],]		
	}

	check <- function(model) {
		### return identifiability check and modified BIC penalty
		if (length(model)==0) return(c(1,log(n)))
		checkpar = rep(ROOT,n)
		for (i in N:1) 
			if (i %in% model) checkpar[which(v[des[i],]==1)] = i
		if (length(levels(as.factor(checkpar))) == length(model)+1) {
			numchange = length(model)
			pen = 3*log(n) + (2*numchange - 1)*log(n)
			for (j in 1:numchange)
				pen = pen + log(length(which(as.factor(checkpar)==levels(as.factor(checkpar))[j])))
			return(c(1,pen)) 
			} else return(c(0,0))
	}	

	score <- function(model,method) {
		if (length(model)==0) {
      output = phylolm(y~1,phy=phy,model="OUfixedRoot")
			logLik = output$logLik
			if (method=="aic") score = -2*logLik + 2*(length(model)+3)
			if (method=="bic") score = -2*logLik + log(n)*(length(model)+3)
			if (method=="mbic") score = -2*logLik + 3*log(n)
			if (method=="saic") score = -2*logLik + 2*(length(model)+3)
			if (method=="sbic") score = -2*logLik + log(n)*(length(model)+3)
			return(list(score = score, output = output))
		}
		X = v[ROOT,]
		for (i in 1:length(model)) X = cbind(X,v[des[model[[i]]],])
		ch = check(model)
		if (ch[1]==0) return(list(score = Inf, output = NULL))
    output = phylolm(y~X-1,phy=phy,model="OUfixedRoot")
		logLik = output$logLik
		if (method=="aic") score = -2*logLik + 2*(length(model)+3)
		if (method=="bic") score = -2*logLik + log(n)*(length(model)+3)
		if (method=="mbic") score = -2*logLik + ch[2]
		if (method=="saic") score = -2*logLik + 2*(2*length(model)+3)
		if (method=="sbic") score = -2*logLik + log(n)*(2*length(model)+3)
    
		return(list(score = score, output = output))
	}

	curmodel = list()	
  current = score(curmodel,method)
	flag = 0
  
	promodel = curmodel
	pro = current
  
	while (flag==0) {		
		for (i in 1:N) {
			if (sum(which(curmodel==i))>0) {	
				pos = which(curmodel==i)
				tempmodel = curmodel[-pos]
        temp = score(tempmodel,method)
        
				if (temp$score < pro$score) {
					promodel = tempmodel
					pro = temp
					flag = flag + 1
				}
			}
      
			if (sum(which(curmodel==i))==0) {
				if (length(curmodel) < nmax) {
					tempmodel = curmodel
					tempmodel[[length(tempmodel)+1]] = i
					temp = score(tempmodel,method)
					
          if (temp$score < pro$score) {
						promodel = tempmodel
						pro = temp
						flag = flag + 1
					}
				}
				if (length(curmodel)>=1)
					for (j in 1:length(curmodel)) 
						if ((anc[i]==des[curmodel[[j]]])||(des[i]==anc[curmodel[[j]]])||(anc[i]==anc[curmodel[[j]]])) {
							tempmodel = curmodel
							tempmodel[j] = i
							temp = score(tempmodel,method)
							
              if (temp$score < pro$score) {
							promodel = tempmodel 
							pro = temp
							flag = flag + 1
							}
						}	
			}		
		}
		if (flag > 0) flag = 0 else flag = 1
		curmodel = promodel
		current = pro
	}
  
  ### return edges from the original tree
	convert = match(phy$edge[, 2], origphy$edge[, 2])
  convertmodel = convert[unlist(curmodel)]
  
	results = list(y = y, phy = origphy, score = current$score,
                 nmax = nmax, nshift = length(convertmodel))
	
  results$alpha = current$output$optpar
	results$sigma2 = current$output$sigma2
	results$mean = as.numeric(current$output$coefficients[1])
	if (results$nshift==0) {
	  results$pshift = NULL
    results$shift = NULL
	} else {
	  results$pshift = convertmodel
    results$shift = current$output$coefficients[-1]
    names(results$shift) = results$pshift
	}
  
	class(results) = "OUshifts"
	return(results)
	}

################################################
plot.OUshifts <-function(x, show.data = TRUE, digits = 3, ...) {
  if (show.data) layout(matrix(1:2,1,2))
  plot(x$phy, ...)  
  edgelabels(edge=x$pshift,pch=8,frame="n",cex=1.1,lwd=2,col="red")
  #x$y = as.vector((x$y-mean(x$y))/sd(x$y))
  #barplot(x$y, horiz=TRUE, xlab="Normalized trait", names.arg="")	
  if (show.data) {
    lastPP = get("last_plot.phylo", envir=.PlotPhyloEnv)
    o = order(lastPP$yy[1:length(x$phy$tip.label)])
    normy = as.vector((x$y-mean(x$y))/sd(x$y))
    barplot(normy[o], horiz=TRUE, xlab="Trait", names.arg="",
          xaxt="n")
    axis(1,at=range(normy),labels=round(range(x$y),digits=digits))
  }
}






