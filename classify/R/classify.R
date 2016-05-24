utils::globalVariables(c("tru","acc","grade","M1","M2","corr","score"))

setClass("classification",
		representation (
				acc = "numeric",
				fp = "numeric",
				fn = "numeric",	
				cand.acc="matrix",
				cand.fn="matrix",
				cand.fp="matrix",
				consist="numeric",
				cand.consist="matrix",
				Xij="matrix",
                kappa = "numeric",
				item.scores = "matrix",
				bnds = "numeric",
				tru.grades = "matrix",
				tru.scores = "matrix",
				acc.by.grade = "matrix",
				fp.by.grade = "matrix",
				fn.by.grade = "matrix",
				consist.by.grade = "matrix",
				labels = "character",
				m.acc="numeric",
				m.consist="numeric",
				m.kappa="numeric",
				m.fp="numeric",
				m.fn="numeric",
				sd.acc="numeric",
				sd.consist="numeric",
				sd.kappa="numeric",
				sd.fp="numeric",
				sd.fn="numeric",
				m.acc.by.grade="numeric",
				m.fp.by.grade="numeric",
				m.fn.by.grade="numeric",
				m.consist.by.grade="numeric",
				max = "numeric"
))

setMethod("summary",
		signature(object="classification"),
		function(object){
			acc.by.grade <- cbind(accuracy=object@m.acc.by.grade,false.positive=object@m.fp.by.grade,false.negative=object@m.fn.by.grade,consistency=object@m.consist.by.grade)
			row.names(acc.by.grade)<- object@labels
			cat("Marginal Classification Accuracy: ",round(object@m.acc,3),"(",round(object@sd.acc,3),")","\n",sep="")
			cat("Marginal Classification Consistency :",round(object@m.consist,3),"(",round(object@sd.consist,3),")","\n",sep="")
			cat("Kappa: ",round(object@m.kappa,3),"(",round(object@sd.kappa,3),")","\n",sep="")
			cat("Marginal False Negative Error Rate: ",round(object@m.fp,3),"(",round(object@sd.fp,3),")","\n",sep="")
			cat("Marginal False Positive Error Rate: ",round(object@m.fn,3),"(",round(object@sd.fn,3),")","\n",sep="")
			cat("Accuracy by grade:","\n")
			print(acc.by.grade)
		})

setGeneric("across.reps",
  function(object)
    standardGeneric("across.reps")
)

setMethod("across.reps",
	signature(object="classification"),
	function(object){
		object@m.acc <- mean(object@acc)
		object@m.consist <- mean(object@consist)
		object@m.kappa <- mean(object@kappa)
		object@m.fp <- mean(object@fp)
		object@m.fn <- mean(object@fn)
		object@sd.acc <- sd(object@acc)
		object@sd.consist <- sd(object@consist)
		object@sd.kappa <- sd(object@kappa)
		object@sd.fp <- sd(object@fp)
		object@sd.fn <- sd(object@fn)
		object@m.acc.by.grade <- rowMeans(object@acc.by.grade,na.rm=T)
		object@m.fp.by.grade <- rowMeans(object@fp.by.grade,na.rm=T)
		object@m.fn.by.grade <- rowMeans(object@fn.by.grade,na.rm=T)
		object@m.consist.by.grade <- rowMeans(object@consist.by.grade,na.rm=T)
		return(object)
	}
)	

plot.classification <- function(x,  type = c("acc", "kappa", "density"), ...) {
if (!inherits(x, "classification"))
        stop("Use only with 'classification' objects.\n")
		type <- match.arg(type)
	if(type=="acc"){
		df <- cbind(melt(slot(x, "tru.grades")),melt(slot(x, "cand.acc")),melt(slot(x, "tru.scores")))
		names(df) <- c("cand","sim","grade","cand2","sim2","acc","cand2","sim2","tru")
		df$grade <- factor(df$grade,labels=slot(x, "labels"))
		p <- ggplot(df,aes(x=factor(tru),y=acc))
		p <- p + geom_boxplot(aes(fill = grade)) 
		p <- p + scale_x_discrete(name="True Score",breaks = NULL)
		p <- p + scale_y_continuous(name="Accuracy",limits=c(0,1)) 
		p <- p + scale_fill_brewer(name="grade",palette="Set1")		
		print(p)
	} else if (type=="kappa"){
		kappa.df <- melt(slot(x, "Xij"))
		names(kappa.df) <- c("M1", "M2", "corr") 
		p <- ggplot(kappa.df, aes(M1, M2, fill=corr))
		p <- p + geom_tile() + theme_bw() 
		p <- p + geom_text(data=kappa.df,aes(label=round(corr,2)),size=3)
		meas <- as.character(unique(kappa.df$M2))
		p <- p + scale_colour_identity() + scale_fill_gradientn(colours= c("grey", "white"), limits=c(0.25,0))
		p <- p + scale_x_discrete(limits=meas) + scale_y_discrete(limits=meas)
		p <- p + xlab(NULL) + ylab(NULL)
            p <- p + theme(axis.ticks=element_blank(),panel.grid.major=element_blank(),legend.position="none")
		print(p)
	} else if (type=="density"){
      	scores <- data.frame(score=rowSums(slot(x, "item.scores")))
	    p <- ggplot(scores,aes(x=score))
      	p <- p + geom_density()
		p <- p + geom_vline(xintercept = slot(x, "bnds"),linetype=2,colour="grey")
		p <- p + scale_x_continuous(limits=c(0,slot(x, "max")))
	    print(p)
	}
  }


classify <- function(cssd,expected,bnds,cats,lbls=NULL){
	#Candidate true scores - most probable
	totals <- round(expected,0)
	n.cands <- ncol(cssd)
	mx.grade <- length(bnds)-1

	cand.grades <- cut(totals, bnds, right=F, labels=FALSE, include.lowest=T)

	cand.accs <- vector(mode = "numeric", length = n.cands)
	cand.fn <- vector(mode = "numeric", length = n.cands)
	cand.fp <- vector(mode = "numeric", length = n.cands)

	p.cumsums <- apply(cssd,2,cumsum)
	
	for(cand in 1:n.cands){
		if(cand.grades[cand]==mx.grade){
			cand.accs[cand] <- 1-p.cumsums[bnds[mx.grade]+1,cand]
			cand.fn[cand] <- 0			
			cand.fp[cand] <- 1- cand.accs[cand]
		} else {
			if(cand.grades[cand]==1){
				cand.accs[cand] <- p.cumsums[bnds[2],cand]
				cand.fp[cand] <- 0
				cand.fn[cand] <- 1- cand.accs[cand]
			} else {
				cand.accs[cand] <- p.cumsums[bnds[cand.grades[cand]+1],cand]-p.cumsums[bnds[cand.grades[cand]],cand]
				#False postive is the cumulative perc of highest mark of grade below
				cand.fp[cand] <-  p.cumsums[bnds[cand.grades[cand]],cand]
				#False negative is 1 minus the cumulative perc of highest mark in grade of candidate
				cand.fn[cand] <- 1 - p.cumsums[bnds[cand.grades[cand]+1],cand]
			}
		}
	}

    #	consistency code
    cands <- length(expected) 
    scores <- 0:(sum(cats-1))
    grades <- data.frame(cbind(cssd,scores,cut(scores,bnds,right=F,include.lowest=T,labels=F)))
    names(grades) <- c(paste("c",1:cands,sep=""),"score","grade")
    grade.p <- ddply(grades, .(grade), function(df) colwise(sum)(df[, 1:cands]))
    grade.p <- grade.p[,-1]
    grade.consistency <- grade.p ^2
    grade.consistency.m <- as.numeric(colSums(grade.consistency))
    sum.consist <- sum(grade.consistency)/cands	

    #kappa code
    Xij <- matrix(0,nrow=length(bnds)-1, ncol = length(bnds)-1)
    if(is.null(lbls)){
    	rownames(Xij) <- paste("C", seq(1, dim(Xij)[1]), sep = "")
    	colnames(Xij) <- paste("C", seq(1, dim(Xij)[2]), sep = "")
    } else {
    	rownames(Xij) <- lbls
    	colnames(Xij) <- lbls
    }
			
    for(i in 1:cands){
		#for Cohen's kappa each row becomes a matrix
		tmp <- grade.p[,i]
		#calculate matrix of probabilities
		outer.prod <- t(t(tmp))%*%t(tmp)
		#add to sum matrix
		Xij <- Xij + outer.prod
	}
	#convert to proportions
	Xij.p <- Xij/sum(Xij)
			
	#sum diagonal
	d <- sum(diag(Xij.p))
	#adjust for chance
	#calculate marginals
	m <- rowSums(Xij.p)
	chance <- sum(m^2)
	kappa <- (d-chance)/(1-chance)

	#Accuracy by grade
	grdes <- length(bnds)-1
	cand.grades <- factor(cand.grades,levels=1:grdes,labels=lbls)
	acc.by.grade <- tapply(cand.accs,cand.grades,mean)
	fp.by.grade <- tapply(cand.fp,cand.grades,mean)
	fn.by.grade <- tapply(cand.fn,cand.grades,mean)
	consist.by.grade <- tapply(grade.consistency.m,cand.grades,mean)

	return(list(cand.accs,cand.fn,cand.fp,sum.consist,Xij.p,kappa,grade.consistency.m, acc.by.grade, fp.by.grade, fn.by.grade, consist.by.grade))

}

classify.bug <- function(sims,scores,bnds,lbls=NULL){

    mdl <-  slot(scores, "model")	
    gibbs <- slot(scores, "estimation")	
    cats <-  slot(scores, "cats")	
    expected <- slot(scores, "expected")	
    #top and tail boundaries
    n.sims <- slot(scores, "sims")	
    bns <- bnds
    bnds <- c(0,bnds,sum(cats-1))
    grades <- length(bnds)-1
	mx <- slot(scores, "max")
	
	if (max(bnds)>mx) stop("Grade boundary higher than maximum score\n")	

    if(is.null(lbls)){
    	for (i in 1:grades){
		if(i<grades){
			lbls <- c(lbls,paste(bnds[i],bnds[i+1]-1,sep=":"))
		} else {
			lbls <- c(lbls,paste(bnds[i],bnds[i+1],sep=":"))
		}
    	}	
    }

	if (length(lbls)!=grades) stop("Mismatch between grade labels and boundaries\n")
	if (length(bns)>1 & bns[1]-bns[2]>=0) stop("Grade boundaries should increase through the vector\n")
	
    n <- slot(scores, "persons")	
    accuracy.m <- matrix(nrow=n,ncol=n.sims)
    consistency.m <- matrix(nrow=n,ncol=n.sims)
    fp.m <- matrix(nrow=n,ncol=n.sims)
    fn.m <- matrix(nrow=n,ncol=n.sims)
    kappa <- vector("numeric",length=n.sims)
    consistency <- vector("numeric",length=n.sims)

    Xij <- matrix(0,nrow=grades, ncol = grades)
    pb <- txtProgressBar(min = 1, max = n.sims, style = 3)
    acc.by.grade <- matrix(nrow=grades,ncol=n.sims)
    fp.by.grade <- matrix(nrow=grades,ncol=n.sims)
    fn.by.grade <- matrix(nrow=grades,ncol=n.sims)
    consist.by.grade <- matrix(nrow=grades,ncol=n.sims)

	for(i in 1:n.sims){
		v <- sims[i,]
		pars <- gpcm.bug(v,cats,mdl,gibbs)
		beta <- pars[[1]]
		theta <- pars[[2]]
            if(mdl=="pcm.bug" || mdl=="rasch.bug"){
			cssd <- gpcm.rc(beta,theta,cats)
		} else {
			alpha <- pars[[3]]
			cssd <- gpcm.rc(beta,theta,cats,alpha)
		}
		out <- classify(cssd,expected[,i],bnds,cats,lbls)
		accuracy.m[,i] <- out[[1]]
		fn.m[,i] <- out[[2]]
		fp.m[,i] <- out[[3]]
		consistency[i] <- out[[4]] 		
		Xij <- Xij + out[[5]]		
		kappa[i] <-  out[[6]]
		consistency.m[,i] <- out[[7]]
		acc.by.grade[,i] <- out[[8]]
		fp.by.grade[,i] <- out[[9]]
		fn.by.grade[,i] <- out[[10]]
		consist.by.grade[,i] <- out[[11]]
		setTxtProgressBar(pb, i)
	}
	m.acc = colMeans(accuracy.m)
	m.fp = colMeans(fp.m)
	m.fn = colMeans(fn.m)
	Xij <- Xij / n.sims
	totals <- round(expected,0)
	cand.grades <- matrix(cut(totals, bnds, right=F, labels=FALSE, include.lowest=T),ncol=n.sims)
	ret <- new("classification", acc = m.acc, fp = m.fp, fn = m.fn, cand.acc= accuracy.m,	cand.fn= fn.m,	cand.fp=fp.m,	consist=consistency,	cand.consist=consistency.m,	Xij=Xij, kappa = kappa, item.scores=slot(scores, "item"), bnds = bns, tru.grades = cand.grades, tru.scores = totals,labels = lbls, acc.by.grade=acc.by.grade,fp.by.grade=fp.by.grade,fn.by.grade=fn.by.grade,consist.by.grade=consist.by.grade,max=mx)
	ret <- across.reps(ret)
	cat("100%")
	close(pb)
	return(ret)
}
