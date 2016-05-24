utils::globalVariables(c("observed","expected","raw","model","alpha","score","frequency","model"))

setClass("scores",
		representation (
				item = "matrix",
				expected="matrix",
				conditional="matrix",
				summed = "numeric",
				persons = "numeric",
				items = "numeric",
				sims = "numeric",
				max = "numeric",
				cats = "numeric",
				model = "character",
				estimation = "character"			
))

plot.scores <- function(x,  type = c("exp","cond","qq"),  alpha = 0.05, ...) {
if (!inherits(x, "scores"))
        stop("Use only with 'scores' objects.\n")	
		type <- match.arg(type)
	if(type=="exp"){
		df <- data.frame(cbind(slot(x, "summed"),slot(x, "expected")))
		names(df)<- c("observed",paste("X",1:slot(x, "sims"),sep=""))
		plot.df <- melt(df,id.vars=c("observed"))
		names(plot.df) <- c("observed","sim","expected")
		p <- ggplot(plot.df,aes(x=observed,y=expected))
		p <- p + geom_point(alpha=alpha)
		p <- p + geom_abline(a=0,b=1)
		p <- p + scale_x_continuous(name="Observed",limits=c(0,slot(x, "max")))
		p <- p + scale_y_continuous(name="Expected",limits=c(0,slot(x, "max")))
		print(p)
	} else if (type=="cond") {
		#total score
		ttl <- slot(x, "max")
		#observed candidates
		cs <- slot(x, "persons")
		u.limit <- apply(slot(x, "conditional"),1,quantile,probs=0.95)
		l.limit <- apply(slot(x, "conditional"),1,quantile,probs=0.05)
	
		u.limit <- u.limit/cs
		l.limit <- l.limit/cs
	
		scres <- slot(x, "summed")
	
		f.table <- as.data.frame(table(factor(scres, levels = 0:ttl)))
	
		obs <- data.frame(score=numeric(ttl+1),frequency=numeric(ttl+1),model=character(ttl+1))
		obs$score<-0:ttl
		obs$frequency <- f.table$Freq/cs
		obs$model <- 1
	
		mdl1 <- data.frame(score=numeric(ttl+1),frequency=numeric(ttl+1),model=character(ttl+1))
		mdl1$score<-0:ttl
		mdl1$model <- 2
		mdl1$frequency <- u.limit
	
		mdl2 <- data.frame(score=numeric(ttl+1),frequency=numeric(ttl+1),model=character(ttl+1))
		mdl2$score<-0:ttl
		mdl2$model <- 3
		mdl2$frequency <- l.limit
	
		df <- rbind(obs,mdl1,mdl2)
		df$model <- factor(df$model,labels=c("observed","upper .95 marginal","lower .95 marginal"))
	
		p <- ggplot(df,aes(x=score,y=frequency,group=model))
		p <- p + geom_line(aes(linetype=model))
		p <- p + scale_linetype_manual(values=c(1,2,2),name="")

		p <- p + scale_x_continuous(name="Summed score distribution")
		p <- p + theme(legend.position='bottom')

		print(p)

	}  else if (type=="qq") {

		#q-q plot according to median values
		ttl <- slot(x, "max")
		scres <- slot(x, "summed")
		cs <-  slot(x, "persons")
		pb <- runif(cs)
		mdn <- apply( slot(x, "conditional"),1,median)/cs
		thresh <- c(cumsum(mdn),1)
		mdl.scres <- sort(laply(pb, function (j) min(which(thresh>=j))-1))
      	df <- data.frame(raw=sort(scres),model=mdl.scres)
		p <- ggplot(df,aes(x=raw,y=model))
		p <- p + geom_jitter(alpha=alpha)
		p <- p + geom_abline(a=0,b=1)
		p <- p + scale_x_continuous(name="Observed Summed Score")
		p <- p + scale_y_continuous(name="Marginal Summed Score")
	      print(p)
	}
  }

scores.gpcm.bug <- function(item.scores,sims, mdl= c("rasch.bug", "pcm.bug", "tpl.bug","gpcm.bug"), gibbs = c("bugs","jags")){
	if(is.data.frame(item.scores))(item.scores <- as.matrix(item.scores))
	if(!is.matrix(item.scores))(stop("Item Scores must be a matrix or data frame"))
	mdl <- match.arg(mdl)
	gibbs <- match.arg(gibbs)
	itms <- ncol(item.scores)
	if(max(item.scores)==1){
		cats <- rep(2,itms)
	} else {
		cats <- as.numeric(apply(item.scores,2,max))
		cats <- cats + 1
	}
	n <- nrow(item.scores) 
	n.sims <- nrow(sims)
	expected.scores <- matrix(0,ncol=n.sims,nrow=sum(cats-1)+1)
      real.exp <- matrix(0,ncol=n.sims,nrow=n)
      pb <- txtProgressBar(min = 1, max = n.sims, style = 3) 
	for (i in 1:n.sims){	
		v <- sims[i,]
		pars <- gpcm.bug(v,cats,mdl,gibbs)
		beta <- pars[[1]]
		theta <- pars[[2]]
            if(mdl=="pcm.bug" || mdl=="rasch.bug"){
			cssd <- gpcm.rc(beta,theta,cats)
			exp.score <- expected.rc(beta,theta,cats)

		} else {
			alpha <- pars[[3]]
			cssd <- gpcm.rc(beta,theta,cats,alpha)
			exp.score <- expected.rc(beta,theta,cats,alpha)
		}
		real.exp[,i] <- exp.score
		expected.scores[,i] <- rowSums(cssd)
		setTxtProgressBar(pb, i)
	}
	cat("100%")
	close(pb)
	ret <- new("scores", item=item.scores, expected=real.exp, conditional=expected.scores, summed = rowSums(item.scores), persons = n,sims = n.sims,max=max(sum(cats-1)),items=itms, cats = cats, model = mdl, estimation = gibbs)
	return(ret)
}
