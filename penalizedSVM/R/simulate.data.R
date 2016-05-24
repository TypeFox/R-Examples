
sim.data <- function(n = 256, ng = 1000, nsg = 100, p.n.ratio = 0.5, 
												sg.pos.factor= 1, sg.neg.factor= -1,
												# correlation info:
												corr=FALSE, corr.factor=0.8,
												# block info:
												blocks=FALSE, n.blocks=6, nsg.block=1, ng.block=5, 
												seed=123, ...){ 
####################################################################################
# simulate microarray data with 
# n - number of samples, logistic regression works well if n>200!
# ng - number of genes
# nsg - number of significant genes 
# p.n.ratio - ratio between positive and negative significant genes (default 0.5)  
# sg.pos.factor - impact factor of _positive_ significant genes on the classifaction (default 1)  
# sg.neg.factor - impact factor of _negative_ significant genes on the classifaction (default -1)  
# all other non-significant genes have in both classes balanced ratios

# for correlated blocks of genes
# n.blocks -  number of correlated blocks of genes
# nsg.block - number of significant genes per block
# ng.block - number of genes per block

# if no blockes (n.blocks=0) are defined and corr=TRUE
# create covarance matrix for all genes! with decrease of correlation 

####################################################################################

	require(MASS)
	set.seed(seed)
	# assumption I: intercept = 0, beta0=0
	b0<- 0
	
	# effect for positive genes = +log2(2) = +1
	#            negative genes =  log2(1/2)= -1 
	m.pos<-m.neg<-m.bal<- 0
	sd.pos<-sd.neg<-sd.bal<- 1
	
	if (!( (p.n.ratio>=0) & (p.n.ratio<=1))) stop("rato between positive and negative significant genes should be in [0;1]")
	pos.nsg<- floor (nsg * p.n.ratio)
	neg.nsg<- nsg - pos.nsg
	
	# covariance matrix
	sigma<- .create.covariance.matrix(sd.pos, sd.neg, sd.bal,
																		ng,nsg, pos.nsg, neg.nsg,
																		corr,corr.factor, 
																		blocks, n.blocks, nsg.block, ng.block)
	
	# better to see
	# sigma.see<-sigma; sigma.see[ sigma.see==0]<- ""; rm(sigma.see) 
			
	#  1 if pos, -1 by neg
	bX<-rep(0,ng); bX[grep("pos",rownames(sigma))] <- sg.pos.factor 
	bX[grep("neg",rownames(sigma))] <- sg.neg.factor
	bX = matrix(bX, ncol=1)
	
	# means 
	means<- rep(0,ng); means[grep("pos",rownames(sigma))] <- m.pos; means[grep("neg",rownames(sigma))] <- m.neg
	
	X <-t(mvrnorm(n, means, sigma))
	colnames(X)<-c(1:n) 
	
	# Outcome 
	L <-  t(X) %*% bX
	# error of the model 
	Y <- ifelse(runif(n) < plogis(L), 1, -1)
		
	return(list("x"=X,"y"=Y[,1], "seed"=seed))
}




#####################################################################################

`.create.covariance.matrix` <-
function(sd.pos, sd.neg, sd.bal,
																			ng,nsg, pos.nsg, neg.nsg,
																			corr,corr.factor, 
																			blocks, n.blocks, nsg.block, ng.block){
	# initialise covariance matrix sigma
	sigma<- diag(ng)
	dimnames(sigma)<- list(c(1:ng), c(1:ng))
	
	# 1.  no correlation: sigma:= I 
	if (!corr){
		# first pos, then neg sig. genes, then balanced
		rownames(sigma)<- c( paste("pos",c(1:pos.nsg),sep=""), 
		paste("neg",c(1:neg.nsg),sep=""),
		paste("bal",c(1:(ng-nsg)),sep="") )
	}
	# 2. correlation & no blocks (sign + non-sig genes)
	if (corr & !blocks){
		# correlation, no blocks with non-sg genes
		# --> first "positive block", then "negative" then "balanced"
		#  each block: diag=1, cov(i,j)=cov(j,i):=cor.factor^|i-j|
		
		# first pos, then neg sig. genes, then balanced
		rownames(sigma)<- c( paste("pos",c(1:pos.nsg),sep=""), 
		paste("neg",c(1:neg.nsg),sep=""),
		paste("bal",c(1:(ng-nsg)),sep="") )
		
		.help.cor.block<- function(len, corr.factor){
			# each block: diag=1, cov(i,j)=cov(j,i):=cor.factor^|i-j|
			tmp.sigma<-diag(len)
			diag(tmp.sigma)[grep("pos",rownames(tmp.sigma))]<- sd.pos^2 		
			diag(tmp.sigma)[grep("neg",rownames(tmp.sigma))]<- sd.neg^2 		
			diag(tmp.sigma)[grep("bal",rownames(tmp.sigma))]<- sd.bal^2 		
			
			for(row.i in 1:(len-1))
			 for(col.j in (row.i+1):len)
				tmp.sigma[row.i,col.j] <- tmp.sigma[col.j,row.i] <- corr.factor^abs(row.i-col.j) * sqrt(tmp.sigma[row.i,row.i])  *  sqrt(tmp.sigma[col.j,col.j])
			return(tmp.sigma)
		}
		
		# pos sig genes
		if (pos.nsg>0){
			sigma[1:pos.nsg, 1:pos.nsg]<- .help.cor.block(len=pos.nsg, corr.factor = corr.factor )
			rownames(sigma)[1:pos.nsg]<- paste("pos",c(1:pos.nsg),sep="")
		}
		# neg sig genes
		if (neg.nsg>0){
			sigma[(pos.nsg+1):(pos.nsg+neg.nsg), (pos.nsg+1):(pos.nsg+neg.nsg)]<-  .help.cor.block(len=neg.nsg, corr.factor = corr.factor )
			rownames(sigma)[(pos.nsg+1):(pos.nsg+neg.nsg) ]<- paste("neg",c(1:neg.nsg),sep="")
		}
		# bal genes
		if ((ng-nsg)>0){
			sigma[((pos.nsg+neg.nsg)+1): ng, (pos.nsg+neg.nsg+1): ng]<-  .help.cor.block(len=(ng-nsg), corr.factor = corr.factor )
			rownames(sigma)[((pos.nsg+neg.nsg)+1): ng ]<- paste("bal",c(1:(ng-nsg) ),sep="")
		}
	}
	# 3. correlation & blocks (sign + non-sig genes)
	if (corr & blocks){
		# ith block: first nsg.block sign. genes then (ng.block-nsg.block) bal genes
		# diag = 1, rest = corr.factor
		# rest = diag 1
		
		# cov(i,j)= cor.factor* sd(i)*sd(j)
		# tricky if sd_pos , sd_neg and sd_bal are different 
		
		# 1. first find position of each gene
		
		for (i in 1:n.blocks){
			bl.start<-(i-1)*ng.block+1
			bl.end<- i*ng.block
			# add rownames
			rownames(sigma)[bl.start : bl.end]<- c( paste("sig",c( ((i-1)*nsg.block + 1) : (i*nsg.block) ) ,sep=""), 
														 paste("bal",c( ((i-1)*(ng.block-nsg.block) + 1) : (i*(ng.block-nsg.block)) ),sep="") )
		}
		# do we have some sg left? if not -> 0.
		sg.rest<- max (nsg - nsg.block * n.blocks, 0 )
		if (sg.rest > 0 ) rownames(sigma)[((ng.block * n.blocks) +1): (ng.block * n.blocks+sg.rest)]<- 
																		paste("sig",c( ( nsg.block * n.blocks  + 1) : (nsg) ) ,sep="")
																			
		# do we have some bal left?
		bal.rest<- max( ( (ng-nsg) -  ((ng.block-nsg.block) * n.blocks)  ), 0 )
		if (bal.rest > 0 ) rownames(sigma)[((ng.block * n.blocks)+sg.rest +1): (ng)]<- paste("bal",c( (((ng.block-nsg.block) * n.blocks) + 1) : (ng-nsg) ) ,sep="")
		
		# correct the names for sig genes: sig --> first pos then neg
		rownames(sigma)[grep("sig",rownames(sigma)) ]<- c( paste("pos",c(1:pos.nsg),sep=""), 
		paste("neg",c(1:neg.nsg),sep="") )
		colnames(sigma)<- rownames(sigma)	
		
		
		# 2. fill the blocks	
			
		for (i in 1:n.blocks){
			bl.start<-(i-1)*ng.block+1
			bl.end<- i*ng.block
			block.i<- sigma[bl.start : bl.end, bl.start : bl.end]
			
			# fill the diagonal with sd_pos, sd_neg, sd_bal
			diag(block.i)[grep("pos",rownames(block.i))]<- sd.pos^2 		
			diag(block.i)[grep("neg",rownames(block.i))]<- sd.neg^2 		
			diag(block.i)[grep("bal",rownames(block.i))]<- sd.bal^2 		
			
			# fill the rest
			# cov(i,j)= cor.factor* sd(i)*sd(j)
			for (row.i in 1:(nrow(block.i)-1))
				for (col.j in (row.i+1):ncol(block.i)){
					block.i[row.i, col.j]<- block.i[col.j,row.i]<- corr.factor * sqrt(block.i[row.i,row.i])  *  sqrt(block.i[col.j,col.j])
				}
			sigma[bl.start : bl.end, bl.start : bl.end]<- block.i
		}
		
		# user-fiendly reading: sigma.read<-sigma; sigma.read[sigma.read==0]<- ""; sigma.read[1:50,1:50];     rm(sigma.read)
	}
	# do we need this simplification?
	# # if cov < 10^-3 set it to 0! (memory space reduction!)
	# sigma[sigma< 10^-3]<- 0
	return(sigma)
}
