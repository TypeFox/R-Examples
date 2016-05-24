#####################################################################
#
# peaksummary.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Karl W Broman wrote the nphe function part of  R/qtl package
# licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: peaksummary
#
######################################################################

######################################################################
#
# peaksummary: Print and draw summary information about QTL.
#
######################################################################

`peaksummary` <-
function(peak.array,cross,exc=data.frame(inf=0,sup=0,chr=NA),graph=FALSE, ...)
{

	require(qtl)	

	if ( ! any(class(peak.array) == "peak.array" ))
		stop('peak.array should have class \"peak.array\". ')
	if ( length(class(cross)) < 2 || class(cross)[2] != "cross")
		stop("Input should have class \"cross\".")


	if( class(exc) != 'data.frame' & any(!(names(exc) %in% c('inf','sup','chr'))) )
		stop("exc should a dataframe with colomns 'inf', 'sup', 'chromosome' ")
	if( nrow(exc) != 1) stop("exc$inf, exc$sup and exc$chr should be vectors of length 1 \n")
	if( !is.numeric(exc$inf) & !is.numeric(exc$sup) ) stop("exc should contain numeric vector >= 0\n")
	if( exc$inf<0 | exc$sup<0 ) stop("exc should contain numeric vector >= 0\n")
	if( !is.logical(graph)) stop("graph should be logical: TRUE or FALSE")

	ntrait_t <- nphe(cross)-1

	n.na <- !is.na(peak.array$mname.peak)

	if( !is.na(exc$chr[1])){
		if( !any(names(peak.array) %in% "peak.bp" ))
			cat("no physical position column defined in peak.array peak.array: parameter exc is not taking account.\n")
		else { for( i in seq(nrow(exc))){
				bool <- peak.array$peak.bp > exc$inf[i] & peak.array$peak.bp < exc$sup[i] & peak.array$chr == exc$chr[i] 
				n.na <- n.na & !bool
			}
		}
	}

	# print(any(n.na))
	

	#NB DE TRAIT AFFECTE PAR UN eQTL
	ntrait_c <- length(unique(peak.array$trait[n.na]))

	#TAUX DE TRAIT AFFECTE PAR UN eQTL
	rt <- ntrait_c*100/ntrait_t

	#NB D'eQTL DETECTE
	nqtl <- as.numeric(summary(n.na)[3])

	#RATIO NB eQTL PAR TRAIT TOTAUX
	rqt <- nqtl*100/ntrait_t

	#RATIO NB eQTL PAR TRAIT AFFECTE
	rqc <- nqtl/ntrait_c

	res_1 <- list(data.frame(	trait=ntrait_t,
					controlled_trait=ntrait_c,
					detected_QTL=nqtl,
					rate_QTL_by_trait=rqt,
					average_number_QTL_by_controlled_trait=rqc
			))
	attributes(res_1)$names <- paste('global_result')

	res <- res_1

	#NB DE TRAIT CONTROLLE PAR 1,2,3 ET PLUS eQTL
	temp <- summary(as.factor(summary(as.factor(peak.array$trait[n.na]),maxsum=nrow(peak.array))))

	more <- NA
	for( i in 1:length(temp)){
		if( i >= 4) more <- c(more,temp[i])
		else next
	}

	if ( length(more) > 1 ) more <- c(temp[1:3],more[-1]) else more<-temp[1:length(temp)]
	res_2 <- list(trait_by_nb_of_eQTL=temp)

	res<- c(res,res_2)

	if(graph && !is.na(more) )
		 pie(	more,
			main='Distribution of trait affected by\none,two,three or more eQTLs',
			labels=c('single eQTL','two eQTLs','three eQTLs','more')
		)

	#NB D'eQTL PAR CHR
	qtlchr <- summary(factor(peak.array$chr[n.na]))

	res_3 <- list( QTL_distribution_by_chr=qtlchr)

	res<- c(res,res_3)

	if(graph) pie(	qtlchr,
			main='Distribution of eQTLs on chromosomes',
			)

	#NB D'eQTL CIS/TRANS
	if( !any(names(peak.array) %in% "type"))
		cat("no eQTL type column defined in peak.array\n")
	else {
		ntype <- summary(factor(peak.array$type[n.na]))
		if(graph) barplot(ntype,main="eQTL type")

		#NB d'EQTL CIS/TRANS/TRANSCHR  PAR CHR
		ntype_by_chr <- table(peak.array$type[n.na],peak.array$chr[n.na])

		if(graph) barplot(ntype_by_chr,legend.text=c("cis","trans"),main="proportion of cis/trans-eqtl by chromosome")

		#RATIO TYPE D'eQTL PAR CHR
		rtype_by_chr <- prop.table(ntype_by_chr,1)

		res_4 <- list(list(	cistrans=ntype,
					ncistrans_by_chr=ntype_by_chr,
					rate_cistrans_by_chr=rtype_by_chr
				))
		attributes(res_4)$names<- paste('types_summary')

		res <- c(res,res_4)

		if(graph) {
			pie(rtype_by_chr[1,],legend.text=rtype_by_chr[1,],main='Distribution of cis-eQTL by chromosome')
			pie(rtype_by_chr[2,],legend.text=rtype_by_chr[2,],main='Distribution of trans-eQTL by chromosome')
		}
	}

	if( !any(names(peak.array) %in% "additive.effect"))
		cat("no eQTL type column defined in peak.array\n")
	else{
		#EFFET ADDITIF (Bay - Sha)
		peak.array$additive.effect <- as.numeric(peak.array$additive.effect)
	
		#TOTAL
		nbAT <- length(peak.array$additive.effect[n.na]) #c debile == nb de qtl
		minAT <- min(peak.array$additive.effect[n.na])
		maxAT <- max(peak.array$additive.effect[n.na])
		medAT <- median(peak.array$additive.effect[n.na])
		meanAT <- mean(peak.array$additive.effect[n.na])

		#Bay-0 DOMINANT
		B <- peak.array$additive.effect > 0
		nbB <- length(peak.array$additive.effect[B & n.na])
		minB <- min(peak.array$additive.effect[B & n.na])
		maxB <- max(peak.array$additive.effect[B & n.na])
		medB <- median(peak.array$additive.effect[B & n.na])
		meanB <- mean(peak.array$additive.effect[B & n.na])
		varB <- var(peak.array$additive.effect[B & n.na])

		#Sha DOMINANT
		A <- peak.array$additive.effect < 0
		nbA <- length(peak.array$additive.effect[A & n.na])
		minA <- min(peak.array$additive.effect[A & n.na])
		maxA <- max(peak.array$additive.effect[A & n.na])
		medA <- median(peak.array$additive.effect[A & n.na])
		meanA <- mean(peak.array$additive.effect[A & n.na])
		varA <- var(peak.array$additive.effect[A & n.na])

		if(graph)
			boxplot(abs(peak.array$additive.effect[B & n.na]),abs(peak.array$additive.effect[A & n.na]),names=c("Allele 2","Allele 1"),col="gray",horizontal=TRUE,main='Additive effect')

		#equal
		N <- peak.array$additive.effect == 0
		nbN <- length(peak.array$additive.effect[N & n.na])

		res_5 <- list(list( 	total=data.frame(min=minAT,max=maxAT,median=medAT,mean=meanAT),
					allele_2=data.frame(min=minB,max=maxB,median=medB,mean=meanB,var=varB),	
					allele_1=data.frame(min=minA,max=maxA,median=medA,mean=meanA,var=varA),
					Null=nbN
			))

		attributes(res_5)$names <- paste('additive.effect')

		res <- c(res,res_5)

		#PROPORTION ADDITIVE_EFFECT
		temp <- c(nbN,nbA,nbB)
		names(temp) <- c("NULL","Allele 1","Allele 2")
		if(graph) pie(temp,main="global proportion of allelic contribution")

		if (all(names(peak.array) %in% c("additive.effect","type"))){

			#EFFET ADDITIF PAR TYPE DE QTL
			#TOTAL
			meanAT_type <- tapply(peak.array$additive.effect[n.na],peak.array$type[n.na],mean)
			minAT_type <- tapply(peak.array$additive.effect[n.na],peak.array$type[n.na],min)
			maxAT_type <- tapply(peak.array$additive.effect[n.na],peak.array$type[n.na],max)
			medAT_type <- tapply(peak.array$additive.effect[n.na],peak.array$type[n.na],median)
			varAT_type <- tapply(peak.array$additive.effect[n.na],peak.array$type[n.na],var)
			#Bay-0 DOMINANT
			nbB_type <- tapply(peak.array$additive.effect[B & n.na],peak.array$type[B & n.na],length)
			meanB_type <- tapply(peak.array$additive.effect[B & n.na],peak.array$type[B & n.na],mean)
			minB_type <- tapply(peak.array$additive.effect[B & n.na],peak.array$type[B & n.na],min)
			maxB_type <- tapply(peak.array$additive.effect[B & n.na],peak.array$type[B & n.na],max)
			medB_type <- tapply(peak.array$additive.effect[B & n.na],peak.array$type[B & n.na],median)
			varB_type <- tapply(peak.array$additive.effect[B & n.na],peak.array$type[B & n.na],var)
			#Sha DOMINANT
			nbA_type <- tapply(peak.array$additive.effect[A & n.na],peak.array$type[A & n.na],length)
			meanA_type <- tapply(peak.array$additive.effect[A & n.na],peak.array$type[A & n.na],mean)
			minA_type <- tapply(peak.array$additive.effect[A & n.na],peak.array$type[A & n.na],min)
			maxA_type <- tapply(peak.array$additive.effect[A & n.na],peak.array$type[A & n.na],max)
			medA_type <- tapply(peak.array$additive.effect[A & n.na],peak.array$type[A & n.na],median)
			varA_type <- tapply(peak.array$additive.effect[A & n.na],peak.array$type[A & n.na],var)

			res_6 <- list(list( 	total=data.frame(min=minAT_type,max=maxAT_type,median=medAT_type,mean=meanAT_type,var=varAT_type),
						Allele_2=data.frame(min=minB_type,max=maxB_type,median=medB_type,mean=meanB_type,var=varB_type),	
						Allele_1=data.frame(min=minA_type,max=maxA_type,median=medA_type,mean=meanA_type,var=varA_type)
				))
			attributes(res_6)$names <- paste('additive_effect_by_eQTL_type')

			res <- c(res,res_6)
		}

		temp <- peak.array$additive_effect[n.na]
		#temp[B & n.na] <- "Bay-0"
		#temp[S & n.na] <- "Sha"
		#temp[N & n.na] <- "Null"
	
		#if(graph)
			#barplot(table(temp,peak.array$type[n.na],1))

		#EFFET ADDITIF PAR CHR
		#TOTAL
		meanAT_chr <- tapply(peak.array$additive.effect[n.na],peak.array$chr[n.na],mean)
		minAT_chr <- tapply(peak.array$additive.effect[n.na],peak.array$chr[n.na],min)
		maxAT_chr <- tapply(peak.array$additive.effect[n.na],peak.array$chr[n.na],max)
		medAT_chr <- tapply(peak.array$additive.effect[n.na],peak.array$chr[n.na],median)
		varAT_chr <- tapply(peak.array$additive.effect[n.na],peak.array$chr[n.na],var)
		#Bay-0 DOMINANT = alleles 2
		nbB_chr <- tapply(peak.array$additive.effect[B & n.na],peak.array$chr[B & n.na],length)
		meanB_chr <- tapply(peak.array$additive.effect[B & n.na],peak.array$chr[B & n.na],mean)
		minB_chr <- tapply(peak.array$additive.effect[B & n.na],peak.array$chr[B & n.na],min)
		maxB_chr <- tapply(peak.array$additive.effect[B & n.na],peak.array$chr[B & n.na],max)
		medB_chr <- tapply(peak.array$additive.effect[B & n.na],peak.array$chr[B & n.na],median)
		varB_chr <- tapply(peak.array$additive.effect[B & n.na],peak.array$chr[B & n.na],var)
		#Sha DOMINANT = alleles 1
		nbA_chr <- tapply(peak.array$additive.effect[A & n.na],peak.array$chr[A & n.na],length)
		meanA_chr <- tapply(peak.array$additive.effect[A & n.na],peak.array$chr[A & n.na],mean)
		minA_chr <- tapply(peak.array$additive.effect[A & n.na],peak.array$chr[A & n.na],min)
		maxA_chr <- tapply(peak.array$additive.effect[A & n.na],peak.array$chr[A & n.na],max)
		medA_chr <- tapply(peak.array$additive.effect[A & n.na],peak.array$chr[A & n.na],median)
		varA_chr <- tapply(peak.array$additive.effect[A & n.na],peak.array$chr[A & n.na],var)

		res_7 <- list(list( 	total=data.frame(min=minAT_chr,max=maxAT_chr,median=medAT_chr,mean=meanAT_chr,var=varAT_chr),
					allele_2=data.frame(min=minB_chr,max=maxB_chr,median=medB_chr,mean=meanB_chr,var=varB_chr),	
					allele_1=data.frame(min=minA_chr,max=maxA_chr,median=medA_chr,mean=meanA_chr,var=varA_chr)
			))
		attributes(res_7)$names <- paste('additive_effect_by_chr')

		res <- c(res,res_7)
	}

	if(graph){
		if( any(names(peak.array) %in% "Rsq")){
			hist(as.numeric(as.vector(peak.array$Rsq)),nclass=50,col='gray',main="R square distribution TOTAL/CIS",xlab="R square value")
			hist(as.numeric(as.vector(peak.array$Rsq[peak.array$type %in% 'cis'])),nclass=50,col='white',add=TRUE)
			hist(as.numeric(as.vector(peak.array$Rsq)),nclass=50,col='gray',main="R square distribution TOTAL/TRANS",xlab="R square value")
			hist(as.numeric(as.vector(peak.array$Rsq[peak.array$type %in% 'trans'])),nclass=50,col='white',add=TRUE)
			hist(as.numeric(as.vector(peak.array$Rsq[peak.array$type %in% 'trans'])),nclass=50,col='gray',main="R square distribution TRANS/CIS",xlab="R square value")
			hist(as.numeric(as.vector(peak.array$Rsq[peak.array$type %in% 'cis'])),nclass=50,col='white',add=TRUE)
		}
	}

	return(as.list(res))
}

