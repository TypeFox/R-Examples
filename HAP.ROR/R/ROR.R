#############################################################################
##  Functions for merging alleles given tag-SNPs for high polymorphic loci ##
##  Xin Huang  xhuang@fhcrc.org                                            ##
##  An R programme                                                         ##
##  Version 0.01: July 21, 2011                                            ##
##  Version 0.02: Oct  06, 2011                                            ##
##                fixed a bug in grp.list                                  ##
##  Version 0.07  Mar. 14, 2012                                            ##
##                packaged the ror procedure and its table/figure output   ##
##                into functions                                           ##
##  Version 0.08  Mar. 14, 2012                                            ##
##                implement function for Amino Acid scan, and forward selec##
##  Version 0.09  Mar. 31, 2012                                            ##
##                can simulate case-control data use control hap only      ##
##  note: the initial grouped alleles without SNP deleted are mostly diff  ##
##        by one based pair, indicates an important mutant, do not combine ##
##        the first groups, try to go back to the library to search for    ##
##        the actual mutant nuclioded                                      ##     
#############################################################################
#rm(list = ls());  ### clear workspace
library(hash);

grp.list <- function(allele, snp.de = -1){
	#Input:
	#	allele: data.frame of all the alleles
	#	snp.de: a list SNP position to delete

	# delete given SNPs
	if(length(snp.de) == ncol(allele)){
	
		grp.info <- cbind(row.names(allele),1);
	
	}else{
		if(sum(snp.de) > 0){
			al <- allele[,-snp.de];
		}else{
			al <- allele;
		}
		if( (ncol(allele) - length(snp.de) ) == 1){
			al <- as.matrix(al);
			single.name <- hash(unique(al),c(1:length(unique(al))));
			grp.info <- cbind(row.names(allele), values(single.name,keys=al));
		}else{
	# index for alleles with identical SNPs
			dup.inx <- which(duplicated(al) | duplicated(al, fromLast=TRUE));
			du <- al[dup.inx,];
	# sort all the duplicated alleles
			du.sort <- du[do.call(order,du),];
	# generate grouping index
			grp <- cumsum(!duplicated(du.sort));
			grp.info <- cbind(row.names(du.sort),grp);
		}
	}
	grp.info;
}

# function for assign group info to samples after collapsing
collapse <- function(case, ctl, lib, names, snp.de = -1){
	#Input:
	#	case:	case samples: 1st_col=haplotype_1, 2nd_col=haplotype_2
	#	ctl:	control samples: 1st_col=haplotype_1, 2nd_col=haplotype_2
	#	lib:	the tag-SNPs library *.4d with the only alleles appear in sample 
	#	names:	corresponding allele names in the same format as appear in sample
	#	snp.de:	the column position of a list of SNPs to be deleted, default no delete

	grouping <- grp.list(lib, snp.de);
	if(nrow(grouping)>0){
	#build grouping dictionary
	grp.hs <- hash(grouping[,1], grouping[,2]);
	sub.hs <- hash(names[,1], names[,2]);
	values(sub.hs, keys=keys(grp.hs)) <- values(grp.hs, keys=keys(grp.hs));
#	dict0.hs <- hash(lib.sub.names[,2], values(sub.hs));
#	tmp <- cbind(keys(dict0.hs),values(dict0.hs));
#	tmp2 <- tmp[-which(duplicated(tmp[,2])),];
#	tmp3 <- hash(tmp2[,2],tmp2[,1]);
#	tmp <- cbind(tmp,values(tmp3,keys=tmp[,2]));
#	dict.hs <- hash(tmp[,1],tmp[,3]);
	dict.hs <- hash(names[,2], values(sub.hs));
	ctl0 <- cbind(0, as.numeric(values(dict.hs, keys=ctl[,1])), as.numeric(values(dict.hs, keys=ctl[,2])));
	case0 <- cbind(1, as.numeric(values(dict.hs, keys=case[,1])), as.numeric(values(dict.hs, keys=case[,2])));
	grp.sample <- rbind(ctl0, case0);

}else{
	grp.sample <- rbind(as.matrix(cbind(0,ctl)),as.matrix(cbind(1,case)));

}
	grp.sample;
} 


### function for AIC/Deviance calculation given index of deleted SNPs
### the default reference level is the one with most common allelels


AIC <- function(case.sub, ctl.sub, lib.sub, lib.sub.names, deleted.snps, ref="NA"){

	if(length(deleted.snps) == ncol(lib.sub)){
		deleted.snps <- -1;
	}	
	
	grp.sample <- collapse(case.sub, ctl.sub, lib.sub, lib.sub.names, deleted.snps);
###################################################################################
#Haploitype-based test
	sample.alleles <- c(grp.sample[,2],grp.sample[,3]);
	cate.alleles <- unique(sample.alleles);


### setup reference
	if (is.na(ref)){
		reference <- names(which.max(table(sample.alleles)));
	}else{
		grp.tmp <- grp.list(lib.sub, deleted.snps);
		if(nrow(grp.tmp)!=0){
			lib.hs <- hash(lib.sub.names[,1], lib.sub.names[,2]);
			lib.hs.inv <- hash(lib.sub.names[,2], lib.sub.names[,1])
			grp.hs <- hash(grp.tmp[,1], grp.tmp[,2]);
			if( ref %in% values(lib.hs, keys=grp.tmp[,1])){
				reference <- values(grp.hs, keys=values(lib.hs.inv, keys=ref));
			}else{
				reference <- ref;
			}
		}else{
			reference <- ref;		
		}
	}

	reference <- as.numeric(reference);
	dummy <- contr.treatment(cate.alleles, base=which(cate.alleles==reference));
	dummy.hs <- hash(row.names(dummy),0)
	for(j in 1:nrow(dummy)){
 		values(dummy.hs, keys=row.names(dummy)[j]) <- dummy[j,];
	}
		
	allele.factor<- t(sapply(c(1:nrow(grp.sample)),function(x) dummy.hs[[as.character(grp.sample[x,2])]] + dummy.hs[[as.character(grp.sample[x,3])]]));
	if(nrow(allele.factor)==1) allele.factor <- t(allele.factor);
	res <- glm(grp.sample[,1] ~ allele.factor, family=binomial);
###################################################################################
# the likelihood
	prob.pred <- fitted(res);
	log.prob.ctl <- log(prob.pred);
	log.prob.case <- log(1-prob.pred);
	logLK <- sum((1-grp.sample[,1])*log.prob.ctl) + sum(grp.sample[,1]*log.prob.case);
	AIC <- res$aic;
	dev <- res$deviance;
	df <- res$df.residual;
	dev.null <- res$null.deviance;
	df.null <- res$df.null;
	list(logLK=logLK, AIC=AIC, res=res, dev=dev, df=df, dev.null=dev.null, df.null=df.null);
}


#### function of searching for the next grouping given deleted SNPs 
	### input:
	### 	dev.now: aic/dev for current model
        ###     df.now: degree of freedom for current model
	###	rank: numbers of pairs with top similarity scores to be investigate
        ###           if FALSE, then deviance is calculated for the step-wise merger, then option "alpha" and "step" is used
	###	cut: cutoff for similarity score to be consider, default is -1, means all scores above 0
        ###     alpha: family-wise error, used for deviance only
        ###     step: index for how many deletions have been carried so far
	###	all other inputs are the same as function AIC

deletion <- function(lib, lib.names, case.sub, ctl.sub, aic.now, dev.now, df.now, rank=FALSE, cut=-1, delete.snp=-1, ref=NA, alpha=0.05, step=0){
	null <- 0;
	null <- null[-1];
	total.snps <- ncol(lib);
	records <- null;
	grp <- grp.list(lib,delete.snp);
	if(nrow(grp) > 0 & !(-1 %in% delete.snp)){
		grp.hs <- hash(grp[,1], grp[,2]);
		name.hs <- hash(row.names(lib),row.names(lib));
		values(name.hs, keys=keys(grp.hs)) <- values(grp.hs, keys=keys(grp.hs));
		name.new <- values(name.hs, keys=row.names(lib));
		if(sum(duplicated(name.new)) == nrow(lib) - 2 | ncol(lib)-length(delete.snp) == 1){
			lib.sub <- as.matrix(lib[!duplicated(name.new), -delete.snp]);
			row.names(lib.sub) <- names(name.new[!duplicated(name.new)]);
		}else{
			lib.sub <- lib[-which(duplicated(name.new)), -delete.snp];
		}
	}else{
		lib.sub <- lib;
	}

	if(-1 %in% delete.snp){
		delete.snp <- null; # null list
	}

# calculate similarity matrix	
	alleles.n <- nrow(lib.sub);
	alleles.n1 <- alleles.n - 1;
	diff <- matrix(-1,alleles.n,alleles.n);
	for(i in 1:alleles.n1){
		ii <- i+1;
		for (j in ii:alleles.n){
			diff[i,j] <- mean(lib.sub[i,]==lib.sub[j,]);
		}
	}

# if given rank number, search smallest deviance change amount top similar pairs
	if(class(rank)=="numeric"){
		found <- 1;
	# select the top ranked similar pairs	
		diff.sort <- sort(diff, decreasing=TRUE);
		diff.rank <- unique(diff.sort[1:rank]);
	# garantee just using the above diagnoal entries
		diff.rank <- diff.rank[diff.rank > cut];
		max.inx <- null;

		for(i.rank in 1:length(diff.rank)){
			max.inx <- rbind(max.inx, which(diff==diff.rank[i.rank],arr.ind=T));
		}

	#	max.inx <- which(diff==max(diff),arr.ind=T)

		if(nrow(max.inx) > 1){
		#	group the pair which has the smallest AIC after grouping
			aic <- df <- dev <- p.aic <- rep(0,nrow(max.inx));
			for(i in 1:nrow(max.inx)){
				T <- which(row.names(lib) %in% row.names(lib.sub[max.inx[i,],]));
				new.del <- rep(0,total.snps);
				if(nrow(grp) == 0){
					del <- which(lib[T[1],]!=lib[T[2],]);
					new.del[del]<-1;
				}else{
					del <- unique(c(delete.snp, which(lib[T[1],]!=lib[T[2],])));
					new.del[del[-which(del %in% delete.snp)]]<-1;
				}
				new.del <- c(T[1],T[2],new.del);
				tmp <- AIC(case.sub, ctl.sub, lib, lib.names, del, ref=ref);
				if(length(del) == total.snps){
					df[i] <- tmp$df.null;
					dev[i] <- tmp$dev.null;
					aic[i] <- aic.now;
					p.aic[i] <- 1 - pchisq((tmp$dev.null-dev.now), (tmp$df.null - df.now));
					new.del <- c(0, p.aic[i], new.del);
				}else{
					aic[i] <- tmp$AIC;
					df[i] <- tmp$df;
					dev[i] <- tmp$dev;
					p.aic[i] <- 1 - pchisq((tmp$dev-dev.now), (tmp$df - df.now));
					new.del <- c(0, p.aic[i], new.del);
				}

				records <- rbind(records, new.del);
			}
#			aic.min <- which.min(abs(aic-aic.now));
			aic.min <- which.min(aic);
			df.min <- df[aic.min];
			dev.min <- dev[aic.min];
			p.aic.min <- p.aic[aic.min];
			T <- which(row.names(lib) %in% row.names(lib.sub[max.inx[aic.min,],]));
			del <- unique(c(delete.snp, which(lib[T[1],]!=lib[T[2],])));
			records[aic.min,1] <- 1;

		}else{
			if (nrow(grp) == 0){
				T <- which(row.names(lib) %in% row.names(lib.sub[max.inx[1,],]));
				del <- which(lib[T[1],]!=lib[T[2],]);
				new.del <- rep(0,total.snps);
				new.del[del]<-1;
				new.del <- c(T[1],T[2],new.del);
				tmp <- AIC(case.sub, ctl.sub, lib, lib.names, del, ref=ref);
				aic <- tmp$AIC;
				df <- tmp$df;
				dev <- tmp$dev;
				p.aic <- 1 - pchisq((tmp$dev-dev.now), (tmp$df - df.now));
				new.del <- c(1, p.aic, new.del);
				records <- rbind(records, new.del);
			}else{
				T <- which(row.names(lib) %in% row.names(lib.sub[max.inx[1,],]));
				del <- unique(c(delete.snp, which(lib[T[1],]!=lib[T[2],])));
				new.del <- rep(0,total.snps);
				new.del[del[-which(del %in% delete.snp)]]<-1;
				new.del <- c(T[1],T[2],new.del);
				tmp <- AIC(case.sub, ctl.sub, lib, lib.names, del, ref=ref);
				if(length(del) == total.snps){
					aic <- aic.now;
					df <- tmp$df.null;
					dev <- tmp$dev.null;
					p.aic <- 1 - pchisq((tmp$dev.null-dev.now), (tmp$df.null - df.now));
					new.del <- c(0, p.aic, new.del);
					records <- rbind(records, new.del);
				}else{
					aic <- tmp$AIC;
					df <- tmp$df;
					dev <- tmp$dev;
					p.aic <- 1 - pchisq((tmp$dev-dev.now), (tmp$df - df.now));
					new.del <- c(1, p.aic, new.del);
					records <- rbind(records, new.del);
				}
			}
		}
		fin.inx <- which.min(aic);
		aic.fin <- aic[fin.inx];
		df.min <- df[fin.inx];
		dev.min <- dev[fin.inx];
		p.aic.min <- p.aic[fin.inx];
		found <- 1;
# otherwise, the default rank=FALSE, perform an analytical search given the significance level
	}else{
		diff.sort <- sort(diff, decreasing=TRUE);
		diff.rank <- unique(diff.sort);
	# garantee just using the above diagnoal entries
		diff.rank <- diff.rank[diff.rank > cut];
		found <- 0;
		for( i.a in 1:length(diff.rank)){
			record.rank <- null;
			max.inx <- which(diff==diff.rank[i.a],arr.ind=T);
			if(nrow(max.inx) > 1){
		#	group the pair which has the smallest AIC after grouping
				aic <- df <- dev <- rep(0,nrow(max.inx));
				for(i.b in 1:nrow(max.inx)){
					T <- which(row.names(lib) %in% row.names(as.matrix(lib.sub[max.inx[i.b,],])));
					del <- unique(c(delete.snp, which(lib[T[1],]!=lib[T[2],])));
					new.del <- rep(0,total.snps);
					new.del[del[-which(del %in% delete.snp)]]<-1;
					new.del <- c(T[1],T[2],new.del);
					tmp <- AIC(case.sub, ctl.sub, lib, lib.names, del, ref=ref);
					
					if(length(del) == total.snps){
						df[i.b] <- tmp$df.null;
						dev[i.b] <- tmp$dev.null;
						aic[i.b] <- 1 - pchisq((tmp$dev.null-dev.now), (tmp$df.null - df.now));
						new.del <- c(0, aic[i.b], new.del);
					}else{
						df[i.b] <- tmp$df;
						dev[i.b] <- tmp$dev;
						aic[i.b] <- 1 - pchisq((tmp$dev-dev.now), (tmp$df - df.now));
						#cat("deviance difference = ", tmp$dev-dev.now, "\n");
						new.del <- c(0, aic[i.b], new.del);
					}
					record.rank <- rbind(record.rank, new.del);
				}
				if(nrow(rbind(1,records)) != 1){
					if(length(del) == total.snps){
#						mul.n <- nrow(records) + nrow(record.rank) + step;
						mul.n <- nrow(records) + nrow(record.rank);	
					}else{
						mul.n <- nrow(records) + nrow(record.rank);
					}
				}else{
					if(length(del) == total.snps){
#						mul.n <- nrow(record.rank) + step;
						mul.n <- nrow(record.rank);	
					}else{
						mul.n <- nrow(record.rank);
					}
				}
#				mul.n <- 1;
#				for(i.aa in 1:i.a){
#					mul.n <- mul.n + nrow(which(diff==diff.rank[i.aa],arr.ind=T));
				#	mul.n <- mul.n + nrow(which(diff==diff.rank[i.aa],arr.ind=T)) + step;
				#	mul.n <- 1 + i.a;
#				}
				alpha.n <- alpha/mul.n;
				if(max(aic) > alpha.n){
					aic.min <- which.max(aic);
					df.min <- df[aic.min];
					dev.min <- dev[aic.min];
					T <- which(row.names(lib) %in% row.names(lib.sub[max.inx[aic.min,],]));
					del <- unique(c(delete.snp, which(lib[T[1],]!=lib[T[2],])));
					record.rank[aic.min,1] <- 1;
					records <- rbind(records, record.rank);
					found <- 1;
					break;
				}else{
					records <- rbind(records, record.rank);
					next;
				}
			}else{
				if (nrow(grp) == 0){  # may need to revise this and the above when nrow(max.inx) > 1 for the first grouping
					T <- which(row.names(lib) %in% row.names(lib.sub[max.inx[1,],]));
					del.snp0 <- which(lib[T[1],]!=lib[T[2],]);
					del <- del.snp0;
					new.del <- rep(0,total.snps);
					new.del[del]<-1;
					new.del <- c(T[1],T[2],new.del);
					tmp <- AIC(case.sub, ctl.sub, lib, lib.names, del, ref=ref);
					aic <- 1 - pchisq((tmp$dev-dev.now), (tmp$df - df.now));
					#cat("deviance difference = ", tmp$dev-dev.now, "\n");
					new.del <- c(0, aic, new.del);
					record.rank <- rbind(record.rank, new.del);
					mul.n <- 1;
#					for(i.aa in 1:i.a){
#						mul.n <- mul.n + nrow(which(diff==diff.rank[i.aa],arr.ind=T));
						#mul.n <- mul.n + nrow(which(diff==diff.rank[i.aa],arr.ind=T)) + step;
						#mul.n <- 1 + i.a;
#					}
					alpha.n <- alpha/mul.n;
					if(aic > alpha.n){
						df.min <- tmp$df;
						dev.min <- tmp$dev;
						record.rank[1,1]<-1;
						records <- rbind(records, record.rank);
						found <- 1;
						break;
					}else{
						records <- rbind(records, record.rank);
						next;
					}

				}else{
					T <- which(row.names(lib) %in% row.names(as.matrix(lib.sub[max.inx[1,],])));
					del <- unique(c(delete.snp, which(lib[T[1],]!=lib[T[2],])));
					tmp <- AIC(case.sub, ctl.sub, lib, lib.names, del, ref=ref);
					new.del <- rep(0,total.snps);
					new.del[del[-which(del %in% delete.snp)]]<-1;
					new.del <- c(T[1],T[2],new.del);

					if(length(del) == total.snps){
						aic <- 1 - pchisq((tmp$dev.null-dev.now), (tmp$df.null - df.now));
						new.del <- c(0, aic, new.del);
						record.rank <- rbind(record.rank, new.del);
					}else{
						aic <- 1 - pchisq((tmp$dev-dev.now), (tmp$df - df.now));
						#cat("deviance difference = ", tmp$dev-dev.now, "\n");
						new.del <- c(0, aic, new.del);
						record.rank <- rbind(record.rank, new.del);
					}
					if(nrow(rbind(1,records)) != 1){
						if(length(del) == total.snps){
							mul.n <- nrow(records) + nrow(record.rank) + step;
						}else{
							mul.n <- nrow(records) + nrow(record.rank);
						}
					}else{
						if(length(del) == total.snps){
							mul.n <- nrow(record.rank) + step;
						}else{
							mul.n <- nrow(record.rank);
						}
					}
#					for(i.aa in 1:i.a){
#						mul.n <- mul.n + nrow(which(diff==diff.rank[i.aa],arr.ind=T));
						#mul.n <- mul.n + nrow(which(diff==diff.rank[i.aa],arr.ind=T)) + step;
						#mul.n <- 1 + i.a;
#					}
					alpha.n <- alpha/mul.n;
					if(aic > alpha.n){
						df.min <- tmp$df;
						dev.min <- tmp$dev;
						record.rank[1,1]<-1;
						records <- rbind(records, record.rank);
						found <- 1;
						break;
					}else{
						records <- rbind(records, record.rank);
						next;
					}
				}

			}
		}
	
		if(found == 0){
			aic <- 0;
			df.min <- 0;
			dev.min <- 0;
		}
		aic.fin <- max(aic);
		p.aic.min <- aic.fin;
	}


	list(del=del, aic=aic.fin, p.aic=p.aic.min, df=df.min, dev=dev.min, stop=found, record=records);
	
}


HAP.ror <- function(case.sub, ctl.sub, lib.sub, lib.sub.names, alpha=0.01, ref.level=NA, display.proc=TRUE){
# case.sub: case subjects, two columns for two haplotypes
# ctl.sub: control subjects, two columns for two haplotypes
# lib.sub: the alleles library contains allele sequences for those only appear in the case and control samples
# the corresponding names of the alleles
# alpha: significance level
# ref.level: name of the reference allele, "NA" use the most common allele as reference, can also specify allele name, for DRB1, it is "101"
  null <- 0;
  null <- null[-1]; # null list
  deleted.snps <- rep(FALSE, ncol(lib.sub));
	deleted.snps.ls <- null
	AIC.list <- null;
	dev.list <- null;
	df.list <- null;
	records <- null;
	del.ls <- -1;
 
	tmp <- AIC(case.sub, ctl.sub, lib.sub, lib.sub.names, -1, ref=ref.level);
	dev.now <- tmp$dev;
	df.now <- tmp$df;
	aic.now <- tmp$AIC;
	stop <- 1;
	step <- 1;
	while ((ncol(lib.sub) - sum(deleted.snps)) > 0 & stop == 1){
		result <- deletion(lib.sub, lib.sub.names, case.sub, ctl.sub, aic.now, dev.now, df.now, rank=FALSE, cut=-1, del.ls, ref=ref.level, alpha=alpha, step=step);
		del.ls <- result$del;
		deleted.snps <- (c(1:ncol(lib.sub)) %in% del.ls);
		deleted.snps.ls <- rbind(deleted.snps.ls, deleted.snps);
		dev.list <- c(dev.list, (result$dev - dev.now));
		df.list <- c(df.list, (result$df - df.now));
		dev.now <- result$dev;
		df.now <- result$df;
		aic.now <- result$aic;
		stop <- result$stop;
		records <- rbind(records, cbind(step, result$record));
		AIC.list <- c(AIC.list, result$aic);
		step <- step + 1;
		if(display.proc==TRUE){
			cat("deleting:", del.ls, "\n");
		}
	}
	if(sum(deleted.snps)==ncol(lib.sub) & stop == 1){
		significant <- 0;
		grp.result <- grp.list(lib.sub, deleted.snps);
		deleted.snps <- c(1:ncol(lib.sub));
		model.summary <- NA;
	}else{
		significant <- 1;
		fin <- length(dev.list) - 1;
		deleted.snps <- c(1:ncol(lib.sub))[deleted.snps.ls[fin,]>0];
		grp.result <- grp.list(lib.sub, deleted.snps);
		result <- AIC(case.sub, ctl.sub, lib.sub, lib.sub.names, deleted.snps, ref=ref.level);
		model.summary <- summary(result$res);
	}

	list(dev.list=dev.list, AIC.list=AIC.list, df.list=df.list, records=records, deleted.snps.ls=deleted.snps.ls, deleted.snps=deleted.snps, grp.result=grp.result, model.summary=model.summary, significant=significant);
}

# function for output tables and figures related to ROR result
ODS.ror <- function(case.sub, ctl.sub, lib.sub, lib.sub.names, records, dev.list, AIC.list, deleted.snps.ls, proteinf, locus="DRB1*", ref.level="101"){
  null <- 0;
  null <- null[-1]; # null list
############################
# structure variants table #
############################
sv <- cbind(records[,c(1:5)],NA);
colnames <- names(lib.sub);
colnames <- gsub("X\\.","-",colnames, perl=TRUE);
colnames <- gsub("X|\\.[1-9]","",colnames, perl=TRUE);

for(i in 1:nrow(records)){

sv[i,4] <- row.names(lib.sub)[as.numeric(sv[i,4])];
sv[i,5] <- row.names(lib.sub)[as.numeric(sv[i,5])];
sv[i,6] <- paste(colnames[records[,-c(1:5)][i,]>0], collapse=",");

}

row.names(sv)<-c(1:nrow(sv));
sv <- as.data.frame(sv);
names(sv) <- c("step","merged","p-value", "grp1", "grp2", "Amino-Acid");

write.table(sv, file="ROR.sv", row.names=FALSE, sep="\t");

################################################
# survival tables for each amino acid position #
################################################
survived <- cbind(as.matrix(colSums(records[records[,2]==0,-c(1:5)])),0,0,0,1,0,0);
row.names(survived)<- colnames;
tSNP.inx <- c(1:nrow(survived));
sv1 <- as.numeric(levels(sv[,1])[sv[,1]])
sv3 <- as.numeric(levels(sv[,3])[sv[,3]])
for(i in 1:nrow(records)){

	removed <- tSNP.inx[records[,-c(1:5)][i,]>0];
	for(j in removed){
		if(survived[j,2]==0){
			survived[j,2] <- sv1[i];
		}
		survived[j,3] <- sv1[i];
		if(sv3[i] < survived[j,5]){
			survived[j,5] <- sv3[i];
			survived[j,4] <- sv1[i];
		}
		if(sv3[i] > survived[j,7]){
			survived[j,7] <- sv3[i];
			survived[j,6] <- sv1[i];
		}
	}

}
write.table(survived, file="aa.sur", row.names=TRUE, sep="\t");

##########################
# grouping history table #
##########################
#library(hash);
history <- matrix(NA, length(AIC.list), nrow(lib.sub.names)+3);
history[,nrow(lib.sub.names)+2] <- dev.list;
history[,nrow(lib.sub.names)+3] <- AIC.list;
for(i in 1:length(AIC.list)){
	deleted.snps <- c(1:ncol(lib.sub))[deleted.snps.ls[i,]>0];
	grp.i <- grp.list(lib.sub, deleted.snps);
	grp.hs <- hash(grp.i[,1], grp.i[,2]);
	sub.hs <- hash(lib.sub.names[,1], lib.sub.names[,2]);
	values(sub.hs, keys=keys(grp.hs)) <- values(grp.hs, keys=keys(grp.hs));
	history[i,c(1:nrow(lib.sub.names))] <- values(sub.hs, keys=keys(sub.hs));
	history[i,nrow(lib.sub.names)+1] <- sum(sv[,1]==i);
}


history <- cbind(c(1:length(AIC.list)),history);
colnames(history) <- c("step",keys(sub.hs),"#trials","dev","P-value");
write.table(history, file="grp.his", row.names=FALSE, sep="\t");

#################################
# grouping history detail table #
#################################
#library(hash);

history <- matrix(".", length(AIC.list), 3*nrow(lib.sub.names)+3);
name.inx <- c(1:nrow(lib.sub.names))*3-2
history[,3*nrow(lib.sub.names)+2] <- dev.list;
history[,3*nrow(lib.sub.names)+3] <- AIC.list;
for(i in 1:length(AIC.list)){
	deleted.snps <- c(1:ncol(lib.sub))[deleted.snps.ls[i,]>0];
  if(ncol(lib.sub)==length(deleted.snps)){
    next
  }
	grp.i <- grp.list(lib.sub, deleted.snps);
	grp.hs <- hash(grp.i[,1], grp.i[,2]);
	sub.hs <- hash(lib.sub.names[,1], lib.sub.names[,2]);
	values(sub.hs, keys=keys(grp.hs)) <- values(grp.hs, keys=keys(grp.hs));
	result <- AIC(case.sub, ctl.sub, lib.sub, lib.sub.names, deleted.snps, ref=ref.level);
	p.val <- summary(result$res)$coef[, "Pr(>|z|)"];
	coef <- summary(result$res)$coef[, "Estimate"];
### new! fix a bug in the following if else
	if(length(coef)>2){
		gname <- gsub("allele.factor","",names(coef), perl=TRUE);
		gname[1] <- unique(values(sub.hs)[-which(values(sub.hs) %in% gname)]);
	}else if(length(coef)==2){
		gname <- c(NA, NA);
		lib.hs.inv <- hash(lib.sub.names[,2], lib.sub.names[,1]);
		ref.grpn <- values(grp.hs, keys=values(lib.hs.inv, keys=ref.level));
		gname[1] <- ref.grpn;
		gname[2] <- unique(values(sub.hs)[which(values(sub.hs)!=ref.grpn)]);
	}else{
		cat("Warning: unknown error!\n");
	}
	p.hs <- hash(gname, p.val);
	coef.hs <- hash(gname, coef);
	history[i, c(1:nrow(lib.sub.names))*3-2] <- values(sub.hs, keys=keys(sub.hs));
	history[i, c(1:nrow(lib.sub.names))*3-1] <- values(coef.hs, keys=values(sub.hs, keys=keys(sub.hs)));
	history[i, c(1:nrow(lib.sub.names))*3] <- values(p.hs, keys=values(sub.hs, keys=keys(sub.hs)));
	history[i,3*nrow(lib.sub.names)+1] <- sum(sv[,1]==i);
}

history <- cbind(c(1:length(AIC.list)),history);
colnames(history) <- c("step",c(sapply(keys(sub.hs),FUN=function(x) c(x,"coef","p-val"))),"#trials","dev","P-value");
write.table(history, file="grp.detail", row.names=FALSE, sep="\t");

###################################################################
# output Newick format of tree representation of grouping history #
###################################################################
#library(hash);
sv.merge <- records[records[,2]==1,];
newick.ls <- rep(NA,nrow(sv.merge));
leaf.hs <- hash(lib.sub.names[,1],0);
prelen.hs <- hash(lib.sub.names[,1],0);
branch <- rep(1,nrow(sv.merge));
for(i in 1:nrow(sv.merge)){
	sv.merge[i,4] <- row.names(lib.sub)[as.numeric(sv.merge[i,4])];
	sv.merge[i,5] <- row.names(lib.sub)[as.numeric(sv.merge[i,5])];
	if(values(leaf.hs,keys=sv.merge[i,4])==0 && values(leaf.hs,keys=sv.merge[i,5])==0){
#		newick.ls[i] <- paste("(",gsub(":","",sv.merge[i,4], perl=TRUE),",",gsub(":","",sv.merge[i,5], perl=TRUE),")",sv.merge[i,1], sep="");
		newick.ls[i] <- paste("(",gsub(":","",sv.merge[i,4], perl=TRUE),":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,4]),",",gsub(":","",sv.merge[i,5], perl=TRUE),":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,5]),")",sv.merge[i,1], sep="");
		values(prelen.hs,keys=sv.merge[i,4]) <- sum(deleted.snps.ls[i,]==1);
		values(prelen.hs,keys=sv.merge[i,5]) <- sum(deleted.snps.ls[i,]==1);
		values(leaf.hs,keys=sv.merge[i,4]) <- i;
		values(leaf.hs,keys=sv.merge[i,5]) <- i;
	}else if(values(leaf.hs,keys=sv.merge[i,4])!=0 && values(leaf.hs,keys=sv.merge[i,5])==0){
#		newick.ls[i] <- paste("(",newick.ls[values(leaf.hs,keys=sv.merge[i,4])],",",gsub(":","",sv.merge[i,5], perl=TRUE),")",sv.merge[i,1], sep="");
		newick.ls[i] <- paste("(",newick.ls[values(leaf.hs,keys=sv.merge[i,4])],":", sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,4]),",",gsub(":","",sv.merge[i,5], perl=TRUE),":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,5]),")",sv.merge[i,1], sep="");
		grp.old <- names(which(values(leaf.hs)==values(leaf.hs,keys=sv.merge[i,4])));
		values(prelen.hs,keys=c(grp.old,sv.merge[i,5])) <- sum(deleted.snps.ls[i,]==1);
		branch[values(leaf.hs,keys=grp.old)] <- 0;
		values(leaf.hs,keys=grp.old) <- i;
		values(leaf.hs,keys=sv.merge[i,5]) <- i;
	}else if(values(leaf.hs,keys=sv.merge[i,4])==0 && values(leaf.hs,keys=sv.merge[i,5])!=0){
#		newick.ls[i] <- paste("(",gsub(":","",sv.merge[i,4], perl=TRUE),",",newick.ls[values(leaf.hs,keys=sv.merge[i,5])],")",sv.merge[i,1], sep="");
		newick.ls[i] <- paste("(",gsub(":","",sv.merge[i,4], perl=TRUE),":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,4]),",",newick.ls[values(leaf.hs,keys=sv.merge[i,5])],":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,5]),")",sv.merge[i,1], sep="");
		grp.old <- names(which(values(leaf.hs)==values(leaf.hs,keys=sv.merge[i,5])));
		values(prelen.hs,keys=c(grp.old,sv.merge[i,4])) <- sum(deleted.snps.ls[i,]==1);
		branch[values(leaf.hs,keys=grp.old)] <- 0;
		values(leaf.hs,keys=sv.merge[i,4]) <- i;
		values(leaf.hs,keys=grp.old) <- i;
	}else if(values(leaf.hs,keys=sv.merge[i,4])!=0 && values(leaf.hs,keys=sv.merge[i,5])!=0){
#		newick.ls[i] <- paste("(",newick.ls[values(leaf.hs,keys=sv.merge[i,4])],",",newick.ls[values(leaf.hs,keys=sv.merge[i,5])],")",sv.merge[i,1], sep="");
		newick.ls[i] <- paste("(",newick.ls[values(leaf.hs,keys=sv.merge[i,4])],":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,4]),",",newick.ls[values(leaf.hs,keys=sv.merge[i,5])],":",sum(deleted.snps.ls[i,]==1)-values(prelen.hs,keys=sv.merge[i,5]),")",sv.merge[i,1], sep="");
		grp1.old <- names(which(values(leaf.hs)==values(leaf.hs,keys=sv.merge[i,4])));
		grp2.old <- names(which(values(leaf.hs)==values(leaf.hs,keys=sv.merge[i,5])));
		values(prelen.hs,keys=c(grp1.old,grp2.old)) <- sum(deleted.snps.ls[i,]==1);
		branch[values(leaf.hs,keys=grp1.old)] <- 0;
		branch[values(leaf.hs,keys=grp2.old)] <- 0;
		values(leaf.hs,keys=grp1.old) <- i;
		values(leaf.hs,keys=grp2.old) <- i;
	}
}


#newick.ls[branch==1];
single <- names(which(values(leaf.hs)==0));
n.remain <- ncol(lib.sub) - sum(deleted.snps.ls[nrow(sv.merge),]==1);
branch.inx <- which(branch==1);
newick.fin <- "(";
for(i in branch.inx){
	newick.fin <- paste(newick.fin, newick.ls[i],":",ncol(lib.sub)-sum(deleted.snps.ls[i,]==1),",", sep="");
}
if(length(c(single,1))!=1){
	for(i in 1:length(single)){
		newick.fin <- paste(newick.fin, gsub(":","",single[i], perl=TRUE),":",ncol(lib.sub),",", sep="");
	}
}
substr(newick.fin,nchar(newick.fin),nchar(newick.fin)) <- ")";
newick.fin <- paste(newick.fin, ";", sep="");
cat(newick.fin, file="newick");



###############################
# draw recursive process tree #
###############################
library(ape);

rec.tr <- read.tree("newick");
p.fin <- history[nrow(sv.merge), c(1:nrow(lib.sub.names))*3]
names(p.fin) <- gsub(":","",names(history[,-1][nrow(sv.merge), c(1:nrow(lib.sub.names))*3-2]), perl=TRUE);
p.hs <- hash(names(p.fin), p.fin);
rgb.palette1 <- colorRampPalette(c("red","orange"), space = "rgb");
rgb.palette2 <- colorRampPalette(c("blue","cyan"), space = "rgb");

p.pos <- sort(as.numeric(unique(p.fin[p.fin>0])),decreasing=TRUE);
p.neg <- sort(as.numeric(unique(p.fin[p.fin<0])),decreasing=FALSE);
p.col.hs <- hash(c(p.pos, p.neg), c(rgb.palette1(length(p.pos)), rgb.palette2(length(p.neg))));
values(p.hs, keys=keys(p.hs)) <- values(p.col.hs, keys=values(p.hs))

plot(rec.tr, type="p",show.node.label=TRUE, root.edge=TRUE, label.offset=1, tip.color=values(p.hs, keys=rec.tr$tip), edge.color="blue", font=1, cex=0.8);
nodelabels(pch=21, frame = "c", bg = "tomato", cex=0.6);


########################################################################
# output significant AA and allele list for searching protein sequence #
########################################################################
  
remain.ls <- survived[survived[,3]==length(AIC.list),];
remaining <- unique(row.names(remain.ls[order(remain.ls[,1],-remain.ls[,7],decreasing=TRUE),]));
names.rv.hs <- hash(names(p.fin),names(history[,-1][nrow(sv.merge), c(1:nrow(lib.sub.names))*3-2]));
tr.ord <- values(names.rv.hs, keys=rec.tr$tip);

#write.table(tr.ord[length(tr.ord):1], file="DRB1.allele.order", row.names=FALSE, col.names=FALSE, sep="\t");
#write.table(remaining, file="DRB1.remain", row.names=FALSE, col.names=FALSE, sep="\t");
seqs <- proteinf
row.names(seqs) <- paste(locus,seqs[,1],sep="");
pro.name <- paste(locus,seqs[,2],sep="");
pro.aa <- names(seqs);
pro.aa <- gsub("X\\.","-",pro.aa, perl=TRUE);
pro.aa <- gsub("X|\\.[1-9]","",pro.aa, perl=TRUE);

remain.aa <- seqs[,-c(1:3)][(pro.name %in% tr.ord[length(tr.ord):1]), (pro.aa[-c(1:3)] %in% remaining)];
names(remain.aa) <- gsub("X","",names(remain.aa), perl=TRUE);
names(remain.aa) <- gsub("\\.","-",names(remain.aa), perl=TRUE);
remain.aa.4digit <- pro.name[pro.name %in% tr.ord[length(tr.ord):1]];
remain.aa <- remain.aa[,sapply(remaining, FUN=function(x) which(names(remain.aa)==x))];
remain.aa <- cbind(row.names(remain.aa),remain.aa.4digit,remain.aa);
names(remain.aa)[1] <- "6digit";
names(remain.aa)[2] <- "4digit";
remain.fin <- null;
for(i in 1:length(tr.ord)){
	remain.fin <- rbind(remain.fin, remain.aa[which(remain.aa.4digit==tr.ord[length(tr.ord):1][i]),]);
}

write.table(remain.fin, file="aa6digit", row.names=FALSE, col.names=TRUE, sep="\t");
write.table(remain.fin[-which(duplicated(remain.fin[,2])),-1], file="aa4digit", row.names=FALSE, col.names=TRUE, sep="\t");

}

## Amino Acid scan, can be used to perform forward amino acid selection
#AAscan <- function(case.sub, ctl.sub, lib.sub, lib.sub.names, ref.level, alpha=0.05, retain=-1, dev.now=0, df.now=0){
# retain: a list of positions to retain/condition on, the position number is the order of the amino acid in the columns of lib.sub, -1 means no retain
# dev.now: deviance for the previous model with selected amino acids, used for conditioning, used only when retain != -1
# df.now: deviance for the previous model
#	nAA <- ncol(lib.sub);
#	p.bon <- -log(alpha/nAA, 10);
#	p.vals <- matrix(NA, nrow=1, ncol=nAA);	
#	devs <- dfs <- aics <- rep(NA, nAA);
#	colnames(p.vals) <- colnames(lib.sub);
#	for(i in 1:nAA){
#		if(-1 %in% retain){
#			deleted.snps <- c(1:nAA)[-i]
#			grp.result <- grp.list(lib.sub, deleted.snps);
#			result <- AIC(case.sub, ctl.sub, lib.sub, lib.sub.names, deleted.snps, ref=ref.level);
#			p <- pchisq((result$dev.null - result$dev), (result$df.null - result$df), lower.tail=FALSE, log.p=TRUE);	
#			p.vals[1,i] <- -log(exp(1),10)*p;
#			devs[i] <- result$dev;
#			aics[i] <- result$AIC;
#			dfs[i] <- result$df;
#			
#		}else{
#			if(!(i %in% retain)){
#				deleted.snps <- c(1:nAA)[-c(i,retain)];
#				grp.result <- grp.list(lib.sub, deleted.snps);
#				result <- AIC(case.sub, ctl.sub, lib.sub, lib.sub.names, deleted.snps, ref=ref.level);
#				if(df.now == result$df){
#					diff <- 1;
#				}else{
#					diff <- df.now - result$df;
#				}
#				p <- pchisq((dev.now - result$dev + 10e-8), diff, lower.tail=FALSE, log.p=TRUE);	
#				p.vals[1,i] <- -log(exp(1),10)*p;
#				devs[i] <- result$dev;
#				aics[i] <- result$AIC;
#				dfs[i] <- result$df;
#			}else{
#				p.vals[1,i] <- NA;
#			}
#		}
#	}
#	if(-1 %in% retain){
#		if(sum(p.vals[1,] >= p.bon)>0){
#			significant <- 1;
#		}else{
#			significant <- 0;
#		}
#	}else{
#		if(sum(p.vals[1,-retain] >= p.bon)>0){
#			significant <- 1;
#		}else{
#			significant <- 0;
#		}
#	}
#	list(p.vals=p.vals, significant=significant, devs=devs, dfs=dfs, aics=aics);
#}

### simulating case-control data given causal amino acids/haplotype alleles
cc.sim <- function(n.ctrl, n.case, beta0, beta1, case.sub, ctl.sub, lib.sub, lib.sub.names, risk.type="AA", risk.inx=2, risk.names=c("301", "302"), min.count=10, ctl.only=FALSE){
##n.ctrl: number of control samples desired to generate
##n.case: number of case samples desired to generate
##beta0: the coefficient of intercept for logistic model
##beta1: the coefficient of the causal SNP for logistic model
##case.sub, ctl.sub: two column matrix of sample HLA alleles
##lib.sub, lib.sub.names: alleles sequence library and its names
##risk.type="AA": simulated from given amino acid position as shown in matrix lib.sub, use risk.inx to input position
##risk.type="allele":simulated from given risk alleles, use risk.names=c("301", "302") to specified those alleles
##min.count: use to calculate the warning if the selected alleles have too small frequencies.


	if(ctl.only==FALSE){
		min.count.freq <- 10/(2*nrow(case.sub)+2*nrow(ctl.sub));
		allele.pool <- c(case.sub[,1], case.sub[,2], ctl.sub[,1], ctl.sub[,2]);
		all.alleles <- table(allele.pool)/(2*nrow(case.sub)+2*nrow(ctl.sub));
		alleles.freq.hs <- hash(names(all.alleles), all.alleles);
		risk.alleles.hs <- hash(names(all.alleles), 0);
	}else{
		min.count.freq <- 10/(2*nrow(ctl.sub));
		allele.pool <- c(ctl.sub[,1], ctl.sub[,2]);
		all.alleles <- table(allele.pool)/(2*nrow(ctl.sub));
		alleles.freq.hs <- hash(names(all.alleles), all.alleles);
		risk.alleles.hs <- hash(names(all.alleles), 0);		
	}
	if(risk.type=="AA"){
		protypes <- table(lib.sub[,risk.inx]);
		risk.type <- names(protypes)[which.min(protypes)]
		risk.alleles.inx <- which(lib.sub[,risk.inx]==risk.type);
		risk.alleles.names <- lib.sub.names[risk.alleles.inx,];
		if(is.matrix(risk.alleles.names)){
			risk.names <- risk.alleles.names[,2];
			cat("risk haplotypes:", risk.names, "\n");
		}else{
			risk.names <- risk.alleles.names[2];
			cat("risk haplotypes:", risk.names, "\n");
		}
		values(risk.alleles.hs, keys=risk.names) <- 1;
		select.freq <- values(alleles.freq.hs, keys=risk.names);
		if(sum(select.freq > min.count.freq) == 0){ 
			cat("Warning: causal alleles have too small allele frequencies, should you change another amino acid position?\n");
		}
	}else if(risk.type=="allele"){
		values(risk.alleles.hs, keys=risk.names) <- 1;
		select.freq <- values(alleles.freq.hs, keys=risk.names);
	}else{
		cat("risk type not accept!\n");
	}
	i.ctrl <- i.case <- 1;
	sim.sample <- NULL;
	phenotype <- NULL;
	while(i.ctrl <= n.ctrl | i.case <= n.case){
		sa.inx <- sample(length(allele.pool),2);
		hapCount <- sum(values(risk.alleles.hs, keys=allele.pool[sa.inx]));
		tmp <- exp(beta0+beta1*hapCount);
		p <- tmp/(1+tmp);
		d <- rbinom(1,1,p);
		if(d==1 & i.case <= n.case){
			sim.sample <- rbind(sim.sample, allele.pool[sa.inx]);
			phenotype <- c(phenotype,d);
			i.case <- i.case + 1;
		}
		if(d==0 & i.ctrl <= n.ctrl){
			sim.sample <- rbind(sim.sample, allele.pool[sa.inx]);
			phenotype <- c(phenotype,d);
			i.ctrl <- i.ctrl + 1;
		}
	}

	list(y=phenotype, x=sim.sample, risk.names=risk.names, select.freq=select.freq);

}


