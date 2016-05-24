ReCiPa <-
function(file, new_name, max_overlap, min_overlap){

no_col	<- max(count.fields(file, sep = "\t"))
dtbse	<- read.table(file, header=FALSE, sep="\t", fill=TRUE, col.names=1:no_col)
dtbse	<- as.matrix(dtbse)
dtbse[,1] <- gsub(" ", "_", dtbse[,1])

########################################################################################
##
## 	NoP: number of pathways in database
## 	sgc: starting (first) gene column
##
## WARNING: NUMBER "3" in sgc MARKS THE STARTING GENE COLUMN, 
## 		MUST BE CHANGED IF THE FILE STRUCTURE IS DIFFERENT
##	
## NUMBER AND % OF GENES FROM PATHWAY 'i' IN COMMON WITH PATHWAY 'j'
## 	Ng		: number of genes of each pathway 
## 	nocg_mat	: number of common genes matrix (pathway 'i' ---> pathway 'j')  
## 	perc_mat	: % of common genes (relative to pathway 'i')
## 	arr		: indexes (pair of pathways) that meet the condition (> coff)
## 	maxNg		: maximum number of genes in some pathway
##	
##

NoP	<- length(dtbse[,1])
sgc	<- 3

nocg_mat <- matrix(NA, NoP, NoP)
for(i in 1:NoP){
	for(j in 1:NoP){
		if(length(which(dtbse[i,]=="")) > 0)
			a <- dtbse[i,sgc:length(dtbse[i,-c(which(dtbse[i,]==""))])]
		else a <- dtbse[i,sgc:length(dtbse[i,])]
		if(length(which(dtbse[j,]=="")) > 0)
			b <- dtbse[j,sgc:length(dtbse[j,-c(which(dtbse[j,]==""))])]
		else b <- dtbse[j,sgc:length(dtbse[j,])]
		nocg <- length(a[a %in% b])
		nocg_mat[i,j] <- nocg
	}
}
Ng		<- diag(nocg_mat)
maxNg		<- max(Ng)
diag(nocg_mat)	<- 0
perc_mat	<- nocg_mat / Ng
arr		<- which(perc_mat > max_overlap, arr.ind=TRUE)
rev.arr 	<- cbind(arr[,2],arr[,1])
arr 		<- arr[which(perc_mat[rev.arr] > min_overlap),]
cat("REPETITIONS:", "\n")

########################################################################################
##
## BIG LOOP, REPEATED WHILE THERE ARE INDEXES THAT MEET THE CONDITION 
## k : number of passings (loops)
##
	
k <- 0
while(length(arr) > 0){

	################################################################################
	##
	## l1  : genes of pathway 'i'
	## l2  : genes of pathway 'j'
	## wm1 : matrix where each row has the indexes of pathways 'i' and 'j', 
	## 		 the number of genes 'i' and 'j', % and number of common genes
	##
	
	arr 	<- check.null.dim(arr)
	wm1 <- cbind(arr,Ng[arr[,1]],Ng[arr[,2]],perc_mat[arr], nocg_mat[arr])

	################################################################################
	##
	## ELIMINATE THE REPEATED PAIRS [i,j] = [j,i] IN wm1
	##

	vec <- NULL
	for(i in 1:length(wm1[,1])){
		r1 <- wm1[i,1]
		c1 <- wm1[i,2]
		for(j in i:length(wm1[,1])){
			r2 <- wm1[j,1]
			c2 <- wm1[j,2]
			if(r1==c2 && c1==r2){
				if(wm1[i,5] < wm1[j,5])
					vec <- c(vec,i)
				else vec <- c(vec,j)
			}
		}
	}
	if(length(vec)>0)
		wm1 <- wm1[-vec,]
	colnames(wm1) <- NULL

	################################################################################
	##
	## VECTOR TO SAVE ROWS (FROM dtbse) THAT SHOULD BE DELETED
	##
			
	row2del <- c()

	################################################################################
	##
	## DIFFERENT CASES:
	## 1) % = 1 AND PATHWAY 'j' IS REPEATED IN wm1
	##
	
	wm1	<- check.null.dim(wm1)
	wm2	<- wm1[which(wm1[,5]==1),]
	wm2	<- check.null.dim(wm2)
	ux	<- unique(wm2[,1])
	u 	<- unique(wm2[,2])
	Lu1	<- length(u)
	temp	<- c()
	if(length(ux) > 1 && Lu1 > 0){
		for(i in 1:Lu1){
			aux.v <- which(wm2[,2] == u[i])
			if(length(aux.v) > 1){
				temp <- unique(c(temp, aux.v))
				var1 <- c(u[i], wm2[aux.v,1])
				n.nam <- dtbse[var1,1]
				for(j in 2:length(n.nam)){
					a3 <- setdiff(unlist(strsplit(n.nam[j],"_+_",fixed=TRUE)),unlist(strsplit(n.nam[1],"_+_",fixed=TRUE)))
					if(length(a3)==0) { newpath <- n.nam[1]
					} else {newpath <- paste(n.nam[1],paste(a3, collapse="_+_"),sep="_+_")}
					n.nam[1] <- newpath
				} 
			row2del <- unique(c(row2del,var1[-1]))
			dtbse[var1[1],1] <- n.nam[1]
			} 
		}
	}

	################################################################################
	##
	## 2) % = 1 AND PATHWAY 'i' IS REPEATED IN wm2
	##

	if(length(ux) > 1 && length(temp) > 0)
		wm2 <- wm2[-temp,] 
	wm2 <- check.null.dim(wm2)
	ux <- unique(wm2[,2])
	u <- unique(wm2[,1])
	Lu2 <- length(u)
	temp <- c()
	if(length(ux) > 1 && Lu2 > 0){
		for(i in 1:Lu2){
			aux.v <- which(wm2[,1] == u[i])
			if(length(aux.v) > 1){
				temp <- unique(c(temp,aux.v))
				var1 <- c(u[i], wm2[aux.v,2])
				var2 <- wm2[aux.v[which(wm2[aux.v,4] == max(wm2[aux.v,4]))],1:2]
				n.nam <- dtbse[var2,1]
				a3 <- setdiff(unlist(strsplit(n.nam[1],"_+_",fixed=TRUE)),unlist(strsplit(n.nam[2],"_+_",fixed=TRUE)))
				if(length(a3)==0){ newpath <- n.nam[2] 
				} else { newpath <- paste(n.nam[2],paste(a3, collapse="_+_"),sep="_+_")}
				n.nam[1] <- newpath	
				row2del <- unique(c(row2del,var1[var1%in%var2][1]))
				dtbse[var2[2],1] <- n.nam[1] 
			}
		}
	}

	################################################################################
	##
	## 3) % = 1 AND NO PATHWAY REPETITIONS IN wm2
	##

	if(length(ux) > 1 && length(temp) > 0)
		wm2 <- wm2[-temp,]
	wm2 <- check.null.dim(wm2)
	u <- unique(wm2[,1])
	Lu3 <- length(u)

	tt <- sum(wm2[,2] %in% wm2[,1])
	if(tt != 0){
		m2 <- wm2
		newwm2 <- matrix(rep(0,6*Lu3),Lu3,6)
		i <- 1
		newwm2[i,] <- m2[i,]
		j <- which(m2[,1] == newwm2[i,2])
		while(i < Lu3){
			i <- i + 1
			if(length(j) > 0){
				newwm2[i,] <- m2[j,]
				m2[j,] <- m2[i,]
				m2[i,] <- rep(NA,6)
			} else { newwm2[i,] <- m2[i,] }
			j <- which(m2[,1] == newwm2[i,2])
		}
		newwm2[which(!is.na(m2[,1])),] <- m2[which(!is.na(m2[,1])),]
		wm2 <- newwm2 
	}
	
	if(Lu3 > 0){
		for(i in 1:Lu3){
			n.nam <- dtbse[wm2[i,1:2],1]
			a3 <- setdiff(unlist(strsplit(n.nam[1],"_+_",fixed=TRUE)),unlist(strsplit(n.nam[2],"_+_",fixed=TRUE)))
			if(length(a3)==0){ newpath <- n.nam[2]
			} else { newpath <- paste(n.nam[2],paste(a3, collapse="_+_"),sep="_+_") }
			n.nam[1] <- newpath	
			row2del <- unique(c(row2del, wm2[i,1]))
			dtbse[wm2[i,1:2],1] <- n.nam[1]
		}
	}

	################################################################################
	##
	## 4) % < 1 AND PATHWAY 'i' IS REPEATED in wm1
	##

	if(length(which(wm1[,5] == 1)) > 0){
		wm3 <- wm1[-which(wm1[,5] == 1),]
	} else wm3 <- wm1
	wm3 <- check.null.dim(wm3)
	newmat <- c()
	u2 <- unique(wm3[,1])
	Lu4 <- length(u2)
	maxNg <- max(Ng)
	if(Lu4 > 0){
		for(i in 1:Lu4){
			aux.v <- which(wm3[,1] == u2[i])
			if(length(aux.v) > 1){
				var1 <- c(u2[i], wm3[aux.v,2])
				var2 <- wm3[aux.v[which(wm3[aux.v,5] == max(wm3[aux.v,5]))],1:2]
				if(!is.null(dim(var2))) var2 <- var2[1,]
				nm4 <- dtbse[var2,1]
				a3 <- setdiff(unlist(strsplit(nm4[1],"_+_",fixed=TRUE)),unlist(strsplit(nm4[2],"_+_",fixed=TRUE)))
				if(length(a3)==0) { newpath <- nm4[2]
				} else { newpath <- paste(nm4[2],paste(a3, collapse="_+_"),sep="_+_")}
				newurl <- dtbse[var2[2],2]
				l1 <- dtbse[var2[1],sgc:(Ng[var2[1]]+2)]
				l2 <- dtbse[var2[2],sgc:(Ng[var2[2]]+2)]
				newgenes <- unique(c(l1,l2))
				if(length(newgenes) < maxNg){
					aux9 <- rep("" ,maxNg - length(newgenes))
					newmat <- rbind(newmat, c(newpath,newurl,newgenes,aux9))
				} else {	lnm <- length(newmat[,1])
						xs <- length(newgenes) - maxNg
						if(xs == 0 || lnm == 0){ 
							newmat <- rbind(newmat, c(newpath,newurl,newgenes))
						} else {	auxmat <- matrix(rep("",lnm*xs),lnm,xs)
								newmat <- cbind(newmat,auxmat)
								newmat <- rbind(newmat, c(newpath,newurl,newgenes))
							}
						maxNg <- length(newgenes)
					}
				row2del <- unique(c(row2del,var2))
			}
		}
		if(length(newmat[1,]) > length(dtbse[1,])){
			xs <- length(newmat[1,]) - length(dtbse[1,])
			ldb <- length(dtbse[,1])
			auxmat <- matrix(rep("",ldb*xs),ldb,xs)
			dtbse <- cbind(dtbse,auxmat)
		}
		dtbse <- rbind(dtbse, newmat)
	}
	
	################################################################################
	##
	## 5) % < 1 AND (PATHWAY 'j' IS REPEATED in wm3 OR NO REPETITIONS)
	##

	if(length(u2[u2 %in% row2del]) > 0)
		wm3 <- wm3[-which(wm3[,1] %in% u2[u2 %in% row2del]),]
	if(is.null(dim(wm3))){ 
		u2 <- wm3[1] 
	} else u2 <- unique(wm3[,1])
	wm3 <- check.null.dim(wm3)
	Lu4 <- length(u2)
	newmat <- c()
	if(Lu4 > 0){
		for(i in 1:length(wm3[,1])){
			nm4 <- dtbse[wm3[i,1:2],1]
			a3 <- setdiff(unlist(strsplit(nm4[1],"_+_",fixed=TRUE)),unlist(strsplit(nm4[2],"_+_",fixed=TRUE)))
			if(length(a3)==0) { newpath <- nm4[2]
			} else {newpath <- paste(nm4[2],paste(a3, collapse="_+_"),sep="_+_")}
			newurl <- dtbse[wm3[i,2],2]
			l1 <- dtbse[wm3[i,1],sgc:(wm3[i,3]+2)]
			l2 <- dtbse[wm3[i,2],sgc:(wm3[i,4]+2)]
			newgenes <- unique(c(l1,l2))
			if(length(newgenes) < maxNg){
				aux9 <- rep("" ,maxNg - length(newgenes))
				newmat <- rbind(newmat, c(newpath,newurl,newgenes,aux9))
			} else {	lnm <- length(newmat[,1])
					xs <- length(newgenes) - maxNg
					if(xs == 0 || lnm == 0){ 
						newmat <- rbind(newmat, c(newpath,newurl,newgenes))
					} else {	auxmat <- matrix(rep("",lnm*xs),lnm,xs)
							newmat <- cbind(newmat,auxmat)
							newmat <- rbind(newmat, c(newpath,newurl,newgenes))
						}
					maxNg <- length(newgenes)
				}
			row2del <- unique(c(row2del,wm3[i,1:2]))
		}
		if(length(newmat[1,]) > length(dtbse[1,])){
			xs <- length(newmat[1,]) - length(dtbse[1,])
			ldb <- length(dtbse[,1])
			auxmat <- matrix(rep("",ldb*xs),ldb,xs)
			dtbse <- cbind(dtbse,auxmat)
		}
		dtbse <- rbind(dtbse, newmat)
	}

	dtbse <- dtbse[-row2del,]

	if(is.null(dim(dtbse))){ cat("You got a single megapathway!. Try with a greater overlap\n") 
				 arr <- NULL
	} else{	NoP <- length(dtbse[,1])
			nocg_mat <- matrix(NA, NoP, NoP)
			for(i in 1:NoP){
				for(j in 1:NoP){
					if(length(which(dtbse[i,]=="")) > 0)
						a <- dtbse[i,sgc:length(dtbse[i,-c(which(dtbse[i,]==""))])]
					else a <- dtbse[i,sgc:length(dtbse[i,])]
					if(length(which(dtbse[j,]=="")) > 0)
						b <- dtbse[j,sgc:length(dtbse[j,-c(which(dtbse[j,]==""))])]
					else b <- dtbse[j,sgc:length(dtbse[j,])]
					nocg <- length(a[a %in% b])
					nocg_mat[i,j] <- nocg
				}
			}
			Ng <- diag(nocg_mat)
			diag(nocg_mat) <- 0
			perc_mat <- nocg_mat / Ng
			arr <- which(perc_mat > max_overlap, arr.ind=TRUE)
			rev.arr <- cbind(arr[,2],arr[,1])
			arr <- arr[which(perc_mat[rev.arr] > min_overlap),]
			maxNg <- max(Ng)
			k <- k + 1
			cat(k, "\n")
		}
}

########################################################################################
##
## END OF BIG LOOP
## GENERATING THE FILES:
##	1) NEW DATABASE
## 	new_name + max_overlap + min_overlap + loops + # of pathways + max # of genes in a pathway
##
##	2) LIST OF SUPERPATHWAYS
##	
##

if(is.null(dim(dtbse))){ cat("No output files are created\n")
} else{ ind <- 1
	plus <- c()
	newname <- c()
	for(i in 1:length(dtbse[,1])){
		lp <- unlist(gregexpr("_+_",dtbse[,1][i],fixed=TRUE))
		if(lp[1] != -1){ 
			plus[i] <- length(lp) + 1 
			nam2 <- paste("Superpathway",ind,max_overlap,min_overlap,plus[i],sep="_")
			newname[i] <- nam2
			ind <- ind + 1
		} else{ 
				plus[i] <- 1 
				newname[i] <- dtbse[,1][i]
			}
	}
	
	super <- cbind(newname,dtbse[,1])
	dtbse[,1] <- newname
	nam0 <- paste(new_name,max_overlap,min_overlap,k,NoP,maxNg,sep="_")
	nam1 <- paste(nam0,"txt",sep=".")
	for(i in 1:length(dtbse[,1])){
		write(dtbse[i,1:(Ng[i]+2)], nam1, sep="\t", ncolumns=(Ng[i]+2),append=TRUE)
	}
	write(t(subset(super, grepl("^Su", newname))), paste("Superpathways",nam1,sep="_"), sep="\t", ncolumns=2)
	cat("DONE!","\n")
    }
}
