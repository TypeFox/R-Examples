#########################
#Load files from GenClone
#########################

convert_GC <- function(data1, num, ele){
	
	if (missing(ele)){
		res <- apply(data1, 1:2, function(x) sort(c(substr(x, 1, num), substr(x, num+1, num*2))))
	} else {
		res <- apply(data1, 1:2, function(x) sort(unlist(strsplit(x, split = ele))))
	}
	mat_all <- as.data.frame(t(apply(res, 2, function(x) rbind(x))), stringsAsFactors = FALSE)

	if (ncol(mat_all) != ncol(data1)*2) {stop("Warning 1")}
	if (nrow(mat_all) != nrow(data1)) {stop("Warning 2")}

	names(mat_all) <- unlist(lapply(names(data1), function(x) paste(x, 1:2, sep = "_")))
mat_all
}


transcript_GC <- function(obj, ele, num1, num2, num3){

	dataGC <- read.csv(obj, sep = ele, header = FALSE)
	N <- dataGC[1,1]
	coord_center <- c(mean(c(0,dataGC[1,2])), mean(c(0,dataGC[1,3])))
	nb_loci <- dataGC[1,4]
	ploid <- dataGC[1,5]
		if (ploid != num1){stop("Ploidy different from indicated")}
	names_loci <- dataGC[1,6:ncol(dataGC)]
	names_loci <- as.vector(apply(names_loci, 1:2, function(x) x <- c(as.vector(x))))
		if (length(names_loci) != nb_loci){stop("Number of loci names different")}
		if (nb_loci != num2){stop("Number of loci different from indicated")}

	dataGC <- dataGC[-1,]
	coord <- dataGC[,2:3]
	rownames(coord) <- dataGC[,1]
	colnames(coord) <- c("x", "y")

	data1b <- dataGC[,4:(4+nb_loci-1)]
	rownames(data1b) <- dataGC[,1]
	colnames(data1b) <- names_loci
		if (num1 == 2){
	data1b <- convert_GC(data1b, num3)
		} else {
	data1b <- data1b
		}	
	list(data_genet = data1b, data_coord = coord, names_loci = names_loci, names_units = rownames(coord))
}

sort_all <- function(data1){

	index_l <- 1:c(ncol(data1)/2)*2-1
	nb_loci <- ncol(data1)/2
	mat <- as.data.frame(matrix(NA, ncol = ncol(data1), nrow = nrow(data1)))
		for (j in index_l){ 
			for (i in 1:nrow(data1)){
				mat[i,c(j, j+1)] <- sort(c(data1[i,j], data1[i,j+1]))
			}
		}
	names(mat) <- names(data1)
	mat_all <- mat
	mat_all
}


###################
#Export to Adegenet
###################

export_genclone_genind <- function(data1, ele){

	index_l <- 1:c(ncol(data1)/2)*2-1 
	nb_loci <- ncol(data1)/2
	N <- nrow(data1)
	tab_genind <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nrow(data1)))

	for (j in index_l){
		if (missing(ele)){
			tab_genind[, c((j+1)/2)] <- paste(data1[,j], data1[,j+1], sep = "")
		} else {
			tab_genind[, c((j+1)/2)] <- paste(data1[,j], data1[,j+1], sep = ele)
		}
	}
	names(tab_genind) <- paste("locus", 1:nb_loci, sep = "_")
tab_genind
}

##################
#Export to Genetix
##################

export_genclone_genetix <- function(data1, haploid = FALSE, ele, name){

	index_l <- 1:c(ncol(data1)/2)*2-1
	nb_loci <- ncol(data1)/2
	N <- nrow(data1)
	tab_genix <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nrow(data1)))

	for (j in index_l){
		if (missing(ele)){
			tab_genix[, c((j+1)/2)] <- paste(data1[,j], data1[,j+1], sep = "")
		} else {
			tab_genix[, c((j+1)/2)] <- paste(data1[,j], data1[,j+1], sep = ele)
		}
	}

	if(haploid){
		tab_genix <- cbind(rep("infile", times = nrow(tab_genix)), rownames(tab_genix), data1)
	} else {
		tab_genix <- cbind(rep("infile", times = nrow(tab_genix)), rownames(tab_genix), tab_genix)
	}
	names(tab_genix) <- c("", "", paste("locus", 1:nb_loci, sep = "_"))
		
	if(missing(name)){
		tab_genix
	} else {
		write.table(tab_genix, name, row.names = FALSE, quote = FALSE, sep = "\t")
	}
}

###################
#Export to Arlequin
###################

export_genclone_arlequin <- function(data1, haploid = FALSE, name){

	mat <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
	mat[1,] <- "[Profile]"
	mat[2,] <- ' Title=" "'
	mat[3,] <- " NbSamples=1"
	mat[4,] <- " DataType=STANDARD"
	mat[5,] <- " GenotypicData=1"
	mat[6,] <- " GameticPhase=0"
	mat[7,] <- " LocusSeparator=WHITESPACE"
	mat[8,] <- " MissingData='0'"
	mat[9,] <- ""
	mat[10,] <- "[Data]"
	mat[11,] <- ""
	mat[12,] <- "[[Samples]]"
	mat[13,] <- '   SampleName="infile"'
	mat[14,] <- paste("   SampleSize=", nrow(data1), sep = "")
	mat[15,] <-"   SampleData={"

	mat_head <- mat

	index_l <- 1:c(ncol(data1)/2)*2-1
	nb_loci <- ncol(data1)/2
	N <- nrow(data1)
	mat <- as.data.frame(matrix(NA, nrow = N*2, ncol = 1))
	index_unit <- 1:c(nrow(mat)/2)*2-1

if (haploid){

	index_l <- ncol(data1)

		for (i in 1:N){
			if (i < 10){
				mat[i,] <- paste(paste("         ", i, " ", 1, " ", sep = ""), data1[i,index_l], sep = "")
			} else {
				if (i < 100){
					mat[i,] <- paste(paste("        ", i, " ", 1, " ", sep = ""), data1[i,index_l], sep = "")
				} else {
					if (i < 1000){
						mat[i,] <- paste(paste("       ", i, " ", 1, " ", sep = ""), data1[i,index_l], sep = "")
					} else {
						if (i < 10000){
							mat[i,] <- paste(paste("      ", i, " ", 1, " ", sep = ""), data1[i,index_l], sep = "")
						} else {
							if (i > 10000) stop("Export to Arlequin is not available for N > 10000")
						}
					}
				}
			}
		}
} else {
		for (i in 1:N){
			if (i < 10){
				mat[(i*2-1),] <- paste(paste("         ", i, " ", 1, " ", sep = ""), paste(data1[i,index_l], collapse = " "), sep = "")
				mat[(i*2),] <- paste("             ", paste(data1[i,index_l+1], collapse = " "), sep = "")
			} else {
				if (i < 100){
					mat[(i*2-1),] <- paste(paste("        ", i, " ", 1, " ", sep = ""), paste(data1[i,index_l], collapse = " "), sep = "")
					mat[(i*2),] <- paste("             ", paste(data1[i,index_l+1], collapse = " "), sep = "")
				} else {
					if (i < 1000){
						mat[(i*2-1),] <- paste(paste("       ", i, " ", 1, " ", sep = ""), paste(data1[i,index_l], collapse = " "), sep = "")
						mat[(i*2),] <- paste("             ", paste(data1[i,index_l+1], collapse = " "), sep = "")
					} else {
						if (i < 10000){
							mat[(i*2-1),] <- paste(paste("      ", i, " ", 1, " ", sep = ""), paste(data1[i,index_l], collapse = " "), sep = "")
							mat[(i*2),] <- paste("             ", paste(data1[i,index_l+1], collapse = " "), sep = "")
						} else {
							if (i > 10000) stop("Export to Arlequin is not available for N > 10000")
						}
					}
				}
			}
		}
}
	mat_genet <- mat

		if (missing(name)){
			tab_arl <- rbind(mat_head, mat_genet, paste("  }"))
			tab_arl
		} else {
			write.table(rbind(mat_head, mat_genet, paste("  }")), name, row.names = FALSE, col.names = FALSE, quote = FALSE)
		}
}


################
#Listing alleles
################


corresp_loci <- function(data1, haploid = FALSE){

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	nb_loci <- ncol(data1)/2
	index_l <- 1:c(ncol(data1)/2)*2-1
		for (i in index_l){
			print(paste(paste("locus",c((i+1)/2), sep = "_"), names(data1)[i], sep = "/"))
		}
}


list_all_obj_core <- function(data1, haploid = FALSE){

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	nb_loci <- ncol(data1)/2
	index_l <- 1:c(ncol(data1)/2)*2-1
		for (i in index_l){
			print(paste(paste("locus",c((i+1)/2), sep = "_"), names(data1)[i], sep = "/"))
			print(c(length(unique(c(data1[,i], data1[,i+1]))), unique(c(data1[,i], data1[,i+1]))))
		}
}


list_all_obj <- function(data1, haploid = FALSE, vecpop = NULL){

		if (length(vecpop) != 0){
			if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
			
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

			datatot <- split(data1, vecpop)
				for (p in 1:length(unique(vecpop))){
					print(unique(vecpop_o)[p])
					list_all_obj_core(datatot[[p]], haploid)
				}
		} else {
			list_all_obj_core(data1, haploid)
		}
}


list_all_tab_core <- function(data1, haploid = FALSE){

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	nb_loci <- ncol(data1)/2
	nb_unit <- nrow(data1)/2
	index_l <- 1:c(ncol(data1)/2)*2-1

	tab <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nb_unit*2))
	colnames(tab) <- paste("locus", 1:nb_loci, sep = "_")
		for (i in index_l){  
			tab[1:length(unique(c(data1[,i], data1[,i+1]))),c((i+1)/2)] <- unique(c(data1[,i], data1[,i+1]))
		}
	tab_list_all <- tab[1:max(apply(tab, 2, function(x) length(which(is.na(x) == "FALSE")))),]
	tab_list_all[is.na(tab_list_all)] <- ""
	tab_list_all
}


list_all_tab <- function(data1, haploid = FALSE, vecpop = NULL){

		if (length(vecpop) != 0){
			if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
			
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

			datatot <- split(data1, vecpop)
			res <- lapply(datatot, function(x) list_all_tab_core(x, haploid))
			names(res) <- unique(vecpop_o)
		} else {
			 res <- list_all_tab_core(data1, haploid)
		}
	res
}


list_all_tab2_core <- function(data1, haploid = FALSE){

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	nb_loci <- ncol(data1)/2
	index_l <- 1:c(ncol(data1)/2)*2-1

	tab <- as.data.frame(matrix(NA, ncol = 2, nrow = 1))
		for (i in index_l){
			tab <- rbind(tab, cbind(paste("locus", c((i+1)/2), sep = "_"), sort(unique(c(data1[,i], data1[,i+1])))))
		}
	tab <- tab[-1,]
	row.names(tab) <- 1:nrow(tab)
	colnames(tab) <- c("locus", "allele")
	mat_list_all2 <- tab
	mat_list_all2
}


list_all_tab2 <- function(data1, haploid = FALSE, vecpop = NULL){

		if (length(vecpop) != 0){			
			if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
			
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

			datatot <- split(data1, vecpop)
			res <- lapply(datatot, function(x) list_all_tab2_core(x, haploid))
			names(res) <- unique(vecpop_o)
		} else {
			 res <- list_all_tab2_core(data1, haploid)
		}
	res
}


############
#Listing MLG
############

MLG_list_core <- function(data1){

	nb_loci <- ncol(data1)/2
	mat_genet2 <- unique(data1)
	list_genet <- list(NULL)

	for (i in 1:nrow(mat_genet2)){
		list_genet[[i]] <- c(which(apply(data1, 1, function(x) which(sum(x == mat_genet2[i,]) == nb_loci*2)) != 0))
			if (length(names(list_genet[[i]])) != 0) {
				names(list_genet[[i]]) <- NULL
			}
	}
	list_genet <- list_genet
}


MLG_list <- function(data1, vecpop = NULL){

		if (length(vecpop) != 0){
			if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
			
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

			datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					res[[p]] <- MLG_list_core(as.data.frame(datatot[[p]]))
				}
			names(res) <- unique(vecpop_o)
		} else {
			 res <- MLG_list_core(data1)
		}
res
}


MLG_tab_core <- function(data1){

	list_genet <- MLG_list(data1)
	data_MLG <- unique(data1)
	nb_loci <- ncol(data1)/2
	G <- nrow(data_MLG)
	list_genet2 <- as.data.frame(matrix(NA, ncol = max(sapply(list_genet, length)), nrow = G))
		for (i in 1:nrow(data_MLG)){
			list_genet2[i, 1:sapply(list_genet, length)[i]] <- c(which(apply(data1, 1, function(x) which(sum(x == data_MLG[i,]) == nb_loci*2)) != 0))
		}
	list_genet2[is.na(list_genet2)] <- ""
	colnames(list_genet2) <- paste("unit", 1:ncol(list_genet2), sep = "_")
	list_genet2
}


MLG_tab <- function(data1, vecpop = NULL){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		res <- lapply(datatot, function(x) MLG_tab_core(x))
		names(res) <- unique(vecpop_o)
	} else {
		 res <- MLG_tab_core(data1)
	}
	res
}



########################################################
#allelic frequencies: with or without Round-Robin method
########################################################


freq_RR_core <- function(data1, haploid = FALSE, genet = FALSE, RR = FALSE){

	if (genet){
		data1 <- unique(data1)
	}

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	index_l <- 1:c(ncol(data1)/2)*2-1
	tab <- as.data.frame(matrix(NA, ncol = 3, nrow = 1)) 
		if (RR){
			for (j in index_l){
				obj <- cbind(paste("locus", c((j+1)/2), sep = "_"), sort(unique(c(data1[,j], data1[,j+1]))), NA)
					for (i in 1:nrow(obj)){
						obj[i,3] <- length(which(obj[i,2] == c(data1[row.names(unique(data1[,-c(j,j+1)])),c(j,j+1)][,1], data1[row.names(unique(data1[,-c(j,j+1)])),c(j,j+1)][,2])))/(2*nrow(data1[row.names(unique(data1[,-c(j,j+1)])),]))
					}

	if (haploid){
		obj[obj==0] <- 1/(nrow(unique(data1)))
	} else {
		obj[obj==0] <- 1/(2*nrow(unique(data1)))
	}
				tab <- rbind(tab, obj)
			}
		} else {
			for (j in index_l){
				obj <- cbind(paste("locus",c((j+1)/2), sep="_"), sort(unique(c(data1[,j], data1[,j+1]))), NA)
					for (i in 1:nrow(obj)){
						obj[i,3] <- length(which(obj[i,2] == data1[,c(j,j+1)]))/(2*nrow(data1))
					}

	if (haploid){
		obj[obj==0] <- 1/(nrow(data1))
	} else {
		obj[obj==0] <- 1/(2*nrow(data1))
	}
				tab <- rbind(tab, obj)
			}
		}
	tab <- tab[-1,]
	row.names(tab) <- 1:nrow(tab)
	colnames(tab) <- c("locus", "allele", "freq")
	mat_freq <- tab
	mat_freq[,3] <- as.numeric(mat_freq[,3])
	mat_freq
}


freq_RR <- function(data1, haploid = FALSE, vecpop = NULL, genet = FALSE, RR = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		res <- lapply(datatot, function(x) freq_RR_core(x, haploid, genet, RR))
		names(res) <- unique(vecpop_o)
	} else {
		res <- freq_RR_core(data1, haploid, genet, RR)
	}
	res
}


###################################
#Probability of a genotype i (Pgen)
###################################

freq_finder <- function(data1, i , j, haploid = FALSE, genet = FALSE, RR = FALSE){

	if (genet){
		data1 <- unique(data1)
	}

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	index_l <- 1:c(ncol(data1)/2)*2-1
	ncol1 <- 2
	ncol2 <- 3
	data2 <- freq_RR(data1, haploid, genet, RR)
		if (sum(j == index_l) != 0){
			freq <- data2[as.numeric(rownames(data2[data2==paste("locus", c((j+1)/2), sep="_"),])[which(data2[data2==paste("locus", c((j+1)/2), sep="_"), ncol1]==data1[i,j])]), ncol2]
		} else {
			j <- j-1
			freq <- data2[as.numeric(rownames(data2[data2==paste("locus", c((j+1)/2), sep="_"),])[which(data2[data2==paste("locus", c((j+1)/2), sep="_"), ncol1]==data1[i,j+1])]), ncol2]
		}
	as.numeric(freq)
}


pgen_core <- function(data1, data2, haploid = FALSE){

	if (haploid){
		index_l <- 1:ncol(data1)
		nb_loci <- ncol(data1)
	} else {
		index_l <- 1:c(ncol(data1)/2)*2-1
		nb_loci <- ncol(data1)/2
	}

	ncol_all <- 2
	ncol_freq <- 3

	recup <- NULL
	for (i in 1:nrow(data1)){
		pgeni <- 1
		h <- nb_loci

	if (haploid){
			for (j in index_l){
				pgeni <- pgeni*
					as.numeric(data2[as.numeric(rownames(data2[data2 == paste("locus", j, sep="_"),])[which(data2[data2 == paste("locus", j, sep="_"), ncol_all] == data1[i,j])]), ncol_freq])
			}
		recup <- c(recup, pgeni)
	} else {
			for (j in index_l){
				if(sum(data1[i,j] == data1[i,j+1]) == length(data1[i,j])) h <- h-1
				pgeni <- pgeni*
					as.numeric(data2[as.numeric(rownames(data2[data2==paste("locus", c((j+1)/2), sep="_"),])[which(data2[data2==paste("locus", c((j+1)/2), sep="_"), ncol_all]==data1[i,j])]), ncol_freq])*
					as.numeric(data2[as.numeric(rownames(data2[data2==paste("locus", c((j+1)/2), sep="_"),])[which(data2[data2==paste("locus", c((j+1)/2), sep="_"), ncol_all]==data1[i,j+1])]), ncol_freq])
			}
		recup <- c(recup, pgeni*2^h)
	}
	}
	res_pgen <- as.data.frame(recup)
	colnames(res_pgen) <- "pgen"
	res_pgen
}


pgen <- function(data1, haploid = FALSE, vecpop = NULL, genet = FALSE, RR = FALSE){

	if (genet & RR){stop("Round-Robin method is not compatible with genet option.")}

	if (genet){
		datafreq <- freq_RR(data1, haploid, vecpop, RR = FALSE, genet = TRUE)
	} else 
	if (RR){
		datafreq <- freq_RR(data1, haploid, vecpop, RR = TRUE, genet = FALSE)
	} else {
		datafreq <- freq_RR(data1, haploid, vecpop, RR = FALSE, genet = FALSE)	
	}

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		res <- list(NULL)
			for (p in 1:length(unique(vecpop))){
				rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
				res[[p]] <- pgen_core(datatot[[p]], datafreq[[p]], haploid)
			}
				names(res) <- unique(vecpop_o)
	} else {
		data2 <- datafreq
		res <- pgen_core(data1, data2, haploid)
	}
	res
}


###################
#Fis, Hatt et Hobs
###################

Fis_core <- function(data1, data2, genet = FALSE, RR = FALSE){

	ncol_freq <- 3
	ncol_all <- 2
	nb_loci <- ncol(data1)/2
	index_l <- 1:c(ncol(data1)/2)*2-1
	tab <- as.data.frame(matrix(0, ncol = 4, nrow = nb_loci))
		if (RR){
		tab2 <- as.data.frame(matrix(NA, ncol = 3, nrow = 1))
			for (i in index_l){
				obj <- cbind(paste("locus", c((i+1)/2), sep="_"), sort(unique(c(data1[,i], data1[,i+1]))), NA)
				sub <- data1[row.names(unique(data1[,-c(i,i+1)])),]
					for (j in 1:nrow(obj)){
						obj[j,3] <- length(which(obj[j,2] == c(sub[,i], sub[,(i+1)])))/(2*nrow(sub))
					} 				
				obj[obj==0] <- 1/(2*nrow(unique(data1)))
				tab2 <- rbind(tab2, obj)
				tab[c((i+1)/2),1] <- paste("locus",c((i+1)/2), sep = "_")
				tab[c((i+1)/2),2] <- sum(sub[,i] != sub[,(i+1)])/nrow(sub)
				tab[c((i+1)/2),3] <- (2*nrow(sub)/(2*nrow(sub)-1))*(1 - sum(as.numeric(obj[,3])^2))
				tab[c((i+1)/2),4] <- (tab[c((i+1)/2),3]-as.numeric(tab[c((i+1)/2),2]))/tab[c((i+1)/2),3]
			}
		tab2 <- tab2[-1,]
		row.names(tab2) <- 1:nrow(tab2)
		colnames(tab2) <- c("locus", "allele", "freq")
		mat_freq <- tab2
		colnames(tab) <- c("locus", "Hobs", "Hatt", "Fis")
		mat_Fis_RR <- tab
		res <- list(mat_freq, mat_Fis_RR)
		} else 
		if (genet){
			for (i in index_l){
				tab[c((i+1)/2),1] <- paste("locus", c((i+1)/2), sep = "_")
				tab[c((i+1)/2),2] <- sum(unique(data1)[,i] != unique(data1)[,i+1])/nrow(unique(data1))
				tab[c((i+1)/2),3] <- (2*nrow(unique(data1))/(2*nrow(unique(data1))-1))*(1 - sum(data2[rownames(data2[data2==paste("locus",c((i+1)/2), sep="_"),]), ncol_freq]^2))
				tab[c((i+1)/2),4] <- (tab[c((i+1)/2),3]-as.numeric(tab[c((i+1)/2),2]))/tab[c((i+1)/2),3]
			}
		colnames(tab) <- c("locus", "Hobs", "Hexp", "Fis")
		tab_Fis_genet <- tab
		res <- tab_Fis_genet
		} else {
			for (i in index_l){
				tab[c((i+1)/2),1] <- paste("locus", c((i+1)/2), sep = "_")
				tab[c((i+1)/2),2] <- sum(data1[,i] != data1[,i+1])/nrow(data1)
				tab[c((i+1)/2),3] <- (2*nrow(data1)/(2*nrow(data1)-1))*(1 - sum(data2[rownames(data2[data2==paste("locus",c((i+1)/2), sep="_"),]), ncol_freq]^2))
				tab[c((i+1)/2),4] <- (tab[c((i+1)/2),3]-as.numeric(tab[c((i+1)/2),2]))/tab[c((i+1)/2),3]
			}
		colnames(tab) <- c("locus", "Hobs", "Hexp", "Fis")
		tab_Fis_ramet <- tab
		res <- tab_Fis_ramet
		}
res
}


Fis <- function(data1, vecpop = NULL, genet = FALSE, RR = FALSE){

	if (genet & RR){stop("Round-Robin method is not compatible with genet option.")}

	if (RR){
		 datafreq <- NULL 
	} else
	if (genet){
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = TRUE, RR = FALSE)
	} else {
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = FALSE, RR = FALSE)
	}

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					res[[p]] <- Fis_core(datatot[[p]], datafreq[[p]], genet, RR)
				}
			names(res) <- unique(vecpop_o)
	} else {
		data2 <- datafreq
		res <- Fis_core(data1, data2, genet, RR)
	}
	res
}


###########################################################
#Probability of a genotype i with H-W deviation : (PgenFis)
###########################################################

pgen_Fis_core <- function(data1, data2, genet = FALSE, RR = FALSE){

	ncol_all <- 2
	ncol_freq <- 3
	recup <- NULL
	nb_loci <- ncol(data1)/2
	index_l <- 1:c(ncol(data1)/2)*2-1
	if (RR){
		res_Fis <- Fis_core(data1, data2, genet, RR)[[2]][,4]
		res_Fis_NA <- which(is.na(res_Fis))
	} else {
		res_Fis <- Fis_core(data1, data2, genet, RR)[,4]
		res_Fis_NA <- which(is.na(res_Fis))
	}

	if (length(res_Fis_NA) != 0){
		nb_loci2 <- c(1:nb_loci)[-res_Fis_NA]
		index_l <- nb_loci2*2-1
		nb_loci <- length(nb_loci2)
	}
		for (i in 1:nrow(data1)){ 
			pgenfisi <- 1
			h <- nb_loci
				for (j in index_l){
					fi <- as.numeric(data2[rownames(data2[data2==paste("locus",c((j+1)/2), sep="_"),])[which(data2[data2==paste("locus",c((j+1)/2), sep="_"), ncol_all]==data1[i,j])], ncol_freq])
					gi <- as.numeric(data2[rownames(data2[data2==paste("locus",c((j+1)/2), sep="_"),])[which(data2[data2==paste("locus",c((j+1)/2), sep="_"), ncol_all]==data1[i,j+1])], ncol_freq])
						if (sum(data1[i,j] == data1[i,j+1]) == length(data1[i,j])) {
							h <- h-1
							pgenfisi <- pgenfisi*(fi*(fi+(1-fi)*res_Fis[c((j+1)/2)]))
						} else {
							pgenfisi <- pgenfisi*(fi*gi*(1+(-1*res_Fis[c((j+1)/2)])))
						}
				}
			recup <- c(recup, pgenfisi*2^h)
		}
	pgenFis <- as.data.frame(recup)
	colnames(pgenFis) <- "pgenFis"
	pgenFis
}


pgen_Fis <- function(data1, vecpop = NULL, genet = FALSE, RR = FALSE){

	if (genet & RR){stop("Round-Robin method is not compatible with genet option.")}

	if (genet){
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = TRUE, RR = FALSE)
	} else 
	if (RR){
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = FALSE, RR = TRUE)
	} else {
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = FALSE, RR = FALSE)
	}

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		res <- list(NULL)
			for (p in 1:length(unique(vecpop))){
				rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
				res[[p]] <- pgen_Fis_core(datatot[[p]], datafreq[[p]], genet, RR)
			}
			names(res) <- unique(vecpop_o)
		} else {
			data2 <- datafreq
			res <- pgen_Fis_core(data1, data2, genet, RR)
		}
	res
}


####################################################
#probability of one genotype with repro events: psex
####################################################

#load("C:/Users/Diane/RClone/data/factoR.rda")

psex_core <- function(data1, data2, haploid = FALSE, MLGsim = FALSE, nbrepeat = NULL, bar = FALSE){

	list_genet <- MLG_list(data1)
	ncol_all <- 2
	ncol_freq <- 3
	factoR <- factoR
	res_pgen <- pgen_core(data1, data2, haploid)
	tab <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(data1)))

	if (length(list_genet[which(sapply(list_genet, length) > 1)]) > 1){

	if (MLGsim){
		for (m in 1:length(list_genet[which(sapply(list_genet, length) > 1)])){
			recup <- NULL
			sub_list <- list_genet[which(sapply(list_genet, length) > 1)][[m]]
			l <- sub_list[1]
			recup <- c(recup, factoR[nrow(data1),2]/(factoR[length(sub_list), 2]*factoR[c(nrow(data1)-length(sub_list)),2])*
				(res_pgen[l,])^(length(sub_list))*
				(1-res_pgen[l,])^(nrow(data1)-length(sub_list)))
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]],2] <- recup
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]],1] <- list_genet[which(sapply(list_genet, length) > 1)][[m]][1]
		}
	} else {
		for (m in 1:length(list_genet[which(sapply(list_genet, length) > 1)])){
			recup <- NULL
				for (l in list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]){
					recup <- c(recup, factoR[nrow(data1),2]/(factoR[which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l),2]*factoR[c(nrow(data1)-which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l)),2])*
					(res_pgen[l,])^(which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l))*
					(1-res_pgen[l,])^(nrow(data1)-which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l)))
				}
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]][-1],2] <- recup
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]][-1],1] <- list_genet[which(sapply(list_genet, length) > 1)][[m]][1]
		}
	}
	names(tab) <- c("genet", "psex")
	psexFR <- tab

	if (length(nbrepeat) != 0){
		psex_recup <- NULL
		nb_loci <- ncol(data1)/2
		index_l <- 1:c(ncol(data1)/2)*2-1
		N <- nrow(data1)

	if (bar){
		total <- nbrepeat
		pb <- txtProgressBar(min = 0, max = total, style = 3)
	}

		for (s in 1:nbrepeat){
			tab2 <- as.data.frame(matrix(NA, nrow = N, ncol = nb_loci*2))
				for (j in 1:N){
					for (i in index_l){
						tab2[j, c(i, i+1)] <- sort(sample(split(data2, data2[,1])[[c((i+1)/2)]][,ncol_all], 2, prob = split(data2, data2[,1])[[c((i+1)/2)]][,ncol_freq], replace = TRUE))
					}
				}
			tab_sim <- tab2
			colnames(tab_sim) <- colnames(data1)
			MLG_sim <- unique(tab_sim)
				if (nrow(MLG_sim) != nrow(tab_sim)){
					pgen_sim <- pgen_core(tab_sim, data2, haploid)
					psex_sim <- NULL
					list_genet_sim <- MLG_list(tab_sim)
	if (MLGsim){
						for (m in 1:length(list_genet_sim[which(sapply(list_genet_sim, length) > 1)])){
							recup <- NULL
							sub_list <- list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]]
							l <- sub_list[1]
							recup <- c(recup, factoR[nrow(tab_sim),2]/(factoR[length(sub_list), 2]*factoR[c(nrow(tab_sim)-length(sub_list)),2])*
								(pgen_sim[l,])^(length(sub_list))*
								(1-pgen_sim[l,])^(nrow(tab_sim)-length(sub_list)))
						}
					psex_sim <- c(psex_sim, recup)
	} else {
						for (m in 1:length(list_genet_sim[which(sapply(list_genet_sim, length) > 1)])){
							recup <- NULL
								for(l in list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]){
									recup <- c(recup, factoR[nrow(tab_sim), 2]/(factoR[which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l), 2]*factoR[c(nrow(tab_sim)-which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l)),2])*
										(pgen_sim[l,])^(which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l))*
										(1-pgen_sim[l,])^(nrow(tab_sim)-which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l)))
								}
							psex_sim <- c(psex_sim, recup)
						}
	}
					psex_recup <- c(psex_recup, psex_sim)
				}
	if (bar){
		setTxtProgressBar(pb, s)
	}
		}

	psexFR_p <- cbind(psexFR, as.data.frame(matrix(NA, ncol = 1, nrow = nrow(data1))))
	for (m in 1:length(list_genet[which(sapply(list_genet, length) > 1)])){
		if (MLGsim){
			for (l in list_genet[which(sapply(list_genet, length) > 1)][[m]]){
				psexFR_p[l,3] <- mean(psexFR_p[l,2] > psex_recup)
			}
	} else {
			for (l in list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]){
				psexFR_p[l,3] <- mean(psexFR_p[l,2] > psex_recup)
			}
	}
	}
	colnames(psexFR_p) <- c("genet", "psex", "pvalue")

	if (bar){
		close(pb)
	}

	psexFR_p[is.na(psexFR_p)] <- ""
	if (length(psex_recup) < 100){print("Warning: Simulated populations contain few repeated genotypes and p-value estimations may be incorrect.")}
	list(psexFR_p, psex_recup)
	} else {
	psexFR[is.na(psexFR)] <- ""
	psexFR
	}
	}
}


psex <- function(data1, haploid = FALSE, vecpop = NULL, genet = FALSE, RR = FALSE, MLGsim = FALSE, nbrepeat = NULL, bar = FALSE){

	if (genet & RR){stop("Round-Robin method is not compatible with genet option.")}

	if (RR){
		datafreq <- freq_RR(data1, haploid, vecpop, genet = FALSE, RR = TRUE)
	} else 
	if (genet){
		datafreq <- freq_RR(data1, haploid, vecpop, genet = TRUE, RR = FALSE)
	} else {
		datafreq <- freq_RR(data1, haploid, vecpop, genet = FALSE, RR = FALSE)
	}

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		res <- list(NULL)
			for(p in 1:length(unique(vecpop))){
				rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
				res[[p]] <- psex_core(data1 = datatot[[p]], data2 = datafreq[[p]], haploid, MLGsim, nbrepeat, bar)
			}
				names(res) <- unique(vecpop_o)
		} else {
			 res <- psex_core(data1, data2 = datafreq, haploid, MLGsim, nbrepeat, bar)
		}
	res
}


#########################################################################
#probability of one genotype with repro events and H-W deviation: psexFis
#########################################################################

psex_Fis_core <- function(data1, data2, MLGsim = FALSE, genet = FALSE, RR = FALSE, nbrepeat = NULL, bar = FALSE){

	ncol_all <- 2
	ncol_freq <- 3
	factoR <- factoR
	pgenFis <- pgen_Fis_core(data1, data2, genet, RR)
	tab <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(data1)))
	list_genet <- MLG_list(data1)

	if (length(list_genet[which(sapply(list_genet, length) > 1)]) > 1){

	if (MLGsim){
		for (m in 1:length(list_genet[which(sapply(list_genet, length) > 1)])){
			recup <- NULL
			sub_list <- list_genet[which(sapply(list_genet, length) > 1)][[m]]
			l <- sub_list[1]
			recup <- c(recup, factoR[nrow(data1),2]/(factoR[length(sub_list), 2]*factoR[c(nrow(data1)-length(sub_list)),2])*
				(pgenFis[l,])^(length(sub_list))*
				(1-pgenFis[l,])^(nrow(data1)-length(sub_list)))
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]],2] <- recup
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]],1] <- list_genet[which(sapply(list_genet, length) > 1)][[m]][1]
		}
	} else {
		for (m in 1:length(list_genet[which(sapply(list_genet, length) > 1)])){
			recup <- NULL
				for (l in list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]){
					recup <- c(recup, factoR[nrow(data1),2]/(factoR[which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l),2]*factoR[c(nrow(data1)-which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l)),2])*
					(pgenFis[l,])^(which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l))*
					(1-pgenFis[l,])^(nrow(data1)-which(list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]==l)))
				}
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]][-1],2] <- recup
			tab[list_genet[which(sapply(list_genet, length) > 1)][[m]][-1],1] <- list_genet[which(sapply(list_genet, length) > 1)][[m]][1]
		}
	}
	names(tab) <- c("genet", "psexFis")
	psexFR_Fis <- tab

	if (length(nbrepeat) != 0){
		psexFis_recup <- NULL
		nb_loci <- ncol(data1)/2
		index_l <- 1:c(ncol(data1)/2)*2-1
		N <- nrow(data1)

	if (bar){
		total <- nbrepeat
		pb <- txtProgressBar(min = 0, max = total, style = 3)
	}

		for (s in 1:nbrepeat){
			tab2 <- as.data.frame(matrix(NA, nrow = N, ncol = nb_loci*2))
				for (j in 1:N){
					for (i in index_l){
						tab2[j, c(i, i+1)] <- sort(sample(split(data2, data2[,1])[[c((i+1)/2)]][,ncol_all], 2, prob = split(data2, data2[,1])[[c((i+1)/2)]][,ncol_freq], replace = TRUE))
					}
				}
			tab_sim <- tab2
			colnames(tab_sim) <- colnames(data1)
			MLG_sim <- unique(tab_sim)
				if (nrow(MLG_sim) != nrow(tab_sim)){
					pgenFis_sim <- pgen_Fis_core(tab_sim, data2)
					psexFis_sim <- NULL
					list_genet_sim <- MLG_list(tab_sim)
	if (MLGsim){
						for (m in 1:length(list_genet_sim[which(sapply(list_genet_sim, length) > 1)])){
							recup <- NULL
							sub_list <- list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]]
							l <- sub_list[1]
							recup <- c(recup, factoR[nrow(tab_sim),2]/(factoR[length(sub_list), 2]*factoR[c(nrow(tab_sim)-length(sub_list)),2])*
								(pgenFis_sim[l,])^(length(sub_list))*
								(1-pgenFis_sim[l,])^(nrow(tab_sim)-length(sub_list)))
						}
					psexFis_sim <- c(psexFis_sim, recup)
	} else {
						for (m in 1:length(list_genet_sim[which(sapply(list_genet_sim, length) > 1)])){
							recup <- NULL
								for (l in list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]){
									recup <- c(recup, factoR[nrow(tab_sim), 2]/(factoR[which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l), 2]*factoR[c(nrow(tab_sim)-which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l)),2])*
										(pgenFis_sim[l,])^(which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l))*
										(1-pgenFis_sim[l,])^(nrow(tab_sim)-which(list_genet_sim[which(sapply(list_genet_sim, length) > 1)][[m]][-1]==l)))
								}
							psexFis_sim <- c(psexFis_sim, recup)
						}
	}
					psexFis_recup <- c(psexFis_recup, psexFis_sim)
				}

	if (bar){
		setTxtProgressBar(pb, s)
	}
		}
	if (length(psexFis_recup) == 0){stop("No clone was found during Simulations.")}

	psexFR_p <- cbind(psexFR_Fis, as.data.frame(matrix(NA, ncol = 1, nrow = nrow(data1))))
		for (m in 1:length(list_genet[which(sapply(list_genet, length) > 1)])){
			if (MLGsim){
				for (l in list_genet[which(sapply(list_genet, length) > 1)][[m]]){
					psexFR_p[l,3] <- mean(psexFR_p[l,2] > psexFis_recup)
				}
			} else {
				for (l in list_genet[which(sapply(list_genet, length) > 1)][[m]][-1]){
					psexFR_p[l,3] <- mean(psexFR_p[l,2] > psexFis_recup)
				}
		}
	}
	colnames(psexFR_p) <- c("genet", "psexFis", "pvalue")

	if (bar){
		close(pb)
	}
	psexFR_p[is.na(psexFR_p)] <- ""
		if (length(psexFis_recup) < 100){print("Warning: Simulated populations contain few repeated genotypes and p-value estimations may be incorrect.")}
	list(psexFR_p, psexFis_recup)
	} else {
	psexFR_Fis[is.na(psexFR_Fis)] <- ""
	psexFR_Fis
	}
	}
}


psex_Fis <- function(data1, vecpop = NULL, genet = FALSE, RR = FALSE, MLGsim = FALSE, nbrepeat = NULL, bar = FALSE){

	if (genet & RR){stop("Round-Robin method is not compatible with genet option.")}

	if (RR){
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = FALSE, RR = TRUE)
	} else 
	if (genet){
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = TRUE, RR = FALSE)
	} else {
		datafreq <- freq_RR(data1, haploid = FALSE, vecpop, genet = FALSE, RR = FALSE)
	}

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		res <- list(NULL)
			for (p in 1:length(unique(vecpop))){
				rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
				res[[p]] <- psex_Fis_core(datatot[[p]], data2 = datafreq[[p]], MLGsim, genet, RR, nbrepeat, bar)
			}
				names(res) <- unique(vecpop_o)
		} else {
			 data2 <- datafreq
			 res <- psex_Fis_core(data1, data2, MLGsim, genet, RR, nbrepeat, bar)
		}
	res
}


##############################################################################
################################Resampling####################################
##############################################################################

###########################
#Resampling l loci X fois :
###########################


sample_loci_core <- function(data1, data2, haploid = FALSE, nbrepeat = 1000, He = FALSE, graph = FALSE, bar = FALSE){

	if (He & haploid){stop("Haploids and He computations are not compatible.")}

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	index_l <- 1:c(ncol(data1)/2)*2-1 
	nb_loci <- ncol(data1)/2
	ncol_freq <- 3

	mat_box_l <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nbrepeat))
	mat_box_l2 <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nbrepeat))
	mat_box_l3 <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nbrepeat)) 
	mat_res_MLG_l <- as.data.frame(matrix(NA, ncol = 5, nrow = nb_loci))
		if (He){
			He_box_l <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nbrepeat))
			He_box_l2 <- as.data.frame(matrix(NA, ncol = nb_loci, nrow = nbrepeat))
			mat_res_all_l <- as.data.frame(matrix(NA, ncol = 7, nrow = nb_loci))
		} else {
			mat_res_all_l <- as.data.frame(matrix(NA, ncol = 5, nrow = nb_loci))
		}	

if (bar){
	total <- max(index_l)
	pb <- txtProgressBar(min = 0, max = total, style = 3)
}


for (j in 1:nbrepeat){
	for (i in index_l){
		indic <- sample(index_l, which(index_l == i))
		mat_test <- data1[,sort(c(indic, indic+1))]
			if (He){
				He_ind <- NULL
					for (l in indic){
						He_ind <- c(He_ind, 2*nrow(mat_test)/(2*nrow(mat_test)-1)*(1 - sum(as.numeric(data2[rownames(data2[data2==paste("locus",c((l+1)/2), sep="_"),]), ncol_freq])^2)))
					}
				He_box_l[j,c((i+1)/2)] <- mean(He_ind)
				He_box_l2[j,c((i+1)/2)] <- var(He_ind)
				colnames(He_box_l) <- colnames(He_box_l2) <- c("1_locus", paste(2:nb_loci, "loci", sep = "_"))
			}
		mat_test <- unique(mat_test)
		mat_box_l[j,c((i+1)/2)] <- nrow(mat_test)
		nb_loc <- NULL
			for (k in index_l[1:c(ncol(mat_test)/2)]){
				nb_loc <- c(nb_loc, length(unique(unlist(mat_test[,c(k,k+1)]))))
			}
		mat_box_l2[j,c((i+1)/2)] <- sum(nb_loc)
		mat_box_l3[j,c((i+1)/2)] <- var(nb_loc)
	}

if (bar){
	setTxtProgressBar(pb, i)
}
}

if (bar){
	close(pb)
}

	colnames(mat_box_l) <- colnames(mat_box_l2) <- colnames(mat_box_l3) <- c("1_locus", paste(2:nb_loci, "loci", sep = "_"))

	mat_res_MLG_l[,1] <- 1:nb_loci
	mat_res_MLG_l[,2] <- apply(mat_box_l, 2, min)
	mat_res_MLG_l[,3] <- apply(mat_box_l, 2, max)
	mat_res_MLG_l[,4] <- apply(mat_box_l, 2, mean)
	mat_res_MLG_l[,5] <- apply(mat_box_l, 2, function(x) sd(x)/sqrt(length(x)))
	colnames(mat_res_MLG_l) <- c("nb_loci", "min", "max", "mean_MLG", "SE")

	mat_res_all_l[,1] <- 1:nb_loci
	mat_res_all_l[,2] <- apply(mat_box_l2, 2, min)
	mat_res_all_l[,3] <- apply(mat_box_l2, 2, max)
	mat_res_all_l[,4] <- apply(mat_box_l2, 2, mean)
	mat_res_all_l[,5] <- sqrt(apply(mat_box_l3, 2, function(x) sum(x^2))/(1:nb_loci))
	mat_res_all_l[c(1, nb_loci),5] <- NA

	if (He){
		mat_res_all_l[,6] <- apply(He_box_l, 2, mean)
		mat_res_all_l[,7] <- sqrt(apply(He_box_l2, 2, function(x) sum(x^2))/(1:nb_loci))
		mat_res_all_l[c(1, nb_loci), 7] <- NA
		colnames(mat_res_all_l) <- c("nb_loci", "min", "max", "mean_all", "SE", "He", "SE")
	} else {
		colnames(mat_res_all_l) <- c("nb_loci", "min", "max", "mean_all", "SE")
	}
	
if (graph){
	boxplot(mat_box_l, ylab = "Number of multilocus genotypes", 
	xlab = "Number of loci sampled")
	title(paste("Genotype accumulation curve"))
}

if (He){
	res <- list("res_MLG" = mat_res_MLG_l, "res_alleles" = mat_res_all_l, "raw_He" = He_box_l, "raw_MLG" = mat_box_l, "raw_all" = mat_box_l2)
	} else {
	res <- list("res_MLG" = mat_res_MLG_l, "res_alleles" = mat_res_all_l, "raw_MLG" = mat_box_l, "raw_all" = mat_box_l2)
	}
res
}


sample_loci <- function(data1, haploid = FALSE, vecpop = NULL, nbrepeat = 1000, He = FALSE, graph = FALSE, export = FALSE, bar = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
			datafreq <- freq_RR(data1, haploid, vecpop, genet = FALSE, RR = FALSE)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					res[[p]] <- sample_loci_core(datatot[[p]], data2 = datafreq[[p]], haploid, nbrepeat, He, graph = FALSE, bar)
					
					if (graph){
						par(ask = TRUE)
						boxplot(res[[p]]$raw_MLG, ylab = "Number of multilocus genotypes", 
							xlab = "Number of loci sampled")
						title(paste("Genotype accumulation curve"))
						par(ask = FALSE)
					}
					if (export){
						postscript(file = paste(paste("sample_loci", unique(vecpop_o)[[p]], sep = "_"), ".eps", sep = ""), onefile = FALSE, paper = "letter")
						boxplot(res[[p]]$raw_MLG, ylab = "Number of multilocus genotypes", 
							xlab = "Number of loci sampled")
						title(paste("Genotype accumulation curve"))
						dev.off()
					}
				}
					names(res) <- unique(vecpop_o)
		} else {
			data2 <- freq_RR_core(data1, haploid, genet = FALSE, RR = FALSE)
			res <- sample_loci_core(data1, data2, haploid, nbrepeat, He, graph, bar)
				if (export){
					postscript(file = "sample_loci.eps", onefile = FALSE, paper = "letter")
					boxplot(res$raw_MLG, ylab = "Number of multilocus genotypes", 
						xlab = "Number of loci sampled")
					title(paste("Genotype accumulation curve"))
					dev.off()
				}
		}
	res
}


################################
#Resampling u individus X fois :
################################

sample_units_core <- function(data1, haploid = FALSE, nbrepeat = 1000, He = FALSE, graph = FALSE, bar = FALSE){

	if (He & haploid){stop("Haploids and He computations are not compatible.")}

	if (haploid){
		data1 <- cbind(data1, data1)
		data1 <- data1[, sort(rep(1:(ncol(data1)/2), 2))]
	}

	index_l <- 1:c(ncol(data1)/2)*2-1 
	nb_loci <- ncol(data1)/2
	L <- nrow(data1)

	mat_box_u <- as.data.frame(matrix(NA, ncol = L, nrow = nbrepeat))
	mat_box_u2 <- as.data.frame(matrix(NA, ncol = L, nrow = nbrepeat))
	mat_box_u3 <- as.data.frame(matrix(NA, ncol = L, nrow = nbrepeat))
	mat_res_MLG_u <- as.data.frame(matrix(NA, ncol = 5, nrow = L))

		if (He){
	He_box_u <- as.data.frame(matrix(NA, ncol = L, nrow = nbrepeat))
	He_box_u2 <- as.data.frame(matrix(NA, ncol = L, nrow = nbrepeat))
	mat_res_all_u <- as.data.frame(matrix(NA, ncol = 7, nrow = L))
		} else {
	mat_res_all_u <- as.data.frame(matrix(NA, ncol = 5, nrow = L))
		}	

if (bar){
	total <- L
	pb <- txtProgressBar(min = 0, max = total, style = 3)
}


for (s in 1:nbrepeat){
	for(u in 1:L){
		indic <- sample(1:L, u, replace = FALSE)
		mat_test <- data1[indic,]

		mat_freq_temp <- freq_RR(mat_test, haploid, genet = FALSE, RR = FALSE)
			
			if (He | haploid){
			} else
			if (He){
				He_unit <- NULL
					for (l in index_l){ 
						He_unit <- c(He_unit, 2*nrow(mat_test)/(2*nrow(mat_test)-1)*(1 - sum(mat_freq_temp[rownames(mat_freq_temp[mat_freq_temp==paste("locus",c((l+1)/2), sep="_"),]),3]^2)))
					}
				He_box_u[s,u] <- mean(He_unit)
				He_box_u2[s,u] <- var(He_unit)
				colnames(He_box_u) <- colnames(He_box_u2) <- c("1_unit", paste(2:L, "units", sep = "_"))
			}

		mat_test <- unique(mat_test) 
		mat_box_u[s,u] <- nrow(mat_test)
		nb_loc <- NULL
			for (k in index_l){
				nb_loc <- c(nb_loc, length(unique(unlist(mat_test[,c(k,k+1)]))))
			}
		mat_box_u2[s,u] <- sum(nb_loc)
		mat_box_u3[s,u] <- var(nb_loc)
	}

if (bar){
	setTxtProgressBar(pb, u)
}
}

if (bar){
	close(pb)
}

	colnames(mat_box_u) <- 1:L
	mat_res_MLG_u[,1] <- 1:L
	mat_res_MLG_u[,2] <- apply(mat_box_u, 2, min)
	mat_res_MLG_u[,3] <- apply(mat_box_u, 2, max)
	mat_res_MLG_u[,4] <- apply(mat_box_u, 2, mean)
	mat_res_MLG_u[,5] <- apply(mat_box_u, 2, function(x) sd(x)/sqrt(length(x)))
	colnames(mat_res_MLG_u) <- c("nb_units", "min", "max", "mean_MLG", "SE")

	mat_res_all_u[,1] <- 1:L
	mat_res_all_u[,2] <- apply(mat_box_u2, 2, min)
	mat_res_all_u[,3] <- apply(mat_box_u2, 2, max)
	mat_res_all_u[,4] <- apply(mat_box_u2, 2, mean)
	mat_res_all_u[,5] <- sqrt(apply(mat_box_u3, 2, function(x) sum(x^2))/(1:L))
	mat_res_all_u[c(1, L),5] <- NA

	colnames(mat_box_u) <- colnames(mat_box_u2) <- colnames(mat_box_u3) <- c("1_unit", paste(2:L, "units", sep = "_"))

	if (He){
		mat_res_all_u[,6] <- apply(He_box_u, 2, mean)
		mat_res_all_u[,7] <- sqrt(apply(He_box_u2, 2, function(x) sum(x^2))/(1:L))
		mat_res_all_u[c(1, L) ,7] <- NA
		colnames(mat_res_all_u) <- c("nb_units", "min", "max", "mean_units", "SE", "He", "SE")
	} else {
		colnames(mat_res_all_u) <- c("nb_units", "min", "max", "mean_units", "SE")
	}
	
if (graph){
	boxplot(mat_box_u, range = 3, ylab = "Number of multilocus genotypes", 
	xlab = "Number of units sampled")
	title(paste("Genotype accumulation curve"))
}

if (He){
	res <- list("res_MLG" = mat_res_MLG_u, "res_alleles" = mat_res_all_u, "raw_He" = He_box_u, "raw_MLG" = mat_box_u, "raw_all" = mat_box_u2)
	} else {
	res <- list("res_MLG" = mat_res_MLG_u, "res_alleles" = mat_res_all_u, "raw_MLG" = mat_box_u, "raw_all" = mat_box_u2)
	}
res
}


sample_units <- function(data1, haploid = FALSE, vecpop = NULL, nbrepeat = 1000, He = FALSE, graph = FALSE, export = FALSE, bar = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					par(ask = TRUE)
					res[[p]] <- sample_units_core(datatot[[p]], haploid, nbrepeat, He, graph, bar)
					par(ask = FALSE)

					if (export){
						postscript(paste(paste("sample_units", unique(vecpop)[[p]], sep = "_"), ".eps", sep = ""), onefile = FALSE, paper = "letter")
						boxplot(res[[p]]$raw_MLG, ylab = "Number of multilocus genotypes", 
							xlab = "Number of units sampled")
						title(paste("Genotype accumulation curve"))
						dev.off()
					}
				}
					names(res) <- unique(vecpop_o)
		} else {
			res <- sample_units_core(data1, haploid, nbrepeat, He, graph, bar)
				if (export){
					postscript("sample_units.eps", onefile = FALSE, paper = "letter")
					boxplot(res[[p]]$raw_MLG, ylab = "Number of multilocus genotypes", 
						xlab = "Number of units sampled")
					title(paste("Genotype accumulation curve"))
					dev.off()
				}
		}
	res
}


##############################################################################
#####################Genetic Distance between units###########################
##############################################################################

##################
#Genetic distances
##################


genet_dist_core <- function(data1, haploid = FALSE, manh = FALSE, manh_w = FALSE, graph = FALSE, breaking = NULL, 
				alpha1 = NULL, alpha2 = NULL){

	if (length(alpha1) != 0 & length(alpha2) != 0){stop("Genetic distances : must choose between percentage of distribution (alpha1) and numeric floor (alpha2).")}
	if (length(alpha1) != 0){
		if (alpha1 < 0 | alpha1 >= 1){stop("Wrong definition of alpha1 (pourcentage between 0 and 1).")}
		if (alpha1 == 0){alpha1 <- NULL}
	}
	if (length(alpha2) != 0){
		if (alpha2 < 0) {stop("Wrong definition of alpha2.")}
		if (alpha2 == 0){alpha2 <- NULL}
	}
	
        tab_MLG <- unique(data1)
        G <- nrow(tab_MLG)
        tab <- matrix(NA, ncol = G, nrow = G)
        colnames(tab) <- rownames(tab) <- 1:G

	if (haploid){
        index_l <- 1:ncol(data1)
        nb_loci <- ncol(data1)
	} else {
        index_l <- 1:c(ncol(data1)/2)*2-1
        nb_loci <- ncol(data1)/2
	}

if (manh){
                for (j in 1:G){
                        for (i in j:G){
                                loc1 <- as.numeric(tab_MLG[i,])
                                loc2 <- as.numeric(tab_MLG[j,])
                                recup <- NULL

if (haploid){
                                        for (k in index_l){
                                                recup <- c(recup, abs(loc1[k]-loc2[k]))
                                        }
} else {
                                        for (k in index_l){
                                                recup <- c(recup, min(abs(loc1[k]-loc2[k])+ abs(loc1[k+1]-loc2[k+1]),
                                                        abs(loc1[k]-loc2[k+1])+ abs(loc1[k+1]-loc2[k])))
                                        }
}
                                tab[i,j] <- sum(recup)
                        }
                }

} else 
if (manh_w){
                for (j in 1:G){
                        for (i in j:G){
                                loc1 <- as.numeric(tab_MLG[i,])
                                loc2 <- as.numeric(tab_MLG[j,])
                                recup <- NULL
if (haploid){
                                        for (k in index_l){
                                                recup <- c(recup, abs(loc1[k]-loc2[k]))
                                        }
} else {
                                        for (k in index_l){
                                                recup <- c(recup, min(abs(loc1[k]-loc2[k])+ abs(loc1[k+1]-loc2[k+1]),
                                                        abs(loc1[k]-loc2[k+1])+ abs(loc1[k+1]-loc2[k])))
                                        }
}
                                tab[i,j] <- sum(recup)/nb_loci
                        }
                }
} else {
                for (j in 1:G){
                        for (i in j:G){
                                tab[i,j] <- length(which(tab_MLG[i,] != tab_MLG[j,]))
                                        }
                        }
}

dist_all <- tab

if (length(alpha1) != 0 | length(alpha2) != 0){
        tab2 <- as.data.frame(matrix(0, ncol = 3, nrow = 1))
                for (l in 1:(G-1)){
                        sub_tab <- cbind((l+1):G, rep(l, times = length((l+1):G)), dist_all[(l+1):G, l])
                        tab2 <- rbind(tab2, sub_tab)
                }
        tab2 <- tab2[-1,]
        names(tab2) <- c("unit_1", "unit_2", "dist")
        tab_sort <- tab2[order(tab2[,3]),]
	if (is.null(alpha2)){
        	alpha <- alpha1
        	indic <- round(alpha*length(as.dist(dist_all)))
		indic_et <- indic
			while(tab_sort[indic_et,3] == tab_sort[indic_et+1,3]){
                        indic_et <- indic_et+1
			}
		tab_sort_red <- tab_sort[1:indic_et,]
	} else {
		beta <- alpha2
		if (length(which(beta == tab_sort[,3])) == 0) {stop("alpha2 is not a genetic distance in matrix distance")}
		indic <- which(beta == tab_sort[,3])[length(which(beta == tab_sort[,3]))]
		tab_sort_red <- tab_sort[1:indic,]
	}
		convert <- cbind(1:G, as.numeric(rownames(tab_MLG)))
		colnames(convert) <- c("genets", "ramets")
			for (i in 1:nrow(tab_sort_red)){
				tab_sort_red[i,1] <- convert[which(tab_sort_red[i,1] == convert[,1]),2]
				tab_sort_red[i,2] <- convert[which(tab_sort_red[i,2] == convert[,1]),2]
			}
}

if (graph){
        fin <- round(max(as.dist(dist_all)))+1
                if (length(breaking) != 0){
				pas <- breaking
                        hist(as.dist(dist_all), breaks = seq(0, fin, pas))
                                if (length(alpha1) != 0 | length(alpha2) != 0){
                                        if(round(tab_sort[indic,3]) == tab_sort[indic,3]){
                                                        if (tab_sort[indic,3] != 0){
                                                                abline(v = tab_sort[indic,3]-pas/2, col = "red", lwd = 2)
                                                        } else {
                                                                abline(v = tab_sort[indic,3]+pas/2, col = "red", lwd = 2)
                                                        }
                                        } else {
                                                abline(v = tab_sort[indic,3], col = "red", lwd = 2)
                                        }
                                }                               
                } else {
                        hist(as.dist(dist_all), breaks = 0:fin)
                                if (length(alpha1) != 0 | length(alpha2) != 0){
                                        if (round(tab_sort[indic,3]) == tab_sort[indic,3]){
                                                        if (tab_sort[indic,3] != 0){
                                                                abline(v = tab_sort[indic,3]-1/2, col = "red", lwd = 2)
                                                        } else {
                                                                abline(v = tab_sort[indic,3]+1/2, col = "red", lwd = 2)
                                                        }
                                        } else {
                                                abline(v = tab_sort[indic,3], col = "red", lwd = 2)
                                        }
                                }
                }
}

if (length(alpha1) != 0 | length(alpha2) != 0){
        list("distance_matrix" = as.dist(dist_all), "potential_clones" = tab_sort_red, "all_pairs" = tab_sort, "sign" = tab_sort[indic,3])
} else {
        list("distance_matrix" = as.dist(dist_all))
}
}


graph_genet_dist <- function(dist_all, breaking, tab_sort, indic, alpha1 = NULL, 
				alpha2 = NULL){

	dist_all <- dist_all
	pas <- breaking 
	alpha1 <- alpha1
	alpha2 <- alpha2

if(length(alpha1) != 0) {if (alpha1 == 0) {alpha1 <- NULL}}
if(length(alpha2) != 0) {if (alpha2 == 0) {alpha2 <- NULL}}

if (length(alpha1) != 0 | length(alpha2) != 0){
	tab_sort <- tab_sort
	indic <- indic
}

        fin <- round(max(as.dist(dist_all)))+1
                if (length(pas) != 0){
                        hist(as.dist(dist_all), breaks = seq(0, fin, pas))
                                if (length(alpha1) != 0 | length(alpha2) != 0){
                                        if(round(tab_sort[indic,3]) == tab_sort[indic,3]){
                                                        if (tab_sort[indic,3] != 0){
                                                                abline(v = tab_sort[indic,3]-pas/2, col = "red", lwd = 2)
                                                        } else {
                                                                abline(v = tab_sort[indic,3]+pas/2, col = "red", lwd = 2)
                                                        }
                                        } else {
                                                abline(v = tab_sort[indic,3], col = "red", lwd = 2)
                                        }
                                }                               
                } else {
                        hist(as.dist(dist_all), breaks = 0:fin)
                                if (length(alpha1) != 0 | length(alpha2) != 0){
                                        if (round(tab_sort[indic,3]) == tab_sort[indic,3]){
                                                        if (tab_sort[indic,3] != 0){
                                                                abline(v = tab_sort[indic,3]-1/2, col = "red", lwd = 2)
                                                        } else {
                                                                abline(v = tab_sort[indic,3]+1/2, col = "red", lwd = 2)
                                                        }
                                        } else {
                                                abline(v = tab_sort[indic,3], col = "red", lwd = 2)
                                        }
                                }
                }
}


genet_dist <- function(data1, haploid = FALSE, vecpop = NULL, manh = FALSE, manh_w = FALSE, graph = FALSE, breaking = NULL, 
			alpha1 = NULL, alpha2 = NULL, export = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		if (length(alpha1) != 0){
			if(length(alpha1) != length(unique(vecpop))) {stop("alpha1 length does not fit with vecpop length.")}
			alpha11 <- split(alpha1, unique(vecpop))
			alpha21 <- NULL
		} else 
		if (length(alpha2) != 0){
			if(length(alpha2) != length(unique(vecpop))) {stop("alpha2 length does not fit with vecpop length.")}
			alpha21 <- split(alpha2, unique(vecpop))
			alpha11 <- NULL
		} else {
			alpha11 <- NULL
			alpha21 <- NULL
		}
		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					par(ask = TRUE)
					res[[p]] <- genet_dist_core(datatot[[p]], haploid, manh, manh_w, graph, breaking, alpha11[[p]], alpha21[[p]])
					par(ask = FALSE)

					if (export){
						postscript(paste(paste("genet_dist", unique(vecpop)[[p]], sep = "_"), ".eps", sep = ""), onefile = FALSE, paper = "letter")
							if (length(alpha1) != 0 | length(alpha2) != 0){
								graph_genet_dist(dist_all = res[[p]]$distance_matrix, breaking, tab_sort = res[[p]]$all_pairs, indic = res[[p]]$sign[[1]], alpha11[[p]], alpha21[[p]])
							} else {
								graph_genet_dist(dist_all = res[[p]]$distance_matrix, breaking)
							}
						dev.off()
					}
				}
			names(res) <- unique(vecpop_o)
		} else {
			res <- genet_dist_core(data1, haploid, manh, manh_w, graph, breaking, alpha1, alpha2)
				if (export){
					postscript("genet_dist.eps", onefile = FALSE, paper = "letter")
						if (length(alpha1) != 0 | length(alpha2) != 0){
							graph_genet_dist(dist_all = res$distance_matrix, breaking, tab_sort = res$all_pairs, indic = res$sign[[1]], alpha1, alpha2)
						} else {
							graph_genet_dist(dist_all = res$distance_matrix, breaking)
						}				
					dev.off()
				}
		}
	res
}


genet_dist_sim_core <- function(data1, haploid = FALSE, nbrepeat = 1000, genet = FALSE, manh = FALSE, manh_w = FALSE, graph = FALSE, 
					breaking = NULL){
        
	tab_MLG <- unique(data1)
	G <- nrow(tab_MLG)
	N <- nrow(data1)
        tab_sim <- matrix(NA, ncol = ncol(data1), nrow = nbrepeat)
        colnames(tab_sim) <- colnames(data1)
        rownames(tab_sim) <- 1:nbrepeat

		if (haploid){
        index_l <- 1:ncol(data1)
        nb_loci <- ncol(data1)
		} else {
        index_l <- 1:c(ncol(data1)/2)*2-1
        nb_loci <- ncol(data1)/2
		}

if (genet){
                for (s in 1:nbrepeat){
				indic <- sample(1:G, 2, replace = FALSE)
                                for (l in index_l){
						if (haploid){
							tab_sim[s,l] <- sample(data1[indic,l], 1)
						} else {
                                        	tab_sim[s,c(l, l+1)] <- sort(c(as.numeric(sample(tab_MLG[indic[[1]],c(l, l+1)], 1)), as.numeric(sample(tab_MLG[indic[[2]],c(l, l+1)], 1))))
						}	
                                }
                }
} else {
                for (s in 1:nbrepeat){
				indic <- sample(1:N, 2, replace = FALSE)
                                for (l in index_l){
							if (haploid){
								tab_sim[s,l] <- sample(data1[indic,l], 1)
							} else {
                                        		tab_sim[s,c(l, l+1)] <- sort(c(as.numeric(sample(data1[indic[[1]],c(l, l+1)], 1)), as.numeric(sample(data1[indic[[2]],c(l, l+1)], 1))))
							}
                                }
                }
}

if (nrow(unique(tab_sim)) != nbrepeat) {print(paste("Number of MLG sim = ", nrow(unique(tab_sim)), sep = ""))}

        tab_sim <- unique(tab_sim)
        genet_dist_core(tab_sim, haploid, manh, manh_w, graph, breaking, alpha1 = NULL, alpha2 = NULL)
}


genet_dist_sim <- function(data1, haploid = FALSE, vecpop = NULL, nbrepeat = 1000, genet = FALSE, manh = FALSE, manh_w = FALSE, graph = FALSE, 
				breaking = NULL, export = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)){stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					par(ask = TRUE)
					res[[p]] <- genet_dist_sim_core(datatot[[p]], haploid, nbrepeat, genet, manh, manh_w, graph, breaking)
					par(ask = FALSE)

					if (export){
						postscript(paste(paste("genet_dist_sim", unique(vecpop)[[p]], sep = "_"), ".eps", sep = ""), onefile = FALSE, paper = "letter")
							graph_genet_dist(dist_all = res[[p]]$distance_matrix, breaking)
						dev.off()
					}
				}
			names(res) <- unique(vecpop_o)
		} else {
			res <- genet_dist_sim_core(data1, haploid, nbrepeat, genet, manh, manh_w, graph, breaking)
				if (export){
					postscript("genet_dist_sim.eps", onefile = FALSE, paper = "letter")
						graph_genet_dist(dist_all = res$distance_matrix, breaking)			
					dev.off()
				}
		}
	res
}


MLL_generator_core <- function(data1, haploid = FALSE, manh = FALSE, manh_w = FALSE, alpha1 = NULL, alpha2 = NULL){

	resmlg <- MLG_list(data1)

	if (length(alpha1) == 0 & length(alpha2) == 0) {resmlg <- resmlg
	} else {
		if (length(alpha1) != 0){
			if (alpha1 != 0){
				resdist <- genet_dist_core(data1, haploid, manh, manh_w, graph = FALSE, alpha1 = alpha1, alpha2 = alpha2)[[2]]

				for (j in 1:nrow(resdist)){
					unit <- resdist[j,1:2]
						for (i in 1:length(resmlg)){
							if (length(which(resmlg[[i]] == unit[[2]])) != 0){
								nb1 <- i
							}
							if (length(which(resmlg[[i]] == unit[[1]])) != 0){
								nb2 <- i
							}
						}

					if (nb1 != nb2){
						resmlg[[nb1]] <- c(resmlg[[nb1]], resmlg[[nb2]])
						resmlg <- resmlg[-nb2]	
					}
				}
			} else {resmlg <- resmlg}
		}

		if (length(alpha2) != 0){
			if (alpha2 != 0){
				resdist <- genet_dist_core(data1, haploid, manh, manh_w, graph = FALSE, alpha1 = alpha1, alpha2 = alpha2)[[2]]

				for (j in 1:nrow(resdist)){
					unit <- resdist[j,1:2]
						for (i in 1:length(resmlg)){
							if (length(which(resmlg[[i]] == unit[[2]])) != 0){
								nb1 <- i
							}
							if (length(which(resmlg[[i]] == unit[[1]])) != 0){
								nb2 <- i
							}
						}

					if (nb1 != nb2){
						resmlg[[nb1]] <- c(resmlg[[nb1]], resmlg[[nb2]])
						resmlg <- resmlg[-nb2]	
					}
				}
			} else {resmlg <- resmlg}
		}
	}
	resmlg
}


MLL_generator <- function(data1, haploid = FALSE, vecpop = NULL, manh = FALSE, manh_w = FALSE, alpha1 = NULL, alpha2 = NULL){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		if (length(alpha1) != 0){
			if(length(alpha1) != length(unique(vecpop))) {stop("alpha1 length does not fit with vecpop length.")}
			alpha11 <- split(alpha1, unique(vecpop))
			alpha21 <- NULL
		} else
		if (length(alpha2) != 0){
			if(length(alpha2) != length(unique(vecpop))) {stop("alpha2 length does not fit with vecpop length.")}
			alpha21 <- split(alpha2, unique(vecpop))
			alpha11 <- NULL
		} else {
			alpha11 <- alpha21 <- NULL
		}
		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					res[[p]] <- MLL_generator_core(datatot[[p]], haploid, manh, manh_w, alpha11[[p]], alpha21[[p]])
				}
			names(res) <- unique(vecpop_o)
		} else {
			res <- MLL_generator_core(data1, haploid, manh, manh_w, alpha1, alpha2)
		}
	res
}


MLL_generator2_core <- function(potential_clones = NULL, res_mlg = NULL){

resdist <- potential_clones
resmlg <- res_mlg

if (length(resdist) != 0){
	for (j in 1:nrow(resdist)){
		unit <- resdist[j,1:2]
			for (i in 1:length(resmlg)){
				if (length(which(resmlg[[i]] == unit[[2]])) != 0){
				nb1 <- i
				}
				if (length(which(resmlg[[i]] == unit[[1]])) != 0){
				nb2 <- i
				}
			}

		if (nb1 != nb2){
			resmlg[[nb1]] <- c(resmlg[[nb1]], resmlg[[nb2]])
			resmlg <- resmlg[-nb2]	
		}
	}
}
resmlg
}


MLL_generator2 <- function(potential_clones = NULL, res_mlg = NULL, vecpop = NULL){

	if (length(vecpop) != 0){			
		if (length(unique(vecpop)) != length(res_mlg)) {stop("vecpop length does not fit with MLG.")}
		if (length(unique(vecpop)) != length(potential_clones)) {stop("vecpop length does not fit with clones.")}
		
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		res <- list(NULL)
			for(p in 1:length(unique(vecpop))){
				res[[p]] <- MLL_generator2_core(potential_clones[[p]], res_mlg[[p]])
			}
			names(res) <- unique(vecpop_o)
	} else {
		res <- MLL_generator2_core(potential_clones, res_mlg)
	}
	res
}

##############################################################################
################Richness/diversity/Evenness Index#############################
##############################################################################

clonal_index_core <- function(data1, listMLL = NULL){

	N <- nrow(data1)
	if (length(listMLL) != 0){
			list_genet <- listMLL
		if (length(unlist(list_genet)) != nrow(data1)){stop("MLL list does not compute")}
			G <- length(listMLL)
		} else {
			G <- nrow(unique(data1))
        		list_genet <- MLG_list(data1)
		}

        if (G != N){
                R <- (G-1)/(N-1) ##Clonal diversity Dorken & Eckert : R

        count_MLG <- unlist(lapply(list_genet, length))
        freq_MLG <- count_MLG/N

G0 <- 1/sum(freq_MLG*freq_MLG) ##Estimate of genotypic diversity/Stoddart
Ge <- 1/(sum(freq_MLG[which(count_MLG > 1)]*freq_MLG[which(count_MLG > 1)])+(sum(freq_MLG[which(count_MLG < 1)])/N)) ##Estimate of expected genotypic diversity under H-W

                Hpp <- -sum((count_MLG/N)*log(count_MLG/N)) ##Shannon-Wiener index estimator : Hpp
                Jp <- Hpp/log(G) ##Pielou evenness: Jp
        Ls <- sum((count_MLG*(count_MLG-1))/(N*(N-1))) ##Simpson unbiased : Ls
                Dp <- 1-Ls ##Simpson complement unbiased : Dp
        Dmin <- (((2*N-G)*(G-1))/(N*N))*(N/(N-1))
        Dmax <- ((G-1)/G)*(N/(N-1))
                V <- (Dp-Dmin)/(Dmax-Dmin) ##Simpson complement index : V
                Hill <- 1/Ls ##Reciprocal of Simpson index unbiased : 1/Ls HILL

        tab <- as.data.frame(matrix(c(G, R, Hpp, Jp, Dp, V, Hill), nrow = 1))

	if(length(listMLL) != 0){
		rownames(tab) <- "MLL"
	} else {
		rownames(tab) <- "MLG"
	}
	names(tab) <- c("G", "R", "H''", "J'", "D", "V", "Hill")
	tab       
}
}


clonal_index <- function(data1, vecpop = NULL, listMLL = NULL){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- as.character(vecpop)
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
			res <- as.data.frame(matrix(NA, ncol = 7, nrow = length(unique(vecpop))))
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					resu <- clonal_index_core(datatot[[p]], listMLL[[p]])
						if(length(resu) != 0){
							res[p,] <- resu
						} else {res[p,] <- rep(0, 7)}
					rownames(res)[p] <- unique(vecpop_o)[p]
				}
			names(res) <- c("G", "R", "H''", "J'", "D", "V", "Hill")
		} else {
			res <- clonal_index_core(data1, listMLL)
		}
	res
}


Pareto_index_core <- function(data1, listMLL = NULL, full = FALSE, graph = FALSE, legends = 1){
	
	N <- nrow(data1)

	if (length(listMLL) != 0){
		G <- length(listMLL)
	} else {
		G <- nrow(unique(data1))
	}

if (G != N){
		if (length(listMLL) != 0){
			list_genet <- listMLL
		} else {
			list_genet <- MLG_list(data1)
		}

if (length(unlist(list_genet)) != nrow(data1)){stop("MLG/MLL list does not compute")}

	count_MLG <- unlist(lapply(list_genet, length))
	freq_MLG <- count_MLG/N
	count_MLL <- as.data.frame(table(count_MLG))[,2]
	xi <- nb_geno <- as.numeric(levels(as.data.frame(table(count_MLG))[,1]))

	cumul_MLL <- NULL
		for (i in 1:length(count_MLL)){
				if (i == 1){cumul_MLL <- count_MLL[i]} else {
					cumul_MLL <- c(cumul_MLL, count_MLL[i]*i+cumul_MLL[i-1])
				}
		}

	coord_Pareto <- NULL
		for (i in 1:length(cumul_MLL)){
			if (i == 1){coord_Pareto <- N} else {
				coord_Pareto <- c(coord_Pareto, coord_Pareto[i-1]-count_MLL[i-1]*nb_geno[i-1])
			}
		}

	yi <- coord_Pareto/N

	Pareto <- -lm(log10(yi)~log10(xi))$coefficients[[2]]
	c_Pareto <- Pareto + 1
	coef <- summary(lm(log10(yi)~log10(xi)))$coefficients
	res_red <- list(Pareto, c_Pareto, coef)
	names(res_red) <- c("Pareto", "c_Pareto", "coefficients")
	res_full <- list(Pareto, c_Pareto, summary(lm(log10(yi)~log10(xi))), cbind(xi, yi))
	names(res_full) <- c("Pareto", "c_Pareto", "regression_results", "coords_Pareto")

if (graph){
	plot(yi~xi, pch = 20, log = "xy", ylab = "Inverse cumulative frequencies", xlab = "Number of replicates", main = "Pareto distribution")
	abline(a = summary(lm(log10(yi)~log10(xi)))$coefficients[1], b = summary(lm(log10(yi)~log10(xi)))$coefficients[2])
		if (legends == 2){
			leg.txt <- c(paste("Pareto = ", round(Pareto, digits=3), sep = ""),
					paste("r2 = ", round(summary(lm(log10(yi)~log10(xi)))$adj.r.squared, digits = 3), sep = ""),
					paste("pvalue = ", round(summary(lm(log10(yi)~log10(xi)))$coefficients[[8]], digits = 3), sep = ""))
			legend("topright", leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white", cex = 0.8)
		} else {
			leg.txt <- c(paste("y=", round(summary(lm(log10(yi)~log10(xi)))$coefficients[2], digits = 3), "x+", round(summary(lm(log10(yi)~log10(xi)))$coefficients[1], digits = 3), sep = ""),
					paste("c_Pareto = ", round(c_Pareto, digits=3), sep = ""))
			legend("topright", leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white", cex = 0.8)
		}
}

if (full){
	res_full
} else {
	res_red
}
}
}


Pareto_index <- function(data1, vecpop = NULL, listMLL = NULL, full = FALSE, graph = FALSE, legends = 1, export = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					par(ask = TRUE)
					res[[p]] <- Pareto_index_core(datatot[[p]], listMLL[[p]], full, graph, legends)
					par(ask = FALSE)

					if (export){
						postscript(paste(paste("Pareto_index", unique(vecpop)[[p]], sep = "_"), ".eps", sep = ""), onefile = FALSE, paper = "letter")
						resu <- Pareto_index_core(datatot[[p]], listMLL[[p]], full = TRUE, graph, legends)

		xi <- resu$coords_Pareto[,1]
		yi <- resu$coords_Pareto[,2]
		Pareto <- resu$Pareto
		c_Pareto <- resu$c_Pareto
	plot(yi~xi, pch = 20, log = "xy", ylab = "Inverse cumulative frequencies", xlab = "Number of replicates", main = "Pareto distribution")
	abline(a = summary(lm(log10(yi)~log10(xi)))$coefficients[1], b = summary(lm(log10(yi)~log10(xi)))$coefficients[2])
		if (legends == 2){
			leg.txt <- c(paste("Pareto = ", round(Pareto, digits = 3), sep = ""),
					paste("r2 = ", round(summary(lm(log10(yi)~log10(xi)))$adj.r.squared, digits = 3), sep = ""),
					paste("pvalue = ", round(summary(lm(log10(yi)~log10(xi)))$coefficients[[8]], digits = 3), sep = ""))
			legend("topright", leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white", cex = 0.8)
		} else {
			leg.txt <- c(paste("y=", round(summary(lm(log10(yi)~log10(xi)))$coefficients[2], digits = 3), "x+", round(summary(lm(log10(yi)~log10(xi)))$coefficients[1], digits = 3), sep = ""),
					paste("c_Pareto = ", round(c_Pareto, digits = 3), sep = ""))
			legend("topright", leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white", cex = 0.8)
		}
						dev.off()
					}
				}
			names(res) <- unique(vecpop_o)
		} else {
			res <- Pareto_index_core(data1, listMLL, full, graph, legends)
				if (export){
					postscript("Pareto_index.eps", onefile = FALSE, paper = "letter")
					resu <- Pareto_index_core(data1, listMLL, full = TRUE, graph, legends)

		xi <- resu$coords_Pareto[,1]
		yi <- resu$coords_Pareto[,2]
		Pareto <- resu$Pareto
		c_Pareto <- resu$c_Pareto
	plot(yi~xi, pch = 20, log = "xy", ylab = "Inverse cumulative frequencies", xlab = "Number of replicates", main = "Pareto distribution")
	abline(a = summary(lm(log10(yi)~log10(xi)))$coefficients[1], b = summary(lm(log10(yi)~log10(xi)))$coefficients[2])
		if (legends == 2){
			leg.txt <- c(paste("Pareto = ", round(Pareto, digits = 3), sep = ""),
					paste("r2 = ", round(summary(lm(log10(yi)~log10(xi)))$adj.r.squared, digits = 3), sep = ""),
					paste("pvalue = ", round(summary(lm(log10(yi)~log10(xi)))$coefficients[[8]], digits = 3), sep = ""))
			legend("topright", leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white", cex = 0.8)
		} else {
			leg.txt <- c(paste("y=", round(summary(lm(log10(yi)~log10(xi)))$coefficients[2], digits = 3), "x+", round(summary(lm(log10(yi)~log10(xi)))$coefficients[1], digits = 3), sep = ""),
					paste("c_Pareto = ", round(c_Pareto, digits = 3), sep = ""))
			legend("topright", leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, bg = "white", cex = 0.8)
		}

					dev.off()
				}
		}
	res
}


#############################################################################
#######################Spatial composant of clonality########################
#############################################################################


######################
#kinship de Loiselle :
######################

kinship_Loiselle_core <- function(data1, haploid = FALSE){

	data2 <- freq_RR_core(data1, haploid, genet = FALSE, RR = FALSE)
	ncol_all <- 2
	ncol_freq <- 3
	N <- nrow(data1)

	if(haploid){
		index_l <- 1:ncol(data1)
	} else {
		index_l <- 1:c(ncol(data1)/2)*2-1
	}

	mat_auto <- matrix(0, ncol = N, nrow = N)
	row.names(mat_auto) <- colnames(mat_auto) <- 1:N

		divideur <- NULL
			for (k in index_l){
					if(haploid){
				freq_temp <- as.numeric(data2[data2 == paste("locus", k, sep = "_"), ncol_freq])
					} else {
				freq_temp <- as.numeric(data2[data2 == paste("locus", c((k+1)/2), sep = "_"), ncol_freq])
					}
				divideur <- c(divideur,sum(freq_temp*(1-freq_temp)))
			}
		div_f <- sum(divideur)
			

for (j in 1:N){
	for (i in j:N){
		sup_tiers <- NULL
		for (k in index_l){

		if (haploid){
			list_all_temp <- data2[data2 == paste("locus", k, sep = "_"), ncol_all]
			sup <- NULL
				nl <- nrow(data1)
			for (l in 1:length(list_all_temp)){
				pila <- length(which(data1[i, k]  == list_all_temp[l]))
				pjla <- length(which(data1[j, k] == list_all_temp[l]))
				pla <- as.numeric(data2[data2 == paste("locus", k, sep = "_"), ncol_freq][l])
				sup <- c(sup,(pila-pla)*(pjla-pla)+(pla*(1-pla))/(nl-1))
			}
		} else {
			list_all_temp <- data2[data2 == paste("locus",c((k+1)/2), sep = "_"), ncol_all]
			sup <- NULL
			nl <- 2*nrow(data1)
			for (l in 1:length(list_all_temp)){
				pila <- length(which(c(data1[i, k], data1[i,k+1]) == list_all_temp[l]))/2
				pjla <- length(which(c(data1[j, k], data1[j,k+1]) == list_all_temp[l]))/2
				pla <- as.numeric(data2[data2 == paste("locus", c((k+1)/2), sep = "_"), ncol_freq][l])
				sup <- c(sup,(pila-pla)*(pjla-pla)+(pla*(1-pla))/(nl-1))
			}	
		}
			sup_tiers <- c(sup_tiers, sum(sup))
		}
	mat_auto[i,j] <- sum(sup_tiers)/div_f
	}
}
as.dist(mat_auto)
}


kinship_Loiselle <- function(data1, haploid = FALSE, vecpop = NULL){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)))

		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					res[[p]] <- kinship_Loiselle_core(datatot[[p]], haploid)
					names(res[p]) <- unique(vecpop_o)[p]
				}
		} else {
			res <- kinship_Loiselle_core(data1, haploid)

		}
	res
}


###################
#kinship de Ritland
###################

kinship_Ritland_core <- function(data1, haploid = FALSE){

	ncol_all <- 2
	ncol_freq <- 3
	N <- nrow(data1)
	data2 <- freq_RR(data1, haploid, genet = FALSE, RR = FALSE)

	if (haploid){
		index_l <- 1:ncol(data1)
	} else {
		index_l <- 1:c(ncol(data1)/2)*2-1
	}

	mat_Rit2 <- matrix(0, ncol = N, nrow = N)
	row.names(mat_Rit2) <- colnames(mat_Rit2) <- 1:N

	divideur <- NULL
		for (k in index_l){
			if (haploid){
				list_all_temp <- data2[data2 == paste("locus", k, sep = "_"), ncol_all]				
			} else {
				list_all_temp <- data2[data2 == paste("locus", c((k+1)/2), sep = "_"), ncol_all]
			}
		divideur <- c(divideur, length(list_all_temp)-1)
		}
	div_f <- sum(divideur)
			

for (j in 1:N){
	for (i in j:N){
		sup_tiers <- NULL
		for (k in index_l){
			if (haploid){
				list_all_temp <- data2[data2 == paste("locus", k, sep = "_"), ncol_all]
				sup <- NULL
					for (l in 1:length(list_all_temp)){
						pila <- as.numeric(length(which(data1[i, k] == list_all_temp[l])))
						pjla <- as.numeric(length(which(data1[j, k] == list_all_temp[l])))
						pla <- as.numeric(data2[data2 == paste("locus", k, sep = "_"), ncol_freq][l])
						sup <- c(sup,(pila*pjla)/pla)
					}
			} else {
				list_all_temp <- data2[data2 == paste("locus", c((k+1)/2), sep = "_"), ncol_all]
				sup <- NULL
					for (l in 1:length(list_all_temp)){
						pila <- as.numeric(length(which(c(data1[i, k], data1[i,k+1]) == list_all_temp[l]))/2)
						pjla <- as.numeric(length(which(c(data1[j, k], data1[j,k+1]) == list_all_temp[l]))/2)
						pla <- as.numeric(data2[data2==paste("locus",c((k+1)/2), sep = "_"), ncol_freq][l])
						sup <- c(sup,(pila*pjla)/pla)
					}
			}	
			sup_tiers <- c(sup_tiers, sum(sup)-1)
		}
	mat_Rit2[i,j] <- sum(sup_tiers)/div_f
	}
}
as.dist(mat_Rit2)
}


kinship_Ritland <- function(data1, haploid = FALSE, vecpop = NULL){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}
		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)))

		datatot <- split(data1, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					res[[p]] <- kinship_Ritland_core(datatot[[p]], haploid)
					names(res[p]) <- unique(vecpop_o)[p]
				}
		} else {
			res <- kinship_Ritland_core(data1, haploid)

		}
	res
}


red_MLL <- function(data1, listMLL){
	N <- nrow(data1)
	mat <- as.data.frame(matrix(0, ncol = N, nrow = 1))
	colnames(mat) <- 1:N
	MLG <- MLG_list(data1)
		for (i in 1:length(MLG)){
			if (length(MLG[[i]]) > 1){
				interm <- MLG[[i]]
					for(j in 1:length(interm)){
						mat[1, which(interm[j] == colnames(mat))] <- length(interm)
					}
			} else {
				mat[1, which(MLG[[i]] == colnames(mat))] <- 1
			}
		}
	poids <- mat

	MLL_red <- listMLL
		for (k in 1:length(listMLL)){
			if (length(listMLL[[k]]) > 1){
				interm <- NULL
				for (l in 1:length(listMLL[[k]])){
					interm <- c(interm, poids[1, which(listMLL[[k]][l] == colnames(poids))])
				}
				interm <- interm/length(interm)
				MLL_red[[k]] <- listMLL[[k]][which(interm == max(interm))][1]
			}
		}

if(length(listMLL) != length(MLL_red)){stop("Problem with MLL weigth generation list")}

	sort(unlist(MLL_red))
}


red_MLL_mod <- function(data1, listMLL){
	N <- nrow(data1)
	mat <- as.data.frame(matrix(0, ncol = N, nrow = 1))
	colnames(mat) <- 1:N
		MLG <- MLG_list(data1)
			for (i in 1:length(MLG)){
				if (length(MLG[[i]]) > 1){
					interm <- MLG[[i]]
						for(j in 1:length(interm)){
							mat[1, which(interm[j] == colnames(mat))] <- length(interm)
						}
				} else {
					mat[1, which(MLG[[i]] == colnames(mat))] <- 1
				}
			}
		poids <- mat

		MLL_red <- listMLL
			for (k in 1:length(listMLL)){
				if (length(listMLL[[k]]) > 1){
					interm <- NULL
					for (l in 1:length(listMLL[[k]])){
						interm <- c(interm, poids[1, which(listMLL[[k]][l] == colnames(poids))])
					}
					interm <- interm/length(interm)
					MLL_red[[k]] <- sample(listMLL[[k]][which(interm == max(interm))], 1)
				}
			}

if(length(listMLL) != length(MLL_red)){stop("Problem with MLL weigth generation list")}

	sort(unlist(MLL_red))
}

####################
#Fonction principale
####################

autocorrelation_core <- function(data1, haploid = FALSE, coords = NULL, listMLL = NULL, 
	Loiselle = FALSE, Ritland = FALSE,
	genet = FALSE, central_coords = FALSE, random_unit = FALSE, weighted = FALSE,
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL,
	graph = FALSE, nbrepeat = NULL){

if (length(listMLL) != 0 & genet){
	if (!central_coords & !random_unit & !weighted){stop("You must choose a MLG genets methods for MLL.")}
}

        ncol_all <- 2
        ncol_freq <- 3
        N <- nrow(data1)
        G <- nrow(unique(data1))
        

	if (length(listMLL) != 0){
                list_genet <- listMLL
			dataMLG <- data1[red_MLL(data1, list_genet),]
	} else {
                list_genet <- MLG_list(data1)
			dataMLG <- unique(data1)
	}

if (length(unlist(list_genet)) != nrow(data1)){stop("MLL list does not compute")}

if (length(listMLL) != 0 & central_coords | length(listMLL) != 0 & random_unit | genet & central_coords | genet & random_unit){
        dataused <- dataMLG
} else {
        dataused <- data1
}

if (Loiselle){
        mat_auto <- as.matrix(kinship_Loiselle_core(dataused, haploid))
} else if (Ritland){
        mat_auto <- as.matrix(kinship_Ritland_core(dataused, haploid))
} else {
stop("You must choose between Loiselle and Ritland kinship methods.")
}

if (genet){
        if (central_coords){
        recup_clo <- NULL
        coord_MLG1 <- cbind(1:N, coords)
        index_clo <- which(lapply(list_genet, length) > 1)
                if (length(index_clo) != 0){
                	for (i in index_clo){
                        recup <- list_genet[[i]]
                        recup_clo <- c(recup_clo, recup)
                        coord_MLG1 <- rbind(coord_MLG1, c(recup[1], mean(coords[recup,1]), 
					mean(coords[recup,2])))
                	}
        coord_MLG1 <- coord_MLG1[-recup_clo,]
                }
        coord_MLG1 <- coord_MLG1[order(coord_MLG1[,1]),]
        rownames(coord_MLG1) <- coord_MLG1[,1]
        coord_MLG1 <- coord_MLG1[,-1]

        dist_geo <- dist(coord_MLG1, method = "euclidean", upper = FALSE, diag = FALSE)
        mat_dist_geo <- as.matrix(dist_geo)
        max_dist <- max(mat_dist_geo)

        } else if (random_unit){
		if (length(listMLL) != 0){
			coord_MLG2 <- coords[red_MLL_mod(data1, list_genet),]
		} else {
        recup_clo <- NULL
        coord_MLG2 <- cbind(1:N, coords)
        index_clo <- which(lapply(list_genet, length) > 1)
                if (length(index_clo) != 0){
                	for (i in index_clo){
                        recup <- list_genet[[i]]
                        recup_clo <- c(recup_clo, recup)
                        shuffle <- sample(recup, 1)
                        coord_MLG2 <- rbind(coord_MLG2, c(shuffle, coords[shuffle,1], 
					coords[shuffle,2]))
                	}
        coord_MLG2 <- coord_MLG2[-recup_clo,]
                }
        coord_MLG2 <- coord_MLG2[order(coord_MLG2[,1]),]
        rownames(coord_MLG2) <- coord_MLG2[,1]
        coord_MLG2 <- coord_MLG2[,-1]
		}

        dist_geo <- dist(coord_MLG2, method = "euclidean", upper = FALSE, diag = FALSE)
        mat_dist_geo <- as.matrix(dist_geo)
        max_dist <- max(mat_dist_geo)

        } else if (weighted){

        mat_pond <- matrix(1, ncol = N, nrow = N)
        index_clo <- which(lapply(list_genet, length) > 1)
                for (i in index_clo){
                        develop <- list_genet[[i]]
                                for (j in develop){
                                        mat_pond[j,] <- mat_pond[j,]*1/length(develop)
                                        mat_pond[,j] <- mat_pond[,j]*1/length(develop)
                                                for (k in develop){
                                                        mat_pond[j,k] <- mat_pond[j,k]*0
                                                }
                                }       
                }
        mat_auto <- mat_auto*mat_pond

        dist_geo <- dist(coords, method = "euclidean", upper = FALSE, diag = FALSE)
        mat_dist_geo <- as.matrix(dist_geo)
        max_dist <- max(mat_dist_geo)
        }
} else {
        dist_geo <- dist(coords, method = "euclidean", upper = FALSE, diag = FALSE)
        mat_dist_geo <- as.matrix(dist_geo)
        max_dist <- max(mat_dist_geo)
}
        
if (class1){
        d <- d
        class_d <- (max_dist/d)*0:d
} else if (class2){
        d <- d
        dif <-  length(dist_geo) - round(length(dist_geo)/d, digits = 0)*d
                if (dif >= 0){
                        class_du <- rep(round(length(dist_geo)/d, digits=0), d) + c(rep(1, dif), rep(0, d-dif))
                        if (sum(class_du) != length(dist_geo)){stop("warning ! Problem with distance class")}
                }
                if (dif < 0){
                        class_du <- rep(round(length(dist_geo)/d, digits=0), d) + c(rep(0, d-abs(dif)), rep(-1, abs(dif)))
                        if (sum(class_du) != length(dist_geo)){stop("warning ! Problem with distance class")}
                }

        L <- ncol(mat_dist_geo)

        tab <- as.data.frame(matrix(0, ncol = 3, nrow = 1))
                for (j in 1:(L-1)){
                        sub_tab <- cbind((j+1):L, rep(j, times = length((j+1):L)), mat_dist_geo[(j+1):L, j])
                        tab <- rbind(tab, sub_tab)
                }
        tab <- tab[-1,]
        names(tab) <- c("unit_1", "unit_2", "dist")
        nrow(tab) == length(dist_geo)
        tab_sort <- tab[order(tab[,3]),]

        class_du <- cumsum(class_du)
        class_du <- c(0, class_du)

} else if (length(vecdist) != 0){
        class_d <- vecdist
} else {
        d <- 10
        class_d <- (max_dist/d)*0:d
} 

if (class2){

	recup_kin_all <- list(NULL)
	recup_dist_all <- list(NULL)
        for (n in 1:(length(class_du)-1)){
                couples_class <- tab_sort[(class_du[n]+1):(class_du[n+1]),]
                recup_kin <- NULL
                        for (i in 1:nrow(couples_class)){
                                recup_kin <- c(recup_kin, mat_auto[couples_class[i,1], couples_class[i,2]])
                        }
                recup_kin_all[[n]] <- recup_kin
                recup_dist_all[[n]] <- couples_class[,3]
        }
	class_d <- c(0, unlist(lapply(recup_dist_all, max)))

} else {

	L <- ncol(mat_dist_geo)
	mat_dist_geo <- round(mat_dist_geo, digits = 5)
	class_d <- round(class_d, digits = 5)
	recup_kin_all <- list(NULL)
	recup_dist_all <- list(NULL)
        for (n in 1:(length(class_d)-1)){
                recup_kin_col <- NULL
                recup_dist_col <- NULL
                        for (j in 1:(L-1)){
                                if (n == length(class_d)-1){
                                        couples_class <- which(mat_dist_geo[(j+1):L, j] <= class_d[n+1] & mat_dist_geo[(j+1):L, j] >= class_d[n])
                                } else {
                                        couples_class <- which(mat_dist_geo[(j+1):L, j] < class_d[n+1] & mat_dist_geo[(j+1):L, j] >= class_d[n])
                                }
                                couples_class <- as.numeric(couples_class)+j
                                        if (length(couples_class) != 0){
                                                recup_kin <- NULL
                                                recup_dist <- NULL
                                                        for (i in 1:length(couples_class)){
                                                                        recup_kin <- c(recup_kin, mat_auto[couples_class[i], j])
                                                                        recup_dist <- c(recup_dist, mat_dist_geo[couples_class[i], j])
                                                        }
                                                recup_kin_col <- c(recup_kin_col, recup_kin)
                                                recup_dist_col <- c(recup_dist_col, recup_dist)
                                        }
                        }
                recup_kin_all[[n]] <- recup_kin_col
                recup_dist_all[[n]] <- recup_dist_col
        }
}

if (sum(unlist(lapply(recup_kin_all, length))) != length(as.dist(mat_dist_geo)) | sum(unlist(lapply(recup_kin_all, length))) != length(as.dist(mat_dist_geo))){stop("Problem with class distances or geographic matrix detected")}

if (length(nbrepeat) != 0){

        if (length(listMLL) != 0 & random_unit | genet & random_unit){
                kin_sim <- matrix(NA, ncol = length(class_d)-1, nrow = 1)
                        for (s in 1:nbrepeat){
					if(length(listMLL) != 0){
						coord_MLG2 <- coords[red_MLL_mod(data1, list_genet),]
					} else {
                                recup_clo <- NULL
                                coord_MLG2 <- cbind(1:N, coords)
                                index_clo <- which(lapply(list_genet, length) > 1)
                                        if (length(index_clo) != 0){
                                        for (i in index_clo){
                                                recup <- list_genet[[i]]
                                                recup_clo <- c(recup_clo, recup)
                                                shuffle <- sample(recup, 1)
                                                coord_MLG2 <- rbind(coord_MLG2, c(shuffle, coords[shuffle,1], coords[shuffle,2]))
                                        }
                                coord_MLG2 <- coord_MLG2[-recup_clo,]
                                        }
                                coord_MLG2 <- coord_MLG2[order(coord_MLG2[,1]),]
                                rownames(coord_MLG2) <- coord_MLG2[,1]
                                coord_MLG2 <- coord_MLG2[,-1]
					}
        
                                dist_geo_sim <- dist(coord_MLG2, method = "euclidean", upper = FALSE, diag = FALSE)
                                mat_dist_geo_sim <- as.matrix(dist_geo_sim)
                                max_dist_sim <- max(mat_dist_geo_sim)

                                recup_coord_all_x <- list(NULL)
                                recup_coord_all_y <- list(NULL)
                                        for (n in 1:(length(class_d)-1)){
                                                recup_coord_x2 <- NULL
                                                recup_coord_y2 <- NULL
                                                        for (j in 1:(L-1)){
                                                                if (n == length(class_d)-1){
                                                                        couples_class <- which(mat_dist_geo_sim[(j+1):L, j] <= class_d[n+1] & mat_dist_geo_sim[(j+1):L, j] >= class_d[n])
                                                                } else {
                                                                couples_class <- which(mat_dist_geo_sim[(j+1):L, j] < class_d[n+1] & mat_dist_geo_sim[(j+1):L, j] >= class_d[n])
                                                                }
                                                                couples_class <- as.numeric(couples_class)+j    
                                                                        if (length(couples_class) != 0){
                                                                                recup_coord_x <- NULL
                                                                                recup_coord_y <- NULL
                                                                                        for (i in 1:length(couples_class)){
                                                                                                recup_coord_x <- c(recup_coord_x, couples_class[i])
                                                                                                recup_coord_y <- c(recup_coord_y, j)
                                                                                        }
                                                                                recup_coord_x2 <- c(recup_coord_x2, recup_coord_x)
                                                                                recup_coord_y2 <- c(recup_coord_y2, recup_coord_y)
                                                                        }
                                                        }
                                                recup_coord_all_x[[n]] <- recup_coord_x2
                                                recup_coord_all_y[[n]] <- recup_coord_y2
                                        }

                                order_sim <- sample(L, L, replace = FALSE)
                                tab_kin_sim <- mat_auto[order_sim, order_sim]
                                recup_kin_sim <- list(NULL)
                                        for (n in 1:length(recup_coord_all_x)){
                                                recup_kin <- NULL
                                                coords_kin_x <- recup_coord_all_x[[n]]
                                                coords_kin_y <- recup_coord_all_y[[n]]
                                                        if (length(coords_kin_x) == length(coords_kin_y)){
                                                                for (i in 1:length(coords_kin_x)){
                                                                        recup_kin <- c(recup_kin, tab_kin_sim[coords_kin_x[i], coords_kin_y[i]])
                                                                }
                                                        }
                                                recup_kin_sim[[n]] <- recup_kin
                                        }
                                kin_sim <- rbind(kin_sim, unlist(lapply(recup_kin_sim, mean)))
                        }

                kin_sim <- kin_sim[-1,]
                colnames(kin_sim) <- paste("class", 1:(length(class_d)-1), sep = "_")
        
        } else {

        recup_coord_all_x <- list(NULL)
        recup_coord_all_y <- list(NULL)
                for(n in 1:(length(class_d)-1)){
                        recup_coord_x2 <- NULL
                        recup_coord_y2 <- NULL
                                for(j in 1:(L-1)){
                                        if(n == length(class_d)-1){
                                                couples_class <- which(mat_dist_geo[(j+1):L, j] <= class_d[n+1] & mat_dist_geo[(j+1):L, j] >= class_d[n])
                                        } else {
                                                couples_class <- which(mat_dist_geo[(j+1):L, j] < class_d[n+1] & mat_dist_geo[(j+1):L, j] >= class_d[n])
                                        }
                                        couples_class <- as.numeric(couples_class)+j
                                                if (length(couples_class) != 0){
                                                        recup_coord_x <- NULL
                                                        recup_coord_y <- NULL
                                                                for(i in 1:length(couples_class)){
                                                                        recup_coord_x <- c(recup_coord_x, couples_class[i])
                                                                        recup_coord_y <- c(recup_coord_y, j)
                                                                }
                                                        recup_coord_x2 <- c(recup_coord_x2, recup_coord_x)
                                                        recup_coord_y2 <- c(recup_coord_y2, recup_coord_y)
                                                }
                                }
                        recup_coord_all_x[[n]] <- recup_coord_x2
                        recup_coord_all_y[[n]] <- recup_coord_y2
                }

        kin_sim <- matrix(NA, ncol = length(recup_coord_all_x), nrow = 1)
                for (s in 1:nbrepeat){
                        order_sim <- sample(L, L, replace = FALSE)
                        tab_kin_sim <- mat_auto[order_sim, order_sim]

                        recup_kin_sim <- list(NULL)
                                for (n in 1:length(recup_coord_all_x)){
                                        recup_kin <- NULL
                                        coords_kin_x <- recup_coord_all_x[[n]]
                                        coords_kin_y <- recup_coord_all_y[[n]]
                                                if (length(coords_kin_x) == length(coords_kin_y)){
                                                        for (i in 1:length(coords_kin_x)){
                                                                recup_kin <- c(recup_kin, tab_kin_sim[coords_kin_x[i], coords_kin_y[i]])
                                                        }
                                                }
                                        recup_kin_sim[[n]] <- recup_kin
                                }
                        kin_sim <- rbind(kin_sim, unlist(lapply(recup_kin_sim, mean)))
                }

        kin_sim <- kin_sim[-1,]
        colnames(kin_sim) <- paste("class", 1:length(recup_coord_all_x), sep = "_")
}

	pval_kin <- NULL
        for (o in 1:length(recup_coord_all_x)){
			pval_kin_sim <- 2*min(mean(kin_sim[,o] >= mean(recup_kin_all[[o]])), mean(kin_sim[,o] <= mean(recup_kin_all[[o]])))
				if (pval_kin_sim > 1) {pval_kin_sim <- 1}
			pval_kin <- c(pval_kin, pval_kin_sim)
        }

	tabsim <- as.data.frame(matrix(NA, ncol = 4, nrow = nbrepeat))
	names(tabsim) <- c("b", "b_log", "Sp", "Sp_log")
	dm <- unlist(lapply(recup_dist_all, mean))
        for (s in 1:nbrepeat){
                sim <- kin_sim[s,]
                tabsim[s, 1] <- lm(sim~dm)$coefficients[[2]]
                tabsim[s, 2] <- lm(sim~log(dm))$coefficients[[2]]
			tabsim[s, 3] <- -tabsim[s, 1]/(1-sim[1])
			tabsim[s, 4] <- -tabsim[s, 2]/(1-sim[1])
        }
}

tabres <- cbind(unlist(lapply(recup_dist_all, min)), 
        unlist(lapply(recup_dist_all, max)), unlist(lapply(recup_dist_all, mean)), 
        log(unlist(lapply(recup_dist_all, mean))), unlist(lapply(recup_dist_all, length)),
        unlist(lapply(recup_kin_all, mean)))
tabres <- as.data.frame(tabres)

if (Loiselle){
        names(tabres) <- c("dist_min", "dist_max", "dist_mean", "ln(dist_mean)", "nb_pairs", "mean_Loiselle")
} else {
        names(tabres) <- c("dist_min", "dist_max", "dist_mean", "ln(dist_mean)", "nb_pairs", "mean_Ritland")
}

	obs <- unlist(lapply(recup_kin_all, mean))
	dm <- unlist(lapply(recup_dist_all, mean))
	b <- lm(obs~dm)$coefficients[[2]]
	b_log <- lm(obs~log(dm))$coefficients[[2]]
	F1 <- obs[1]
	Sp <- -b/(1-F1)
	Sp_log <- -b_log/(1-F1)

	tabres2 <- as.data.frame(cbind(b, b_log, Sp, Sp_log))

if (length(nbrepeat) != 0){
	tabres <- cbind(tabres, pval_kin)

	tabres3 <- as.data.frame(cbind(b, b_log, Sp, Sp_log))
	tabres3 <- rbind(tabres3, apply(tabsim, 2, mean), apply(tabsim, 2, sd),
	c(sort(tabsim[,1])[nbrepeat*0.025], sort(tabsim[,2])[nbrepeat*0.025], sort(tabsim[,3])[nbrepeat*0.025], sort(tabsim[,4])[nbrepeat*0.025]),
	c(sort(tabsim[,1])[nbrepeat*0.975], sort(tabsim[,2])[nbrepeat*0.975], sort(tabsim[,3])[nbrepeat*0.975], sort(tabsim[,4])[nbrepeat*0.975]),
	c(sort(tabsim[,1])[nbrepeat*0.05], sort(tabsim[,2])[nbrepeat*0.05], sort(tabsim[,3])[nbrepeat*0.05], sort(tabsim[,4])[nbrepeat*0.05]),
	c(sort(tabsim[,1])[nbrepeat*0.95], sort(tabsim[,2])[nbrepeat*0.95], sort(tabsim[,3])[nbrepeat*0.95], sort(tabsim[,4])[nbrepeat*0.95]),
	c(mean(tabsim[,1] <= b), mean(tabsim[,2] <= b_log), mean(tabsim[,3] <= Sp), mean(tabsim[,4] <= Sp_log)),
	c(mean(tabsim[,1] >= b), mean(tabsim[,2] >= b_log), mean(tabsim[,3] >= Sp), mean(tabsim[,4] >= Sp_log)),
	c(if (2*min(mean(tabsim[,1] <= b), mean(tabsim[,1] >= b)) > 1) {1} else {2*min(mean(tabsim[,1] <= b), mean(tabsim[,1] >= b))},
	if (2*min(mean(tabsim[,2] <= b_log), mean(tabsim[,2] >= b_log)) > 1) {1} else {2*min(mean(tabsim[,2] <= b_log), mean(tabsim[,2] >= b_log))},
	if (2*min(mean(tabsim[,3] <= Sp), mean(tabsim[,3] >= Sp)) > 1) {1} else {2*min(mean(tabsim[,3] <= Sp), mean(tabsim[,3] >= Sp))},
	if (2*min(mean(tabsim[,4] <= Sp_log), mean(tabsim[,4] >= Sp_log)) > 1) {1} else {2*min(mean(tabsim[,4] <= Sp_log), mean(tabsim[,4] >= Sp_log))}))

	rownames(tabres3) <- c("obs_value", "mean_sim", "sd_sim", "0.95_inf", "0.95_sup", "0.9_inf", "0.9_sup", "pval_upper", "pval_lower", "pval_2sides")
}

if (graph){
	plot(unlist(lapply(recup_dist_all, mean)), unlist(lapply(recup_kin_all, mean)), type = "b", pch = 20, 
        ylim = c(min(unlist(lapply(recup_kin_all, min))), max(unlist(lapply(recup_kin_all, max)))),
        main = "Spatial autocorrelation analysis", xlab = "Spatial distance", ylab = "Coancestry (Fij)")
	abline(h = 0, lty = 3)
}

if (length(nbrepeat) != 0){
        list(Main_results = tabres, Slope_and_Sp_index = tabres3, Slope_resample = tabsim, Kinship_resample = kin_sim, Matrix_kinship_results = as.dist(mat_auto), Class_kinship_results = recup_kin_all, Class_distance_results = recup_dist_all)
} else {
        list(Main_results = tabres, Slope_and_Sp_index = tabres2, Matrix_kinship_results = as.dist(mat_auto), Class_kinship_results = recup_kin_all, Class_distance_results = recup_dist_all)
}
}


graph_autocorrelation <- function(Class_distance_results, Class_kinship_results){

	recup_dist_all <- Class_distance_results
	recup_kin_all <- Class_kinship_results

	plot(unlist(lapply(recup_dist_all, mean)), unlist(lapply(recup_kin_all, mean)), type = "b", pch = 20, 
        	ylim = c(min(unlist(lapply(recup_kin_all, min))), max(unlist(lapply(recup_kin_all, max)))),
       	main = "Spatial autocorrelation analysis", xlab = "Spatial distance", ylab = "Coancestry (Fij)")
	abline(h = 0, lty = 3)
}


autocorrelation <- function(data1, haploid = FALSE, coords = NULL, vecpop = NULL, listMLL = NULL, 
	Loiselle = FALSE, Ritland = FALSE,
	genet = FALSE, central_coords = FALSE, random_unit = FALSE, weighted = FALSE,
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL, 
	graph = FALSE, nbrepeat = NULL, export = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		coordtot <- split(coords, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					par(ask = TRUE)
					res[[p]] <- autocorrelation_core(datatot[[p]], haploid, coordtot[[p]], listMLL[[p]], Loiselle, 
								Ritland, genet, central_coords, random_unit, weighted,
								class1, class2, d, vecdist, graph, nbrepeat)
					par(ask = FALSE)

					if (export){
						postscript(paste(paste("autocorrelation", unique(vecpop)[[p]], sep = "_"), ".eps", sep = ""), onefile = FALSE, paper = "letter")
							graph_autocorrelation(res[[p]]$Class_distance_results, res[[p]]$Class_kinship_results)
						dev.off()
					}
				}
					names(res) <- unique(vecpop_o)
		} else {
			res <- autocorrelation_core(data1, haploid, coords, listMLL, Loiselle, 
				Ritland, genet, central_coords, random_unit, weighted,
				class1, class2, d, vecdist, graph, nbrepeat)
				if (export){
					postscript("autocorrelation.eps", onefile = FALSE, paper = "letter")
						graph_autocorrelation(res$Class_distance_results, res$Class_kinship_results)				
					dev.off()
				}
		}
	res
}


#############
#Subfunctions
#############

id_clonal <- function(data1){

	N <- nrow(data1)
	nb_all <- ncol(data1)
	res_id_clonal <- matrix(NA, ncol = N, nrow = N)
		for (j in 1:N){
			for (i in j:N){
				if (sum(data1[i,] == data1[j,]) == nb_all){
				res_id_clonal[i,j] <- 0
				} else res_id_clonal[i,j] <- 1
			}
		}
	res_id_clonal
}


id_clonal_mod <- function(data1){

	N <- nrow(data1)
	nb_all <- ncol(data1)
	res_id_clonal_mod <- matrix(NA, ncol = N, nrow = N)
		for (j in 1:N){
			for (i in j:N){
				if (sum(data1[i,] == data1[j,]) == nb_all){
				res_id_clonal_mod[i,j] <- 1
				} else res_id_clonal_mod[i,j] <- 0
			}
		}
	res_id_clonal_mod
}


id_clonal_MLL <- function(data1, listMLL = NULL){

        N <- nrow(data1)
        list_genet <- listMLL
if (length(unlist(list_genet)) != nrow(data1)){stop("MLL list does not compute")}
        indic <- which(lapply(list_genet, length) > 1)
        res_id_clonal <- matrix(1, ncol = N, nrow = N)
                for (k in indic){
                        liste <- as.numeric(list_genet[[k]])
                                for (i in liste){
                                        for (j in liste){
                                                res_id_clonal[i,j] <- 0
                                        }
                                }
                }
rownames(res_id_clonal) <- colnames(res_id_clonal) <- 1:N
res_id_clonal <- res_id_clonal[as.numeric(rownames(data1)), as.numeric(rownames(data1))]
res_id_clonal
}


id_clonal_mod_MLL <- function(data1, listMLL = NULL){

	N <- nrow(data1)
	list_genet <- listMLL
if (length(unlist(list_genet)) != nrow(data1)){stop("MLL list does not compute")}
	indic <- which(lapply(list_genet, length) > 1)
	res_id_clonal_mod <- matrix(0, ncol = N, nrow = N)
		for (k in indic){
			liste <- as.numeric(list_genet[[k]])
				for (i in liste){
					for (j in liste){
						res_id_clonal_mod[i,j] <- 1
					}
				}
		}
res_id_clonal_mod
}


nearest_unit <- function(coords){

	nb1 <- nrow(coords)
	dist_geo <- dist(coords, method = "euclidean", upper = FALSE, diag = FALSE)
	obj1 <- as.matrix(dist_geo)
	coord_min <- remp <- as.data.frame(matrix(0, ncol = 4, nrow = 1))
	listing <- 1:nb1
		for (i in 1:nb1){
			index <- rownames(obj1)[which(obj1[i,] == min(obj1[i,listing[-i]]))]
				for (j in 1:length(index)){
					remp[,1] <- 1/length(index)
					remp[,2] <- sort(as.numeric(c(index[j],i)))[1]
					remp[,3] <- sort(as.numeric(c(index[j],i)))[2]
					remp[,4] <- min(obj1[i,listing[-i]])
				}
			coord_min <- rbind(coord_min, remp)
		}

	coord_min <- coord_min[-1,]
	names(coord_min) <- c("Freq", "Unit", "Unit", "Dist_min")
	row.names(coord_min) <- 1:nrow(coord_min)

		for (i in 1:nrow(coord_min)){
			index <- which(coord_min[,2] == coord_min[i,2] & coord_min[,3] == coord_min[i,3])
				if (sum(coord_min[index,1] != 1) != 0){
				coord_min[index,1] <- mean(coord_min[index,1])
				}
			}
	coord_min_red <- coord_min[rownames(unique(coord_min[,2:3])),]
	coord_min_red <- apply(coord_min_red, 2, as.numeric)
	coord_min_red
}


################
#Clonal Subrange
################

clonal_sub_core <- function(data1, coords = NULL, listMLL = NULL, 
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL){
	
	N <- nrow(data1)
	mat_dist_geo <- as.matrix(dist(coords, method = "euclidean", upper = FALSE, diag = FALSE))
	dist_geo <- dist(coords, method = "euclidean", upper = FALSE, diag = FALSE)
	max_dist <- max(as.dist(mat_dist_geo))

if (length(vecdist) != 0){
	class_d <- vecdist
} else if (class1){
	d <- d
	class_d <- (max_dist/d)*0:d
} else if (class2){
	d <- d
	dif <- length(dist_geo) - round(length(dist_geo)/d, digits = 0)*d
		if (dif >= 0){
			class_du <- rep(round(length(dist_geo)/d, digits = 0), d) + c(rep(1, dif), rep(0, d-dif))
			if (sum(class_du) != length(dist_geo)){stop("warning !")}
		}
		if (dif < 0){
			class_du <- rep(round(length(dist_geo)/d, digits = 0), d) + c(rep(0, d-abs(dif)), rep(-1, abs(dif)))
			if (sum(class_du) != length(dist_geo)){stop("warning !")}
		}

	L <- ncol(mat_dist_geo)

	tab <- as.data.frame(matrix(0, ncol = 3, nrow = 1))
		for (j in 1:(L-1)){
			sub_tab <- cbind((j+1):L, rep(j, times = length((j+1):L)), mat_dist_geo[(j+1):L, j])
			tab <- rbind(tab, sub_tab)
		}
	tab <- tab[-1,]
	names(tab) <- c("unit_1", "unit_2", "dist")
	if (nrow(tab) != length(dist_geo)){stop("Problem with geographic matrix")}
	tab_sort <- tab[order(tab[,3]),]

	class_du <- cumsum(class_du)
	class_du <- c(0, class_du)
} else {
	d <- 10
	class_d <- (max_dist/d)*0:d
}

if (length(listMLL) != 0){
	res_id_clonal_mod <- id_clonal_mod_MLL(data1, listMLL)
} else {
	res_id_clonal_mod <- id_clonal_mod(data1)
}

if (sum(as.dist(res_id_clonal_mod)) != 0){

mat_dist_geo <- round(mat_dist_geo, digits = 7)

if (class2){

recup_sub_all <- list(NULL)
recup_sub_all2 <- list(NULL)
	for (n in 1:(length(class_du)-1)){
		recup_sub_col <- NULL
		couples_class <- tab_sort[(class_du[n]+1):(class_du[n+1]),]
			for (i in 1:nrow(couples_class)){
				recup_sub_col <- c(recup_sub_col, (res_id_clonal_mod*mat_dist_geo)[couples_class[i,1], couples_class[i,2]])
			}
		recup_sub_all[[n]] <- recup_sub_col
		recup_sub_all2[[n]] <- couples_class[,3]
	}
class_d <- round(c(0, unlist(lapply(recup_sub_all2, max))), digits = 7)

} else {
	
class_d <- round(class_d, digits = 7)

recup_sub_all <- list(NULL)
recup_sub_all2 <- list(NULL)
	for (n in 1:(length(class_d)-1)){
		recup_sub_col2 <- NULL
		recup_sub_col3 <- NULL
			for (j in 1:(N-1)){
				if (n == length(class_d)-1){
					couples_class <- which(mat_dist_geo[(j+1):N, j] <= class_d[n+1] & mat_dist_geo[(j+1):N, j] >= class_d[n])
				} else {
					couples_class <- which(mat_dist_geo[(j+1):N, j] < class_d[n+1] & mat_dist_geo[(j+1):N, j] >= class_d[n])
				}
				couples_class <- as.numeric(couples_class)+j
					if (length(couples_class) != 0){
						recup_sub <- NULL
						recup_sub2 <- NULL
						for (i in 1:length(couples_class)){
							recup_sub <- c(recup_sub, (res_id_clonal_mod*mat_dist_geo)[couples_class[i], j])
							recup_sub2 <- c(recup_sub2, mat_dist_geo[couples_class[i], j])
						}
					recup_sub_col2 <- c(recup_sub_col2, recup_sub)
					recup_sub_col3 <- c(recup_sub_col3, recup_sub2)
				}
		}
	recup_sub_all[[n]] <- recup_sub_col2
	recup_sub_all2[[n]] <- recup_sub_col3
}
}

if (sum(unlist(lapply(recup_sub_all2, length))) != length(as.dist(mat_dist_geo)) | sum(unlist(lapply(recup_sub_all, length))) != length(as.dist(mat_dist_geo))){stop("Problem with class distances or geographic matrix detected")}

recup_sub_all_1 <- lapply(recup_sub_all, function(x) x[x > 0])
clonal_sub_res <- max(unlist(recup_sub_all_1))

tabres <- cbind(lapply(recup_sub_all2, length), lapply(recup_sub_all2, min), 
	lapply(recup_sub_all2, max), lapply(recup_sub_all2, mean),
	unlist(lapply(recup_sub_all_1, length))/unlist(lapply(recup_sub_all, length)), 
	log10(unlist(lapply(recup_sub_all_1, length))/unlist(lapply(recup_sub_all, length)))
)

colnames(tabres) <- c("nb_pairs", "dist_min", "dist_max", "dist_mean", "Fr", "log(Fr)")
rownames(tabres) <- 1:nrow(tabres)

list(clonal_sub_res, clonal_sub_tab = tabres)
}
}


clonal_sub <- function(data1, coords = NULL, vecpop = NULL, listMLL = NULL, class1 = FALSE, class2 = FALSE, 
	d = NULL, vecdist = NULL){
	
	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		coordtot <- split(coords, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					rownames(coordtot[[p]]) <- 1:nrow(coordtot[[p]])
					res[[p]] <- clonal_sub_core(datatot[[p]], coordtot[[p]], listMLL[[p]], class1, class2, d, vecdist)
				}
					names(res) <- unique(vecpop_o)
		} else {
			res <- clonal_sub_core(data1, coords, listMLL, class1, class2, d, vecdist)
		}
	res
}


#####################
#Indice d'Aggregation
#####################


agg_index_core <- function(data1, coords = NULL, nbrepeat = 1, bar = FALSE, listMLL = NULL){
        
        N <- nrow(data1)
		rownames(data1) <- rownames(coords) <- 1:N
        dist_geo <- dist(coords, method = "euclidean", upper = FALSE, diag = FALSE)
        mat_dist_geo <- as.matrix(dist_geo)

if (length(listMLL) != 0){
        res_id_clonal <- id_clonal_MLL(data1, listMLL)
} else {
        res_id_clonal <- id_clonal(data1)
}

if (sum(as.dist(res_id_clonal)) != length(as.dist(mat_dist_geo))){

        coord_min_red <- nearest_unit(coords)

        recup <- NULL
                for (i in 1:nrow(coord_min_red)){
                        recup <- c(recup, res_id_clonal[coord_min_red[i,3], coord_min_red[i,2]])
                }

        Psg <- mean(as.dist(res_id_clonal))
        Psp <- mean(coord_min_red[,1]*recup)
        Ac <- (Psg-Psp)/Psg

if (bar){
        total <- nbrepeat
        pb <- txtProgressBar(min = 0, max = total, style = 3)
}
        recup_p <- NULL
                for (s in 1:nbrepeat){
                        order_sim <- sample(1:N, N, replace = FALSE)
                        mat_new <- data1[order_sim,]
if (length(listMLL) != 0){
        res_id_clonal_new <- id_clonal_MLL(mat_new, listMLL)
} else {
        res_id_clonal_new <- id_clonal(mat_new)
}

                        Psg_sim <- mean(as.dist(res_id_clonal_new))

                        recup <- NULL
                                for (i in 1:nrow(coord_min_red)){
                                        recup <- c(recup, res_id_clonal_new[coord_min_red[i,3], coord_min_red[i,2]])
                                }
                        Psp_sim <- mean(coord_min_red[,1]*recup)

                        Ac_sim <- (Psg_sim-Psp_sim)/Psg_sim
                        recup_p <- c(recup_p, Ac_sim)
if (bar){
                        setTxtProgressBar(pb, s)
}
                }

if (bar){
        close(pb)
}
		pval <- 2*min(mean(recup_p <= Ac), mean(recup_p >= Ac))
		if (pval > 1){pval <- 1}
		
list("results" = cbind(Ac, pval, nbrepeat), "simulation" = recup_p)
}
}


agg_index <- function(data1, coords = NULL, vecpop = NULL, nbrepeat = 1, bar = FALSE, listMLL = NULL){
	
	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		coordtot <- split(coords, vecpop)
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					rownames(coordtot[[p]]) <- 1:nrow(coordtot[[p]])
					res[[p]] <- agg_index_core(datatot[[p]], coordtot[[p]], nbrepeat, bar, listMLL[[p]])
				}
					names(res) <- unique(vecpop_o)
		} else {
			res <- agg_index_core(data1, coords, nbrepeat, bar, listMLL)
		}
	res
}


##############
#Effet de Bord
##############

edge_effect_core <- function(data1, coords = NULL, center = NULL, nbrepeat = 1, bar = FALSE, listMLL = NULL){
        
        N <- nrow(data1)
		rownames(data1) <- rownames(coords) <- 1:N

if (length(listMLL) != 0){
        list_genet <- listMLL
	  data_MLG <- data1[red_MLL(data1, listMLL),]
} else {
        list_genet <- MLG_list(data1)
        data_MLG <- unique(data1)
        G <- nrow(data_MLG)
}

if (length(unlist(list_genet)) != nrow(data1)){stop("Problem with MLL list detected")}

		count_MLG <- unlist(lapply(list_genet, length))

if (length(count_MLG) != N){

		rownames(data1) <- 1:N

        dist_geo_EE <- dist(rbind(center, coords), method = "euclidean", upper = FALSE, diag = FALSE)
        mat_dist_geo_EE <- as.matrix(dist_geo_EE)
        mat_dist_geo_EE <- mat_dist_geo_EE[1,][-1]
        mat_dist_geo_EE <- as.data.frame(mat_dist_geo_EE)
        rownames(mat_dist_geo_EE) <- 1:N

        mat_un_MLL <- cbind(count_MLG, data_MLG)
        mat_dist_recup <- mat_dist_geo_EE[rownames(mat_un_MLL[mat_un_MLL[,1] == 1,]),]
        mat_dist_recup <- as.data.frame(mat_dist_recup)
        row.names(mat_dist_recup) <- rownames(mat_un_MLL[mat_un_MLL[,1] == 1,])

	Du <- mean(mat_dist_recup[,1])
	Da <- mean(mat_dist_geo_EE[,1])

	Ee <- (Du-Da)/Da

if (bar){
        total <- nbrepeat
        pb <- txtProgressBar(min = 0, max = total, style = 3)
}

	recup_p_Ee <- NULL
for (s in 1:nbrepeat){
        order_sim <- sample(1:N, N, replace = FALSE)
        mat_sim <- data1[order_sim,]
        data_sim <- unique(mat_sim)

        	list_genet_sim <- lapply(list_genet, function(x) rownames(mat_sim)[x])
		data_sim <- mat_sim[unlist(lapply(X = list_genet_sim, FUN = `[[`, 1)),]

        count_MLG_sim <- unlist(lapply(list_genet_sim, length))
        mat_un_MLL_sim <- cbind(count_MLG_sim, data_sim)
        mat_dist_recup_sim <- mat_dist_geo_EE[rownames(mat_un_MLL_sim[mat_un_MLL_sim[,1] == 1,]),]
        mat_dist_recup_sim <- as.data.frame(mat_dist_recup_sim)
        row.names(mat_dist_recup_sim) <- rownames(mat_un_MLL_sim[mat_un_MLL_sim[,1] == 1,])
        
        Du_sim <- mean(mat_dist_recup_sim[,1])
        Ee_sim <- (Du_sim-Da)/Da
        
        recup_p_Ee <- c(recup_p_Ee, Ee_sim)

if (bar){
                        setTxtProgressBar(pb, s)
}

}

if (bar){
        close(pb)
}

pval_Ee <- 2*min(mean(recup_p_Ee <= Ee), mean(recup_p_Ee >= Ee))
	if (pval_Ee > 1){pval_Ee <- 1}

list("results" = as.data.frame(cbind(Ee, pval_Ee, nbrepeat)), "simulations" = recup_p_Ee)
}
}


edge_effect <- function(data1, coords = NULL, center = NULL, vecpop = NULL, nbrepeat = 1, bar = FALSE, listMLL = NULL){
	
	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datatot <- split(data1, vecpop)
		coordtot <- split(coords, vecpop)
		centertot <- center
			res <- list(NULL)
				for (p in 1:length(unique(vecpop))){
					rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
					rownames(coordtot[[p]]) <- 1:nrow(coordtot[[p]])
					res[[p]] <- edge_effect_core(datatot[[p]], coordtot[[p]], centertot[[p]], nbrepeat, bar, listMLL[[p]])
				}
			names(res) <- unique(vecpop_o)
		} else {
			res <- edge_effect_core(data1, coords, center, nbrepeat, bar, listMLL)
		}
	res
}


#################
#Fonction bonus :
#################

popsimgen <- function(data1, haploid = FALSE){

	N <- nrow(data1)

	if (haploid){
		index_l <- 1:ncol(data1)
	} else {
		index_l <- 1:c(ncol(data1)/2)*2-1
	}

	tab_sim <- as.data.frame(matrix(0, ncol = ncol(data1), nrow = N))

	for(n in 1:N){	
		indic <- sample(1:N, 2, replace = FALSE)
			for (l in index_l){
 				if (haploid){
					tab_sim[n,l] <- sample(data1[indic,l], 1)
				} else {
					tab_sim[n,c(l, l+1)] <- sort(c(as.numeric(sample(data1[indic[[1]],c(l, l+1)], 1)), as.numeric(sample(data1[indic[[2]],c(l, l+1)], 1))))
				}
			}
	}
	colnames(tab_sim) <- colnames(data1)
	rownames(tab_sim) <- rownames(data1)
tab_sim
}


genclone_core <- function(data1, data2, haploid = FALSE, coords = NULL, listMLL = NULL, nbrepeat = NULL, bar = FALSE){

	N <- nrow(data1)
	if (length(listMLL) != 0){
			list_genet <- listMLL
if (length(unlist(list_genet)) != nrow(data1)){stop("MLL list does not compute")}
			G <- length(listMLL)
			if(G == nrow(unique(data1))){
				L <- "MLG"
			} else {
				L <- "MLL"
			}
		} else {
			G <- nrow(unique(data1))
        		list_genet <- MLG_list(data1)
			L <- "MLG"
		}

	if(haploid){
		nb_all <- mean(sapply(apply(data1, 2, unique), length))
		SE <- sd(sapply(apply(data1, 2, unique), length))/sqrt(length(sapply(apply(data1, 2, unique), length)))
		res_Fis <- NA
		res_Fis_WR <- NA
		pval_Fis <- NA
		pval_Fis_WR <- NA
	} else {
		index_l <- 1:c(ncol(data1)/2)*2-1
		recup <- NULL
		for(i in index_l){ 
			recup <- c(recup, length(unique(c(data1[,i], data1[,i+1]))))
		}
		nb_all <- mean(recup)
		SE <- sd(recup)/sqrt(length(recup))

	res_Fis <- mean(Fis_core(data1, data2, genet = FALSE, RR = FALSE)[,4])
	res_Fis_WR <- mean(Fis_core(data1, data2, genet = TRUE, RR = FALSE)[,4])

		if(length(nbrepeat) != 0){
			Fis_sim <- NULL
			FisWR_sim <- NULL
			
			if (bar){
				total <- nbrepeat
				pb <- txtProgressBar(min = 0, max = total, style = 3)
			}
			
				for(s in 1:nbrepeat){
					tab_sim <- popsimgen(data1, haploid = FALSE)
					Fis_sim <- c(Fis_sim, mean(Fis(tab_sim, RR = FALSE, genet = FALSE)[,4]))
					FisWR_sim <- c(FisWR_sim, mean(Fis(tab_sim, RR = FALSE, genet = TRUE)[,4]))
					if (bar){
						setTxtProgressBar(pb, s)
					}
				}
				
				if (bar){
					close(pb)
				}
				
			pval_Fis <- 2*min(mean(Fis_sim <= res_Fis), mean(Fis_sim >= res_Fis))
			pval_Fis_WR <- 2*min(mean(FisWR_sim <= res_Fis_WR), mean(FisWR_sim >= res_Fis_WR))
		} else {
			pval_Fis <- NA
			pval_Fis_WR <- NA
		}
	}
	
if (G != N){

	Beta_P <- Pareto_index_core(data1, listMLL)[1]

	res1 <- autocorrelation_core(data1, haploid, coords, listMLL, 
	Loiselle = TRUE, Ritland = FALSE,
	genet = FALSE, central_coords = FALSE, random_unit = FALSE, weighted = FALSE,
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL,
	graph = FALSE, nbrepeat)

	Sp_L <- res1$Slope_and_Sp_index[1,3]
	pval_SpL <- res1$Slope_and_Sp_index[10,3]

	res2 <- autocorrelation_core(data1, haploid, coords, listMLL, 
	Loiselle = TRUE, Ritland = FALSE,
	genet = TRUE, central_coords = TRUE, random_unit = FALSE, weighted = FALSE,
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL,
	graph = FALSE, nbrepeat)

	Sp_L_WR <- res2$Slope_and_Sp_index[1,3]
	pval_SpLWR <- res1$Slope_and_Sp_index[10,3]

	res3 <- autocorrelation_core(data1, haploid, coords, listMLL, 
	Loiselle = FALSE, Ritland = TRUE,
	genet = FALSE, central_coords = FALSE, random_unit = FALSE, weighted = FALSE,
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL,
	graph = FALSE, nbrepeat)

	Sp_R <- res3$Slope_and_Sp_index[1,3]
	pval_SpR <- res1$Slope_and_Sp_index[10,3]

	res4 <- autocorrelation_core(data1, haploid, coords, listMLL, 
	Loiselle = FALSE, Ritland = TRUE,
	genet = TRUE, central_coords = TRUE, random_unit = FALSE, weighted = FALSE,
	class1 = FALSE, class2 = FALSE, d = NULL, vecdist = NULL,
	graph = FALSE, nbrepeat)

	Sp_R_WR <- res4$Slope_and_Sp_index[1,3]
	pval_SpRWR <- res1$Slope_and_Sp_index[10,3]


	R <- (G-1)/(N-1)
	count_MLG <- unlist(lapply(list_genet, length))
	freq_MLG <- count_MLG/N
 	Hpp <- -sum((count_MLG/N)*log(count_MLG/N))
	Jp <- Hpp/log(G)
	Ls <- sum((count_MLG*(count_MLG-1))/(N*(N-1)))
	Dp <- 1-Ls
	Dmin <- (((2*N-G)*(G-1))/(N*N))*(N/(N-1))
	Dmax <- ((G-1)/G)*(N/(N-1))
	V <- (Dp-Dmin)/(Dmax-Dmin)
	Hill <- 1/Ls
      
} else {
	L <- "no_clone"
	R <- Beta_P <- Sp_L <- pval_SpL <- Sp_L_WR <- pval_SpLWR <- Sp_R <- pval_SpR <-  
		Sp_R_WR <- pval_SpRWR <- Hpp <- Jp <- Dp <- V <- Hill <- NA
}

	tab <- as.data.frame(matrix(c(N, L, G, nb_all, SE, res_Fis, pval_Fis, res_Fis_WR, 
		pval_Fis_WR, R, Beta_P, Sp_L, pval_SpL, Sp_L_WR, pval_SpLWR, Sp_R, pval_SpR, 
		Sp_R_WR, pval_SpRWR, Hpp, Jp, Dp, V, Hill), nrow = 1))

	names(tab) <- c("N", "Lineage", "nb_L", "nb_all", "SE", "Fis", "pval_2sides", "Fis_WR", 
		"pval_2sides", "R", "Pareto_index", "Sp_Loiselle", "pval_2sides", "Sp_L_WR", 
		"pval_2sides", "Sp_Ritland", "pval_2sides", "Sp_R_WR", "pval_2sides", "H''", "J'", 
		"D", "V", "Hill")
tab
}


genclone <- function(data1, haploid = FALSE, coords = NULL, vecpop = NULL, listMLL = NULL, nbrepeat = NULL, bar = FALSE){

	if (length(vecpop) != 0){			
		if (length(vecpop) != nrow(data1)) {stop("vecpop length is not equal to the number of rows of your dataset.")}

		vecpop_o <- vecpop
		vecpop <- c(rep(1:length(unique(vecpop)), times = table(vecpop)[unique(vecpop)]))

		datafreq <- freq_RR(data1, haploid, vecpop, genet = FALSE, RR = FALSE)
		datatot <- split(data1, vecpop)
		coordtot <- split(coords, vecpop)
		res <- as.data.frame(matrix(NA, ncol = 24, nrow = 1))
		colnames(res) <- c("N", "Lineage", "nb_L", "nb_all", "SE", "Fis", "pval_2sides", "Fis_WR", 
			"pval_2sides", "R", "Pareto_index", "Sp_Loiselle", "pval_2sides", "Sp_L_WR", 
			"pval_2sides", "Sp_Ritland", "pval_2sides", "Sp_R_WR", "pval_2sides", "H''", "J'", 
			"D", "V", "Hill")
			for (p in 1:length(unique(vecpop))){
				rownames(datatot[[p]]) <- 1:nrow(datatot[[p]])
				rownames(coordtot[[p]]) <- 1:nrow(coordtot[[p]])
				res <- rbind(res, genclone_core(datatot[[p]], datafreq[[p]], haploid, coordtot[[p]], listMLL[[p]], nbrepeat, bar))
			}
		res <- res[-1,]
		rownames(res) <- unique(vecpop_o)
	} else {
		datafreq <- freq_RR(data1, haploid, vecpop, genet = FALSE, RR = FALSE)
		res <- genclone_core(data1, datafreq, haploid, coords, listMLL, nbrepeat, bar)
	}
	res
}
