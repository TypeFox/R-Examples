###create function to create zcor from kinship matrix
mat2zcor <- function(kmat, id) {
    cls <- table(id)
    ncls <- length(unique(id))
    nzcor <- sum(cls*(cls-1)/2)
    zcor <- vector(mode="numeric", length=nzcor)
    zcor.idx <- 1
    mat.idx <- 1
    for(ncl in cls) {
        mat.start <- mat.idx
        mat.end <- mat.start - 1 + ncl

        if(ncl!=1) {
            kmat.tmp <- kmat[mat.start:mat.end,mat.start:mat.end]

            suppressMessages(zcor.tmp <- kmat.tmp[lower.tri(kmat.tmp)])
            zcor.start <- zcor.idx
            zcor.end <- zcor.start - 1 + length(zcor.tmp)
            zcor[zcor.start:zcor.end] <- zcor.tmp
            zcor.idx <- zcor.end + 1
        }
        mat.idx <- mat.end + 1
    }
    return(zcor)
}

#####reconditioning PCA matrix so that it will have a full rank
recond <- function(pca_mat) {
    pca_mat <- as.matrix(pca_mat)
    npc <- dim(pca_mat)[2]
    for(i in npc:1) {
        if(rankMatrix(as.matrix(pca_mat[,1:i])) == i) {
            break
        }
    }
    return(pca_mat[,1:i])
}

#maximum identity length contrast distance between string a and b start at position between k and k + 1
MILC <- function(a, b, k) {
   if(nchar(a)!=nchar(b)) {return("error")}
   if(k<1 || k >= nchar(a)) {return("index k need to be smaller than length of the strings")}
    str_a <- unlist(strsplit(a,""))
    str_b <- unlist(strsplit(b,""))
    len <- length(str_a)
    dist <- 0

    if(str_a[k]==str_b[k]) dist <- dist + 1

    for(i in (k+1):len) {
        if(str_a[i]==str_b[i]) dist <- dist + 1
        else break
    }

    for(i in (k-1):1) {
        if(str_a[i]==str_b[i]) dist <- dist + 1
        else break
    }
    return(dist)
}

#adding two haplotypes to make a genotype
mk.genotype <- function(hap1,hap2) {
    h1 <- unlist(strsplit(hap1,""))
    h2 <- unlist(strsplit(hap2,""))
    g <- vector("character",length=length(h1))
    for(i in 1:length(h1)) {
        g[i] <- as.character(as.numeric(h1[i]) + as.numeric(h2[i]))
    }
    return(paste(g,collapse=""))
}

#make genotype vector from haplotype vector
mk.genotype.vec <- function(dat) {
    ninds <- length(dat)/2
    g.vec <- vector("character", length=ninds)
    for(i in 1:ninds) {
        h1.idx <- i*2-1
        h2.idx <- i*2
        g.vec[i] <- mk.genotype(dat[h1.idx],dat[h2.idx])
    }
    return(g.vec)
}

#make diplotype vector where each element is
# an ordered concatenation of two haplotypes id separated by "_"
mk.diplotype.vec <- function(dat) {
    ninds <- length(dat) / 2
    diplo.vec <- vector("character", length=ninds)
    unique.haps <- unique(sort(dat))
    for(i in 1:ninds) {
        h1.idx <- i*2-1
        h2.idx <- i*2
        hap1 <- which(unique.haps==dat[h1.idx])
        hap2 <- which(unique.haps==dat[h2.idx])
        if(hap1 > hap2) {
            tmp <- hap1
            hap1 <- hap2
            hap2 <- tmp    
        } 
        diplo.vec[i] <- paste0(hap1,"_",hap2)
    }
    return(diplo.vec)
}

#make a map from genotype table to diplotype table
#the map contains all possible diplotypes for each genotype separated by comma
mk.geno.diplo.map <- function(geno.vec, diplo.vec, geno.table, diplo.table) {
    #make map from geno.table to diplo.table
    geno.diplo.map <- vector("character", length(geno.table))
    for(i in 1:length(geno.table)) {
        idx <- which(geno.vec == names(geno.table)[i])
        diplo <- unique(diplo.vec[idx])
        geno.diplo.map[i] <- paste0(which(names(diplo.table) %in% diplo),collapse=",")
    }
    return(geno.diplo.map)
}

####optimized version of X_hap
####pay attention to unique haplotype
####need to handle homozygous
mk.design.mat <- function(dat, geno.diplo.map, geno.table, geno.vec, diplo.table, diplo.vec) {
    ninds <- length(dat)/2
    unique.haps <- unique(sort(dat))
    nhaps <- length(unique.haps)
    X_hap <- matrix(0,ncol=nhaps,nrow=ninds)
    hap.length <- nchar(unique.haps[1])
    for(i in 1:ninds) {
        geno.table.idx <- which(names(geno.table)==geno.vec[i])
        diplo.table.idx <- geno.diplo.map[geno.table.idx]
        diplo.table.idx.arr <- as.numeric(unlist(strsplit(diplo.table.idx,",")))
        ngeno <- unname(geno.table[geno.table.idx])
        
        ####finding and making hap idx
        hap.idx <- as.numeric(unlist(strsplit(names(diplo.table[diplo.table.idx.arr]),"_")))
        ndiplo <- unname(rep(diplo.table[diplo.table.idx.arr],each=2))
        if(length(hap.idx) != length(unique(hap.idx))) {
            X_hap[i,hap.idx] <- ndiplo / (ngeno)
        } else {
            X_hap[i,hap.idx] <- ndiplo / (2*ngeno)
        }
    }
    return(X_hap)
}

####S is the similarity matrix of all distinct haplotypes in the data
mk.S <- function(dat, k) {
    unique.haps <- unique(sort(dat))
    nhaps <- length(unique.haps)
    hap.length <- nchar(unique.haps[1])
    S <- diag(rep(hap.length,nhaps))

    for(i in 1:(nhaps-1)) {
        for(j in (i+1):nhaps) {
            #cat(paste0("i:",i,",j:",j,"\n"))
            S[i,j] <- MILC(unique.haps[i], unique.haps[j], k)
            S[j,i] <- S[i,j]
        }
    }
    return(S)
}

####X_sum is simply the matrix product of X_hap and S
mk.X_sum <- function(X_hap,S) {
    return(X_hap %*% S)
}

####X_max[i,j] is the maximum of hadamard product between
####the i-th row of X_hap and j-th col of S
####introduce apply to speed up the function
mk.X_max <- function(X_hap,S) {
    X_max <- matrix(nrow=dim(X_hap)[1], ncol=dim(X_hap)[2])
    max.hadam <- function(x,hap) {
        return(max(hap * x))    
    }
    for(i in 1:dim(X_hap)[1]) {
        hap <- X_hap[i,]
        X_max[i,] <- apply(S,2,max.hadam, hap=hap)
    }
    return(X_max)
}

####Transform the design matrix into PCA matrix
####the pca.thres determines how many Principal Components to retain
####as a function of the percentage of explained variability
####I should add the loading matrix here 11/11/2015
mk.X_pca <- function(X_mat, pca.thres=0.95) {
    ###return null if X_mat is empty
    if(dim(X_mat)[2]==0) return(NULL)

    X_pca <- PCA(X_mat, ncp = rankMatrix(X_mat)[1], graph=F)
    #some eigenvalues are NaN

    #####check if pca failed
    if(all(is.na(X_pca$eig$eigenvalue))) {
        return(NULL)    
    }

    #if the conditional number of X_mat is too large
    if(any(is.nan(X_pca$eig$eigenvalue))) {
        nan.idx <- is.nan(X_pca$ind$coord[1,])
        rank_pca <- dim(X_pca$ind$coord)[2]
        X_pca_eig <- X_pca$eig[1:rank_pca,]
        X_pca_eig <- X_pca_eig[!nan.idx,]
        X_pca_eig$"percentage of variance" <- X_pca_eig$eigenvalue / sum(X_pca_eig$eigenvalue)
        X_pca_eig$"cumulative percentage of variance" <- cumsum(X_pca_eig[,2])*100
        npca <- sum(X_pca_eig[,3] <= (pca.thres*100))
    } else {
        npca <- sum(X_pca$eig[,3] <= (pca.thres*100))
    }

    #prevent the case when only one or no PC left
    #if(npca==1 || npca==0) npca <- length(X_pca$eig$eigenvalue)
    if(npca==0) npca <- 1	

    if(any(is.nan(X_pca$ind$coord[1,]))) {
        nan.idx <- is.nan(X_pca$ind$coord[1,])
        X_pca_dat <- X_pca$ind$coord[,!nan.idx]
        X_pca_dat <- as.matrix(X_pca_dat[,1:npca])
    } else {
        X_pca_dat <- X_pca$ind$coord[,1:npca]
        X_pca_dat <- as.matrix(recond(X_pca_dat))
    }
    npca <- dim(X_pca_dat)[2]
    colnames(X_pca_dat) <- paste0("PC",1:npca)
    return(X_pca_dat)
}

#####make kinship matrix
#####require kinship2 package
#####make sure phen is sorted by famid
mk.kinship.matrix <- function(phen, famid="famid", iid="iid", pid="pid", mid="mid", gender="sex", missid="") {
    phen <- phen[order(phen[,famid]),]
    missid <- as.character(missid)
    fam_id <- phen[,famid]
    ind_id <- paste0(fam_id,"_",phen[,iid])
    dad_id <- ifelse(phen[,pid]==missid, missid, paste0(fam_id,"_",phen[,pid]))
    mom_id <- ifelse(phen[,mid]==missid, missid, paste0(fam_id,"_",phen[,mid]))

    ####adding ungenotyped id
    allid <- unique(c(ind_id,dad_id[dad_id!=missid],mom_id[mom_id!=missid]))
    additional_id <- allid[which(!(allid %in% ind_id))]
    additional_famid <- gsub("^(.*)_.*", "\\1", additional_id)
    dad_idx <- which(additional_id %in% dad_id)
    mom_idx <- which(additional_id %in% mom_id)
    sex <- rep(0, length(additional_id))
    sex[dad_idx] <- 1
    sex[mom_idx] <- 2

    additional_ped <- data.frame(
                famid = as.character(additional_famid),
                id = as.character(additional_id),
                pid = as.character(missid),
                mid = as.character(missid),
                sex = sex,
                stringsAsFactors=F)

    ori.ped <- data.frame(famid = as.character(fam_id),
                        id = as.character(ind_id),
                        pid = as.character(dad_id),
                        mid = as.character(mom_id),
                        sex = phen[,gender],
                        stringsAsFactors=F)


    new.ped <- rbind(ori.ped, additional_ped)

    ####constructing kinship matrix
    pedAll <- pedigree(id=new.ped$id,
                    dadid=new.ped$pid,
                    momid=new.ped$mid,
                    sex=new.ped$sex,
                    famid=new.ped$famid,
                    missid=as.character(missid))

    kmat <- kinship(pedAll)
    ori.id <- which(kmat@Dimnames[[1]] %in% ind_id)
    kmat.ori <- kmat[ori.id,ori.id]

    return(kmat.ori)
}

#make contrast matrix for the general test
#nvar : number of independent variables
#ncov : number of covariates
#e.g. : Testing : B1=B2=B3 with 4 covariates
# 0 1 0 0 0 0 0 0
# 0 0 1 0 0 0 0 0
# 0 0 0 1 0 0 0 0
mk.contrast.matrix <- function(nvar, ncov) {
    C <- matrix(c(rep(0,nvar), diag(nvar)), ncol=nvar+1) 
    if(ncov > 0) {
        C <- cbind(C, matrix(rep(0, nvar * ncov), ncol=ncov, nrow=nvar))  
    }   
    return(C)
}

#####return pvalue from the test
#####X_pca is the design matrix
SHLR <- function(X_pca, phen, outcomes, out_dist="binomial", cov="", test="Rao", is.pval=T) {

    ####if X_pca has no variables
    if(dim(X_pca)[2]==0) return(1)
   
    ####make formula
    colnames(X_pca) <- paste0("PC",1:dim(X_pca)[2])

    ####if X_pca only has one variables
    if(dim(X_pca)[2]==1) {
    	formula_pca <- (paste0(outcomes," ~ PC1 "))
    } else {
    	formula_pca <- (paste0(outcomes," ~ PC1 ", paste0("+ PC", 2:dim(X_pca)[2], collapse=" ")))
    }
    if(cov[1]!="") formula_pca <- paste0(formula_pca, " + ",paste0(cov,collapse=" + "))
    formula_pca <- as.formula(formula_pca)

    dat_pca <- cbind(phen, X_pca)
    formula_null <- outcomes
    if(cov[1]!="") {
        formula_null <- paste0(formula_null, " ~ ",paste0(cov,collapse=" + "))
    } else {
        formula_null <- paste0(formula_null, " ~ 1")
    }
    formula_null <- as.formula(formula_null)
    glm.fit <- glm(formula_pca, family=out_dist, data=dat_pca)
    glm.fit.null <- glm(formula_null, family=out_dist, data=dat_pca)
    pval <- anova(glm.fit,glm.fit.null,test=test)[2,"Pr(>Chi)"]
 
    if(is.pval) {
        return(pval)
    }
    else {
        return(glm.fit)
    }
}

#####return pvalue from the test, depend on X_pca
SHLR.fam <- function(X_pca, phen, outcomes, out_dist="gaussian", cov="", fid="famid", iid="iid", pid="pid", mid="mid", gender="sex", missid="0", corstr="exchangeable", std.err="naive") {

    ####make formula
    colnames(X_pca) <- paste0("PC",1:dim(X_pca)[2])
    formula_pca <- (paste0(outcomes," ~ ", paste0("PC", 1:dim(X_pca)[2], collapse=" + ")))
    if(cov[1]!="") formula_pca <- paste0(formula_pca, " + ",paste0(cov,collapse=" + "))
    formula_pca <- as.formula(formula_pca)

    ####make kinship matrix
    ####WARNING: no "_" character can be used as identifier in id columns
    if(corstr=="kinship") {
        kmat <- mk.kinship.matrix(phen, famid=fid, iid=iid, pid=pid, mid=mid, gender=gender, missid=missid)
    }

    ####remove missing data
    if(cov[1] != "") {
       #na.idx <- unique(which(is.na(phen[,c(outcomes,cov)])) %% dim(phen)[1])
       na.idx <- !complete.cases(phen[,c(outcomes,cov)])
    } else {
       na.idx <- !complete.cases(phen[,outcomes])
       #na.idx <- unique(which(is.na(phen[,outcomes])) %% dim(phen)[1])
    }
    if(sum(na.idx)!=0) {
        phen <- phen[!(na.idx),]
        X_pca <- as.matrix(X_pca[!(na.idx),])
        colnames(X_pca) <- paste0("PC",1:dim(X_pca)[2])
        if(corstr=="kinship")
            kmat <- kmat[!na.idx,!na.idx]
    }

    dat_pca <- cbind(phen, X_pca)

    ####get gee fit
    gee.fit <- tryCatch({
        if(corstr=="kinship") {
            zcor <- mat2zcor(kmat, phen[,fid])
            geeglm(formula_pca, family=out_dist, data=dat_pca, id=dat_pca[,fid], corstr="fixed", zcor=zcor)
        } else if(corstr=="exchangeable") {
            geeglm(formula_pca, family=out_dist, data=dat_pca, id=dat_pca[,fid], corstr="exchangeable")
        } else {
            geeglm(formula_pca, family=out_dist, data=dat_pca, id=dat_pca[,fid], corstr="independence")
        }
    }, warning=function(w) {
        print("warning")
    }, error=function(e) {
        print("error")
    }, finally={
    })
    ###This part is the test statistic part
    ###Prereq : gee.fit, cov, X_pca, std.err
    ###Contrast matrix
    npc <- dim(X_pca)[2]
    if(cov[1]!="") {
        C_pca <- mk.contrast.matrix(npc, length(cov))
    } else {
        C_pca <- mk.contrast.matrix(npc,0)
    }
    q <- rankMatrix(C_pca)[1]

    ##check if error
    if(is.character(gee.fit)) {
        pval <- 2
    } else {
        B <- as.vector(gee.fit$coefficients)
        if(std.err=="naive") {
            S <- gee.fit$geese$vbeta.naiv
        } else {
            S <- gee.fit$geese$vbeta
        }
        test <- t(C_pca%*%B) %*% ginv(C_pca %*% S %*% t(C_pca)) %*% C_pca %*% B
        pval <- pchisq(test, df=q, lower.tail=F)
    }
    return(pval)
}


gee.pca.null <- function(X_pca, phen, outcomes, out_dist="gaussian", cov="", famid="famid", iid="iid", pid="pid", mid="mid", gender="sex", missid="0", corstr="exchangeable", std.err="naive") {

    ####make formula

    if(cov[1]=="") {
        formula_pca <- paste0(outcomes, " ~ 1")
    } else {
        formula_pca <- paste0(outcomes, " ~ ",paste0(cov,collapse=" + "))
    }
    formula_pca <- as.formula(formula_pca)
    
    ####make kinship matrix
    ####WARNING: no "_" character can be used as identifier in id columns
    if(corstr=="kinship") {
        kmat <- mk.kinship.matrix(phen, famid=famid, iid=iid, pid=pid, mid=mid, gender=gender, missid=missid)
    }

    ####remove missing data
    if(cov[1] != "") {
       #na.idx <- unique(which(is.na(phen[,c(outcomes,cov)])) %% dim(phen)[1])
       na.idx <- !complete.cases(phen[,c(outcomes,cov)])
    } else {
       na.idx <- !complete.cases(phen[,outcomes])
       #na.idx <- unique(which(is.na(phen[,outcomes])) %% dim(phen)[1])
    }
    if(sum(na.idx)!=0) {
        phen <- phen[!(na.idx),]
        X_pca <- as.matrix(X_pca[!(na.idx),])
        colnames(X_pca) <- paste0("PC",1:dim(X_pca)[2])
        if(corstr=="kinship")
            kmat <- kmat[!na.idx,!na.idx]
    }

    dat_pca <- cbind(phen, X_pca)

    ####get gee fit
    gee.fit <- tryCatch({
        if(corstr=="kinship") {
            zcor <- mat2zcor(kmat, phen[,famid])
            geeglm(formula_pca, family=out_dist, data=dat_pca, id=famid, corstr="fixed", zcor=zcor)
        } else if(corstr=="exchangeable") {
            geeglm(formula_pca, family=out_dist, data=dat_pca, id=famid, corstr="exchangeable")
        } else {
            geeglm(formula_pca, family=out_dist, data=dat_pca, id=famid, corstr="independence")
        }
    }, warning=function(w) {
        print("warning")
    }, error=function(e) {
        print("error")
    }, finally={
    })

    return(gee.fit)
}

####read data from shapeit phased haplotypes format
####return array of string
####the successive row pair (2k-1,2k) corresponds to the first haplotype and
####the second haplotype of the k-th individual.
read.shapeit.haps <- function(infile)
{
	shapeit.haps <- read.table(infile, header=F, stringsAsFactors=F)
	out.array <- list()
	out.array$marker.map <- shapeit.haps[,2:3]
	shapeit.haps <- shapeit.haps[,-(1:5)]

	out.dat <- t(shapeit.haps)
	out.array$haps <- unname(apply(out.dat,1,paste,collapse=""))
	out.array
}

########
#run.gee.pca <- function(X_pca, phen, outcomes, out_dist="gaussian", cov="", fid="famid", iid="iid", pid="pid", mid="mid", gender="sex", missid="0", corstr="exchangeable", std.err="naive") {
#run.glm.pca <- function(X_pca, phen, outcomes, out_dist="binomial", cov="", missid="0", test="Rao", is.pval=T) {
#The list of general parameters. These parameters will be used for both SHLR and SHLR-fam:
#haps : Haplotypes array
#phen : Phenotype table
#marker.map : The filename that contain the map of marker name to marker position
#outcome : The name of the outcome variable
#out.dist : The assumed distribution of the outcome variable
#cov : The name of the covariates
#missing.code : The character used to represent missing data
#model : SHLR/SHLR-fam
#pca.thres : The proportion of the variance threshold. This parameter will determine how many principal components will be used in the model. Default : NULL
#hap.freq.thres : The haplotype frequency threshold. This parameter will determine how many haplotypes will be used in the model. Default : 0.005
#window.size : The number of markers to be included in the haplotype model. Default : 100
#threads : The number of threads to be used to run this function in parallel. Default : NULL
#nSplits : The number of chunks. The output will be written incrementally by chunk. Chunk size will be determined by the number of splits.

#Parameters specific to SHLR-fam
#corstr : Correlation structure (exchangeable, kinship, independent, or unstructured)
#std.err : Variance estimator for SHLR-fam ("naive","sandwich","j1k"). Default : "naive"

#These additional parameters will be required if corstr=="kinship"
#gender : The name of the gender variable in the phenotype table
#fid : The name of the family id in the phenotype table
#iid : The name of the individual id in the phenotype table
#pid : The name of the paternal id in the phenotype table
#mid : The name of the maternal id in the phenotype table

#output 
#outfile : The name of the file to store the output table
#output format
#MarkerName
#Chr
#Position
#SHLR_Sum_Pval
#SHLR_Max_Pval


run.SHLR.scan <- function(haps, phen, marker.map, outfile, outcome, out.dist="Gaussian", cov="", missing.id.code="0", method, nSplits=1, pca.thres=NULL, hap.freq.thres=0.005, window.size=100, threads=NULL, corstr="kinship", std.err="naive", gender="", fid="", iid="", pid="", mid="") {

	#check whether haps and phen have the same number of individuals.
	if(length(haps)/2 != dim(phen)[1]) stop("This function requires the same number of individuals in the haplotype array and the phenotype table")

	#check outcome and cov
	colvar <- c(outcome,cov)
	if(!(outcome %in% colnames(phen))) stop("Outcome variable is not found in the phenotype table")
	if(cov[1] != "") {
		idx.cov <- cov %in% colnames(phen)
		if(!all(idx.cov)) {
			stop(paste0("The following variables are not found in the phenotype table: ", paste0(cov[!idx.cov], collapse=", ")))
		}
	}

	#check outcome distribution
	if(!(out.dist %in% c("gaussian", "binomial"))) stop("Currently the assumed distribution of the outcome variable is restricted to either 'gaussian' or 'binomial'")	
	
	#check window.size
	if(!(window.size >=2 && window.size <= nchar(haps[1]))) stop("Window size needs to be set between 2 and the length of the haplotype")

	#check hap.freq.thres
	if(!(hap.freq.thres >=0 && hap.freq.thres <=1)) stop("hap.freq.thres argument needs to be between 0 and 1") 

	#check pca.thres
	if(!(is.null(pca.thres)) && !(pca.thres >=0 && pca.thres <=1)) stop("pca.thres argument needs to be between 0 and 1")

	#check threads : need to be positive integer or NULL for single thread
	if(!(is.null(threads)) && (round(threads) < 1)) stop("threads argument need to be larger than 1")
	
	#determine method type: SHLR/SHLR-fam
	if(!(method %in% c("SHLR", "SHLR-fam"))) stop("The method argument needs to be one of 'SHLR' or 'SHLR-fam'")

	#check SHLR-fam specific parameters: corstr and std.err
	if(method=="SHLR-fam") {
		if(!(corstr %in% c("unstructured", "kinship", "independence", "exchangeable"))) stop("The corstr argument need to be one of 'unstructured', 'kinship', 'independence', or 'exchangeable'")
		
		if(!(std.err %in% c("naive", "sandwich"))) stop("The std.err argument need to be one of 'naive' or 'sandwich'")

		#check 'kinship corstr' specific parameters : gender, fid, iid, pid, mid
		if(corstr=="kinship") {
			fam.var <- c(gender, fid, iid, pid, mid)
			idx.fam.var <- fam.var %in% colnames(phen)
			
			if(any(fam.var=="")) stop("The following arguments need to be initialized : gender, fid, iid, pid, mid")

			if(!all(idx.fam.var)) {
				stop(paste0("The following variables are not found in the phenotype table: ", paste0(fam.var[!idx.fam.var], collapse=", ")))
			}
		}
	}

	#create output dataframe
	max.hap.length <- nchar(haps[1])

	#k is the middle position of the window
	k <- floor(window.size / 2)	

	#prepare threads
	registerDoParallel(cores=threads)

	#the actual computation starts here
	###need to split the data for efficiency
	split.size <- ceiling(max.hap.length / nSplits)
	if(split.size <= window.size) stop("The chunk size is too small. Please decrease the number of splits")

	last.split.size <- max.hap.length - ((nSplits-1) * split.size)
	split.start <- vector(length=nSplits)
	split.end <- vector(length=nSplits)	
	for(s in 1:nSplits) {
		if(s==1) split.start[s] <- 1
		else split.start[s] <- (s-1)*split.size + 2 - window.size
		
		if(s==nSplits) split.end[s] <- split.size * (s-1) + last.split.size 
		else split.end[s] <- split.size * s
	}	

	####write the header to outfile
	header <- c("Marker.Name","Position","SHLR_sum_pval", "SHLR_max_pval")
	write(paste0(header,collapse=","), file=outfile, sep=",")
	
	print("SHLR scan started....")
	for(s in 1:nSplits) {
		
		splitted.haps <- substr(haps,split.start[s],split.end[s])
		nTests <- split.end[s] - (split.start[s]-1) - window.size + 1
		out.df <- data.frame(matrix(nrow=nTests, ncol=4))
		marker.idx <- 0
		pval.df <- foreach(marker.idx=1:nTests, .combine=rbind) %dopar% {
			windowed.haps <- substr(splitted.haps, marker.idx, marker.idx+window.size-1)
			
			hap.table <- table(windowed.haps)
			if(hap.freq.thres==0) {
				unique.haps.idx <- which(hap.table > 1)
			} else {
				hap.thres = hap.freq.thres * sum(hap.table)
				unique.haps.idx <- which(hap.table >= hap.thres)
				if(length(unique.haps.idx)==0) stop(paste0("There are no haplotypes with haplotype frequency >= ", hap.freq.thres*100, "%"))
			}

			#genotype data structure	
			geno.vec <- mk.genotype.vec(windowed.haps)
			geno.table <- table(geno.vec)

			#diplotype data structure
			diplo.vec <- mk.diplotype.vec(windowed.haps)
			diplo.table <- table(diplo.vec)

			#geno-diplo map structure
			geno.diplo.map <- mk.geno.diplo.map(geno.vec, diplo.vec, geno.table, diplo.table)

			#similarity matrix
			S <- mk.S(windowed.haps, k)

			#Defining design matrices
			X_hap <- mk.design.mat(windowed.haps, geno.diplo.map, geno.table, geno.vec, diplo.table, diplo.vec)
			X_sum <- mk.X_sum(X_hap, S)
			X_max <- mk.X_max(X_hap, S)

			X_sum <- X_sum[,unique.haps.idx]
			X_max <- X_max[,unique.haps.idx]

			#compute PCA
			if(!is.null(pca.thres)) {
				X_sum <- mk.X_pca(X_sum, pca.thres)
				X_max <- mk.X_pca(X_max, pca.thres)
			}

			###depend on method
			if(method=="SHLR") {
				sum.pval <- SHLR(X_sum, phen, outcome, out.dist, cov)
				max.pval <- SHLR(X_max, phen, outcome, out.dist, cov)
			} else if(method=="SHLR-fam") {
				sum.pval <- SHLR.fam(X_sum, phen, outcome, out.dist, cov, fid, iid, pid, mid, gender, missing.id.code, corstr, std.err) 
				max.pval <- SHLR.fam(X_max, phen, outcome, out.dist, cov, fid, iid, pid, mid, gender, missing.id.code, corstr, std.err) 
			}
			pval <- c(sum.pval, max.pval)
		}
		
		out.df <- cbind(marker.map[(split.start[s]+k-1):(split.start[s]+k+nTests-2),1],marker.map[(split.start[s]+k-1):(split.start[s]+k+nTests-2),2],pval.df)
		write.table(out.df, file=outfile, row.names=F, quote=F, sep=",", append=T, col.names=F)
		print(paste0(s/nSplits * 100, "% ..."))
	}
	print("SHLR scan has been completed...")
}

