#
# cSSBR.R
# Claas Heuer, September 2014
#
# Copyright (C)  2014 Claas Heuer
#
# This file is part of cpgen.
#
# cpgen is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# cpgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# in file.path(R.home("share"), "licenses").  If not, see
# <http://www.gnu.org/licenses/>.
#

# cSSBR

cSSBR <- function(data, M, M.id, X=NULL, par_random=NULL, scale_e=0, df_e=0, niter=5000, burnin=2500, seed=NULL, verbose=TRUE) {


# double check some stuff to prevent clmm from failing 
# after having already done some heavy computations

if(!is.null(X)) {
  if(!class(X)%in%c("matrix","dgCMatrix")) stop("'X' has to be of type 'matrix' or 'dgCMatrix'") 
  if(nrow(X) != nrow(data)) stop("'X' must have as many rows as 'data'")
  if(anyNA(X)) stop("No NAs allowed in 'X'")
} else {
    X = array(1,dim=c(nrow(data),1))
  }

if(!is.null(par_random)) {
  if(length(par_random)!=2) stop("'par_random' has to be of length 2")
  for(i in 1:length(par_random)) {
    allowed_methods = c("fixed","ridge","BayesA")
    if(is.null(par_random[[i]]$method)) stop(paste("Define a method for random-effect: ",i,sep=""))
    if(par_random[[i]]$method %in% allowed_methods == FALSE) stop(paste("Method must be one of: ",paste(allowed_methods,collapse=" , "),sep=""))
  }
}

# obtain the model terms
model_terms <- cSSBR.setup(data, M, M.id, verbose)
X <- X[match(model_terms$ids,data$id),,drop=FALSE]


#################
### Run Model ###
#################

if(verbose) cat(" Running Model\n")
mod <- clmm(y=model_terms$y,
	    X=X,
	    Z=list(model_terms$Marker_Matrix,model_terms$Z_residual),
# update 14.09.2015: pass ginverse
            ginverse = list(NULL, model_terms$ginverse_residual),
	    par_random=par_random,
	    scale_e=scale_e,
	    df_e=df_e,
	    niter=niter,
	    burnin=burnin,
	    seed=seed,
	    verbose=verbose) 

# obtain breeding values
bv <- mod$Predicted - (X %c% mod[[3]]$posterior$estimates_mean)

# add animals that were not in model to the bv-vector
gt_bv <- (M %c% mod[[4]]$posterior$estimates_mean)[,1]
gt_not_in_model <- M.id[!M.id %in% model_terms$ids]
bv <- c(bv,gt_bv[!M.id %in% model_terms$ids])
names(bv) <- c(model_terms$ids,gt_not_in_model)

# add some important information concerning SSBR
mod$SSBR <- list(ids = model_terms$ids, y=model_terms$y, X=X, Marker_Matrix=model_terms$Marker_Matrix, 
                 Z_residual = model_terms$Z_residual, ginverse_residual = model_terms$ginverse_residual,  
                 Breeding_Values = bv)

return(mod)

}



# cSSBR.setup

cSSBR.setup <- function(data, M, M.id, verbose=TRUE) {


# double check some stuff to prevent clmm from failing 
# after having already done some heavy computations
if(class(data)!="data.frame") stop("'data' has to be of type 'data.frame'") 

if(is.null(data$y)) stop("'data' has to include a phenotype-vector 'y'")
if(is.null(data$id)) stop("'data' has to include an id-vector 'id'")
if(is.null(data$sire)) stop("'data' has to include a sire-vector 'sire'")
if(is.null(data$dam)) stop("'data' has to include a dam-vector 'dam'")

if(class(M)!="matrix") stop("'M' has to be of type 'matrix'") 
if(is.null(M.id)) stop("Please provide rownames for 'M' as 'M.id'")
if(length(M.id)!=nrow(M)) stop("'M.id' must have as many elements as nrow(M)")
if(anyNA(M)) stop("No NAs allowed in marker matrix")

# FIXME repeated measures are not implemented yet - find a way not to copy too much
if(length(unique(data$id))<nrow(data)) stop("Repeated Measures not implmented yet")
 
if(verbose) cat("\n Gathering information and processing pedigree\n")


# Constructing pedigree: adding missings and ordering by generation.
# Individuals with markers might be missing in the data set, so we add those
# as unrelated founders

data$id <- as.character(data$id)
data$sire <- as.character(data$sire)
data$dam <- as.character(data$dam)
M.id <- as.character(M.id)

marker_without_ped <- M.id[!M.id %in% c(data$id,data$sire,data$dam)]
if(length(marker_without_ped)>0) {

  temp <- data.frame(id = marker_without_ped, sire = NA, dam = NA, y = NA)
  data <- cbind(data,temp)

}

# R-package: pedigreemm
# as of now (Sep. 2014) 'editPed' is very slow on large pedigrees.
# For a faster version install this modified package:
# library(devtools)
# install_github("cheuerde/pedigreemm", ref = "master", build_vignettes=FALSE)

# Update: 14.09.2015 - Dont force people to change pedigreem version, simply
# call an internally provided version of 'editPed' = 'editPed_fast'. See below
ped <- editPed_fast(label=data$id,sire=data$sire,dam=data$dam,verbose=FALSE)
ped <- pedigreemm::pedigree(label=ped$label,sire=ped$sire,dam=ped$dam)

# construct A-inverse
# Ainv <- as(getAInv(ped),"dgCMatrix")

T_Inv <- as(ped, "sparseMatrix")
D_Inv <- Diagonal(x=1/Dmat(ped))
Ainv<-t(T_Inv) %*% D_Inv %*% T_Inv
dimnames(Ainv)[[1]]<-dimnames(Ainv)[[2]] <-ped@label
Ainv <- as(Ainv, "dgCMatrix")

# now check which individuals have phenotypes and provide pedigree and/or marker information
DAT <- data.frame(id=ped@label, y = as.numeric(NA), has_y = 0, has_gt = 0, has_ped = 0, stringsAsFactors=FALSE)
n = nrow(DAT)

DAT[match(data$id,DAT$id),"y"] <- data[match(DAT[match(data$id,DAT$id),"id"],data$id),"y"]
DAT[!is.na(DAT$y),"has_y"] <- 1
DAT[DAT$id %in% M.id,"has_gt"] <- 1
DAT[DAT$id %in% ped@label,"has_ped"] <- 1

ids <- 1:n

genotyped <- ids[match(M.id,DAT$id)]

# from the non genotyped we only need those with phenotype
non_genotyped <- ids[DAT$has_y & !DAT$id %in% DAT[genotyped,"id"]]

# if non-genotyped dont provide phenotypes then SSBR is not a good model choice
if(length(non_genotyped) == 0) stop("There are no non-genotyped individuals that provide phenotypes")

# allocate the combined marker matrix for the model (only individuals with phenotype)
nrow_gt = sum(DAT[genotyped,"has_y"])
nrow_non_gt = length(non_genotyped)
tmp <- DAT[genotyped,]
index_gt <- seq(1,nrow(M))[match(tmp[tmp[,"has_y"]==1,"id"] ,M.id)]
# set up model ids
model_ids <- c(M.id[index_gt],DAT$id[non_genotyped])

#######################
### Imputation Step ###
#######################

### this part is recplaced now by the very memory efficient c++-function 'cSSBR_impute'
# M_combined <- matrix(double(),nrow = nrow_gt + nrow_non_gt, ncol = ncol(M))
# if(nrow_gt > 0){
#  M_combined[1:nrow_gt,] <- M[index_gt,]
# }

if(verbose) cat(" Allocating combined Marker matrix ( n =",nrow_gt + nrow_non_gt,", p =",ncol(M),")\n")
if(verbose) cat(" Imputing non-genotyped individuals\n")
# 'csolve' uses Eigen's 'SimplicalLLT-Solver' which is very fast for this purpose.
# Now impute the non_genotyped directly into the model matrix
#
# M_combined[(nrow_gt+1):nrow(M_combined),] <- csolve(Ainv[non_genotyped,non_genotyped],
# (-Ainv[non_genotyped,genotyped]) %c% M)

# This function runs in parallel with an absolute minimum of temporary memory allocation
M_combined <- .Call("cSSBR_impute",Ainv[non_genotyped,non_genotyped],
                   (-Ainv[non_genotyped,genotyped]), M, as.integer(index_gt), options()$cpgen.threads)

# obtain the cholesky factor for residual error of non_genotyped
# using pedigreemm::relfactor
#L <- t(as(relfactor(ped),"dgCMatrix"))

# set zeros for genotyped individuals
#if(nrow_gt > 0){
#  L[match(M.id[index_gt],DAT$id),] <- 0 
#}

# remove individuals we dont need and order
#L <- L[match(model_ids,DAT$id),non_genotyped]

# Update 14.09.2015: Insted of the cholesky we now use the
# submatrix of Ainv directly and pass it as 'ginverse' argument
# to clmm

Ainv11 <- Ainv[non_genotyped, non_genotyped]

Z <- sparse.model.matrix( ~ -1 + id, data=DAT, drop.unused.levels=FALSE)

if(nrow_gt > 0){
  Z[match(M.id[index_gt],DAT$id),] <- 0 
}

Z <- Z[match(model_ids,DAT$id),non_genotyped]


# phenotype vector
pheno <- DAT$y[match(model_ids,DAT$id)]

return(list(ids = model_ids, y=pheno, Marker_Matrix=M_combined, Z_residual = Z, ginverse_residual = Ainv11))

}




# this is an internal function that improves 'editPed' from
# 'pedigreemm' speed-wise.
# main body is copied and modified from the pedigreemm function 'editPed'

editPed_fast <- function(sire, dam, label, verbose = FALSE)
{
    nped <- length(sire)
    if (nped != length(dam))  stop("sire and dam have to be of the same length")
    if (nped != length(label)) stop("label has to be of the same length than sire and dam")
    tmp <- unique(sort(c(as.character(sire), as.character(dam))))
    
    missingP <-NULL
    if(any(completeId <- ! tmp %in% as.character(label))) missingP <- tmp[completeId]
    labelOl <- c(as.character(missingP),as.character(label))
    sireOl <- c(rep(NA, times=length(missingP)),as.character(sire))
    damOl  <- c(rep(NA, times=length(missingP)),as.character(dam))
    sire <- as.integer(factor(sireOl, levels = labelOl))
    dam  <- as.integer(factor(damOl, levels = labelOl))
    nped <-length(labelOl)
    label <-1:nped
    sire[!is.na(sire) & (sire<1 | sire>nped)] <- NA
    dam[!is.na(dam) & (dam < 1 | dam > nped)] <- NA
    pede <- data.frame(id= label, sire= sire, dam= dam, gene=rep(NA, times=nped))
    noParents <- (is.na(pede$sire) & is.na(pede$dam))
    pede$gene[noParents] <- 0
    pede$gene<-as.integer(pede$gene)

# new version that calls the C-function "get_generation"
# Note: there is no return value, the R-objects are directly being changed

    .Call("get_generation", pede$sire, pede$dam, pede$id, pede$gene, as.integer(verbose), package="cpgen")    

    ord<- order(pede$gene)
    ans<-data.frame(label=labelOl, sire=sireOl, dam=damOl, gene=pede$gene,
                    stringsAsFactors =F)
    ans[ord,]
}























