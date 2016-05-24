f.make.design <- function(maternal, response, info, test = F, ret.characteristics = F){
#
# THE PROGRAM ESTIMATES EFFECTS OF SEVERAL ALLELES IN A CASE-TRIAD, CASE-CONTROL-TRIAD OR CASE-CONTROL DESIGN
# CREATES AN APPROPRIATE DESIGN MATRIX FOR USE IN f.tri.glm
#
# THE MODEL ASSUMES HARDY-WEINBERG EQUILIBRIUM (AND RARE DISEASE, IF NECESSARY, IN CC)
# IMPORTANT: OBSERVED FREQUENCIES (observed) MUST COME FROM A COMPLETE
# GRID IN APPROPRIATE ORDER! THE ORDER IS:
# TRIAD: BY MOTHER'S FIRST ALLELE AND SECOND ALLELE,
# THEN THE FATHER'S FIRST AND SECOND
# CASE-CONTROL-TRIAD: BY MOTHER'S FIRST ALLELE AND SECOND ALLELE,
# THEN THE FATHER'S FIRST AND SECOND, AND FINALLY CONTROL (1) AND CASE (2), 
# SO THAT TOP HALF OF THE DATA ARE CONTROLS
# CASE-CONTROL: BY FIRST ALLELE, SECOND ALLELE
# AND FINALLY CONTROL (1) AND CASE (2), SO THAT TOP HALF OF THE DATA ARE
# CONTROLS
#
# maternal DETERMINES WHETHER OR NOT MATERNAL EFFECTS SHOULD BE ESTIMATED
# (NOT AVAILABLE FOR THE CASE-CONTROL DESIGN)
#
#
## EXTRACT VARIABLES
.ref.cat <- info$haplos$ref.cat
.xchrom <- info$model$xchrom
.design <- info$model$design
.poo <- info$model$poo
.sel.sex <- info$variables$sel.sex
.sel.sex.yes <- !is.null(.sel.sex)
.comb.sex <- info$model$comb.sex
.covar.yes <- !is.null(info$variables$covar)
.A <- sum(info$haplos$selected.haplotypes) ## NUMBER OF ALLELES (HAPLOTYPES)
#
###########
# CREATE "BASE" GRID OF PARENTAL GENOTYPES (AND POSSIBLY SEX AND/OR CASE-CONTROL STATUS)
###########
.tmp <- factor(1:.A)
.mf.tmp <- list(m1 = .tmp, m2 = .tmp, f1 = .tmp, f2 = .tmp, sex = c(0,1), cc = c(0,1))
## HERE, 0 AND 1 (MALES AND FEMALES, RESP.) ARE USED AS DUMMY FOR sex. ORIGINAL CODING (OUTSIDE f.tri.glm AND f.make.design) SHOULD BE 1 (MALES) AND 2 (FEMALES)
## HERE, 0 AND 1 (CONTROL AND CASE, RESP.) ARE USED AS DUMMY FOR cc. ORIGINAL CODING (OUTSIDE f.tri.glm AND f.make.design) SHOULD BE 1 (CONTROL) AND 2 (CASE). IN INPUT FILE, THE LARGEST IS ALWAYS THE CASE
.cond.triad <- (.design %in% c("triad", "cc.triad"))
.cond.f1 <- !.xchrom & (.design != "cc")
.cond.sex <- .xchrom & !.sel.sex.yes
.cond.cc <- (.design %in% c("cc", "cc.triad"))
#
## m1
if(!.cond.triad) .mf.tmp[["m1"]] <- NULL
## m2
## f1
if(!.cond.f1)  .mf.tmp[["f1"]] <- NULL
## f2
## sex
if(!.cond.sex) .mf.tmp[["sex"]] <- NULL
## cc
if(!.cond.cc) .mf.tmp[["cc"]] <- NULL
#
##
#if((.design == "cc") & .xchrom) stop("Not implemented")## MERK: DE OVENSTAAENDE SELEKSJONENE SKAL FUNGERE FOR DENNE SITUASJONEN OGSAA!
if(.poo & (.design == "cc"))stop("Not implemented")## DETTE BLIR DET OGSAA SJEKKET FOR I f.check.pars

if(.covar.yes) .mf.tmp$covar <- factor(seq(along = info$variables$covar.codes))

if(ret.characteristics){
	.char <- sapply(.mf.tmp, length)
	## LITT AD-HOC, HAAPER DET IKKE TRENGS ANDRE STEDER:
	if(.design == "cc"){
		names(.char)[names(.char) == "m2"] <- "c1"
		names(.char)[names(.char) == "f2"] <- "c2"
	}
	return(.char)
}
#
## EXPAND GRID
.mf <- do.call("expand.grid", .mf.tmp)
#
###########
# CREATE DUMMY VARIABLES FOR PARENTAL GENOTYPES
###########
if(.cond.triad){
	.m1.dum <- model.matrix( ~ -1 + m1, data = .mf)
}
.m2.dum <- model.matrix( ~ -1 + m2, data = .mf)
if(.cond.f1){
	.f1.dum <- model.matrix( ~ -1 + f1, data = .mf)
}
.f2.dum <- model.matrix( ~ -1 + f2, data = .mf)
#
###########
# CREATE COUNTING VARIABLES FOR PARENTS AND CHILDREN
###########
## MOTHER
if(.cond.triad){
	.m.dumsum <- .m1.dum + .m2.dum
	dimnames(.m.dumsum)[[2]] <- paste("m", 1:.A, sep = "")
}else{
	.m.dumsum <- .m2.dum
}
## FATHER
if(.cond.f1){
	.f.dumsum <- .f1.dum + .f2.dum
}else{
	.f.dumsum <- .f2.dum
}
#
## PARENTS COMBINED, FOR ESTIMATING ALLELE FREQUENCIES
.parents.dumsum <- .m.dumsum + .f.dumsum
dimnames(.parents.dumsum)[[2]] <- paste("mf", 1:.A, sep = "")
#
## CHILD, DEFINE .c.dumsum:
if(.xchrom){
	## X-INACTIVATION(comb.sex == "double"): BOYS GET TWICE THE EFFECT. GIRLS GET ONE OR TWO.
	## MEANS THAT SINGLE DOSE IN BOYS CORRESPONDS TO DOUBLE DOSE IN GIRLS,
	## SINGLE DOSE IN GIRLS IS SQUARE ROOT OF DOUBLE DOSE (WHEN response = "mult")
	## BOYS: 1, RR^2, GIRLS: 1, RR, RR^2 (response = "mult")
	## BOYS: 1, RR1^2*RR2, GIRLS: 1, RR1, RR1^2*RR2 (response = "free")
	#
	## DOSE-RESPONSE(comb.sex == "single"): BOYS GET SINGLE DOSE EFFECT. GIRLS GET ONE OR TWO.
	## MEANS THAT SINGLE DOSE IN BOYS CORRESPONDS TO SINGLE DOSE IN GIRLS,
	## DOUBLE DOSE IN GIRLS IS SQUARE OF SINGLE DOSE. (WHEN response = "mult")
	## BOYS: 1, RR, GIRLS: 1, RR, RR^2 (response = "mult")
	## BOYS: 1, RR1, GIRLS: 1, RR1, RR1^2*RR2 (response = "free")
	#
	## RECALL THAT sex IS NOT PART OF .mf IF ONLY ONE SEX IS SELECTED
	if(.comb.sex %in% c("single", "double")){
		.sex <- .mf$sex
	}
	if(.comb.sex == "males"){
		.sex <- rep(0, nrow(.mf))
	}
	if(.comb.sex == "females"){
		.sex <- rep(1, nrow(.mf))
	}
	.m2.girls.dum <- .m2.dum * .sex # ALLELE FROM MOTHER, IN GIRLS
	.f2.girls.dum <- .f2.dum * .sex # ALLELE FROM FATHER, IN GIRLS
	.m2.boys.dum <- .m2.dum * (1 - .sex) # ALLELE FROM MOTHER, IN BOYS
	#
	## DEFINE COMBINATION
	if(.comb.sex == "double"){
		.c.dumsum <- .m2.girls.dum + .f2.girls.dum + 2 * .m2.boys.dum
	}else{
		.c.dumsum <- .m2.girls.dum + .f2.girls.dum + 1 * .m2.boys.dum
	}
}else{
	## CHILD, STANDARD AUTOSOMAL MODEL
	.c.dumsum <- .m2.dum + .f2.dum
}
dimnames(.c.dumsum)[[2]] <- paste("c", 1:.A, sep = "")
#
## WHEN RESPONSE IS FREE, ADD STAR TERM
if(response == "free"){
	## SPECIAL CODING FOR HOMOZYGOTES:
	.c.dd <- (.c.dumsum == 2) + 0	# MODEL: R^2*Rstar
	# SET APPROPRIATE COLUMN NAMES:
	dimnames(.c.dd)[[2]] <- paste("cdd", 1:.A, sep = "")	#
	#
	if(maternal){
		.m.dd <- (.m.dumsum == 2) + 0	# MODEL: R^2*Rstar FOR MAT. EFF
		# SET APPROPRIATE COLUMN NAMES:
		dimnames(.m.dd)[[2]] <- paste("mdd", 1:.A, sep = "")	
	}# END if(maternal)
}# END if(response == "free")
#
##

if(.covar.yes){
	.tmpmat <- cbind(.parents.dumsum, .mf[, c("covar"), drop = F])
	.form <- paste("mf", 1:.A, sep = "", collapse = " + ")
	.form <- paste("~ -1 + (", .form, "):covar", sep = "")
	.form <- formula(.form)
	.parents.dumsum <- model.matrix(.form, data = .tmpmat)
	cat("kontroller at sortering etc. blir rett!\n")
	.navn <- colnames(.parents.dumsum)
	.navn <- gsub(":", "_", .navn)

#	.navn <- gsub("_covar1", "", .navn)
#	cat("ad hoc!\n")

	colnames(.parents.dumsum) <- .navn

}



# SET UP FINAL DESIGN MATRIX, REMOVING REDUNDANT .ref.cat-ROW.
# NOTE: DIALLELIC SITUATION REQUIRES THE REMOVAL OF ONE DOUBLE DOSE PARAMETER:
#
## START WITH sex AND cc (WHICH ARE NULL IF NOT RELEVANT)
.design.matrix0 <- cbind(sex = .mf$sex, cc = .mf$cc)
## ADD .parents.dumsum, FOR ESTIMATING HAPLOTYPE FREQUENCIES
.design.matrix0 <- cbind(.design.matrix0, .parents.dumsum)
## KEEP THE EFFECTS PART SEPARATE TO BEGIN WITH, IN CASE DESIGN INCLUDES cc
.design.matrix1 <- NULL
#
## IF response IS MORE THAN simple, ADD THE NECESSARY DESIGN VARIABLES
if(response != "simple"){
	if(.poo){
		if(.xchrom){
			if(.comb.sex == "single"){
				## A (SINGLE) DOSE IN MALES CORRESPONDS TO A SINGLE, MATERNALLY INHERITED, DOSE IN FEMALES
				.mat <- .m2.girls.dum + 1 * .m2.boys.dum
				.pat <- .f2.girls.dum
			}
			if(.comb.sex == "double"){
				## A (SINGLE) DOSE IN MALES CORRESPONDS TO A DOUBLE DOSE IN FEMALES, I.E. A SINGLE MATERNALLY INHERITED + A SINGLE PATERNALLY INHERITED + THE Rstar DOUBLE DOSE
				## SO THAT A (SINGLE) DOSE IN MALES CORRESPONDS TO THE "FULL PACKAGE" IN FEMALES
				.mat <- .m2.girls.dum + 1 * .m2.boys.dum
				.pat <- .f2.girls.dum + 1 * .m2.boys.dum
			}
			if(.comb.sex == "males"){
				.mat <- .m2.boys.dum
				.pat <- NULL
			}
			if(.comb.sex == "females"){
				.mat <- .m2.girls.dum
				.pat <- .f2.girls.dum
			}
			## RENAME COLUMNS AND DROP REF CAT
			colnames(.mat) <- sub("m2", "cm", colnames(.mat))
			.mat <- .mat[,  - .ref.cat, drop = F]
			if(.comb.sex != "males"){
				colnames(.pat) <- sub("f2", "cf", colnames(.pat))
				.pat <- .pat[,  - .ref.cat, drop = F]
			}
			.design.matrix1 <- cbind(.design.matrix1, .mat, .pat)
		}else{
			## RENAME COLUMNS
			colnames(.m2.dum) <- sub("m2", "cm", colnames(.m2.dum))
			colnames(.f2.dum) <- sub("f2", "cf", colnames(.f2.dum))
			## DROP REF CAT
			.m2.dum <- .m2.dum[,  - .ref.cat, drop = F]
			.f2.dum <- .f2.dum[,  - .ref.cat, drop = F]
			.design.matrix1 <- cbind(.design.matrix1, .m2.dum, .f2.dum)
		}
	}else{
		.c.dumsum <- .c.dumsum[,  - .ref.cat, drop = F]
		.design.matrix1 <- cbind(.design.matrix1, .c.dumsum)
	}



	if(response == "free"){
		if(.A == 2) .c.dd <- .c.dd[,  - .ref.cat, drop = F]
		.design.matrix1 <- cbind(.design.matrix1, .c.dd)
	}
	if(maternal){
		.m.dumsum <- .m.dumsum[,  - .ref.cat, drop = F]
		.design.matrix1 <- cbind(.design.matrix1, .m.dumsum)
		if(response == "free"){
			if(.A == 2) .m.dd <- .m.dd[,  - .ref.cat, drop = F]
			.design.matrix1 <- cbind(.design.matrix1, .m.dd)
		}
	}
}
#
## LET THE FORMULA PART APPLY ONLY FOR CASES.
## (WARNING! RARE DISEASE ASSUMPTION HERE!)
if(.cond.cc){
	.design.matrix1 <- .design.matrix1 * .mf$cc
}
#
## JOIN AND CONVERT TO DATA FRAME
.design.matrix <- as.dframe(cbind(.design.matrix0, .design.matrix1))
#
##
return(.design.matrix)
}
