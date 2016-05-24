f.pos.to.haplocomb <- function(pos, A, comb, fam = "mf"){
#
#
# pos = POSITION IN GRID (numeric vector), A = vector WITH NUMBER OF ALLELES
# AT EACH MARKER
# comb = ALLELE COMBINATION (numeric matrix with four X length(A) columns)
#
# FROM ALLELE COMBINATIONS AT EACH MARKER, COMPUTES POSITION (ROW NUMBER) IN FULL 
# HAPLOTYPE GRID WITH FOUR X length(A) COLUMNS, OR INVERSELY,
# COMPUTES GRID "COORDINATES" (ALLELE VALUES IN A ROW) FOR A GIVEN POSITION (ROW NUMBER).
#
# fam SAYS WHETHER THIS IS A MOTHER-FATHER ALLELE COMBINATION, OR JUST CHILD ALLELES
#
	.nmarkers <- length(A)
	.nhaplo <- prod(A)
#
	if(fam == "mf"){
		.navn0 <- c("m1", "m2", "f1", "f2")
	}
	if(fam == "mfx"){
		.navn0 <- c("m1", "m2", "f2")
	}
	if(fam == "c"){
		.navn0 <- c("c1", "c2")	
	}
	.rep <- length(.navn0)

	if(T){
		if(!missing(pos)){
		# FROM POSITION TO COMBINATION
			.haplocomb <- f.pos.in.grid(pos = pos, A = rep(A, .rep))
			.navn <- paste("l", 1:.nmarkers, sep = "")
			.navn <- as.character(outer(.navn, .navn0, function(x,y) paste(x,y,sep=".")))
			dimnames(.haplocomb) <- list(pos, .navn)
		return(.haplocomb)	
		}

		if(!missing(comb)){
		# FROM COMBINATION TO POSITION, CONVERTS SINGLE NUMERIC TO MATRIX IF NECESS.
			if(!is.matrix(comb)){
				if(is.numeric(comb) && (length(comb) == .rep*.nmarkers)){
					comb <- matrix(comb, ncol = .rep*.nmarkers)
				}else stop("Input is not a numeric matrix!")
			}
			.pos <- f.pos.in.grid(comb = comb, A = rep(A, .rep))
		return(.pos)
			}
		}
#
#
}
