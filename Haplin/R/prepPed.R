prepPed <- function(pedfile, outdir, create.map = F, ask = T){
##
## EXTRACT PED-INFORMATION, CREATE INDEX FILE AND PHENOTYPE FILE
##
#
## PREPARE FILE NAMES
.basename <- basename(pedfile)
.basename <- file_path_sans_ext(.basename)
.index.name <- paste(outdir, "/", .basename, ".pedIndex", sep = "")
.pheno.name <- paste(outdir, "/", .basename, ".ph", sep = "")
.map.name <- paste(outdir, "/", .basename, ".map", sep = "")
#
##
if(ask){
	## PREVENT FILES FROM BEING UNINTENTIONALLY OVERWRITTEN
	if(create.map){
		.filnavn <- c(.index.name, .pheno.name, .map.name)
	}else{
		.filnavn <- c(.index.name, .pheno.name)
	}
	.test <- do.call("file.exists", as.list(.filnavn))
	if(any(.test)){
		cat("The following file(s) already exist(s):", .filnavn[.test], sep = "\n")
		.answer <- readline(paste('Overwrite file(s)? (y/n)', sep = ""))
		if(.answer != "y"){
			cat("Stopped without overwriting\n")
			return(invisible())
		}
	}
}
#
## EXTRACT FAMILY INFORMATION FROM PED FILE
cat("Extracting family, sex and case/control information from ped file...\n")
.fam <- as.dframe(extract_text_file_columns(pedfile, 1:6))
cat("Done\n")
gc()
names(.fam) <- c("family", "id", "father", "mother", "sex", "cc")
.pheno <- .fam[, c("id", "sex", "cc")]
.fam <- .fam[, c("family", "id", "father", "mother")]
#
##
.sex.u <- sort(unique(.pheno$sex))
if(length(.sex.u) > 2) stop("More than 2 different codes in the sex column (column 5)", call. = F)
if(any(!is.element(.sex.u, c("0", "1")))){
	if(all(is.element(.sex.u, c("1", "2")))){
		.pheno$sex[.pheno$sex == "2"] <- "0" ## RECODE FEMALES
	}else stop("Invalid codes in the sex column (column 5)", call. = F)
}
#
## CREATE INDEXING TO BE USED LATER WHEN CONVERTING TO HAPLIN FORMAT
.pedIndex <- f.make.index(.fam, output = "ids")
write.table(.pedIndex, file = .index.name, col.names = T, row.names = F, quote = T)
#
## WRITE A PHENOTYPE FILE TO BE USED IF NO OTHER IS PROVIDED
write.table(.pheno, file = .pheno.name, col.names = T, row.names = F, quote = T)
#
## IF REQUESTED, WRITE AN ARTIFICIAL MAP FILE
if(create.map){
	.nsnps <- (length(scan(pedfile, what = "", nlines = 1)) - 6)/2
	.map <- cbind(chrom  = 0, name = 1:.nsnps, position = 1:.nsnps)
	write.table(.map, file = .map.name, col.names = T, row.names = F, quote = T)
}
#
## FEEDBACK
cat("\nFiles written:\n")
cat(.index.name, "\n", .pheno.name, "\n", sep = "")
if(create.map) cat(.map.name, "\n", sep = "")
cat("\n")
#
return(invisible())
}
