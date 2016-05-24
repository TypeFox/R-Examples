TPL <-
function(splist, genus=NULL, species=NULL, infrasp=NULL, infra=TRUE, abbrev=TRUE, corr=TRUE, diffchar=2, max.distance=1, version="1.1", encoding="UTF-8", file="") {
splist2 <- NULL
try(splist2 <- splist, silent=TRUE)
if(!is.null(splist2) && (!is.null(genus) || !is.null(species) || !is.null(infrasp))) {
stop("argument 'splist' incompatible with arguments 'genus' and 'species'")
} else
if(is.null(splist2) && ((is.null(genus) && !is.null(species)) || (!is.null(genus) && is.null(species)))) {
stop("arguments 'genus' and 'species' must be provided")
} else
if(is.null(splist2) && !is.null(genus) && !is.null(species)) {
if(infra==TRUE && !is.null(infrasp)) {
splist <- paste(genus, species, infrasp)
} else
if(infra==FALSE || is.null(infrasp)) {
splist <- paste(genus, species)
}
}
#if(abbrev==TRUE) {
#abbr <- rep(NA, length(splist))
#splist0 <- splist

#vec0 <- c("nothossp. ", "nothossp.", " nothossp ", "nothosubsp. ", "nothosubsp.", " nothosubsp ", "cultivar. ", "cultivar.", " cultivar ", "subfo. ", " subfo ", "subf. ", "subf.", " subf ", " subproles ", "cf. ", "cf.", " cf ", "aff. ", "aff.", " aff ", "s.l. ", "s.l.", "s.l ", "s.str. ", "s.str.", "s.str ", "\u00D7", "x. ", "x.", " x ", "X. ", "X.", " X ", "X ", "f. ", "f.", " f ", "fo. ", "fo.", " fo ", " forma ", "subvar.", " subvar ", "var. ", "var.", " var ", "subsp. ", "subsp.", " subsp ", "ssp. ", "ssp.", " ssp ", " gama ", " grex ", "lus. ", "lus.", " lus ", " lusus ", "monstr. ", " monstr ", "nm. ", "nm.", " nm ", "prol. ", "prol.", " prol ", " proles ", " race ", "subvar. ", "cv. ", "cv.", " cv ")
#for(j in 1:length(vec0)) {
#for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep(vec0[j], splist[i], fixed=TRUE))>0, vec0[j], abbr[i])}
#splist <- sub(vec0[j], " ", splist, fixed=TRUE)
#}
#abbr <- gsub(" ", "", abbr, fixed=TRUE)
#ns <- sum(!is.na(abbr))
#print(paste(ns, "substitutions of standard annotations in specific and infraspecific epithets"))
#}

#splist <- gsub("      ", " ", splist, fixed=TRUE)
#splist <- gsub("     ", " ", splist, fixed=TRUE)
#splist <- gsub("    ", " ", splist, fixed=TRUE)
#splist <- gsub("   ", " ", splist, fixed=TRUE)
#splist <- gsub("  ", " ", splist, fixed=TRUE)
#splist <- gsub("_", "", splist, fixed=TRUE)
#for(i in 1:length(splist)) {splist[i] <- ifelse(substr(splist[i], 1, 1)==" ", substr(splist[i], 2, nchar(splist[i])), splist[i])} 


TPLck2 <- function(d) {
TPLck(sp=d, corr=corr, diffchar=diffchar, max.distance=max.distance, infra=infra, abbrev=abbrev, version=version, encoding=encoding)
}
results <- do.call("rbind", lapply(splist, TPLck2))
#results$Infraspecific <- ifelse(results$Infraspecific=="NA", "", as.character(results$Infraspecific))
#results <- data.frame(results[,1:2], Abbrev=as.character(abbr), results[,-c(1:2)])
if(infra==FALSE) {
results <- results[,-c(4,13)]
}
if(nchar(file)>0) {
write.csv(results, file=file, row.names=F)
}
return(results)
}
