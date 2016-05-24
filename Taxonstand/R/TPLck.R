TPLck <-
function(sp, corr=TRUE, diffchar=2, max.distance=1, infra=TRUE, abbrev=TRUE, version="1.1", encoding="UTF-8") {
sp <- as.character(sp)

abbr <- NA
if(abbrev==TRUE) {
vec0 <- c("nothossp. ", "nothossp.", " nothossp ", "nothosubsp. ", "nothosubsp.", " nothosubsp ", "cultivar. ", "cultivar.", " cultivar ", "subfo. ", " subfo ", "subf. ", "subf.", " subf ", " subproles ", "cf. ", "cf.", " cf ", "aff. ", "aff.", " aff ", "s.l. ", "s.l.", "s.l ", "s.str. ", "s.str.", "s.str ", "\u00D7", "x. ", "x.", " x ", "X. ", "X.", " X ", "f. ", "f.", " f ", "fo. ", "fo.", " fo ", " forma ", "subvar.", " subvar ", "var. ", "var.", " var ", "subsp. ", "subsp.", " subsp ", "ssp. ", "ssp.", " ssp ", " gama ", " grex ", "lus. ", "lus.", " lus ", " lusus ", "monstr. ", " monstr ", "nm. ", "nm.", " nm ", "prol. ", "prol.", " prol ", " proles ", " race ", "subvar. ", "cv. ", "cv.", " cv ")

for(j in 1:length(vec0)) {
abbr <- ifelse(length(grep(vec0[j], sp, fixed=TRUE))>0, vec0[j], abbr)
sp <- sub(vec0[j], " ", sp, fixed=TRUE)
}
abbr <- gsub(" ", "", abbr, fixed=TRUE)
sp <- gsub("      ", " ", sp, fixed=TRUE)
sp <- gsub("     ", " ", sp, fixed=TRUE)
sp <- gsub("    ", " ", sp, fixed=TRUE)
sp <- gsub("   ", " ", sp, fixed=TRUE)
sp <- gsub("  ", " ", sp, fixed=TRUE)
sp <- gsub("_", "", sp, fixed=TRUE)
sp <- ifelse(substr(sp, 1, 1)==" ", substr(sp, 2, nchar(sp)), sp)
}
Abbrev <- abbr

genus <- unlist(strsplit(sp," "))[1]
species <- unlist(strsplit(sp," "))[2]
if(corr==TRUE) {species <- tolower(species)}
infrasp <- ifelse(length(unlist(strsplit(sp," "))) > 2, unlist(strsplit(sp," "))[3], "")
if(corr==TRUE) {infrasp <- tolower(infrasp)}
vv <- ifelse(version == "1.0", "", version)

searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",genus,"+",species, "&csv=true", sep="")
table.sp <- NULL
try(table.sp <- read.csv(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character", as.is = TRUE, encoding=encoding),silent=TRUE)
Genus <- genus
Species <- species
Infraspecific <- infrasp
New.Hybrid.marker <- ""
marker <- FALSE
marker.infra <- FALSE

# Loop 1
if(is.null(table.sp)) {
Family <- ""
Taxonomic.status <- ""
Plant.Name.Index <- FALSE
New.Genus <- Genus
New.Species <- Species
New.Infraspecific <- Infraspecific
Authority <- ""
Typo <- FALSE
WFormat <- TRUE
ID <- ""
New.ID <- ""
} else
# Loop 2
if(!is.null(table.sp)) {
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]
# Loop 2.1
if(k == 1) {
Family <- ""
Plant.Name.Index <- FALSE
Taxonomic.status <- ""
New.Genus <- Genus
New.Species <- Species
New.Infraspecific <- Infraspecific
Authority <- ""
Typo <- FALSE
WFormat <- FALSE
ID <- ""
New.ID <- ""
} else
# Loop 2.2
if(k > 1) {
# Loop 2.2.1
if(z == 0) {
Family <- ""
Plant.Name.Index <- FALSE
Taxonomic.status <- ""
New.Genus <- Genus
New.Species <- Species
New.Infraspecific <- Infraspecific
Authority <- ""
Typo <- FALSE
WFormat <- FALSE
ID <- ""
New.ID <- ""
} else
# Loop 2.2.2
if(z > 1) {
# Loop 2.2.2.1
if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 1 && corr==TRUE) {
spx <- length(agrep(species, "sp", max.distance=0)) + length(agrep(species, "sp.", max.distance=0))
mf <- c(as.character(1:1000))
is.mf <- agrep(species, mf, max.distance=list(deletions=1), value=TRUE)
cck <- agrep(species, table.sp$Species, value=TRUE, max.distance=max.distance)
ddf <- abs(nchar(cck) - nchar(species))
if(length(cck) > 0) {
cck <- cck[ddf==min(ddf)]
ddf <- abs(nchar(cck) - nchar(species))
}
levs <- length(unique(cck))
if(length(is.mf) == 0 && length(cck) > 0 && ddf <= diffchar && levs == 1 && spx == 0) {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",genus,"+",cck[1], "&csv=true", sep="")
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE,colClasses = "character", as.is = TRUE, encoding=encoding), silent=T)
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]
marker <- TRUE
}
}
# Loop 2.2.2.2
if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 1) {
Family <- ""
Plant.Name.Index <- FALSE
Taxonomic.status <- ""
New.Genus <- Genus
New.Species <- Species
New.Infraspecific <- Infraspecific
Authority <- ""
Typo <- FALSE
WFormat <- FALSE
ID <- ""
New.ID <- ""
}
# Loop 2.2.2.3
if (length(unique(paste(table.sp$Genus, table.sp$Species))) == 1) {
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
if ((length(ngrep) == 0 || abs(ngrep-nchar(infrasp)) > 0) && corr==TRUE && !is.na(infrasp) && nchar(infrasp)>0) {
cck.infra <- agrep(infrasp, table.sp$Infraspecific.epithet, value=TRUE, max.distance=max.distance)
ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
if(length(cck.infra) > 0) {
cck.infra <- cck.infra[ddf.infra==min(ddf.infra)]
ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
}
levs <- length(unique(cck.infra))
if(length(cck.infra) > 0 && levs == 1) {
infrasp <- unique(cck.infra)
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
marker.infra <- TRUE
}
}
# Loop 2.2.2.3.A
if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
#table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
if (dim(table.sp)[1]>1 && sum(!is.na(grep(Abbrev, table.sp$Infraspecific.rank, ignore.case=TRUE)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.rank==Abbrev,]
}
} else
# Loop 2.2.2.3.B
if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
# Loop 2.2.2.3.B.1
if((table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet))==FALSE && length(grep(Species, table.sp$Infraspecific.epithet))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
} else
# Loop 2.2.2.3.B.2
if((table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet))==FALSE && length(grep(Species, table.sp$Infraspecific.epithet))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==Species, ]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
}
}
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]
# Loop 2.2.2.3.1
if (z == 0) {
Family <- ""
Plant.Name.Index <- FALSE
Taxonomic.status <- ""
New.Genus <- Genus
New.Species <- Species
New.Infraspecific <- Infraspecific
Authority <- ""
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- ""
New.ID <- ""
} else
# Loop 2.2.2.3.2 (Search for accepted names)
if (any(table.sp$Taxonomic.status.in.TPL=="Accepted")) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted", ]
Plant.Name.Index <- TRUE
Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species[1]
if (infra == T && length(grep(infrasp, table.sp$Infraspecific.epithet))>0) {
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
} else
if (infra == F || length(grep(infrasp, table.sp$Infraspecific.epithet))==0) {
New.Infraspecific <- ""
}
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- table.sp[1,1]
New.ID <- ID
} else
# Loop 2.2.2.3.3 (Search for synonyms or missapplied names)
if (sum(table.sp$Taxonomic.status.in.TPL=="Accepted")==0 && (sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0||sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0)) {
if(sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Synonym", ]
} else
if(sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Misapplied", ]
}
if(sum(table.sp$Confidence.level=="H")>0) {
table.sp <- table.sp[table.sp$Confidence.level=="H", ]
} else
if(sum(table.sp$Confidence.level=="M")>0) {
table.sp <- table.sp[table.sp$Confidence.level=="M", ]
}
if(nrow(table.sp)>1) {
warning(paste(sp, "has more than one valid synonym"))
}
table.sp.id <- table.sp[1,1]
ID <- table.sp.id
at <- readLines(paste("http://www.theplantlist.org/tpl", vv, "/record/", table.sp.id, sep=""), encoding=encoding)
if (sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0) {
if (version=="1.1") {
az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
} else
if (version=="1.0") {
az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
}
} else 
if (sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0) {
if (version=="1.1") {
az <- "<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to"
} else 
if (version=="1.0") {
az <- "<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to"
}
}
n <- pmatch(az, at)
nsen <- at[n]
nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
if (version=="1.1") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[13],"+",ifelse(nsen[17]=="\u00D7", nsen[21], nsen[17]), "&csv=true", sep="")
} else
if (version=="1.0") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[11],"+",ifelse(nsen[15]=="\u00D7", nsen[19], nsen[15]), "&csv=true", sep="")
}
kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
if(infra==T && kup > 0) {
if (version=="1.1") {
infrasp <- ifelse(nsen[17]=="\u00D7", nsen[29], nsen[25])
} else
if (version=="1.0") {
infrasp <- ifelse(nsen[15]=="\u00D7", nsen[27], nsen[23])
}
} else
if(kup == 0) {
infrasp <- ""
}
if (sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0) {
Taxonomic.status <- "Synonym"
} else 
if (sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0) {
Taxonomic.status <- "Misapplied"
}
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character", as.is = TRUE, encoding=encoding),silent=TRUE)
if (version=="1.1" && nsen[17]=="\u00D7") {
table.sp <- table.sp[table.sp$Species.hybrid.marker=="\u00D7", ]
} else
if (version=="1.0" && nsen[15]=="\u00D7") {
table.sp <- table.sp[table.sp$Species.hybrid.marker=="\u00D7", ]
}
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
# Loop 2.2.2.3.3.A
if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
} else
# Loop 2.2.2.3.3.B
if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
}
Plant.Name.Index <- TRUE
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species[1]
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
New.ID <- table.sp[1,1]
} else
# Loop 2.2.2.3.4 (Search for unresolved names)
if (sum(table.sp$Taxonomic.status.in.TPL=="Accepted")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Synonym")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Unresolved")>0) {
Plant.Name.Index <- TRUE
Taxonomic.status <- "Unresolved"
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species[1]
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- table.sp[1,1]
New.ID <- ID
}}} else
# Loop 2.2.3
if(z == 1) {
# Loop 2.2.3.1
if (is.na(table.sp$Taxonomic.status.in.TPL)) {
Taxonomic.status <- ""
Plant.Name.Index <- FALSE
Family <- ""
New.Genus <- Genus
New.Species <- Species
New.Infraspecific <- Infraspecific
Authority <- ""
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- ""
New.ID <- ""
} else
# Loop 2.2.3.2
if (table.sp$Taxonomic.status.in.TPL=="Synonym"||table.sp$Taxonomic.status.in.TPL=="Misapplied") {
ID <- table.sp[1,1]
at <- readLines(paste("http://www.theplantlist.org/tpl", vv, "/search?q=", genus,"+", species, sep=""), encoding=encoding)
if (table.sp$Taxonomic.status.in.TPL=="Synonym") {
if (version=="1.1") {
az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
} else
if (version=="1.0") {
az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
}
} else 
if (table.sp$Taxonomic.status.in.TPL=="Misapplied") {
if (version=="1.1") {
az <- "<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to"
} else 
if (version=="1.0") {
az <- "<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to"
}
}
n <- pmatch(az, at)
nsen <- at[n]
nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
if (version=="1.1") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[13],"+",ifelse(nsen[17]=="\u00D7", nsen[21], nsen[17]), "&csv=true", sep="")
} else
if (version=="1.0") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[11],"+",ifelse(nsen[15]=="\u00D7", nsen[19], nsen[15]), "&csv=true", sep="")
}
kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
if(infra==T && kup > 0) {
if (version=="1.1") {
infrasp <- ifelse(nsen[17]=="\u00D7", nsen[29], nsen[25])
} else
if (version=="1.0") {
infrasp <- ifelse(nsen[15]=="\u00D7", nsen[27], nsen[23])
}
} else
if(kup == 0) {
infrasp <- ""
}
if (table.sp$Taxonomic.status.in.TPL=="Synonym") {
Taxonomic.status <- "Synonym"
} else 
if (table.sp$Taxonomic.status.in.TPL=="Misapplied") {
Taxonomic.status <- "Misapplied"
}
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character", encoding=encoding),silent=TRUE)
if (version=="1.1" && nsen[17]=="\u00D7") {
table.sp <- table.sp[table.sp$Species.hybrid.marker=="\u00D7", ]
New.Hybrid.marker <- "\u00D7"
} else
if (version=="1.0" && nsen[15]=="\u00D7") {
table.sp <- table.sp[table.sp$Species.hybrid.marker=="\u00D7", ]
New.Hybrid.marker <- "\u00D7"
}
Plant.Name.Index <- TRUE
Family <- table.sp$Family[1]
if (version=="1.1") {
New.Genus <- nsen[13]
New.Species <- ifelse(nsen[17]=="\u00D7", nsen[21], nsen[17])
Authority <- ifelse(nsen[17]=="\u00D7", nsen[25], nsen[21])
} else
if (version=="1.0") {
New.Genus <- nsen[11]
New.Species <- ifelse(nsen[15]=="\u00D7", nsen[19], nsen[15])
Authority <- ifelse(nsen[15]=="\u00D7", nsen[23], nsen[19])
}
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
if (length(ngrep) == 0 && corr==TRUE && !is.na(infrasp) && nchar(infrasp)>0) {
cck.infra <- agrep(infrasp, table.sp$Infraspecific.epithet, value=TRUE, max.distance=max.distance)
ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
if(length(cck.infra) > 0) {
cck.infra <- cck.infra[ddf.infra==min(ddf.infra)]
ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
}
levs <- length(unique(cck.infra))
if(length(cck.infra) > 0 && ddf.infra <= diffchar && levs == 1) {
infrasp <- cck.infra
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
marker.infra <- TRUE
}
}
# Loop 2.2.2.3.3.A
if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
if (dim(table.sp)[1]>1 && sum(!is.na(grep(Abbrev, table.sp$Infraspecific.rank, ignore.case=TRUE)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.rank==Abbrev,]
}
} else
# Loop 2.2.2.3.3.B
if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
}
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
New.ID <- table.sp[1,1]
} else
# Loop 2.2.3.3
if (table.sp$Taxonomic.status.in.TPL=="Accepted"||table.sp$Taxonomic.status.in.TPL=="Unresolved") {
Plant.Name.Index <- TRUE
Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
Family <- table.sp$Family
New.Genus <- table.sp$Genus
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species
New.Infraspecific <- table.sp$Infraspecific.epithet
Authority <- table.sp$Authorship
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- table.sp[1,1]
New.ID <- ID
}
} # End of loop 2.2.3
} # End of loop 2.1
} # End of loop 2
results <- data.frame(Genus, Species, Abbrev, Infraspecific, ID, Plant.Name.Index, TPL_version=version, Taxonomic.status, Family, New.Genus, New.Hybrid.marker, New.Species, New.Infraspecific, Authority, New.ID, Typo, WFormat, stringsAsFactors = FALSE)
}
