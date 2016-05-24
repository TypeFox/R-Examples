pedToHaplin <-
  function(indata, outdata, merge = F, na.strings = "0", sep, colnames.out = F)
{
  ## READS AND CONVERTS DATA FROM PEDIGREE FORMAT (LINKAGE FORMAT, WITH 6 LEADING COLUMNS) TO TRIAD FORMAT.
  ## INDATA IS THE FILE CONTAINING PEDIGREE DATA. OUTDATA IS THE 
  ## ASCII FILE OF THE CONVERTED DATA.
                                        #
  ## NOTE THAT GENETIC DATA (EVERYTHING AFTER COLUMN 6) COULD HAVE ALLELES IN
  ## SEPARATE COLUMNS OR IN SAME COLUMN, ALWAYS TREATED THE SAME
                                        #
                                        # (mapfile-ARGUMENT NOT IMPLEMENTED)
                                        #
  ##
#
## INFO TO USER
cat("#### CONVERTING PED FILE TO HAPLIN FORMAT ####\n")
cat("\nNote: First 6 columns should always be:\nfamily id, individual id, father's id, mother's id, sex and casetype\nThere should be no column names in input file.\n") 
#
if(missing(indata)) stop('File name argument "indata" must be specified!')
if(missing(outdata)) stop('File name argument "outdata" must be specified!')
#
## READ DATA
cat("\nReading data...")
if(!missing(sep)){
	.indata <- read.table(file = indata, sep = sep, as.is = T, na.strings = na.strings, colClasses = "character", stringsAsFactors = F)
}else{# TAKE ADVANTAGE OF read.table'S INTERPRETATION OF WHITE SPACE, ETC.
	.indata <- read.table(file = indata, as.is = T, na.strings = na.strings, colClasses = "character", stringsAsFactors = F)
}
cat("Done")
                                        #
  ##
  .nlines <- dim(.indata)[1]
  .nvars <- dim(.indata)[2]
  .nmarkers <- .nvars - 6
  if(.nmarkers <= 0) stop("Not enough columns in file! Did you specify the correct separator?")
  .vardata <- .indata[, 1:6]
  .gendata <- as.matrix(.indata[, -(1:6), drop = F])
                                        #
  if(merge){## GENETIC DATA MUST BE PASTED WITH A ";" TO CONFORM
    cat("\nPasting data columns pairwise...")
    if((.nmarkers %% 2) != 0) stop("Number of columns is not even")
    ## PICK FIRST AND SECOND ALLELES FOR EACH MARKER
    .mat1 <- .gendata[, seq(1, .nmarkers, by = 2)]
    .mat2 <- .gendata[, seq(2, .nmarkers, by = 2)]
    .mat1 <- as.character(.mat1)
    .mat2 <- as.character(.mat2)
    ## PASTE FIRST AND SECOND
    .gendata <- paste(.mat1, .mat2, sep = ";")
    .gendata <- matrix(.gendata, ncol = .nmarkers/2)
    .nmarkers <- .nmarkers/2
    cat("Done")
  }
                                        #
  ## NAME COLUMNS
  if(F && !missing(mapfile)){
	mapfile <- NULL # BARE TULL, MEN MAA UNNGAA FEILMELDING I R CMD check
    .map <- read.table(mapfile, stringsAsFactors = F)
    names(.map)[c(1,2,4)] <- c("chrom", "marker", "pos")
    dimnames(.gendata)[[2]] <- .map$marker
  } else dimnames(.gendata)[[2]] <- paste("mark", 1:.nmarkers, sep = "")
                                        #
  ## MAKE SURE ALL MISSING ARE NA, NOT "NA;NA"
  .gendata[.gendata == "NA;NA"] <- NA
                                        #
  ##
  names(.vardata) <- c("family", "id", "father", "mother", "sex", "casetype")
                                        #
  ## CHECK family VARIABLE
  cat("\nCheck family and id variables...")
  if(any(table(.vardata$family, exclude = NULL) > 3)) cat("(Found family size larger than 3!  Will extract trios from general pedigree.)")
  if(any(.vardata$id == "0")) stop('Cannot use "0" in id code!')
  if(any(is.na(.vardata$id))) stop('id variable cannot contain missing')
  if(any(is.na(.vardata$family))) stop('family variable cannot contain missing')
  .famcodes <- unique(.vardata$family)
                                        #.nfams <- length(.famcodes)
                                        #
  ## CREATE A VARIABLE THAT UNIQUELY IDENTIFIES INDIVIDUALS
  .tagit <- "-xXx-" # I BET THAT ONE'S UNIQUE!
  .tag <- paste(.vardata$family, .vardata$id, sep = .tagit)
  .dupl <- duplicated(.tag)
  if(any(.dupl)) {
	cat("\n")
	.mess <- paste("Individual id appears several times within one family!\nFor instance, family ", .vardata$family[.dupl][1], " contains more than one of individual ", .vardata$id[.dupl][1], ".", sep = "")	
	stop(.mess)
  }
  ## CREATE THE SAME FOR MOTHERS AND FATHERS
  .mother.tag <- paste(.vardata$family, .vardata$mother, sep = .tagit)
  .father.tag <- paste(.vardata$family, .vardata$father, sep = .tagit)
                                        #
  ## IDENTIFY MOTHERS, FATHERS AND CHILDREN
  .mother.pos <- is.element(.tag, .mother.tag)
  .father.pos <- is.element(.tag, .father.tag)

  if (FALSE) { ### PREVIOUS CODE ###
    .child.pos <- !.mother.pos & !.father.pos
  }
  else { ### NEW CODE ###
    ## 1) Child = person having a mother or father
    .child.pos.1 <- !is.na(.vardata$mother) | !is.na(.vardata$father)
    ## 2)    or = person not having mother and father and does not have any children
    .child.pos.2 <- !.child.pos.1 & (!.mother.pos & !.father.pos)

    .child.pos <- .child.pos.1 | .child.pos.2
  }
  

  ## SOME CHECKING
  if(any(.father.pos & .mother.pos)){
    stop("Sorry, father and mother cannot have the same id!")
  }
  if(F){
    if(any(!is.element(.mother.tag[!is.na(.vardata$mother)], .tag)))stop("Mother not found in file")
    if(any(!is.element(.father.tag[!is.na(.vardata$father)], .tag)))stop("Father not found in file")
  }# RAPPORTER DISSE SOM MULIGE FEIL!?
  ## ALLE DISSE KUNNE OPPGI FAMILIER SPESIFIKT!



  if (FALSE) { ### PREVIOUS CODE ###

    
    ## MAKES A NEW id-VARIABLE, WITH CODES 1, 2 OG 3 FOR MOTHER, FATHER, CHILD, RESP.
    cat("Done")
    cat("\nExpand data to complete families...")
    .id.ny <- (.mother.pos * 1) + (.father.pos * 2) + (.child.pos * 3)
    .tag.ny <- paste(.vardata$family, .id.ny, sep = .tagit)
    if(any(duplicated(.tag.ny))) stop("Problem!")

    ## CREATE COMPLETE FILE WITH ALL MOTHERS, FATHERS, CHILDREN
    .full <- expand.grid(ind = 1:3, fam = .famcodes)
    .full$fam <- as.character(.full$fam) # DUSTEGREIE
    .full.tag <- paste(.full$fam, .full$ind, sep = .tagit)
                                        #
    ## MATCH OBSERVED DATA TO FULL FILE
    .match <- match(.full.tag, .tag.ny)
    .utdata <- cbind(as.matrix(.vardata), .gendata)
    .utdata <- .utdata[.match,]

    ## IMPUTE FAMILY ID FOR THE MISSING
    .utdata[is.na(.match), "family"] <- .full[is.na(.match), "fam"] 
  }
  else { ### NEW CODE ###
    ## Find index position of child and of their parents
    .index.child <- which(.child.pos)
    .index.mother.of.child <- match(.mother.tag, .tag)[.child.pos]
    .index.father.of.child <- match(.father.tag, .tag)[.child.pos]

    ## Exclude trios with only missing genotypes
    .na.geno <- apply(is.na(.gendata), 1, prod)
    .na.geno.child <- .na.geno[.index.child]
    .na.geno.mother <- .na.geno[.index.mother.of.child]
    .na.geno.father <- .na.geno[.index.father.of.child]
    .na.geno.mother[is.na(.na.geno.mother)] <- 1 ## Set missing mothers = 0
    .na.geno.father[is.na(.na.geno.father)] <- 1 ## Set missing fathers = 0
    .na.all <- .na.geno.child & .na.geno.mother & .na.geno.father

    .index.child <- .index.child[!.na.all] # exclude missing
    .index.mother.of.child <- .index.mother.of.child[!.na.all]
    .index.father.of.child <- .index.father.of.child[!.na.all]
    
    .index <- as.vector(rbind(.index.mother.of.child,
                              .index.father.of.child,
                              .index.child))

    .utdata <- cbind(as.matrix(.vardata), .gendata)
    .family <- rep(.utdata[.index.child, "family"], each=3)
    
    .utdata <- .utdata[.index,]

    ## IMPUTE FAMILY ID FOR THE MISSING
    .utdata[is.na(.index), "family"] <- .family[is.na(.index)]
  }
  
  
                                        #
### REMOVE REDUNDANT COLS (KEEP casetype BY DEFAULT)
  ## vars SPECIFIES WHAT VARIABLES (COVARIATES) TO INCLUDE AFTER CONVERSION
  vars <- c("family", "sex", "casetype") # KEEP BY DEFAULT
  .keep <- match(vars, c("family", "id", "father", "mother", "sex", "casetype"))
  .rem <- setdiff(1:6, .keep)
  .utdata <- .utdata[, -.rem]
  cat("Done")
                                        #
  ## RESHAPE DATA
  cat("\nReshaping data...")
  .navn <- dimnames(.utdata)[[2]]
  .navn <- as.vector(t(outer(.navn, c("m", "f", "c"), paste, sep = ".")))
  .utdata <- lapply(seq(dim(.utdata)[2]), function(i) .utdata[, i])# CONV. TO LIST
  ###.utdata <- f.matrix.to.list(.utdata)
  .utdata <- lapply(.utdata, function(x) matrix(x, ncol = 3, byrow = T))
  .utdata <- do.call("cbind", .utdata)
  dimnames(.utdata) <- list(NULL, .navn)
  cat("Done")
                                        #
  ##
  cat("\nRemove redundant columns...")
  if(is.element("family", vars)){
    ## CHECK THAT FAMILY ID'S ARE THE SAME FOR MOTHERS, FATHERS AND CHILDREN (SHOULD BE!)
    if(any(.utdata[, "family.m"] != .utdata[, "family.f"]) | any(.utdata[, "family.m"] != .utdata[, "family.c"])) stop("Something's wrong with the family variable")
    ## KEEP ONLY ONE COLUMN
    .navn <- dimnames(.utdata)[[2]]
    .rem <- match(c("family.f", "family.c"), .navn)
    dimnames(.utdata)[[2]][.navn == "family.m"] <- "family"
    .utdata <- .utdata[, -.rem]
  }
  if(is.element("sex", vars)){
    ## CHECK THAT MOTHER'S AND FATHER'S GENDER ARE CORRECT (ALLOW MISSING)
    if(any(!is.element(.utdata[, "sex.m"], c(NA, 2)))) stop("There are mothers with sex variable not equal to 2")
    if(any(!is.element(.utdata[, "sex.f"], c(NA, 1)))) stop("There are fathers with sex variable not equal to 1")
    ## KEEP ONLY CHILD
    .navn <- dimnames(.utdata)[[2]]
    .rem <- match(c("sex.m", "sex.f"), .navn)
    dimnames(.utdata)[[2]][.navn == "sex.c"] <- "sex"
    .utdata <- .utdata[, -.rem]
  }
  if(is.element("casetype", vars)){
    ## REMOVE CASETYPE FOR MOTHER AND FATHER
                                        #if(any(!is.element(.utdata[, "sex.m"], c(NA, 2)))) stop("There are mothers with sex variable not equal to 2")
                                        #if(any(!is.element(.utdata[, "sex.f"], c(NA, 1)))) stop("There are fathers with sex variable not equal to 1")
    ## KEEP ONLY CHILD
    .navn <- dimnames(.utdata)[[2]]
    .rem <- match(c("casetype.m", "casetype.f"), .navn)
    dimnames(.utdata)[[2]][.navn == "casetype.c"] <- "casetype"
    .utdata <- .utdata[, -.rem]
  }
  cat("Done")
                                        #
  ## WRITE DATA TO FILE, "STANDARD" HAPLIN FORMAT
  cat("\nWrite data to file...")
  write.table(.utdata, file = outdata, sep = "\t", append = FALSE, quote = FALSE,row.names = FALSE, col.names = colnames.out, na = "NA")
  cat("Done\n")
                                        #
  ##
  return(invisible(.utdata))
}
