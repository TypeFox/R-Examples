f.make.index <- function(vardata, output = "line numbers"){
##
## vardata ARE THE 4 FIRST COLUMNS OF A STANDARD PED FILE.
## FROM vardata AND INDEX MATRIX IS CONSTRUCTED, WITH ONE 
## ROW FOR EACH FAMILY IDENTIFIED IN vardata.
## COLUMNS IN OUTPUT MATRIX REPRESENT INDIVIDUAL ID'S FOR MOTHER, FATHER AND CHILD
## IF standard.output = T, THE COLUMNS GIVE ROW NUMBERS RATHER THAN ID VALUES
#									#

## HUSK: HVIS DENNE SKAL BRUKES HELLER ENN pedToHaplin ER DET ET PROBLEM MED
## KRAV OM UNIK ID!


## CHECK family VARIABLE
cat("\nChecking family and id variables...\n")
if(any(table(vardata$family, exclude = NULL) > 3)) warning("Found family size larger than 3!  Will extract trios from general pedigree.", call. = F)
if(.test <- (any(vardata$family == "0") | any(vardata$id == "0"))) stop(paste('Cannot use "0" in family or id code!\n', 'Found on line(s): ', paste(which(.test), collapse = " "), sep = ""), call. = F)
if(.test <- (any(is.na(vardata$family)) | any(is.na(vardata$id)))) stop(paste('Cannot have missing values in family or id variable\n', 'Found on line(s): ', paste(which(.test), collapse = " "), sep = ""), call. = F)
#
## RECODE ZEROES TO MISSING
vardata$mother[vardata$mother == "0"] <- NA
vardata$father[vardata$father == "0"] <- NA
#
## CREATE A VARIABLE THAT UNIQUELY IDENTIFIES INDIVIDUALS
.tagit <- "<>" # I BET THAT ONE'S UNIQUE!
.tag <- paste(vardata$family, vardata$id, sep = .tagit)
.tag.mother <- paste(vardata$family, vardata$mother, sep = .tagit)
.tag.father <- paste(vardata$family, vardata$father, sep = .tagit)
## RETAIN MISSING
.tag.mother[is.na(vardata$mother)] <- NA
.tag.father[is.na(vardata$father)] <- NA
#
## CHECK FOR DUPLICATES
.dupl <- duplicated(.tag)
if(any(.dupl)) {
	cat("\n")
	.mess <- paste("Individual id appears several times within one family!\nFor instance, family ", vardata$family[.dupl][1], " contains more than one of individual ", vardata$id[.dupl][1], ".", sep = "")
	stop(.mess, call. = F)
}
#
## IDENTIFY MOTHERS AND FATHERS
.is.mother <- is.element(.tag, .tag.mother)
.is.father <- is.element(.tag, .tag.father)
#
## ELIMINATE PARENT CODINGS THAT REFER TO NON-EXISTENT INDIVIDUALS,
## SAY, IF LINES OF THE FILE HAVE BEEN REMOVED
.rem.mother <- !is.element(.tag.mother, c(NA, .tag))
.rem.father <- !is.element(.tag.father, c(NA, .tag))
.sm <- sum(.rem.mother)
.sf <- sum(.rem.father)
if(.sm + .sf > 0){
	.tag.mother[.rem.mother] <- .tag.father[.rem.father] <- NA
	.mess <- paste(.sm, " mother code(s) and ", .sf, " father code(s) refer to non-existing individuals and have been set to missing", sep = "")
	warning(.mess, call. = F)
}
#
## DEFINE AND IDENTIFY "CHILD"
## 1) Child = person having a mother or father
.is.child.1 <- !is.na(.tag.mother) | !is.na(.tag.father)
## 2) or = person not having mother nor father but also does not have any children, i.e. "all alone"
.is.child.2 <- !.is.child.1 & (!.is.mother & !.is.father)
#
.is.child <- .is.child.1 | .is.child.2
#  
## SOME CHECKING
if(any(.test <- .is.father & .is.mother)){
	stop(paste("Sorry, the same individual cannot be both father and mother!\n", "Found on line(s): ", paste(which(.test), collapse = " "), sep = ""), call. = F)
}
#
## FIND LINE NUMBER OF CHILD AND OF ITS PARENTS
.line.child <- which(.is.child)
.line.mother <- match(.tag.mother, .tag)[.is.child]
.line.father <- match(.tag.father, .tag)[.is.child]
#
##
if(output == "line numbers"){
	## OUTPUT GIVES LINE NUMBERS, NOT IDs.
	## THIS COULD BE USEFUL FOR GenABEL SINCE FAMILY IDS ARE NOT RETAINED
	## (BUT SHOULD MAKE SURE SUBSETTING IN GenABEL IS TAKEN INTO ACCOUNT WHEN CONVERTING TO HAPLIN)
	#
	## JOIN LINE NUMBERS INTO MATRIX
	.ut <- cbind(line.child = .line.child, line.mother = .line.mother, line.father = .line.father)
}else if(output == "ids"){
	## USES SAME IDS AS IN FILE. REQUIRES IDS TO BE UNIQUE!
	## CAN BE USED WITH GenABEL SINCE THEY INSIST THAT IDS SHOULD BE UNIQUE.
	if(any(duplicated(vardata$id))) stop("Duplicated individual id", call. = F)
	.fam.child <- vardata$family[.line.child]
	.fam.mother <- vardata$family[.line.mother]
	.fam.father <- vardata$family[.line.father]
	#
	.id.child <- vardata$id[.line.child]
	.id.mother <- vardata$id[.line.mother]
	.id.father <- vardata$id[.line.father]
	#
	if(any(.fam.child != .fam.mother, na.rm = T) | any(.fam.child != .fam.father, na.rm = T)) stop("Problem with family identification!", call. = F)
	#
	## JOIN FAMILY ID, CHILD ID, MOTHER ID, FATHER ID INTO MATRIX
	.ut <- cbind(family = .fam.child, id.child = .id.child, id.mother = .id.mother, id.father = .id.father)
}else if(output == "tags"){
	## USES TAGS CREATED FROM FAMILY ID COMBINED WITH IDS AS IN FILE.
	## CAN BE USED IN GENERAL HAPLIN CONVERSIONS WHICH DO NOT REQUIRE ID ITSELF TO BE UNIQUE
	## BUT WILL NOT WORK WITH GenABEL SINCE LATTER DOES NOT RETAIN FAM IDS
	.fam.child <- vardata$family[.line.child]
	.fam.mother <- vardata$family[.line.mother]
	.fam.father <- vardata$family[.line.father]
	#
	.id.child <- .tag[.line.child]
	.id.mother <- .tag[.line.mother]
	.id.father <- .tag[.line.father]
	#
	if(any(.fam.child != .fam.mother, na.rm = T) | any(.fam.child != .fam.father, na.rm = T)) stop("Problem with family identification!", call. = F)
	#
	## JOIN FAMILY ID, CHILD ID, MOTHER ID, FATHER ID INTO MATRIX
	.ut <- cbind(family = .fam.child, id.child = .id.child, id.mother = .id.mother, id.father = .id.father)
}
#
##
return(.ut)
}
