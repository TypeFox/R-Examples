
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.bc <- function(nb)
#TITLE  sequence for a loop
#DESCRIPTION
# This function returns \code{1:nb} when \code{nb > 0} and
#         \code{numeric(0)} otherwise.\cr
# Quite useful to prevent starting
# a loop of length nought
#DETAILS
#PKEYWORDS helpful
#KEYWORDS iteration
#INPUTS
#{nb}    <<length of the loop>>
#[INPUTS]
#VALUE
# \code{1:nb} if \code{nb > 0}
# else \code{numeric(0)}.
#EXAMPLE
# r.bc(0);
# r.bc(5);
#REFERENCE
#SEE ALSO r.bf
#CALLING
#COMMENT
# r.bc for BouCle (loop in French)
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_18
#REVISED 09_04_05
#--------------------------------------------
{
if (is.null(nb)) {return(numeric(0));}
if (length(nb)!=1) {
    r.erreur(nb,"r.bc deals only with scalar nb");
}
if (nb > 0) {return(1:max(1,round(nb)));}
numeric(0);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.bd <- function(n1,n2)
#TITLE  sequence for insertions
#DESCRIPTION
# This function returns \code{n1:n2} when \code{n1<=n2} and
#         \code{numeric(0)} otherwise.\cr
# Quite useful when some insertion must be done within
# a sequence
#DETAILS
#PKEYWORDS helpful
#KEYWORDS iteration
#INPUTS
#{n1}    <<first element>>
#{n2}    <<second element>>
#[INPUTS]
#VALUE
# \code{n1:n2} if \code{n1<n2}
# else \code{numeric(0)}.
#EXAMPLE
# xx <- 1:5;
# for (ii in 1:6) { print(c(xx[r.bd(1,ii-1)],10,xx[r.bd(ii,5)]));}
#REFERENCE
#SEE ALSO r.bf
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_01_12
#REVISED 11_01_12
#--------------------------------------------
{
if (n1 <= n2) {return(n1:n2);}
numeric(0);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.bf <- function(x)
#TITLE  sequence for looping over an object
#DESCRIPTION
# This function returns \code{1:length(x)} when \code{length(x) > 0} and
#         \code{numeric(0)} otherwise.
# Quite useful to prevent starting
# a loop of length nought
#DETAILS
#PKEYWORDS helpful
#KEYWORDS iteration
#INPUTS
#{x}    <<vector>>
#[INPUTS]
#VALUE
# \code{1:length(x)} if \code{length(x) > 0}
# else \code{numeric(0)}.
#EXAMPLE
# r.bf(0);
# r.bf(5);
# r.bf(character(0));
# r.bf(letters);
#REFERENCE
#SEE ALSO r.bc
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 08_09_05
#REVISED 08_09_05
#--------------------------------------------
{
if (length(x) > 0) { return(1:length(x));}
numeric(0);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.erreur <- function(x,...,w=FALSE)
#TITLE  issues an error message and concludes accordingly
#DESCRIPTION
# when called this function prints x, then displays a message before stopping 
# the process.
#DETAILS
#PKEYWORDS
#KEYWORDS error
#INPUTS
#{x} <<object to be printed before the message. When \code{is.null(x)}
#      nothing is printed. When it is a list, all components
#      of the list are successively printed.>>
#{\dots}<<pieces of message to display after pasting>>
#[INPUTS]
#{w} << Indicates if it is a simple warning>> 
#VALUE
# nothing is returned
#EXAMPLE
# r.erreur(matrix(1:4,2),"This matrix is not symmetric",w=TRUE)
# r.erreur("Are you sure of it?",w=TRUE);
#
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_03
#REVISED 13_04_29
#--------------------------------------------
{
r.form3repeat("~",60,TRUE);
print(x);
message <- paste(...);
cat("<<<<< MESSAGE >>>>>\n");
print(message);
if (w) {
    cat("SIMPLE WARNING:\n");
} else {
    on.exit(traceback());
}
r.form3repeat("~",60,TRUE);
invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.check4tyle <- function(x,typ,len=-1,con=numeric(0),
                       message=NULL,fatal=TRUE,na.rm=FALSE)
#TITLE  checks the type, the length and the content
# of some standard object
#DESCRIPTION
#  checks the type, the length and the content
# of some standard object.
# If not correct, a message and a fatal error are issued.
# NA are detected and considered as wrong.
#DETAILS
# 'integer' has not got the meaning in \code{is.integer} R
# standard function. 'null' must be understood as
# resulting TRUE with function \code{is.null}.\cr
# 'Named character' is different from  
# 'character': \code{rbsb.chara} can be used for this purpose.
#PKEYWORDS
#KEYWORDS error
#INPUTS
#{x} <<object to be checked.>>
#{typ} <<The list of correct types, among
# 'null', 'integer', 'numeric', 'character',
# 'logical', 'list', 'any', 'function', 'data.frame',
# 'matrix'. Also are 'nlogical', 'ninteger', 'nnumeric' and 
# 'ncharacter' for named structures.
# As understood, 'any' implies that no
# check of the type is performed.>>
#[INPUTS]
#{len} <<If \code{length(len)==1}, the exact length
# to be checked, if not must be of length two for
# the possible range of \code{length(x)}. When -1,
# no check on the length is performed.\cr
# For data.frame, it is the number of columns.
# When \code{na.rm} the check is made after
# removing the NA values.>>
#{con} << If \code{length(con)>0}, some check about
# the content of \code{x} is done for some of the 
# types. More precisely for (integer, numeric): all
# values must belong to the interval \code{[con[1],con[2]]}
# and for (character), the set of possible \code{character(1)}
# is described by \code{con}.>>
#{message} << Some additional message to be
#            issued before stopping in case of error.>>
#{fatal} << what to do when discovering
#           an inconsistency ? TRUE: this function prints the message
#           and stops; FALSE: this function returns
#           the message as a character.>>
#{na.rm} << Must the \code{NA} values be removed before checking?
#           This is implemented only for the types integer, numeric,
#           character and logical.>>
#VALUE
# a character with the message(s), when everything is OK
# returns a \code{character(0)}.
# But when \code{fatal} if not OK prints the
# message and stops.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE improves the treatment of 'NA's
#AUTHOR J.-B. Denis
#CREATED 09_05_02
#REVISED 13_04_29
#--------------------------------------------
{
  # checking the arguments and preparation
  # accepted types by r.check4tyle
  # ....
  invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.form3justifie <- function(chaine,
                          nbc=8,
                          format=3,
                          tronc=TRUE,
                          carac=" ")
#TITLE  formats a character string
#DESCRIPTION
# Formats character string(s).
# The main use of this function is to produce
# aligned columns lists while printing
# the rows not at the same time.
#DETAILS
#PKEYWORDS
#KEYWORDS print
#INPUTS
#{chaine}<<the character string to be printed, can be a vector.>>
#[INPUTS]
#{nbc} << Desired number of characters for the result; when
#         \code{chain} is a vector can be a vector of the same length>>
#{format} << Indicates the type of alignment:\cr
#   0 no aligment (no supplementary character added to reach \code{nbc})\cr
#   1 to the left side\cr
#   2 centered\cr
#   3 to the right side>>
#{tronc} << If true, no more than
#     \code{nbc} characters are returned and
# possibly the string is truncated. In that
# case, \code{$} is added for indicating the fact.>>
#{carac} << Character to use for enlarging the string>>
#VALUE a character string
#EXAMPLE
# r.form3justifie("vers")
# r.form3justifie("versification",5)
# r.form3justifie(letters[1:5],format=2,carac="+");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 1999_05_25
#REVISED   10_07_22
#--------------------------------------------
{
# checking
# the null case
if (length(chaine)==0) { return(character(0));}
# preparing
nbc[nbc < 3] <- 8;
if (length(nbc) < length(chaine)) {
    nnbc <- rep(nbc,length(chaine));
} else {
    nnbc <- nbc;
}
#
itr <- "$"; # truncation indicator
rres <- cchaine <- chaine;
for (rr in r.bf(rres)) {
    res <- rres[rr];
    nbc <- nnbc[rr];
    chaine <- cchaine[rr];
    # truncation
    if ( (nchar(res) > nbc) & tronc ) {
     if (format <= 1) {
      res <- substring(chaine,1,nbc);
      res <- paste(res,itr,sep="");
      }
     else {
      if (format == 2) {
       otg <- (nchar(chaine) - nbc) %/% 2;
       res <- substring(chaine,1+otg);
       res <- substring(res,1,nbc);
       res <- paste(itr,res,itr,sep="");
       }
      else {
       res <- substring(chaine,1+nchar(chaine)-nbc,
			nchar(chaine));
       res <- paste(itr,res,sep="");
       }
      }
     }
    if ((nchar(res) < nbc) & (format != 0)) {
     if (format == 1) {
      while (nchar(res) < nbc) res <-
	     paste(res,"",collapse="",sep=carac);
      }
     else {
      if (format == 2) {
       raj <- (nbc - nchar(res)) %/% 2;
       if (raj > 0) {
	for (jbd in 1:raj) res <-
	 paste(res,"",collapse="",sep=carac);
	}
       while (nchar(res) < nbc) {
	res <- paste("",res,collapse="",sep=carac);
	}
       }
      else {
       while (nchar(res) < nbc) {
	res <- paste("",res,collapse="",sep=carac);
	}
       }
      }
     }
    rres[rr] <- res;
}
# returning
rres;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.form3names <- function(nbn,nom=character(0),prefix="",
                           upca=TRUE,nume=14)
#TITLE  provides systematic names for items
#DESCRIPTION
# Provides systematic names for a series of items according their 
# number taking care of previous names.
#DETAILS (see the code)
#PKEYWORDS
#KEYWORDS utilities
#INPUTS
#{nbn} <<Number of new item names to generate>>
#[INPUTS]
#{nom} << Already present names (to avoid identical names).>>
#{prefix} << Systematic prefix of all names to generate. Must
#                 comprise the dot, if one wants such a separator
#                 between it and the specific part of the name. 
#                 Of course can be 'underscore' or whatever else.>>
#{upca} << Indicates whether the letters constituting the new
#          names must be uppercase or not.>>
#{nume} << Its absolute value gives the number of the letter to use
#          when the alphabet is not sufficient. When negative, alphabet
#          is not considered as a first possibility.>>
#VALUE
# vector with \code{nbn} different strings associated to new names
#EXAMPLE
# r.form3names(2);
# r.form3names(2,nume=-3);
# r.form3names(2,prefix="rbsb.");
# r.form3names(2,upca=FALSE);
# r.form3names(5,"D");
# r.form3names(5,"Y");
# r.form3names(30);
#REFERENCE
#SEE ALSO r.form3nume
#CALLING
#COMMENT
#FUTURE
# Monitor the number of digits to obtain series of "Z01", "Z02",... ,"Z79", "Z80".
#AUTHOR J.-B. Denis
#CREATED 07_10_19
#REVISED 10_02_15
#--------------------------------------------
{
# checking
#
if (upca) { Letters <- LETTERS;
} else { Letters <- letters;}
#
if (is.null(nbn)) { return(character(0));}
if (nbn < 1) { return(character(0));
} else {
    if (prefix != "") {
        # keeping only the names having got the prefix then
        # removing the prefixes from them (to be added further)
        decom <- sapply(strsplit(nom,prefix),function(ll){length(ll);});
        nom <- nom[decom == 2];
        nom <- sapply(strsplit(nom,prefix),function(ll){ll[2];});
    }
    # looking for the maximum letter in noms
    if ( length(nom) == 0 ) { mama <- 0;
    } else { mama <- max(1*outer(nom,Letters,"==") %*% matrix(1:26,ncol=1));}
    if ((nbn < (27-mama)) & (nume>0)) {
        # adding letters 
        res <- Letters[mama+(1:nbn)];
    } else {
        # adding numbered nodes
        ajou <- 0; nu <- 1; res <- character(0);
        while ( ajou < nbn ) {
          nono <- paste(Letters[abs(nume)],nu,sep="");
          if (all(nono != nom)) {
              ajou <- ajou + 1;
              res <- c(res,nono);
          }
          nu <- nu+1;
        }
    }
}
# adding the prefix
res <- paste(prefix,res,sep="");
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.numero <- function(qui,dans,strict=TRUE,repe=TRUE)
#TITLE  returns the positions of identifiers
#DESCRIPTION
#  Returns the positions of identifiers.
# From a reference set of
# identifiers returns the positions of some of them.
# If they don't belong to the list, an error is issued
# according to \code{strict}. If there are repetitions
# an error is issued or not according to \code{repe}.
#DETAILS
#PKEYWORDS helpful
#KEYWORDS
#INPUTS
#{qui}    <<Identifiers the position of them is looked for.
# Must be a character or a numeric.>>
#{dans}    <<The ordered vector of all identifiers.
# Must be a character or a numeric.>>
#[INPUTS]
#{strict} <<Must \code{qui} be a subset of \code{dans}.>>
#{repe} <<Are repeatitions allowed?>>
#VALUE
# A numeric vector of the positions. It is named with the identifiers.
# When the identifier does not belong to the list and \code{strict} is
# FALSE, NA values are returned.
#EXAMPLE
# r.numero(10:1,1:10);
# r.numero(c("Z","E","P","L","I","N"),LETTERS);
# r.numero(c("B","o","F","!"),LETTERS,strict=FALSE);
# r.numero(c("B","E","B","E"),LETTERS);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 10_03_09
#REVISED 10_10_12
#--------------------------------------------
{
# checking
qui  <- as.character(qui);
dans <- as.character(dans);
if (strict) {
    uu <- union(qui,dans);
    if (length(uu)!=length(dans)) {
        r.erreur(list(setdiff(uu,dans),dans),"Some identifiers does not belong to the list");
    }
}
# computing
res <- as.numeric(outer(qui,dans,"==")%*%r.bf(dans));
names(res) <- qui;
res[res==0] <- NA;
# repetitions
if (!repe) {
    rr <- res[!is.na(res)];
    if (length(unique(rr)) < length(rr)) {
        r.erreur(rr,"Some repetitions were found while forbidden!");
    }
}  
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.expr3present <- function(sch,ch,exact=FALSE,how="a")
#TITLE  indicates inclusion of character string
#DESCRIPTION
# Checks if some string(s) belong(s) to a series of strings.
# When \code{exact} is FALSE, returns TRUE if the 
# character string \code{sch} is included at least once
# into the character string \code{ch}.\cr
# \code{sch} can be a vector, in that case, the check is made for each
# of its components. According to \code{how} the returned
# value is vectorial or scalar. When \code{sch} is zero length then
# \code{TRUE} is returned but if not and \code{ch} is zero length
# then \code{FALSE} is returned.
#DETAILS
#PKEYWORDS
#KEYWORDS utilities
#INPUTS
#{sch} <<(\code{character}) the character string(s) to be found.>>
#{ch}  <<(\code{character}) the character string(s) to investigate.>>
#[INPUTS]
#{exact} << When exact, one component must
# be strictly identical, if not a subtring is sufficient.>>
#{how} << Indicates what to do when \code{length(sch)>1}. The choice are 
# \code{'v'}: a logical vector gives back each check independently;
# \code{'1'}: returns \code{TRUE} when at least one of the component belongs
# to the series \code{ch} and \code{'a'} when all components must comply to get TRUE.>>
#VALUE A logical vector with length of \code{sch}.
#EXAMPLE
# r.expr3present('a','non');
# r.expr3present('o',c('non','oui'));
# r.expr3present(c("o","oui"),c('non','oui'));
# r.expr3present(c("o","oui"),c('non','oui'),how="v");
# r.expr3present(c("A[SEX]","A[AGE]"),c("A[AGE]","A[SEX]"),how="a")
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 08_08_28
#REVISED 12_05_16
#--------------------------------------------
{
#
# degenerate cases
if (length(sch)==0) { return(TRUE);}
if (length( ch)==0) { return(FALSE);}
#
#
res <- logical(length(sch));
#
for (ii in r.bf(sch)) {
    if (exact) {
        res[ii] <- sum(sch[ii]==ch) > 0;
    } else {
        res[ii] <- length(grep(sch[ii],ch)) > 0;
    }
}
# dealing with scalar cases
if (how == "1") { res <- any(res);}
if (how == "a") { res <- all(res);}
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
r.form3repeat <- function(cha="-",nb=10,imp=FALSE,cr=imp)
#TITLE  prints a repeated given string
#DESCRIPTION
# Without adding breaking line characters, prints
# \code{nb} times a given string of characters
#DETAILS
#PKEYWORDS
#KEYWORDS print
#INPUTS
#[INPUTS]
#{cha} << The string to repeat>> 
#{nb} << Number of repetitions>> 
#{imp} << Printing when TRUE or returning (default)>>
#{cr} << Must a line feed be added?>>
#VALUE
# character string or printing according to \code{imp}
#EXAMPLE
# r.form3repeat('-+',20,TRUE)
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_24
#REVISED 08_09_25
#--------------------------------------------
{
nb <- max(0,round(nb));
res <- "";
for (jbd in r.bc(nb)) { res <- paste(res,cha,sep="");}
if (cr) { res <- paste(res,"\n",sep="");}
if(!imp) { return(res);}
cat(res);
invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
