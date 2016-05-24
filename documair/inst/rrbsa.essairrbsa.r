
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrplaces4text6tags <- function (text,tags,check=FALSE)
#TITLE returns the places of a series of tags within a /text/
#DESCRIPTION From a text, 
# independently for a series of tags, gives back the positions 
# where they are in the text. A position is given by four
# values \samp{(li1,co1,li2,co2)}, respectively the numbers of
# \samp{li}ne and \samp{co}lumn of the first character and last
# character of the tag. \cr
# Non-overlapping between tags can be checked.
#DETAILS
# One tags cannot be upon two successive lines. The same line
# can have more than one tag. This is a clear generalization of the 
# function \samp{rrrtext3places8word} which is called to elaborate the result.\cr
# If some tags have an intersection in the text, the case is considered 
# as no acceptable and a fatal error is issued. For instance in
# "The computation is performed" the two tags "perf" and "formed" are
# not compatible; that is\cr
#  \samp{rrrplaces4text6tags("is performed",c("perf","formed"))}
# stops with a fatal error.
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text.>>
#{tags}<< \samp{character()} the series of tags to be considered.>>
#[INPUTS]
#{check} << Must overlapping be checked?>>
#VALUE
# A five column data frame: a line for each found tag indicated by a factor;
# four columns giving respectively the \samp{line},
# the \samp{beg}inning and \samp{end}ing positions plus the fifth column
# to indicate the corresponding tag value. When they
# exist \samp{names(tags)} are used for the row dimnames.
#EXAMPLE
# rrrplaces4text6tags("Bonjour Monsieur","on");
# rrrplaces4text6tags("Bonjour Monsieur",c("on","ons","mon"));
# \dontrun{rrrplaces4text6tags("Bonjour Monsieur",c("on","ons","mon"),check=TRUE)};
# rrrplaces4text6tags(rrrrbsa0$text1$v,c("1","2","7"));
# rrrplaces4text6tags(rrrrbsa0$text2$v[1:3],"t");
# rrrplaces4text6tags(rrrrbsa0$text3$v,"uu");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_04
#REVISED 14_08_07
#--------------------------------------------
{
  # getting all positions of the tags
  posi <- matrix(NA,0,5);
  dimnames(posi)[[2]] <- c("tag","li1","co1","li2","co2");
  ftag <- character(0);
  for (ww in rrrbf(tags)) {
    rr <- rrrtext3places8word(text,tags[ww]);
    if (!all(rr[,1]==rr[,3])) {
      stop("Unexpected result found in 'rrrplaces4text6tags'");
    }
    posi <- rbind(posi,cbind(rep(ww,nrow(rr)),rr));
    fftag <- rep(names(tags[ww]),nrow(rr));
    # adding numbering when more than one occurrence
    if (length(fftag)>1) {
      fftag <- paste(rrrbf(fftag),fftag,sep=".");
    }
    ftag <- c(ftag,fftag);
  }
  # processing the results
  if (nrow(posi)>1) {
    # ordering the result
    ooo <- order(posi[,"li1",drop=FALSE],posi[,"co1",drop=FALSE]);
    posi <- posi[ooo,,drop=FALSE];
    # data.framing
    res <- as.data.frame(posi[,c("li1","co1","li2","co2"),drop=FALSE]);
    res$tag <- as.factor(tags[posi[,"tag",drop=FALSE]]);
    if (!is.null(names(tags))) {
      dimnames(res)[[1]] <- ftag[ooo];
    }
    # checking possible non overlapping
    if (check) {
      mama <- max(posi[,"co1"],posi[,"co2"])+2;
      debu <- posi[,"li1"] * mama + posi[,"co1"];
      fini <- posi[,"li1"] * mama + posi[,"co2"];
      crit <- c((debu[-1] <= fini[-length(fini)]),FALSE);
      if (!all(!crit)) {
        ouou <- which(crit);
        ouou <- sort(unique(pmin(pmax(1,c(ouou,ouou-1,ouou+1)),nrow(posi))));
        print(res[ouou,,drop=FALSE]);
        rrrerreur(NULL,"Overlapping discovered...");
      }
    }
  } else {
    res <- as.data.frame(cbind(posi[,c("li1","co1","li2","co2"),drop=FALSE],tag=character(0)));
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrexplore8list <- function(lili,monitor=rrrrbsa0$monitor$v)
#TITLE  returns the structure of a list in matrix
#DESCRIPTION
# Recursively explores the branches of a list, returning them into a 
# character matrix: a row for each branch, the columns describing
# some characteristics.\cr
# The columns are 
# 'number' (the number of this component from its branch),
# 'numbers' (the succesions of all numbers leading to 
#            this component separated with spaces),
# 'name' (<NA> if does not exist),
# 'names' (succession of names by default separated with spaces),
# 'level' (the branch level, i.e. a numbering of the branches),
# 'depth' (the branch depth from the root, implicitely 0 for the root),
# 'class' (the classes of the component separated with spaces).
#DETAILS
# \samp{rrrrbsa0$sep1$v} is used to join the names, as a consequence it must
# not be present into the names. In case an error is issued.
#KEYWORDS IO
#INPUTS
#{lili} << The list structure to be explored.>>
#[INPUTS]
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# The resulting character matrix (see the description section).
#EXAMPLE
# rrrexplore8list(rrrrbsa0$lis1$v);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_04_01
#REVISED 10_06_23
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(lili,"list",-1,mensaje="lili must be a list");
}
#
# starting the table
nta <- c("numbers","number","names","name","depth","level","classes");
tab <- matrix(NA,0,length(nta));
dimnames(tab) <- list(NULL,nta);
if (length(lili) == 0) { return(tab);
} else {
    for (ii in rrrbf(lili)) {
        nom <- names(lili)[ii];
        if (is.null(nom)) { nom <- "<NA>";}
        if (nom=="") { nom <- "<NA>";}
        if (rrrbelong9(rrrrbsa0$sep1$v,nom)) {
            rrrerreur(list(rrrrbsa0$sep1$v,nom),
                   "This name for the list comprises 'rrrrbsa0$sep1$v' which is not accepted by 'rrrexplore8list'");
        }
        tab <- rbind(tab,c(ii,ii,nom,nom,1,1,NA));
    }
}
#
qq <- which(is.na(tab[,"classes"]));
if (length(qq) > 0) { qq <- qq[1];
} else { qq <- 0; }
# filling the table
niv <- 1;
while (qq > 0) {
    # the qq-th component must be explored
    coco <- paste("lili[[",
                  paste(strsplit(tab[qq,"numbers"]," ")[[1]],
                        collapse="]][["),
                  "]]",sep="");
    coco <- paste("coco <-",coco);
    eval(parse(text=coco));
    # completing the class
    tab[qq,"classes"] <- paste(class(coco),collapse=" ");
    # if a list adding more components
    if (tab[qq,"classes"] == "list") {
        pro <- length(strsplit(tab[qq,"numbers"]," ")[[1]])+1;
        niv <- niv + 1;
	for (ii in rrrbf(coco)) {
	    nom <- names(coco)[ii];
	    if (is.null(nom)) { nom <- "<NA>";}
	    if (nom=="") { nom <- "<NA>";}
            noms <- paste(tab[qq,"names"],nom,sep=rrrrbsa0$sep1$v);
            iii <- paste(tab[qq,"numbers"],ii);
	    tab <- rbind(tab,c(iii,ii,noms,nom,pro,niv,NA));
	}
    }
    # something more to explore?
    qq <- which(is.na(tab[,"classes"]));
    if (length(qq) > 0) { qq <- qq[1];
    } else { qq <- 0; }
}
# returning
tab;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrget8comp7list <- function(lili,tata,
                       monitor=rrrrbsa0$monitor$v)
#TITLE  returns components from a list structure
#DESCRIPTION
# This is not a end-user function at all! To understand it, one must
# have in mind how works the function \samp{rrrexplore8list}.\cr
# Returns components from a list structure as a one level list.
# The list \samp{lili} must have been explored with \samp{rrrexplore8list}
# and the branch(es) to return are indicated through their line numbers
#  (\samp{tata}) in the table it generates.
#DETAILS
# Names of the produced list are the stacked names of the initial list 
# \samp{lili}.
#KEYWORDS IO
#INPUTS
#{lili} << The list structure components of which have to be extracted.>>
#{tata} << The lines of the table provided by \samp{rrrexplore8list}.>>
#[INPUTS]
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# The resulting list with as many component as indicated rows.
#EXAMPLE
# uu <- list(A=1:3,
#            B=matrix(letters[1:20],nrow=5,dimnames=list(1:5,c("on","tw","th","fo"))),
#            C=list(a=1:3,b=letters,c=list(final=array(1:8,c(2,2,2)))));
# vv <- rrrexplore8list(uu);
# rrrget8comp7list(uu,vv[7,])[[1]];
##
# uu <- list(1:3,
#            matrix(letters[1:20],nrow=5,dimnames=list(1:5,c("on","tw","th","fo"))),
#            list(a=1:3,b=letters,c=list(final=array(1:8,c(2,2,2)))));
# vv <- rrrexplore8list(uu);
# rrrget8comp7list(uu,vv[7,])[[1]];
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_04_01
#REVISED 10_04_20
#--------------------------------------------
{
# in case of
if (length(tata)==0) {
    return(vector("list",0));
}
nta <- c("numbers","number","names","name","depth","level","classes");
# checking
if (monitor$chk$v) {
    rrrobject9(lili,"list",-1,mensaje="lili must be a list");
    #
    rrrobject9(tata,"character",-1,
               mensaje="'tata' must be a 'character'");
    if (!is.matrix(tata)) {
        rrrobject9(tata,"character",length(nta),
                   mensaje="When not a matrix, 'tata' must have the same length that 'nta'");
    } else {
        if (ncol(tata) != length(nta)) {
            rrrerreur(tata,"The 'tata' matrix must have got as many columns as the length of nta");
        }
        if (is.null(dimnames(tata)[[2]])) {
            rrrerreur(tata,"The 'tata' matrix must have named columns");
        }
        if (length(union(dimnames(tata)[[2]],nta))>length(nta)) {
            rrrerreur(list(tata,nta),"The 'tata' matrix must have named columns with 'nta'");
        }
    }
}
# preparation
if (!is.matrix(tata)) {
   tata <- matrix(tata,1,length(nta));
   dimnames(tata) <- list(NULL,nta);
}
#
# starting the resulting list
res <- vector("list",nrow(tata));
nono <- tata[,"names"];
nunu <- tata[,"numbers"];
if (length(unique(nono))<length(nono)) {
    nono <- nunu;
}
names(res) <- nono;
for (ii in rrrbc(nrow(tata))) {
    coco <- paste("lili[[",
                  paste(strsplit(nunu[ii],rrrrbsa0$sep0$v)[[1]],
                        collapse="]][["),
                  "]]",sep="");
    coco <- paste("coco <-",coco);
    eval(parse(text=coco));
    res[[ii]] <- coco;
}
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrlist2text <- function(lili,
                      tags=rrrrbsa0$tag1$v,
                      stag=c("/",";","/"),
                      comment="#",
                      comments=character(0),
                       monitor=rrrrbsa0$monitor$v
                     )
#TITLE  transforms a list into a character
#DESCRIPTION
# The reverse operation of \samp{rrrtext2list}. The list must be a
# rsbsa-list, that is complying
# some properties : all components of the [sub-]lists must be either
# a list or a [named] character vector/matrix/array. The number of nested list
# must not be greater than the number of rows in matrix \samp{tags}.
# Every list component must be named.\cr
# The idea is to get a character compatible with \samp{rrrtext2list} to produce back
# the object \samp{lili}.\cr
# Some comments are added to the content of the list by the function itself and/or
# according to the whish of the user, this is why an escaping character
# for comments is asked.
#DETAILS
# The character strings of the structure must not comprise
# the \samp{rrrrbsa0$sep0$v} constant but this global constant can be
# conveniently modified. The same for \samp{rrrrbsa0$sep1$v}.\cr
# Use is made of the general constant \samp{monitor$chk$v} for the checking.
#KEYWORDS IO
#INPUTS
#{lili} << The list structure to be transformed into a 
#          \samp{character}.>>
#[INPUTS]
#{tags} << Character matrix with two columns indicating the opening and closing
#          tags of the [sub]lists at different levels. Its row numbers gives the maximum
#          accepted number of levels. Opening tags must be different.>>
#{stag} << Three character strings indicating the tagging to define the separator for each
#          character vector \samp{stag[2]} between \samp{stag[1]} and \samp{stag[3]}.>>
#{comment} <<At the beginning of a line, it indicates that this line must not be
#          considered.>>
#{comments} <<Comments that the user want to be added at the beginning of the file.>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# The resulting character.
#EXAMPLE
# rrrlist2text(rrrrbsa0$lis1$v);
# rrrtext2list(rrrlist2text(rrrrbsa0$lis1$v));
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_03_25
#REVISED 10_09_13
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(lili,"list",-1,mensaje="lili must be a list");
    #
    rrrobject9(tags,"character",c(2,Inf),mensaje="tags must be a CHARACTER matrix with two columns");
    if (!is.matrix(tags)) {rrrerreur(tags,"tags must be a character MATRIX of two columns");}
    if (ncol(tags)!=2) {rrrerreur(tags,"tags must be a character matrix OF TWO COLUMNS");}
    if (length(unique(tags[,1]))!= nrow(tags)) { rrrerreur(tags,"Opening tags are not unique!");}
    #
    rrrobject9(stag,"character",3,mensaje="stag must be a character of length three");
    if (any(nchar(stag) == 0)) { rrrerreur(stag,"stag must not have zero length components");}
    if (length(grep(" ",stag[c(1,3)])) > 0) { rrrerreur(stag,"stag must not comprise spaces in first and/or third position");}
    #
    rrrobject9(comment,"character",1,mensaje="'comment' must be a character of length one");
    #
    rrrobject9(comments,"character",-1,mensaje="'comment' must be a character of any length");
}
# 
# preparing
res <- character(0);
#
# issuing the starting comments
res <- c(res,comment);
res <- c(res,paste(comment,"  This character was created on",rrrnow("d"),"by list2file of rbsa package"));
res <- c(res,paste(comment,"  from the object :",deparse(substitute(x))));
res <- c(res,comment);
res <- c(res,comment);
for (ii in rrrbf(comments)) {
    res <- c(res,paste(comment,comment,rrrrbsa0$sep0$v,rrrrbsa0$sep0$v,comments[ii]));
}
res <- c(res,comment);
res <- c(res,comment);
#
# processing the list to translate it
if (length(lili) == 0) {
    res <- c(res,paste(comment,"The proposed list was of length zero!"));
    res <- c(res,comment);
} else {
    # exploring the list
    liliex <- rrrexplore8list(lili);
    # checking the list-vma nature 
    if(!all(liliex[,"classes"] %in% c("list",
                                      "logical","integer","numeric","character",
                                      "matrix","array"))) {
        rrrerreur(liliex,"The 'lili' list is not composed of list/vector/matrix/array");
    }
    # checking the existence of names for all levels
    if (any(is.na(liliex[,"name"]))) {
        rrrerreur(liliex,"Not all 'lili' components have got a name");
    }
    # checking the levels
    xlev <- max(as.numeric(liliex[,"depth"]));
    if (xlev > nrow(tags)) {
        rrrerreur(list(liliex,tags),"'lili' has got too many levels for the proposed 'tags'");
    }
    #
    # ordering the table
    uv <- order(liliex[,"numbers"]);
    # writing down the file
    for (ii in uv) {
        niv <- as.numeric(liliex[ii,"depth"]);
        if (niv == 1) { res <- c(res,paste0(comment,rrrform3repeat("=",20)));}
        if (niv == 2) { res <- c(res,comment);}
        rrr <- paste0(tags[niv,1],liliex[ii,"name"],tags[niv,2]);
        if (liliex[ii,"classes"] != "list") {
            coco <- rrrget8comp7list(lili,liliex[ii,],monitor=monitor)[[1]];
            caca <- rrrvma2text(coco,rrrrbsa0$sep1$v);
            rrr <- paste0(rrr,stag[1],stag[2],stag[3],caca$type);
            res <- c(res,rrr);
            res <- c(res,paste(caca$character,collapse=stag[2]));
        } else {
            res <- c(res,rrr);
        }
    }
}
#
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrvma2text <- function(x,xsep=rrrrbsa0$sep1$v,
                       monitor=rrrrbsa0$monitor$v)
#TITLE  transforms a vector (or matrix, or array) into a character
#DESCRIPTION
# from a vector, or a matrix, or
# an array, builds a \samp{character} vector. More or less the 
# inverse function of \samp{rrrtext2vma}. This vector is the first
# component of a returned list, the second component of the list 
# gives the type (\samp{vector}, \samp{matrix} or \samp{array}) of \samp{x},
# the converted object.
#DETAILS
# When some dimnames exist, the possible missing
# ones will be added.
#KEYWORDS IO
#INPUTS
#{x} << The object to transform.>>
#[INPUTS]
#{xsep} << \samp{character(1)} to be use for the separations.>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE a list with two components: \samp{[[1]]} the coded character vector and
# \samp{[[2]]} the type according to \samp{rrrtext2vma}.
#EXAMPLE
#####
## vectors
# rrrvma2text(letters);
# x <- letters; names(x) <- LETTERS;
# xx <- rrrvma2text(x);
# rrrtext2vma(xx[[1]],xx[[2]]);
# rrrvma2text(character(0));
#####
## matrices
# x <- matrix(1:20,4);
# rrrvma2text(x);
# dimnames(x) <- list(letters[1:4],LETTERS[1:5]);
# rrrvma2text(x);
# x1 <- matrix(NA,3,0);
# xx1 <- rrrvma2text(x1);
# rrrtext2vma(xx1[[1]],xx1[[2]]);
# dimnames(x1) <- list(c("i","ii","iii"),NULL);
# xx1 <- rrrvma2text(x1);
# rrrtext2vma(xx1[[1]],xx1[[2]]);
#####
## arrays
# x <- array(1:24,2:4);
# rrrvma2text(x);
# dimnames(x) <- list(1:2,c("i","ii","iii"),c("I","II","III","IV"));
# rrrvma2text(x,xsep="|||");
# x0 <- array(NA,c(3,0,2));
# dimnames(x0) <- list(1:3,NULL,c("i","ii"));
# xx0 <- rrrvma2text(x0);
# rrrtext2vma(xx0[[1]],xx0[[2]]);
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_03_29
#REVISED 10_06_29
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    #
    rrrobject9(xsep,"character",1,mensaje="must indicate the character string of separation");
    #
    if (!is.vector(x) &
        !is.matrix(x) &
        !is.array(x)) {
        rrrerreur(class(x),"'x' must be a vector or a matrix or an array!");
    }
}
#
res <- vector("list",2);
names(res) <- c("character","type");
#
if (is.array(x)&(!is.matrix(x))) {
    # dealing with an array
    res[[1]] <- c(as.character(dim(x)),xsep);
    if (is.null(dimnames(x))) {
        res[[2]] <- rrrrbsa0$vma$v["a"];
    } else {
        nna <- dimnames(x);
        for (hh in rrrbf(nna)) {
            if (is.null(nna[[hh]])) {
                nna[[hh]] <- rrrbc(dim(x)[hh]);
            }
        }
        if (is.null(names(nna))) {
            for (ii in rrrbf(nna)) {
                res[[1]] <- c(res[[1]],nna[[ii]],xsep);
            }
        res[[2]] <- rrrrbsa0$vma$v["A"];
        } else {
            res[[1]] <- c(res[[1]],names(nna),xsep);
            for (ii in rrrbf(nna)) {
                res[[1]] <- c(res[[1]],nna[[ii]],xsep);
            }
        res[[2]] <- rrrrbsa0$vma$v["B"];
        }
    }
    res[[1]] <- c(res[[1]],as.character(x));
} else {
    if (is.matrix(x)) {
        # dealing with a matrix
        if (is.null(dimnames(x))) {
            res[[1]] <- character(0);
            for (ii in rrrbc(nrow(x))) {
                res[[1]] <- c(res[[1]],x[ii,],xsep);
            }
            res[[1]] <- res[[1]][-length(res[[1]])];
            res[[2]] <- rrrrbsa0$vma$v["m"];
        } else {
            nna <- dimnames(x);
            for (hh in rrrbf(nna)) {
                if (is.null(nna[[hh]])) {
                    nna[[hh]] <- rrrbc(dim(x)[hh]);
                }
            }
            res[[1]] <- c(as.character(nna[[2]]),xsep);
            for (ii in rrrbc(nrow(x))) {
                res[[1]] <- c(res[[1]],nna[[1]][ii],as.character(x[ii,]),xsep);
            }
            res[[1]] <- res[[1]][-length(res[[1]])];
            res[[2]] <- rrrrbsa0$vma$v["p"];
        }
    } else {
        # dealing with a simple vector
        if (is.null(names(x))) {
            res[[1]] <- as.character(x);
            res[[2]] <- rrrrbsa0$vma$v["v"];
        } else {
            res[[1]] <- c(names(x),as.character(x));
            res[[2]] <- rrrrbsa0$vma$v["V"];
        }
    }
}
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrrbsa0 <- 
#TITLE list of the /rbsa/ constants
#DESCRIPTION
# Just a list of constants. Each constant is defined with a name, a definition
# and a value. The initial values can be modified by the user. Some
# specifies the way the \pkg{rbsa} calculations and displays are made. Values
# definitions and names can be obtained with the function \samp{rrrdisplay8k}.\cr
# An important role of the constant is to drive the behavior of many functions
# through \samp{rrrrbsa0$monitor$v} which is itself a \samp{rbsa} constant.
#DETAILS
# It is a named list, one component for each constant.
# A named sublist is associated to each constant with two components: \samp{$d}
# for the definition and \samp{$v} for the value. Be aware that the value
# can be any object (vector, list, matrix, function,...)
#KEYWORDS misc helpful
#INPUTS
#[INPUTS]
#VALUE
# A list see the details for the description.
#EXAMPLE
# print(rrrdisplay8k(rrrrbsa0,"names"));
# print(rrrdisplay8k(rrrrbsa0));
# print(rrrdisplay8k(rrrrbsa0$monitor$v,"names"));
# print(rrrdisplay8k(rrrrbsa0$monitor$v));
# str(rrrrbsa0$monitor$v);
#REFERENCE
#SEE ALSO rrrrbsa7list9 rrrdisplay8k
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 09_09_16
#REVISED 14_09_01
#--------------------------------------------
{
  # additional variables
  axfpx <- "graph";
  vma <- c("c","C","v","V","u","U","m","n","o","p","M","N","O","P","a","A","b","B");
  names(vma) <- vma;
  # definition of the different constants
  list(
       mwi=list(v=70,
                d="The width (nbr of characters) for printing paragraphs"),
       mfi=list(v=TRUE,
                d="Must the results be directed to files (if not to the screen)?"),
       mgr=list(v="pdf",
                d="Type of graphics files"),
       mnd=list(v=3,
                d="Number of decimals when printing"),
       mep=list(v=1,
                d="When printing an object: *e*m*p*hasize level"),
       fpx=list(v=axfpx,
                d="Prefix for the resulting files"),
       fou=list(v=paste(axfpx,"txt",sep="."),
                d="Standard file for text outputs"),
       fin=list(v=paste(axfpx,"int.txt",sep="."),
                d="Standard intermediary file"),
       tyl=list(v=c("A","a","n"),
                d="Types of bullets for standard lists"),
       l00=list(v="~",
                d="Symbol for the null bullets (see read8list)"),
                # due to the use of easyp3cut, no nesting parenthesis are allowed!
                # also opening and closing must be different
       cpt=list(v= matrix(c("{{","}}",
                            "(|","|)",
                            "[" ,"]" ,
                            "<<",">>"),ncol=2,
                          byrow=TRUE,dimnames=list(c("nodes","rounding","variables","vectors"),
                                                   c("opening","closing"))),
                d="Different closing parentheses as a dimnamed matrix"),
       nch1=list(v=structure("1", .Names = "A"),
                 d="Example 1 of a named character objects"),
       nnu1=list(v=structure(1, .Names = "A"),
                d="Example 1 of a named numeric objects"),
       lis1=list(v=list(A=1:3,
                        B=matrix(letters[1:20],nrow=5,dimnames=list(1:5,c("on","tw","th","fo"))),
                        C=list(a=1:3,b=letters,c=list(final=array(1:8,c(2,2,2))))),
                 d="Example1 of a list objects"),
       dfr1=list(v=data.frame(F1=factor(rep(1:3,each=4)),F2=factor(rep(1:4,3))),
                 d="Example1 of data.frame object"),
       dfr2=list(v=data.frame(SEX=factor(c("M","M","F","F","F","M")),
                              AGE=20+5*(1:6),
                              HGT=185-1:6),
                 d="Example2 of data.frame object"),
       dfr3=list(v=data.frame(SEX=factor(c(rep("M",50),rep("F",50))),
                              AGE=round(40+sin(1:100)*20),
                              HGT=c(round(170+cos(1:50)*17),round(155+cos(1:50)*16)),
                              WGT=c(round(75+cos(50:1)*15),round(60+cos(50:1)*20))),
                d="Example3 of data.frame object"),
       #
       vma=list(v=vma,
                d="Different types of vma (see rrrtext2vma for details)"),
       sep0=list(v=" ",
                 d="Basic string separator"),
       sep1=list(v="//",
                 d="Second string separator"),
       sep2=list(v=";",
                 d="Third string separator"),
       tag1=list(v=matrix(c("<<",">>",
                            "[[","]]",
                            "((","))"),
                          ncol=2,byrow=TRUE),
                 d="Standard tagging for text file associated to list"),
       tag2=list(v=matrix(c("[[","]]",
                            "((","))"),
                          ncol=2,byrow=TRUE),
                d="Alternative tagging for text file associated to list"),
       # 
       monitor=list(v=list(chk=list(v=TRUE, d="must arguments be checked?"),
                           pau=list(v=FALSE, d="must call to rrrpause be effective in some functions?"),
                           war=list(v=FALSE,d="must warnings be considered as fatal errors"),
                           fat=list(v=TRUE, d="must fatal errors be considered as fatal errors"),
                           ind=list(v=2,    d="size for indentation")),
                    d="global variable to monitor the behavior of some /rbsa/ functions"),
       tgs=list(v=list(
                  deb=list(v="<<<",t="Start",d="Starting point"),
                  mil=list(v="---",t="Middle",d="Intermediary point"),
                  fin=list(v=">>>",t="End",d="Ending point")
                 ),
                d="Example of tagging for text"),
       text1=list(v=c(paste0(rep(c(1:9,0),1),collapse=""),
                      paste0(rep(c(1:9,0),2),collapse=""),
                      paste0(rep(c(1:9,0),3),collapse=""),
                      paste0(rep(c(1:9,0),4),collapse=""),
                      paste0(rep(c(1:9,0),3),collapse=""),
                      paste0(rep(c(1:9,0),2),collapse=""),
                      paste0(rep(c(1:9,0),1),collapse="")),
                  d="text example number 1"),
       text2=list(v=c("Once upon a time there was a sweet little girl.",
                      "Everyone who saw her liked her, but most of all her grandmother",
                      "Once she gave her a little cap made of red velvet.",
                      "Because it suited her so well, and she wanted to wear it all the time,",
                      "she came to be known as Little Red Riding Hood.",
                      "One day her mother said to her: ",
                      "Come Little Red Riding Hood. Here is a piece of cake and a bottle of wine.",
                      "Take them to your grandmother.",
                      "She is sick and weak, and they will do her well.",
                      "Mind your manners and give her my greetings.",
                      "Behave yourself on the way, and do not leave the path,",
                      "or you might fall down and break the glass,",
                      "and then there will be nothing for your sick grandmother."),
                  d="text example number 2"),
       text3=list(v=c("<<< A gentle introduction!",
                      "---",
                      "A quite interesting first part",
                      "---",
                      "The core of the story (1)",
                      "The core of the story (2)",
                      "The core of the story (3)",
                      "---",
                      "Some concluding remarks.",
                      ">>>"),
                  d="text example number 3 (tagged)"),
       text4=list(v=c("Every paragraph must be followed ",
                      "with an empty line",
                      "",
                      "A quite interesting first part",
                      "---",
                      "---",
                      "",
                      "{#} << First item>>",
                      "{#} <<Second item>>",
                      "{#} << Third item",
                      "       can comprise several",
                      "       lines!>>",
                      "{#} <<Fourth item>>",
                      "",
                      "Some concluding remarks."),
                  d="text example number 4 (tagged)"),
       text5=list(v=c("Every paragraph must be followed ",
                      ")BULLETS( * + -",
                      "* FIRST LEVEL 1",
                      " + Second Level 1",
                      "   - third level 1",
                      "   - third level 2",
                      " + Second Level 2",
                      " + Second Level 3",
                      "* FIRST LEVEL 2",
                      " + Second Level 1",
                      " + Second Level 2",
                      "   - third level 1",
                      "   - third level 2",
                      "* FIRST LEVEL 3"),
                  d="text example number 5 (tagged)")
    );
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3places8word <- function(text,word,
                            column=c(1,Inf),
                            which=c(1,Inf)
                            )
#TITLE returns the places of one word
#DESCRIPTION from a text,
# gives back the positions of a word indicated by
# the associated intervals as defined in the function \samp{rrrtext3interval}. 
# The search can be constrained to specific columns. The output can
# also be restricted.
#DETAILS
# The word cannot be extended upon two successive lines but the same line
# can have more than one word. Be aware that overlapping patterns are not 
# all detected (see one of the examples).
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#{word}<< \samp{character(1)} the word to be found.>>
#[INPUTS]
#{column} << The columns where the first character of the word must found.
# \samp{c(1,1)} means that it must be at the very start of a line. \samp{c(10,12)}
# means that it must start on the 10th, 11th or 12th column of a line.>>
#{which} <<Which occurences of \samp{word} (not the line numbers)
# must be returned defined by
# the number of the first one and the number of the last one.\cr
# So \samp{c(2,2)} will designate the second and only the second;
# \samp{c(1,5)} will ask for the first five. When the components are
# both negative, the numbering is done from the end, so \samp{c(-1,-1)}
# means the last one and \samp{c(-1,-3)} asks for the last three ones
# given starting from the last. >>
#VALUE
# A four column matrix, each row corresponding to a word place with the help of
# an interval.\cr
# For negative values of \samp{which}, the order of occurences is reversed: the last
# found will be in the first row of the matrix output.
#EXAMPLE
# rrrtext3places8word(letters,"j");
# rrrtext3places8word(letters,"J");
# rrrtext3places8word(c("Il etait une fois","un petit et rouge chaperon"),"et");
# rrrtext3places8word(c("Il etait une fois","un petit et rouge chaperon"),"et",which=c(2,3));
# rrrtext3places8word(c("Il etait une fois","un petit et rouge chaperon"),"et",which=-c(1,3));
# rrrtext3places8word(c("# Il etait une fois"," #un petit et rouge chaperon"),"#");
# rrrtext3places8word(c("# Il etait une fois"," #un petit et rouge chaperon"),"#",column=c(1,2));
# rrrtext3places8word(c("# Il etait une fois"," #un petit et rouge chaperon"),"#",column=c(2,2));
# # overlapping pattern  
# rrrtext3places8word("aaaa","aa");
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE Think of a way to introduce "end of line" as a possible word.
#       Improve the case of overlapping patterns.
#AUTHOR J.-B. Denis
#CREATED 13_07_08
#REVISED 13_10_30
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  rrrobject9(word,"character",len=1);
  rrrobject9(which,"numeric",len=2);
  rrrobject9(column,"numeric",len=2);
  # initializing
  res <- matrix(NA,0,4,dimnames=list(NULL,c("i1","p1","i2","p2")));
  # degenerate cases
  if (length(text)==0) { return(res);}
  if (length(word)==0) { return(res);}
  # reversing if negative numbering
  if (which[1] < 0) {
    reve <- TRUE;
    text <- rev(text);
    which <- -which;
  } else { reve <- FALSE;}
  # checking
  if (diff(which) < 0) {
    stop("'which' not consistent");
  }
  # looking for the word
  ouhou <- gregexpr(word,text,fixed=TRUE);
  # looping onto the lines
  nbw <- 0; nlw <- nchar(word); nbtl <- length(text);
  for (nbl in rrrbf(text)) {
    # looking in the line
    if (ouhou[[nbl]][1]>0) {
      # at least an occurrence
      for (oo in rrrbc(length(ouhou[[nbl]]))) {
        nbw <- nbw + 1;
        if (nbw >= which[1]) {
          if ((ouhou[[nbl]][oo] >= column[1]) &
              (ouhou[[nbl]][oo] <= column[2])) {
            res <- rbind(res,c(nbl,ouhou[[nbl]][oo],
                               nbl,ouhou[[nbl]][oo]+nlw-1));
          }
        }
      }
    }
    if (nbl==which[2]) { return(res);}
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrobject9 <- function(x,typ,len=-1,con=numeric(0),
                       mensaje=NULL,fatal=TRUE,na.rm=FALSE,
                       cla=NULL,speci=NULL)
#TITLE  checks the type, the length and the content of some standard object
#DESCRIPTION
#  checks the type, the length and the content
# of some standard object.
# If not correct, a message and a fatal error are issued.
# \samp{NA} are detected and considered as wrong or not.\cr
# When \samp{typ} is 'directory' or 'file', \samp{x} must be
# a character containing path(s) of directories or files.
#DETAILS
# 'integer' has not got the meaning in \samp{is.integer} R
# standard function. 'null' must be understood as
# resulting TRUE with function \samp{rrrvoid9}.
#KEYWORDS error
#INPUTS
#{x} <<object to be checked.>>
#{typ} << Character vector indicating the accepted types
# for the object to be checked.
# See the code for the list of accepted types, among
# them 'null', 'integer', 'numeric', 'character',
# 'logical', 'list', 'any', 'function', 'data.frame',
# 'matrix','file','directory',...
# As indicated, 'any' implies that no
# check of the type is performed.>>
#[INPUTS]
#{len} <<If \samp{length(len)==1}, the exact length
# to be checked, if not must be of length two for
# the possible range of \samp{length(x)}. When -1,
# no check on the length is performed.\cr
# For data.frame, it is the number of columns.
# When \samp{na.rm} the check is made after
# removing the NA values.>>
#{con} << If \samp{length(con)>0}, some check about
# the content of \samp{x} is done for some of the 
# types. More precisely for (integer, numeric): all
# values must belong to the interval \samp{[con[1],con[2]]}
# and for (character), the set of possible \samp{character(1)}
# is described by \samp{con}.>>
#{mensaje} << Some additional message to be
#            issued before stopping in case of error.>>
#{fatal} << what to do when discovering
#           an inconsistency ? TRUE: this function prints the message
#           and stops; FALSE: this function returns
#           the message as a character.>>
#{na.rm} << Must the \samp{NA} values be removed before checking?
#           This is implemented only for the types integer, numeric,
#           character and logical.>>
#{cla} << The class to which the object is supposed to belong. When \samp{NULL}
# no check is made with respect to the class.>>
#{speci} << Some additional specific check for some object types.
#           The following possibilities are implemented:\cr
#           - \samp{NULL} when not used,\cr
#           - \samp{"named"} to indicate that a vector must be named,\cr
#           - a \samp{matrix} with two rows and two column for matrix objects
#             to indicate intervals for the row number (first
#             row) and the column number (second row); \samp{NA} values means
#             no check,\cr
#           - a \samp{matrix} with \samp{ndim+1} rows and two columns for
#             array objects to indicate intervals for the number of
#             dimensions and the size of each of them; \samp{NA} values means
#             no check.
#         >>
#VALUE
# When everything is OK returns a \samp{character(0)}.
# If not according to the argument \samp{fatal} print an error message
# and stops, or just returns the message as a warning.
#EXAMPLE
# rrrobject9(letters,c("numeric","character"),c(20,30),mensaje="!OK");
# rrrobject9(structure("A",.Names="a"),"ncharacter",-1,LETTERS);
# rrrobject9("A","ncharacter",-1,mensaje="BAD",fatal=FALSE);
# rrrobject9(structure("A",.Names="a"),"character",-1,letters,"bad",fatal=FALSE);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE (i) improves the treatment of 'NA's; (ii) allows the check of forbidden values
#AUTHOR J.-B. Denis
#CREATED 09_05_02
#REVISED 14_08_16
#--------------------------------------------
{
# checking the arguments and preparation
# accepted types by rrrobject9
typena <- c("integer","numeric","character",
           "logical","ninteger","nnumeric",
           "ncharacter","nlogical",
            "file","directory");
# In typena are those types from which NA values 
# can be removed
types <- c(typena,"list","function","any",
           "data.frame","matrix","array","null");
type <- match(typ,types);
if (all(is.na(type))) {
    rrrerreur(list(types,typ),
           "OBJECT9: the proposed type was not in the list!"
          );
}
res <- character(0);
# dealing with possible NAs
narm <- 0;
if (na.rm) { if (all(!is.na(match(typ,typena)))) { 
    nu <- !is.na(x);
    narm <- length(x) - length(nu);
    x <- x[nu];
}}
#
# the null case as 'x'
if (is.null(x)) {
    if ("null" %in% typ) { return(character(0));
    } else { 
        res <- paste("The proposed object is NULL not in '",
                     paste(typ,collapse="/"),"'");
	if (!is.null(mensaje)) { res <- c(res,mensaje);}
        cat(paste0(rep("~",60),collapse=""),"\n");
	cat("From rrrobject9: in fact 'x' is '",
	    deparse(substitute(x)),"'\n");
	print.default(x);
	if (fatal) {rrrerreur("Message:",paste(res));}
    }
}
if (all(!is.na(match(typ,typena)))) { if (any(is.na(x))) { 
    rrrerreur(typ,"OBJECT9: 'x' is or comprises 'NA' which were not removed!");
}}
# possible type for the check
type <- types[type[!is.na(type)]];
if (!is.numeric(len)) {
    rrrerreur(len,"OBJECT9: 'len' must be NUMERIC of size one or two");
}
if ((length(len)==0)|(length(len)>2)) {
    rrrerreur(len,"OBJECT9: 'len' must be numeric of size ONE or TWO");
}
# class check
if (!is.null(cla)) {
  if (!(class(x)%in%cla)) {
    res <- c(res,paste("The class is ",class(x)," not in /",
                       paste(cla,collapse="/"),"/",sep=""));
  }
}
# processing for the length
if (length(len)==1) {
    if (len>=0) { if (length(x)!=len) {
        if (narm>0) { res <- c(res,paste(narm,"components were removed"));}
        res <- c(res,paste("OBJECT9: 'x' is not of length ",len));
    }}
} else {
    if ((length(x)<len[1])|(length(x)>len[2])) {
        if (narm>0) { res <- c(res,paste(narm,"components were removed"));}
        res <- c(res,paste("OBJECT9: 'x' has got a length outside:",paste(len,collapse=":")));
    }
}
# processing for the type
if (!("any" %in% type)) {
    ty <- FALSE;
    for (tt in type) {
        #
        if (tt=="integer")   { if (is.numeric(x)) {
            if (all(x==round(x))) {
                ty <- TRUE;
            }
        }}
        if (tt=="ninteger")   { if (is.numeric(x)) {
            if (all(x==round(x))) { if (length(names(x))>0) {
                ty <- TRUE;
            }}
        }}
        #
        if (tt=="numeric")   { if (is.numeric(x))   {
            ty <- TRUE;}
        }
        if (tt=="nnumeric")   { if (is.numeric(x)) {
            if (is.character(attributes(x)$names)) { ty <- TRUE;}
        }}
        #
        if (tt=="character") { if (is.character(x)) {
            ty <- TRUE;}
        }
        if (tt=="ncharacter")   { if (is.character(x)) {
            if (is.character(attributes(x)$names)) { ty <- TRUE;}
        }}
        #
        if (tt=="logical")   { if (is.logical(x))   {
            ty <- TRUE;}
        }
        if (tt=="nlogical")   { if (is.logical(x)) {
            if (is.character(attributes(x)$names)) { ty <- TRUE;}
        }}
        #
        if (tt=="function")  { if (is.function(x))  { ty <- TRUE;}}
        if (tt=="list")      { if (is.list(x))      { ty <- TRUE;}}
        if (tt=="data.frame"){ if (is.data.frame(x)){ ty <- TRUE;}}
        if (tt=="matrix")    { if (is.matrix(x))    { ty <- TRUE;}}
        if (tt=="array")     { if (is.array(x))     { ty <- TRUE;}}
        if (tt=="null")      { if (rrrvoid9(x))        { ty <- TRUE;}}
        #
        if (tt=="file") {
          if (all(rrrfidi9(x) == "f")) { ty= TRUE;
          } else { res <- "Not all components were detected as file(s)";}
        }
        if (tt=="directory") {
          if (all(rrrfidi9(x) == "d")) { ty= TRUE;
          } else { res <- "Not all components were detected as directory(ies)";}
        }
    }
    if (!ty) {
        res <- c(res,paste("Among type = ",paste(type,collapse="/")));
        rr <- "The class of 'x' is '";
        if (is.character(attributes(x)$names)) {
            rr <- c(rr,"named ");
        }
        rr <- c(rr,paste0(class(x),"'!"));
        res <- c(res,paste(rr,collapse=""));
        res <- c(res,"OBJECT9: 'x' does not belong to any of these!");
    }
}
#
# proceding for the content
if (identical(res,character(0))) { if (length(con)>0) {
    if (is.numeric(x)) {
        if (length(con) >= 2) {
            cond <- (x >= con[1]) & (x <= con[2]);
            if (!all(cond)) { res <- c(res,"Does not belong to the prescribed interval"); }
        } else {
            rrrerreur(list(typ,con),
                   paste("OBJECT9:",
                         "For this 'typ', non empty 'con'",
                         "must define an interval with at",
                         "least two numerical values"));
        }
    }
    #
    if (is.character(x)) {
        cond <- (x %in% con);
        if (!all(cond)) { res <- c(res,"Does not belong to the prescribed values"); }
    }
}}
#
# checking the specific requirement
if (!is.null(speci)) {
  # for vectors
  if (is.vector(x)) {
    if (speci=="named") {
      if (is.null(names(x))) {
        res <- c(res,"This vector was required to be named and isn't.");
      }
    } else {
      rrrerreur(speci,"For vectors, the only valid value for 'speci' is 'named'");
    }
  }
  # for matrices
  if (is.matrix(x)) {
    if (is.matrix(speci)) {
      if (all(dim(speci)==2)) {
        if (!any(is.na(speci[1,]))) {
          if (abs(rrrinterv7belonging(nrow(x),speci[1,]))>1) {
            inter <- paste(speci[1,],collapse="-");
            res <- c(res,paste("The number of rows is not in the specified",inter,"interval"));
          }
        }
        if (!any(is.na(speci[2,]))) {
          if (abs(rrrinterv7belonging(ncol(x),speci[2,]))>1) {
            inter <- paste(speci[2,],collapse="-");
            res <- c(res,paste("The number of columns is not in the specified",inter,"interval"));
          }
        }
      } else {
        rrrerreur(speci,"For matrices, the only valid value for speci is a matrix of dimension '2x2'");
      }
    } else {
      rrrerreur(speci,"For matrices, the only valid value for speci is a matrix of dimension '2x2'");
    }
  }
  # for arrays
  if (is.array(x) & !(is.matrix(x))) {
    if (is.matrix(speci)) {
      if (ncol(speci)==2) {
        if (!any(is.na(speci[1,]))) {
          if (abs(rrrinterv7belonging(length(dim(x)),speci[1,]))>1) {
            inter <- paste(speci[1,],collapse="-");
            res <- c(res,paste("The number of dimensions is not in the specified",inter,"interval"));
          }
        }
        for (dd in rrrbf(dim(x))) { if (nrow(speci)>dd) {
          if (!any(is.na(speci[dd+1,]))) {
            if (abs(rrrinterv7belonging(dim(x)[dd],speci[dd+1,]))>1) {
              res <- c(res,paste("The dimension size number",dd,"is not in the specified",
                                 paste(speci[dd+1,],collapse="-"),"interval"));
            }
          }
        }}
      } else {
        rrrerreur(speci,"For arrays, the only valid value for speci is a matrix with '2' columns");
      }
    } else {
      rrrerreur(speci,"For arrays, the only valid value for speci is a matrix with '2' columns");
    }
  }
}
#
# returning
if (!identical(res,character(0))) {
    if (!is.null(mensaje)) { res <- c(res,mensaje);}
    if (fatal) {
        cat(paste(rep("~",60),collapse=""),"\n");
	cat("From rrrobject9: in fact 'x' is '",
	    deparse(substitute(x)),"'\n",sep="");
	print.default(x);
        print(res);
        rrrerreur(NULL,"FATAL ERROR from rrrobject9");
    }
}
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3crop <- function(chaine,bef=rrrrbsa0$sep0$v,aft=bef,
                      mxm=Inf,empty=FALSE,
                      monitor=rrrrbsa0$monitor$v)
#TITLE  removes framing characters from a character string
#DESCRIPTION
# removes \samp{bef}s before and \samp{aft}s after a character string.
#DETAILS
#KEYWORDS utilities
#INPUTS
#{chaine} <<The character string to refine. 
#           Can be a vector.>>
#[INPUTS]
#{bef} << What to repeatedly remove at the beginning.>> 
#{aft} << What to repeatedly remove at the end.>>
#{mxm} << Maximum number of tags to remove.>>
#{empty} << Must remaining empty lines be removed?>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# \samp{character} after removings
#EXAMPLE
# rrrform3crop('IMPORTANT','IM',' ANT');
# rrrform3crop(c('   OUF ',' FOU ',''),' ','',1);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 08_12_10
#REVISED 14_01_22
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(chaine,"character",-1,mensaje="rrrform3crop: 'chaine' must be a character");
    rrrobject9(bef,"character",1,mensaje="rrrform3crop: Vector are not accepted for 'bef'");
    rrrobject9(aft,"character",1,mensaje="rrrform3crop: Vector are not accepted for 'aft'");
    rrrobject9(mxm,"numeric",1,mensaje="rrrform3crop: mxm must be numeric(1)");
}
# null case
if (length(chaine) == 0) { return(chaine);}
lb <- nchar(bef);
la <- nchar(aft);
for (ich in rrrbf(chaine)) {
    cha <- chaine[ich];
    # removing at the beginning of the string
    if (lb>0) {
	nbr <- 0;
	repeat {
	    deb <- substr(cha,1,lb);
	    if ((deb == bef) & (nbr < mxm)) {
		cha <- substring(cha,lb+1);
		nbr <- nbr+1;
	    } else { break;}
	}
    }
    # removing at the end of the string
    if (la>0) {
	nbr <- 0;
	repeat {
	    lc <- nchar(cha);
	    fin <- substr(cha,lc-la+1,lc);
	    if ((fin == aft) & (nbr < mxm)) {
		cha <- substring(cha,1,lc-la);
		nbr <- nbr+1;
	    } else { break;}
	}
    }
    chaine[ich] <- cha;
}
# removing empty lines
if (empty) {
  nbli <- length(chaine);
  for (ii in rev(rrrbc(nbli))) {
    if (chaine[ii] == "") { chaine <- chaine[-ii];}
  }
}
# returning
chaine;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrfile2text <- function(file,path="",
                      clean=TRUE,
                      ended=")STOP(",
                      comment="#",
                      skip=matrix(c(")START_SKIPPING(",")END_SKIPPING(",
                                    "/*","*/"),ncol=2,byrow=TRUE),
                      include=c(")FILE("),
                       monitor=rrrrbsa0$monitor$v
                     )
#TITLE  reads a file and transforms it in a single character
#DESCRIPTION
# reads a conveniently tagged file to produce a vector of characters
# where all non used lines are eliminated. Each component of the resulting
# vector is associated to an original line.
#DETAILS
# All tags are supposed to be in the first position of the line after
# cleaning when asked.\cr
# Successive performed actions are : (i) cleaning the lines, i.e. removing starting and
# ending spaces, (ii) eliminating commented lines, (iii) eliminating
# lines after a 'stop', (iv) including indicated files and (v) skipping 
# sequences of lines.
#KEYWORDS IO
#ALIAS inputting
#INPUTS
#{file} << file which have to be read and transformed into a list.>>
#[INPUTS]
#{path} << Directory containing the file.>>
#{clean} <<Indicates if starting and ending spaces must be eliminated at first.>>
#{ended} << The tag indicating the line from which to stop the reading.>>
#{comment} <<At the beginning of a line, it indicates that this line must not be
#          considered. More than one commenting character can be considered when
#          it is a vector. For instance \samp{c("#","\%")} means that
#          lines starting with an hash or a percent are comment lines.
#          If no comment line must be filtered, just give \samp{comment} the
#          value of \samp{character(0)}.>>
#{skip} << To indicate set(s) of lines to be skipped. Must be a character matrix
#          where the two columns correspond respectively to the opening and 
#          closing tags, and where each row is associate to a couple of tags.
#          Tags are considered successively following the order of these matrix rows;
#          that is skipping with the first row is performed, then with the remaining
#          lines, skipping witht the second row is performed, and so on.>>
#{include} << Tags to indicate a file (including possible path) by
#             a \samp{character(1)} to include at this point its contents
#             as a text file with the same tags specifications.
#             Including files can be recursive. >>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# a character of length, the number of retained lines.
#EXAMPLE
# sink("rbsa.text.txt")
# cat("# comments can be included as well\n")
# cat(" something\n");
# cat("/* skipping from here\n");
# cat("blablabla\n");
# cat("  */ until here (this line is ALSO eliminated\n");
# cat(" interesting:\n");
# cat("un dos tres\n");
# cat(")STOP(\n");
# cat(" It is finished!\n");
# cat(" Don't insist!\n");
# sink();
# rrrfile2text("rbsa.text.txt");
# unlink("rbsa.text.txt");
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 08_12_10
#REVISED 14_08_04
#--------------------------------------------
{
#
# checking
if (monitor$chk$v) {
    rrrobject9(file,"character",1,mensaje="file must indicate the name of one file");
    #
    rrrobject9(path,"character",1,mensaje="path must indicate the name of one directory");
    #
    if (path!="") { fifi <- paste(path,file,sep="/");} else { fifi <- file;}
    if (is.na(file.info(fifi)$size)) { rrrerreur(fifi,"This file seems not to exist.");}
    #
    rrrobject9(clean,"logical",1,mensaje="clean must be a logical(1)");
    #
    rrrobject9(ended,"character",1,mensaje="ended must be a character of length one");
    #
    rrrobject9(comment,"character",-1,mensaje="comment must be a character");
    #
    if (length(skip) == 0) {
      skip <- matrix("",0,2);
    } else {
      rrrobject9(skip,"matrix",speci=matrix(c(1,2,Inf,2),2),
              mensaje="skip must be a character MATRIX with two columns");
      rrrobject9(skip,"character",
              mensaje="skip must be a CHARACTER matrix with two columns");
    }
    #
    rrrobject9(include,"character",1,mensaje="include must be a character of length one");
    #
}
#
# reading the proposed file
if (path!="") {file <- paste(path,file,sep="/");}
lu <- readLines(file);
#
# cleaning
if (clean) { for (ii in rrrbf(lu)) {
    # removing framing spaces
    lu[ii] <- rrrform3crop(lu[ii]," "," ");
}}
#
# removing the lines after a possible stop
sto <- which(substr(lu,1,nchar(ended))==ended);
if (length(sto)>0) {
    lu <- lu[1:(sto[1]-1)];
}
#
# removing the empty lines
lu <- lu[nchar(lu)>0];
#
# removing the commented lines
if (!rrrvoid9(comment)) {
  lu <- rrrfilter8text(lu,comment,remove=TRUE,
                         exact=8,lower=FALSE,
                         monitor=monitor);
}
#
# removing the skipped lines
for (nn in rrrbc(nrow(skip))) {
    deb <- skip[nn,1]; fin <- skip[nn,2];
    debu <- substr(lu,1,nchar(deb)) == deb;
    fini <- substr(lu,1,nchar(fin)) == fin;
    sk1 <- apply(outer(which(debu)  ,rrrbf(lu),"<="),2,sum);
    sk2 <- apply(outer(which(fini)+1,rrrbf(lu),"<="),2,sum);
    nsk <- ((sk1-sk2) < 1);
    lu <- lu[nsk];
}
#
# including the indicated files
inclus <- which(substr(lu,1,nchar(include))==include);
for (ii in rrrbf(inclus)) {
    jj <- length(inclus) + 1 - ii;
    kk <- inclus[jj];
    fi <- strsplit(lu[kk]," ",fixed=TRUE)[[1]][2];
    plus <- Recall(fi,path=path,
                      clean=clean,
                      ended=ended,
                      comment=comment,
                      skip=skip,
                      include=include);
    plus <- c(lu[rrrbc(kk-1)],plus);
    if (kk < length(lu)) { plus <- c(plus,lu[(kk+1):length(lu)]);}
    lu <- plus;
}
#
# returning
lu;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrfidi9 <- function(chemins,retour=c("-","f","d"))
#TITLE  checks files and directories
#DESCRIPTION
#  checks if \samp{chemin} is a file or a directory
#DETAILS
#  the distinction between files and directories is
# based on \samp{file.info} function.
#KEYWORDS file directory
#INPUTS
#{chemins} << \samp{character}; for each element of it, the
#       existence as file or directory must be checked.>>
#[INPUTS]
#{retour} << \samp{character(3)} indicating what to return
#            in case of (non existence), (existence of a file)
#            and (existence of a directory).>>
#VALUE
# a character of same length as \samp{chemins} with one of the
# components of \samp{retour}.
#EXAMPLE
# rrrfidi9(c("/","~","~/.bashrc",getwd(),"toto.txt"));
#REFERENCE
#SEE ALSO rrrdipa
#CALLING
#COMMENT
#FUTURE introduce the notion of accessibility
#AUTHOR J.-B. Denis
#CREATED 14_01_22
#REVISED 14_06_25
#--------------------------------------------
{
  # checking the arguments
  rrrobject9(chemins,"character");
  rrrobject9(retour,"character",3);
  # initializing
  res <- rep(retour[1],length(chemins));
  # proceeding
  for (ii in rrrbf(chemins)) {
    if (nchar(chemins[ii])==0) {
      stop(paste0(" The ",ii,"th path is void!"));
    } else {
      fidi <- chemins[ii];
      if (file.exists(fidi)){
        fidir <- file.info(fidi)$isdir;
        if (!rrrvoid9(fidir)) {
          res[ii] <- retour[2+fidir];
        }
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrnow <- function(what="dh",format="red",seps=c("_","@",":"))
#TITLE  returns a character giving the present moment
#DESCRIPTION
# Returns a character giving the present moment
# with different components
#DETAILS
# based on Sys.time R function
#KEYWORDS utilities
#INPUTS
#[INPUTS]
#{what} <<\samp{character(1)} indicating which components to include
#  among:
#  \cr\samp{d} for the day (including year and month),
#  \cr\samp{h} for the hour,
#  \cr\samp{m} for the minute (including hour),
#  \cr\samp{s} for the second (including hour and minutes.>>
#{format} <<\samp{r} or \samp{"red"} for a reduced output.>>
#{seps} <<\samp{character(3)} providing the separations for the
#         different components of the reduced returned output.>>
#VALUE
# a character
#EXAMPLE
# cat("Now is",rrrnow(),"\n");
# cat(rrrnow("dhms","verbose"),"\n");  
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_04_19
#REVISED 14_06_04
#--------------------------------------------
{
# checking
rrrobject9(seps,"character",c(1,Inf));
if (length(seps) < 3) { seps <- c(seps,seps,seps);}
# adding the consequences to what
if (rrrbelong9("s",what)) { what <- paste0(what,"m");}
if (rrrbelong9("m",what)) { what <- paste0(what,"h");}
# getting the moment
mnt <- Sys.time();
maintenant <- as.character(mnt);
# getting each piece of the moment
an <- strsplit(maintenant,"-")[[1]][1];
mois <- strsplit(maintenant,"-")[[1]][2];
jour <-  strsplit(maintenant,"-")[[1]][3];
heure <- strsplit(jour," ")[[1]][2];
jour <- strsplit(jour," ")[[1]][1];
minute <- strsplit(heure,":")[[1]][2];
seconde <- strsplit(heure,":")[[1]][3];
heure <- strsplit(heure,":")[[1]][1];
# building the result
res <- "";
if (format[1]=="r") {
  seps <- rep(seps[1],3);
  format <- "red";
}
if (format[1]=="red") {
    if (rrrbelong9("d",what)) {
        res <- paste(res,paste(substr(an,3,100),mois,jour,sep=seps[1]),sep="");
        if (rrrbelong9("h",what)) { res <- paste0(res,seps[2]);}
    }
    if (rrrbelong9("h",what)) {
        res <- paste0(res,heure);
        if (rrrbelong9("m",what)) { res <- paste0(res,seps[3]);}
    }
    if (rrrbelong9("m",what)) {
        res <- paste0(res,minute);
        if (rrrbelong9("s",what)) { res <- paste0(res,seps[3]);}
    }
    if (rrrbelong9("s",what)) {
        res <- paste0(res,seconde);
    }
} else {
    if (rrrbelong9("d",what)) {
        res <- paste(res,"le",weekdays(mnt),jour,months(mnt),an);
        if (rrrbelong9("h",what)) { res <- paste0(res," `a ");}
    }
    if (rrrbelong9("h",what)) {
        res <- paste0(res,heure,"H");
        if (rrrbelong9("m",what)) { res <- paste0(res,"");}
    }
    if (rrrbelong9("m",what)) {
        res <- paste0(res,minute,"m");
        if (rrrbelong9("s",what)) { res <- paste0(res,"");}
    }
    if (rrrbelong9("s",what)) {
        res <- paste0(res,seconde,'s');
    }
}
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrbc <- function(nb)
#TITLE  sequences for easy looping
#DESCRIPTION
# This function returns \samp{1:nb} when \samp{nb > 0} and
#         \samp{numeric(0)} otherwise.
# Quite useful to prevent starting
# a loop of length nought
#DETAILS
# nb is rounded before use
#ALIAS looping
#KEYWORDS iteration helpful
#INPUTS
#{nb}    <<length of the loop>>
#[INPUTS]
#VALUE
# \samp{1:nb} if \samp{nb > 0}
# else \samp{numeric(0)}.
#EXAMPLE
# rrrbc(0);
# rrrbc(5);
# rrrbc(pi);
# rrrbc(4*pi);
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
# rrrbc for BouCle (loop in French)
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_18
#REVISED 09_04_05
#--------------------------------------------
{
if (is.null(nb)) {return(numeric(0));}
if (length(nb)!=1) {
    rrrerreur(nb,"rrrbc deals only with scalar nb");
}
if (nb > 0) {return(1:max(1,round(nb)));}
numeric(0);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrbf <- function(x)
#TITLE  sequence for looping over an object
#DESCRIPTION
# This function returns \samp{1:length(x)} when \samp{length(x) > 0} and
#         \samp{numeric(0)} otherwise.
# Quite useful to prevent starting
# a loop of length nought
#DETAILS
#KEYWORDS iteration helpful
#ALIAS
#INPUTS
#{x}    <<vector>>
#[INPUTS]
#VALUE
# \samp{1:length(x)} if \samp{length(x) > 0}
# else \samp{numeric(0)}.
#EXAMPLE
# rrrbf(0);
# rrrbf(5);
# rrrbf(character(0));
# rrrbf(letters);
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
# rrrbf for Boucle For the elements of an object
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
rrrbelong9 <- function(sch,ch,exact=FALSE,how="a",lower=FALSE,
                       monitor=rrrrbsa0$monitor$v)
#TITLE  indicates inclusion of character string
#DESCRIPTION
# Checks if some string(s) belong(s) to a series of strings.
# When \samp{exact} is FALSE, returns TRUE if the 
# character string \samp{sch} is included at least once
# into the character string \samp{ch}.\cr
# \samp{sch} can be a vector, in that case, the check is made for each
# of its components. According to \samp{how} the returned
# value is vectorial or scalar. When \samp{sch} is zero length then
# \samp{TRUE} is returned. When \samp{ch} is zero length
# then \samp{FALSE} is returned except when \samp{length(sch)} is zero.
#DETAILS
# More possibilities exists when \samp{exact} is \samp{6} or \samp{8}.
# Namely that the substring must be at the end or at the beginning of the
# reference chain.
#KEYWORDS utilities
#INPUTS
#{sch} <<The character string(s) to be found.>>
#{ch}  <<The character string(s) to investigate.>>
#[INPUTS]
#{exact} << When exact, one component must
# be strictly identical, if not a subtring is sufficient.
# See also the \bold{Details} section for two other possibilities.>>
#{how} << Indicates what to do when \samp{length(sch)>1}. The choice are 
# \samp{v}: a logical vector gives back each check independently;
# \samp{1}: returns \samp{TRUE} when at least one of the component belongs
# to the series \samp{ch} and \samp{a} when all components must comply to get TRUE.>>
#{lower} << Must the comparisons being done after case lowering?>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE A logical vector with length of \samp{sch}; or \samp{TRUE} or \samp{FALSE} when
#      \samp{how} is \samp{1} or \samp{a}.
#EXAMPLE
# rrrbelong9('a','non');
# rrrbelong9('o',c('non','oui'));
# rrrbelong9('o',c('non','oui'),6);
# rrrbelong9('o',c('non','oui'),8);
# rrrbelong9('O',c("oui"));
# rrrbelong9('O',c("oui"),lower=TRUE);
# rrrbelong9(c('o','n'),c('non','oui'),6,how='v');
# rrrbelong9(c('o','n'),c('non','oui'),8,how='v');
# rrrbelong9(c("o","oui"),c('non','oui'));
# rrrbelong9(c("o","oui"),c('non','oui'),how="v");
# rrrbelong9(c("A[SEX]","A[AGE]"),c("A[AGE]","A[SEX]"),how="a")
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 08_08_28
#REVISED 14_07_31
#--------------------------------------------
{
  # degenerate cases
  if (length(sch)==0) { return(TRUE);}
  if (length( ch)==0) { return(FALSE);}
  #
  if (monitor$chk$v) {
      rrrobject9(  sch,"character",c(1,Inf));
      rrrobject9(   ch,"character",c(1,Inf));
      if (is.logical(exact)) {exact <- exact*10;}
      rrrobject9(as.character(exact),"character",1,c(0,6,8,10));
      rrrobject9(  how,"character",1);
      # because 'rrrobject9' calls 'rrrbelong9' for this type of check!
      if (!(how %in% c("v","1","a"))) {
          rrrerreur(how,"Not accepted value for 'how'");
      }
    rrrobject9(lower,"logical",1);
  }
  if (lower) {
    sch <- tolower(sch);
    ch <- tolower(ch);
  }
  #
  res <- logical(length(sch));
  #
  #
  for (ii in rrrbf(sch)) {
    chaine <- sch[ii];
    if (exact== 6) {
      ch <- substr(ch,pmax(1,nchar(ch)-(nchar(chaine)-1)),nchar(ch));
    }
    if (exact== 8) {
      ch <- substr(ch,1,nchar(chaine));
    }
    if (exact== 0) {
      res[ii] <- (length(grep(chaine,ch,fixed=TRUE)) > 0);
    } else {
      res[ii] <- (sum(chaine==ch) > 0);
    }
  }
  # dealing with a synthetic question
  if (how == "1") { res <- any(res);}
  if (how == "a") { res <- all(res);}
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3replace <- function(text,translate,lowering=FALSE,
                           b.accepted=NULL,a.accepted=NULL,
                           b.rejected=NULL,a.rejected=NULL,
                           b.line=TRUE,e.line=TRUE)
#TITLE translating some words
#DESCRIPTION from a text,
# gives it back after replacing some words (in fact sequences of characters)
# by others. The translation can be conditionned by the surrounding
# characters of the sequence.
#DETAILS
# The words to translate cannot be extended upon two successive lines.\cr
# The translations are done sequentially along the rows of \samp{translate}
# so modifying their order can produce different results.\cr
# For obvious reason, \samp{b.accepted} and \samp{b.rejected} cannot be
# non \samp{NULL} simultaneously; the same for \samp{a.accepted} and
# \samp{a.rejected}.\cr
# Resulting empty lines are not removed neither already present empty lines.
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#{translate}<< A two-columns character matrix with, row by row, a word to
# translate and its translations. The second column can comprise \samp{""}
# not the first.>>
#[INPUTS]
#{lowering} << Must a lowering being applied before the selection
#              of the words to translate?>>
#{b.accepted} << A \samp{character} providing the surrounding characters 
#                before the words to determine the acceptance of the
#                translation.>>
#{a.accepted} << A \samp{character} providing the surrounding characters 
#                after the words to determine the acceptance of the
#                translation.>>
#{b.rejected} << A \samp{character} providing the surrounding characters 
#                before the words to determine the rejection of the
#                translation.>>
#{a.rejected} << A \samp{character} providing the surrounding characters 
#                after the words to determine the rejection of the
#                translation.>>
#{b.line} << Must words at the beginning of a line be translated?>>
#{e.line} << Must words at the end of a line be translated?>>
#VALUE
# The translated \samp{text}
#EXAMPLE
# let <- c("a","e","i","o","u","y");
# tra <- matrix(c(let,tolower(let)),ncol=2);
# rrrtext3replace(letters,tra);
# rrrtext3replace(c("Il etait une fois","un petit et rouge chaperon"),tra);
# rrrtext3replace(c("Il etait une fois","un petit et rouge chaperon"),tra,b.accepted=" ");
#REFERENCE
#SEE ALSO rrrtext3acceptance
#CALLING
#COMMENT This function was written to change the names of defined objects
#        within an R script
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 14_09_03
#REVISED 14_09_03
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  if (length(translate)==0) { return(text);}
  rrrobject9(translate,"character");
  if (length(translate)==2) { translate <- matrix(translate,1);}
  rrrobject9(translate,"matrix",speci=matrix(c(NA,2,NA,2),2));
  rrrobject9(b.accepted,c("null","character"));
  rrrobject9(b.rejected,c("null","character"));
  rrrobject9(a.accepted,c("null","character"));
  rrrobject9(a.rejected,c("null","character"));
  if ((!is.null(b.accepted)) & (!is.null(b.rejected))) {
    rrrerreur(list(b.accepted,b.rejected),"'b.accepted' and 'b.rejected' cannot be active simultaneously");
  }
  if ((!is.null(a.accepted)) & (!is.null(a.rejected))) {
    rrrerreur(list(a.accepted,a.rejected),"'a.accepted' and 'a.rejected' cannot be active simultaneously");
  }
  rrrobject9(b.line,"logical");
  rrrobject9(e.line,"logical");
  # looping on each word
  for (nn in rrrbc(nrow(translate))) {
    amo <- translate[nn,1];
    nmo <- translate[nn,2];
    # looking for the presence of the word to translate
    if (lowering) {
      ouou <- rrrtext3places8word(tolower(text),amo);
    } else {
      ouou <- rrrtext3places8word(text,amo);
    }
    # dealing with each occurrence
    for (oo in rev(rrrbc(nrow(ouou)))) {
      lign <- text[ouou[oo,1]];
      # determining the acceptance
      acce <- rrrtext3acceptance(lign,c(ouou[oo,2],ouou[oo,4]),
                              b.accepted=b.accepted,a.accepted=a.accepted,
                              b.rejected=b.rejected,a.rejected=a.rejected,
                              b.line=b.line,e.line);
      if (acce) {
        # making the translation
        nlig <- "";
        if (ouou[oo,2]>1) { nlig <- paste0(nlig,substr(lign,1,ouou[oo,2]-1));}
        nlig <- paste0(nlig,nmo);
        if (ouou[oo,4] < nchar(lign)) {
          nlig <- paste0(nlig,substr(lign,ouou[oo,4]+1,nchar(lign)));
        }
        text[ouou[oo,1]] <- nlig;
      }
    }
  }
  # returning
  text;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3brackets <- function(text,bra=c("{","}"),
                          col1=c(1,Inf),col2=c(1,Inf),
                          rembra=TRUE,
                          which=c(1,Inf))
#TITLE returns the content of pairs of brackets
#DESCRIPTION from a text,
# gives back the contents of a couple of tags (opening and
# closing tags) indicated by \samp{bra} under the constraint
# that the first character of each delimiter be in the column
# interval of each line indicated with \samp{col1} (respectively
# \samp{col2}).
#DETAILS
# Tags cannot be upon two successive lines. The same line
# can have more than one tag. The result is built by an analysis
# of the result provided by \samp{rrrtext3places8word}.\cr
# When an opening tag is discovered, the following opening
# tag is not considered before a closing tag has been encountered 
# for the previous one.
# So the sequence \samp{"\{ toto \{tut\} bof\} \{deux\}"} will gives two
# contents, respectively \samp{" toto \{tut"} and \samp{"deux"}.\cr
# When an opening tags remains not consistently closed a fatal
# error is issued with some indication.
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#[INPUTS]
#{bra}<< \samp{character(2)} the pair of tags to use.>>
#{col1} <<Positions within a line where the opening bracket has to be found.>>
#{col2} <<Positions within a line where the closing bracket has to be found.>>
#{rembra} <<Must the brackets be removed?>>
#{which} <<Which occurence of \samp{tag} (not the line numbers)
# must be returned defined by the
# the number of the first one and the number of the last one.>>
#VALUE
# A \samp{list} of texts. Portion within two braces are eliminated.
#EXAMPLE
# rrrtext3brackets(paste(letters,collapse=""),c("j","u"));
# rrrtext3brackets(c(" juste {un","deux ou trois} suffiront !"));
#REFERENCE
#SEE ALSO
#CALLING rrrtext3places8brackets
#COMMENT
#FUTURE Think to a way to introduce "end of line" as a possible tag.
#AUTHOR J.-B. Denis
#CREATED 13_10_28
#REVISED 13_10_31
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  rrrobject9(bra,"character",len=2);
  rrrobject9(which,"numeric",len=2);
  rrrobject9(col1,"numeric",len=2);
  rrrobject9(col2,"numeric",len=2);
  # looking for bracket places
  uuu <- rrrtext3places8brackets(text,bra,col1,col2,which);
  # getting the different components
  res <- vector("list",0);
  for (con in rrrbc(dim(uuu)[3])) {
    # getting the intervals
    res[[con]] <- rrrtext3stext(text,as.numeric((uuu[,,con])));
    # removing the closing parenthese
    if (rembra) {
      if (length(res[[con]]) > 0) {
        # closing one
        dli <- length(res[[con]]); dca <- nchar(res[[con]][dli]);
        res[[con]][dli] <- substr(res[[con]][dli],1,dca-nchar(bra[2]));
        # opening one
        dca <- nchar(res[[con]][1]);
        res[[con]][1] <- substr(res[[con]][1],nchar(bra[1])+1,dca);
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3preparation <- function(text,preparation="rbR",
                         content=c(1,length(text)),
                         ccom="#",
                         llength=50)
#TITLE Analyses and prepares a /text/ 
#DESCRIPTION from a text returns another text after some
# transformations.\cr
# Not intended for the standard user.
#DETAILS
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# to prepare.>>
#[INPUTS]
#{preparation} <<A \samp{character(1)} whose characters indicates the actions
# to perform (in this order):
#     \cr\samp{r} to remove the starting characters of \samp{ccom}.
#     \cr\samp{b} to remove 'blank' characters at both ends of lines.
#     \cr\samp{B} to remove repeated 'blank' characters within lines.
#     \cr\samp{R} to remove empty lines.
#     \cr\samp{w} to return the first word of the first line delimited with either blank, \samp{=}
#                 of \samp{<-}.
#     \cr \cr The following options are hierarchized (for instance \samp{j} is equivalent to \samp{cvVSj}).
#     \cr\samp{c} to concatanate lines (between lists if any) into a unique line;
#                 one space is added between two initial lines.
#     \cr\samp{v} to return the vector with one word one component from the concatenated line of \samp{c}.
#     \cr\samp{V} the same as \samp{v} but eliminated repeated words.
#     \cr\samp{S} the same as \samp{v} but sorting the found words.
#     \cr\samp{j} to gather the words obtained after \samp{v} into text components having less that
#                 \samp{llength} characters or only one word.>>
#{content} << Indicates which component of \samp{text} to prepare.
#             Usually \samp{numeric(2)} to indicate the interval of lines to consider.
#             when \samp{numeric(1)}, only this line. When \samp{0} or \samp{diff(content)<0} 
#             \samp{character(0)} is returned.>>
#{ccom} <<A \samp{character(1)} indicating which character(s) at the beginning
#         of lines must possibly be removed (the character(s), not the complete line). >>
#{llength} <<Maximum number of characters for a line (except when it comprises only one word).>>
#VALUE
# The transformed text, that is a \samp{character}.
#EXAMPLE
# rrrtext3preparation(rrrrbsa0$text2$v[1:3],preparation="j",llength=10)
# uu <- c("Il etait une fois un petit et rouge chaperon",
#         "qui voulait aller voir sa mere-grand");
# rrrtext3preparation(uu,"j",llength=20);
# rrrtext3preparation(uu,"j",llength=80);
# rrrtext3preparation( c(" Je veux   voir  "," et re-voir  "),"rbBc")
# rrrtext3preparation(c("# Je veux   voir  "," et re-voir  "),"rbBc")
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE Other functionalities could be added.
#AUTHOR J.-B. Denis
#CREATED 13_10_14
#REVISED 14_07_10
#--------------------------------------------
{
  # initializing
  res <- character(0);
  if (length(text)==0) { return(res);}
  # checking
  rrrobject9(text,
             "character",
             mensaje="Bad 'text' argument");
  rrrobject9(content,
             "numeric",len=1:2,con=c(0,Inf),
             mensaje="Bad 'content' argument (1)");
  if (content[1]==0)   { return(res);}
  if (diff(content)<0) { return(res);}
  if (length(content)!=2) { content <- rep(content[1],2);}
  rrrobject9(content,"numeric",len=2,con=c(1,length(text)),
               mensaje="Bad 'content' argument (2)");
  #
  rrrobject9(preparation,
             "character",len=1,
             mensaje="Bad 'preparation' argument");
  rrrobject9(ccom,
             "character",len=1,
             mensaje="Bad 'ccom' argument");
  rrrobject9(llength,
             "numeric",len=1,con=c(1,Inf),
             mensaje="Bad 'llength' argument");
  # restricting the treatment
  text <- text[content[1]:content[2]];
  ## removing starting tags
  if (rrrbelong9("r",preparation)) {
    nbc <- nchar(ccom);
    if (nbc>0) {
      for (tt in rrrbf(text)) {
        tex <- text[tt];
        if (substr(tex,1,nbc)==ccom) {
          text[tt] <- substr(tex,nbc+1,nchar(tex));
        }
      }
    }
  }
  ## removing ending blanks
  if (rrrbelong9("b",preparation)) {
    for (tt in rrrbf(text)) {
      tex <- text[tt];
      while (substr(tex,1,1)==" ") {
        tex <- substr(tex,2,nchar(tex));
      }
      while (substr(tex,nchar(tex),nchar(tex))==" ") {
        tex <- substr(tex,1,nchar(tex)-1);
      }
      text[tt] <- tex;
    }
  }
  ## removing repeated blanks
  if (rrrbelong9("B",preparation)) {
    for (tt in rrrbf(text)) {
      tex <- text[tt];
      while (rrrbelong9("  ",tex)) {
        tex <- gsub("  "," ",tex);
      }
      text[tt] <- tex;
    }
  }
  res <- text;
  ## removing empty lines
  if (rrrbelong9("R",preparation)) {
    for (tt in rev(rrrbf(text))) {
      if (nchar(text[tt])==0) {
        text <- text[-tt];
      }
    }
  }
  res <- text;
  ## looking for the first word
  if (rrrbelong9("w",preparation)) {
    if (length(text)== 0) {
      rrrerreur(NULL,"No line found to get a first word");
    }
    tex <- text[1];
    n1 <- nchar(strsplit(tex," ",fixed=TRUE)[[1]][1]);
    n2 <- nchar(strsplit(tex,"<-",fixed=TRUE)[[1]][1]);
    n3 <- nchar(strsplit(tex,"=",fixed=TRUE)[[1]][1]);
    ftr <- min(n1,n2,n3);
    return(substr(tex,1,ftr));
  }
  ## concatenation
  if (rrrbelong9("c",preparation)|
      rrrbelong9("v",preparation)|
      rrrbelong9("V",preparation)|
      rrrbelong9("S",preparation)|
      rrrbelong9("j",preparation)) {
    res <- paste(res,collapse=" ");
  }
  ## vectorisation
  if (rrrbelong9("v",preparation)|
      rrrbelong9("V",preparation)|
      rrrbelong9("S",preparation)|
      rrrbelong9("j",preparation)) {
    res <- strsplit(res," ",fixed=TRUE)[[1]];
  }
  ## eliminating redundancy
  if (rrrbelong9("V",preparation)|
      rrrbelong9("S",preparation)|
      rrrbelong9("j",preparation)) {
    res <- unique(res);
  }
  ## sorting
  if (rrrbelong9("S",preparation)|
      rrrbelong9("j",preparation)) {
    res <- sort(res);
  }
  ## justification
  if (rrrbelong9("j",preparation)) {
    if (length(res)>0) {
      rrr <- res[-1]; res <- res[1];
      for (cc in rrrbf(rrr)) {
        if (nchar(res[length(res)])+1+nchar(rrr[cc]) <= llength) {
          res[length(res)] <- paste0(res[length(res)]," ",rrr[cc]);
        } else {
          res <- c(res,rrr[cc]);
        }
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrdipa <- function(chemins,sep=.Platform$file.sep)
#TITLE  normalize directory paths
#DESCRIPTION
#  returns the paths contained into \samp{chemins} in
#  a normalized way. No check of existence of the path(s)
#  is performed, for that the function \samp{rrrfidi9(rrrdipa(...))} could be used.
#DETAILS
#  Multiple \samp{sep} are removed as well as any final \samp{sep}
#  with the exception of \samp{sep} being the unique path, that is
#  the root path in Linux operating system.
#KEYWORDS directory
#INPUTS
#{chemins} << \samp{character} supposed to contain directory paths.>>
#[INPUTS]
#{sep} << The separator to use.>>
#VALUE
# a character with same length as \samp{chemins} 
#EXAMPLE
# rrrdipa(c("/","~","~/.bashrc",getwd(),"toto.txt"));
#REFERENCE
#SEE ALSO rrrfidi9
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_01_23
#REVISED 14_06_25
#--------------------------------------------
{
  # checking the arguments
  rrrobject9(chemins,"character");
  rrrobject9(sep,"character",1);
  # initializing
  se2 <- paste0(sep,sep); lse <- nchar(sep);
  res <- chemins;
  # proceeding
  for (ii in rrrbf(res)) {
    # removing multiple "sep"
    while (length(grep(se2,res[ii],fixed=TRUE))>0) {
      res[ii] <- gsub(se2,sep,res[ii],fixed=TRUE);
    }
    # dealing with the last character
    if (nchar(res[ii])>lse) {
      laca <- substr(res[ii],nchar(res[ii])-lse+1,nchar(res[ii]));
      if (laca == sep) { res[ii] <- substr(res[ii],1,nchar(res[ii])-lse);}
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrerreur <- function(x,...,w=FALSE,monitor=rrrrbsa0$monitor$v)
#TITLE  issues an error message and concludes accordingly
#DESCRIPTION
# when called this function prints x, then displays a message before stopping 
# the process except if it is a warning or if the general constant
# \samp{monitor$fat$v} is true.
#DETAILS
#KEYWORDS error
#INPUTS
#{x} <<object to be printed before the message. When \samp{rrrvoid9(x)}
#      nothing is printed. When it is a list, all components
#      of the list are successively printed.>>
#{\dots}<<pieces of message to display after pasting>>
#[INPUTS]
#{w} << Indicates if it is a simple warning>> 
#{monitor} << List of monitoring constants, see \samp{rrrrbsa0$monitor} to
#             know its structure.>>
#VALUE
# nothing is returned
#EXAMPLE
# rrrerreur(matrix(1:4,2),"This matrix is not symmetric",w=TRUE)
# rrrerreur("Are you sure of it?",w=TRUE);
#
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_03
#REVISED 13_10_01
#--------------------------------------------
{
 cat(paste(rep("~",60),collapse=""),"\n");
if (!rrrvoid9(x)) {
    if (is.list(x)) {
        for (i in 1:length(x)) {
             cat(paste(rep("~",40),collapse=""),"\n");
            cat("<< Displaying ",deparse(substitute(x[[i]]))," >>\n")
            print(x[[i]]);
        }
    } else {
        cat("<< Displaying ",deparse(substitute(x))," >>\n")
        print(x);
    }
}
mensaje <- paste(...);
cat("<<<<< MESSAGE >>>>>\n");
print(mensaje);
if (w) {
    cat("/rbsa/","SIMPLE WARNING:\n");
} else {
    on.exit(traceback());
    cat("/rbsa/","ERREUR FATALE\n");
    cat(paste0(rep("~",60),collapse=""),"\n");
    if (monitor$fat$v) { stop("stopped by rbsa");}
}
invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3places8brackets <- function(text,bra=c("{","}"),
                                 col1=c(1,Inf),col2=c(1,Inf),
                                 which=c(1,Inf))
#TITLE returns the places of some pairs of brackets
#DESCRIPTION from a text,
# gives back the positions of a couple of brackets (opening and
# closing tags) indicated by \samp{bra} under the constraint
# that the first character of each delimiter be in the column
# interval of each line indicated with \samp{col1} (respectively
# \samp{col2}).
#DETAILS
# The tags cannot be upon two successive lines. The same line
# can have more than one tag. The result is built by an analysis
# of the result provided by \samp{rrrtext3places8word}.\cr
# When an opening tag is discovered, the following opening
# tag is considered only when a closing tag has been encountered.
# so the sequence \samp{"\{ toto \{tut\} bof\} \{deux\}"} will gives two
# contents, respectively \samp{" toto \{tut"} and \samp{"deux"}.\cr
# When an opening tags remains not consistently closed a fatal
# error is issued with some indications.
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#[INPUTS]
#{bra}<< \samp{character(2)} the pair of tags to use.>>
#{col1} <<Positions within a line where the opening bracket has to be found.>>
#{col2} <<Positions within a line where the closing bracket has to be found.>>
#{which} <<Which occurences of \samp{tag} (not the line numbers)
# must be returned defined by the
# the number of the first one and the number of the last one.>>
#VALUE
# A three dimensions array. First dimension for the 
# interval associated to the bracket (see function \samp{rrrtext3interval} for
# an explanation; second dimension for the first and last
# position of the content of each bracket; third dimension for the number of
# of discovered brackets.\cr
# Be aware that the bracket tags are included (for a technical reason).
#EXAMPLE
# rrrtext3places8brackets(paste(letters,collapse=""),c("j","u"));
# rrrtext3places8brackets(paste(letters,collapse=""),c("ab","xyz"));
# rrrtext3places8brackets(c(" juste {un","deux ou trois} suffiront !"));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE Think of a way to introduce "end of line" as a possible tag.
#AUTHOR J.-B. Denis
#CREATED 13_10_15
#REVISED 13_10_30
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  rrrobject9(bra,"character",len=2);
  rrrobject9(which,"numeric",len=2);
  rrrobject9( col1,"numeric",len=2);
  rrrobject9( col2,"numeric",len=2);
  # looking for opening positions
  ooo <- rrrtext3places8word(text,bra[1],which=which,column=col1);
  # looking for closing positions
  ccc <- rrrtext3places8word(text,bra[2],which=which,column=col2);
  # translating them into absolute rankings
  ooa <- matrix(NA,nrow(ooo),2);
  cca <- matrix(NA,nrow(ccc),2);
  for (ii in rrrbc(nrow(ooo))) {
    ooa[ii,] <- c(rrrtext3ij2n(ooo[ii,1:2],text),
                  rrrtext3ij2n(ooo[ii,3:4],text));
  }
  for (ii in rrrbc(nrow(ccc))) {
    cca[ii,] <- c(rrrtext3ij2n(ccc[ii,1:2],text),
                  rrrtext3ij2n(ccc[ii,3:4],text));
  }
  # initializing
  res <- array(NA,dim=c(2,2,nrow(ooa)),dimnames=list(c("li","co"),bra,NULL));
  # matching opening and closing
  nbb <- 0;
  while(nrow(ooa)>=1) {
    nbb <- nbb+1;
    # getting the next opening
    res[,1,nbb] <- ooo[1,1:2];
    # looking for the first associated closing
    fipo <- 1 + rrrtext3ij2n(ooo[1,3:4],text);
    bons <- which(cca[,1] >= fipo);
    cca <- cca[bons,,drop=FALSE];
    if (nrow(cca) == 0) {
      print(text);
      stop("'rrrtext3places8brackets' found an opening bracket without a closing one!");
    }
    ccc <- ccc[bons,,drop=FALSE];
    res[,2,nbb] <- ccc[1,3:4];
    # removing the closed opening brackets
    fipo <- 1 + rrrtext3ij2n(ccc[1,3:4],text);
    bons <- which(ooa[,1] >= fipo);
    ooa <- ooa[bons,,drop=FALSE];
    ooo <- ooo[bons,,drop=FALSE];
  }
  # removing non used brackets
  keep <- which(!is.na(res[1,1,]));
  res <- res[,,keep,drop=FALSE];
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrvoid9 <- function(x,vector=FALSE) 
#TITLE  tests the nullness of objects
#DESCRIPTION
# Returns TRUE is the structure is empty.\cr
# This function was proposed because the \samp{NULL} cannot replace
#  any kind of objects.\cr
#  Notice that \samp{rrrvoid9("")} is \samp{TRUE}.
#DETAILS
# Non-existing objects are not detected as void!
# Functions are considered as not void.
# Mainly \samp{NULL}, \samp{NA} and length-zero object are void
# then a \samp{data.frame} with at least one variable is declared as not void.
# Character \samp{""} is considered as void.
#KEYWORDS programming
#INPUTS
#{x}    <<object to be scrutinazed>>
#[INPUTS]
#{vector} <<This argument is considered only when \samp{x} is
#           a \samp{character} or a (\samp{numeric} or \samp{logical}). When \samp{FALSE} it is considered void
#           when all components are \samp{""} or \samp{NA}; when \samp{TRUE}, a vector
#           with same length as \samp{x} is returned indicating which
#           components are \samp{""} or \samp{NA}.>>
#VALUE
# TRUE when the object is considered as void
# FALSE if not. Can also be a vector for characters when \samp{vector}.
# See the proposed examples.
#EXAMPLE
# \dontrun{rrrvoid9(toto);}
# rrrvoid9(numeric(0));
# rrrvoid9(NULL);
# rrrvoid9("");
# rrrvoid9(c("","A"));
# rrrvoid9(c("","A"),vector=TRUE);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE Find a trick for non-existing objects: difficulty of the kind is that \samp{TRUE} 
#       is detected as void and more
#AUTHOR J.-B. Denis
#CREATED 07_10_15
#REVISED 14_08_28
#--------------------------------------------
{
  # non existing object
  #if (!exists(deparse(substitute(x)))) { return(TRUE);}
  # covering some special non null cases
  if (is.function(x))   { return(FALSE);}
  # standard cases
  if (is.null(x))                      { return(TRUE);}
  if (length(x)==0)                    { return(TRUE);}
  # looking when character
  if (is.character(x)) {
    if (vector) {
      return(x=="");
    } else {
      return(all(x==""));
    }
  }
  # looking when numeric
  if (is.numeric(x)|is.logical(x)) {
    if (vector) {
      return(is.na(x));
    } else {
      return(all(is.na(x)));
    }
  }
  # returning
  FALSE;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3title <- function(tit,empha=3,indent=2+2*empha,imp=TRUE)
#TITLE  prints or prepares a title
#DESCRIPTION
# prints or prepares the character string \samp{tit}
# with more or less emphasis.
# This function is a shortcut of the hidden function \samp{rrrform3titre},
# with some specialized calls.
#DETAILS
#KEYWORDS print
#INPUTS
#{tit}<<the title to print (just one line)>>
#[INPUTS]
#{empha} << Level of emphasizing.\cr
#          (0): single line without carriage return\cr
#          (1): single line\cr
#          (2): underlined\cr
#          (3): underlined and overlined\cr
#          (4): (2) + 1 line before\cr
#          (5): (3) + 1 line after\cr
#          (6): (2) + 2 lines before and after\cr
#          (7): corners + 1 line before and after (plus surrounding)\cr
#          (8): box + 1 lines before and after (plus surrounding)\cr>>
#{indent} << Number of spaces to introduce before the title>>
#{imp} << Printing is performed and nothing is returned.
#                If FALSE, the character string is returned 
#                (including possible new lines)>>
#VALUE
# either nothing or a character string according to imp
#EXAMPLE
# for (ii in 0:8) {rrrform3title("Some Title",ii,imp=TRUE)};
#REFERENCE
#SEE ALSO
#CALLING {rrrform3repeat}
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_24
#REVISED 10_09_15
#--------------------------------------------
{
# adjusting
empha <- round(max(0,min(8,empha)));
# preparing the argument for rrrform3titre
if (length(tit)>1) { tit <- paste(tit,collapse=" ");}
if (empha == 0) { tit <- paste0("<",tit,">")}
if (empha == 1) { tit <- paste0("(*)",tit,"(*)")}
sbef <- round(max(indent,0));
caret <- (empha != 0);
saft=""; lbef <- 0; laft <- 0;
box <- "no";
if (empha == 2) { box <- "un";}
if (empha == 3) { box <- "unov";}
if (empha == 4) { box <- "unov"; lbef <- 1;}
if (empha == 5) { box <- "unov"; laft <- 1; lbef <- 1;}
if (empha == 6) { box <- "unov"; laft <- 2; lbef <- 2;}
if (empha == 7) {
  box <- "cor" ; laft <- 1; lbef <- 1;
  tit<- paste0(" ",tit," ");
}
if (empha == 8) {
  box <- "box" ; laft <- 1; lbef <- 1;
  tit<- paste0(" ",tit," ");
}
# calling rrrform3titre
res <- rrrform3titre(tit,box=box,
                  lbef=lbef,sbef=sbef,
                  saft=saft,laft=laft,
                  charbox=c("+","|","+",
                            "-","+","|",
                            "+","-"," "),
                  alig=2,caret=caret,imp=imp);
# returning
if (imp) { return(invisible());}
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrparse8text <- function(text,item1=c("{","}"),item2=c("<<",">>"),
                       numb="#",bull="*",lsep="-")
#TITLE returns the parsed content of a text
#DESCRIPTION from a text comprising paragraphs and items
# finds the different components and returns them by means of a list.
# When the component is a paragraph, it is a \samp{character},
# When the component is an item list, it is a named list.\cr
# This function is not intended for standard users.
#DETAILS
# Each item of a list must comprise two values, framed with \samp{item1} and
# \samp{item2}. When the first value is \samp{numb}, it is a numbered
# item; when the first value is \samp{bull}, it is a bullet item;
# if not it is a labelled item.\cr
# Each component, paragraphs and items are supposed to be proposed on
# non overlapping lines.\cr
# Successive items are considered to belong to the list of items knowing 
# that empty lines (comprising zero characters) are first eliminated 
# (a line with a blank is not empty and will be considered as a paragraph).
# Also are eliminated lines starting with \samp{lsep}, their role is to
# separate distinct paragraphs and lists.\cr
# When the braces for items are not consistent, no error is reported but the
# staff is interpreted as part of a paragraph.\cr
# When two list items have got identical labels, an error is reported.
#KEYWORDS 
#INPUTS
#{text}<< The \samp{text} to be parsed. For the moment just a \samp{character}
# vector.>>
#[INPUTS]
#{item1}<< \samp{character(2)} the pair of tags to use to define the
# first value of an item. When \samp{numb} interpreted as an enumeration,
# when \samp{*} interpreted as an itemized list, if not a description
# list. The first character of \samp{item1[1]} must start at the beginning
# a line and the two braces must be on the same line.>>
#{item2}<< \samp{character(2)} the pair of tags to use to define the
# second value of an item.>>
#{numb} << \samp{character(1)} code to indicate automatically numbered items.>>
#{bull} << \samp{character(1)} code to indicate bullet items.>>
#{lsep} << \samp{character(1)} Each line starting with \samp{lsep} is
# considered a tagging line to separate two paragraphs or two item lists.
# They can be used to separate a paragraph and an item list but are useless.
# Separating lines within list items are not considered as separating.
# Successive separating lines are considered as a unique separating line.
# They are eliminated in the resulting list.>>
#VALUE
# A named \samp{list}. The names for paragraphs start with \samp{P}, those
# for item lists with \samp{L}.
#EXAMPLE
# rrrparse8text(c("{a}","<<","pour voir",">>"));
# uu <- c("1rst paragraph","","2d paragraph","",
#         "{#} <<un>>","{#}","<<deux>>","","3rd and last paragraph");
# rrrparse8text(uu);
# vv <- c("1rst paragraph","","2d paragraph","",
#         "{AA} <<un>>","{BBB}","<<deux>>","","3rd and last paragraph");
# rrrparse8text(vv);
# rrrparse8text(rrrrbsa0$text4$v);
#REFERENCE
#SEE ALSO
#CALLING 
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_12_30
#REVISED 14_08_07
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  rrrobject9(item1,"character",len=2);
  if (any(item1=="")) {
    rrrerreur(item1,mensaje="(1) A closing cannot be ''");
  }
  rrrobject9(item2,"character",len=2);
  if (any(item2=="")) {
    rrrerreur(item2,mensaje="(2) A closing cannot be ''");
  }
  rrrobject9(numb,"character",len=1);
  rrrobject9(bull,"character",len=1);
  rrrobject9(lsep,"character",len=1);
  # initialization
  res <- ouca <- vector("list",0);
  nt <- np <- nn <- ni <- nd <- nitem <- 0;
  lili <- rrrbf(text);
  # removing empty lines
  vides <- which(text=="");
  if (length(vides)>0) { text <- text[-vides];}
  # corner case
  if (length(text) == 0) { return(vector("list",0));}
  # sticking to its previous line
  # every line starting with an \samp{item2[1]}
  # (the used algorithm fails in that circumstance)
  colle <- (substr(text,1,nchar(item2[2]))==item2[2]);
  colle[1] <- FALSE;
  for (ii in rev(which(colle))) {
    text[ii-1] <- paste0(text[ii-1],text[ii]);
    text <- text[-ii];
  }
  # looking for the positions of valid list items
  b1 <- rrrtext3places8word(text,item1[1],column=c(1,1));
  if (nrow(b1)>0) {
    # possibly some lists are present
    eli <- TRUE;
    nb1 <- rrrtext3ij2n(b1[,1:2],text);
    b2 <- rrrtext3places8word(text,item1[2]);
    if (nrow(b2)>0) {
      nb2 <- rrrtext3ij2n(b2[,3:4],text);
      b3 <- rrrtext3places8word(text,item2[1]);
      if (nrow(b3)>0) {
        nb3 <- rrrtext3ij2n(b3[,1:2],text);
        b4 <- rrrtext3places8word(text,item2[2]);
          if (nrow(b4)>0) {
            nb4 <- rrrtext3ij2n(b4[,3:4],text);
            # finding which braces must be considered
            bra <- matrix(NA,0,4);
            while (length(nb1) > 0) {
              bb1 <- nb1[1]; nb1 <- nb1[-1];
              nb2 <- nb2[which(nb2>bb1)];
              if (length(nb2)>0) {
                bb2 <- nb2[1]; nb2 <- nb2[nb2>bb2];
                nb3 <- nb3[which(nb3>bb2)];
                if (length(nb3)>0) {
                  bb3 <- nb3[1]; nb3 <- nb3[nb3>bb3];
                  nb4 <- nb4[which(nb4>bb3)];
                  if (length(nb4)>0) {
                    bb4 <- nb4[1]; nb4 <- nb4[nb4>bb4];
                    bra <- rbind(bra,c(bb1,bb2,bb3,bb4));
                  }
                }
              }
            }
          } else {
            eli <- FALSE;
          }
      } else {
        eli <- FALSE;
      }
    } else { 
      eli <- FALSE;
    }
    
  } else {
    # no list
    eli <- FALSE;
  }
  ## identifying the role of the different lines
  # The default paragraphs
  lig <- rep("P",length(text));
  # The separating lines
  if (nchar(lsep)>0) {
    lig[substr(text,1,nchar(lsep))==lsep] <- "S";
  }
  # The item lists
  if (eli) {
    # adjusting for the brace widths
    bra[,1] <- bra[,1]+nchar(item1[1]);
    bra[,2] <- bra[,2]-nchar(item1[2]);
    bra[,3] <- bra[,3]+nchar(item2[1]);
    bra[,4] <- bra[,4]-nchar(item2[2]);
    # going back to (line,column) positions
    pbra <- array(NA,dim=c(4,2,nrow(bra)),
                  dimnames=list(c("s1","e1","s2","e2"),
                                c("row","col"),rrrbc(nrow(bra))));
    pbra["s1",,] <- t(rrrtext3n2ij(bra[,1],text));
    pbra["e1",,] <- t(rrrtext3n2ij(bra[,2],text));
    pbra["s2",,] <- t(rrrtext3n2ij(bra[,3],text));
    pbra["e2",,] <- t(rrrtext3n2ij(bra[,4],text));
    #
    for (li in rrrbc(dim(pbra)[3])) {
      lig[pbra["s1","row",li]:pbra["e2","row",li]] <- paste0("L",li,"-");
    }
    # checking that the item names are on the same line
    egal <- (pbra["s1","row",]!=pbra["e1","row",]);
    if (any(egal)) {
      print(text[sort(unique(c(pbra["s1","row",],pbra["e1","row",])))]);
      rrrerreur(pbra[,,egal,drop=FALSE],"Some item name not on the same line");
    }
  }
  ## numbering the different roles
  ## and initialization of the resulting list
  npa <- nli <- 0; courant <- "X";
  res <- vector("list",0); loli <- 0;
  for (lili in rrrbf(lig)) {
    if (lig[lili] == "P") {
      if (courant != "P") {
        npa <- npa+1;
        loli <- loli+1;
        res[[loli]] <- character(0);
        names(res)[loli] <- paste0("P",npa);
      }
      lig[lili] <- paste0(lig[lili],npa);
      courant <- "P";
    }
    if (substr(lig[lili],1,1) == "L") {
      if (courant != "L") {
        nli <- nli+1;
        loli <- loli+1;
        res[[loli]] <- vector("list",0);
        names(res)[loli] <- paste0("L",nli);
      }
      lig[lili] <- paste0(lig[lili],nli);
      courant <- "L";
    }
    if (lig[lili] == "S") {
      courant <- "S";
    }
  }
  ## decoding the lig coding
  llig <- lig[substr(lig,1,1)=="L"];
  llig <- substr(llig,2,nchar(llig));
  llig <- strsplit(llig,"-");
  lli1 <- as.numeric(sapply(llig,function(a)a[1]));
  lli2 <- as.numeric(sapply(llig,function(a)a[2]));
  ## filling the different components
  for (coco in rrrbf(res)) { 
    quoi <- names(res)[coco];
    if (substr(quoi,1,1)=="P") {
      # a paragraph
      res[[coco]] <- text[lig==quoi];
    } else {
      # a list
      # looking for the corresponding items
      nuli <- as.numeric(substr(quoi,2,nchar(quoi)));
      quelles <- (lli2==nuli);
      quelles <- unique(lli1[quelles]);
      labt <- character(0);
      for (ququ in rrrbf(quelles)) {
        labnum <- as.numeric(t(pbra[c("s1","e1"),,quelles[ququ]]));
        labite <- rrrtext3stext(text,labnum);
        if (labite==numb) { labite <- as.character(ququ);}
        if (labite==bull) { labite <- paste0("*",ququ);}
        if (labite %in% labt) {
          rrrerreur(labite,"This label was not unique as list item name");
        } else {
          labt <- c(labt,labite);
        }
        valnum <- as.numeric(t(pbra[c("s2","e2"),,quelles[ququ]]));
        valite <- rrrtext3stext(text,valnum);
        res[[coco]][[labite]] <- valite;
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrfile2list <- function(file,path="",
                      clean=TRUE,ended=")STOP(",comment="#",
                      skip=matrix(c(")SKIPPING(",")READING("),1),
                      include=c(")FILE("),
                      tags=rrrrbsa0$tag1$v,
                      sep=rrrrbsa0$sep0$v,
                      rsep="no_action",
                      stag=c("/","/")
                     )
#TITLE  reads a file and transforms it in a list (of lists) of characters
#DESCRIPTION
# reads a conveniently tagged text file into nested lists.
# It is just the linking of the two functions \samp{rrrfile2text}
# and \samp{rrrtext2list}, see their comments for the description of the arguments.
#DETAILS
#KEYWORDS IO
#INPUTS
#{file} << file to be read and transformed into a list.>>
#[INPUTS]
#{path} << Directory containing the file.>>
#{clean} <<Indicates if starting and ending spaces must be eliminated at first.>>
#{ended} << To indicate the line at which to stop the reading.>>
#{comment} <<At the beginning of a line, it indicates that this line must not be
#          considered. More than one commenting character can be considered when
#          it is a vector. For instance \samp{c("#","-")} means that
#          lines starting with an hash or a hyphen are comment lines.>>
#{skip} << To indicate set(s) of lines to be skipped. Must be a character matrix
#          where the two columns correspond respectively to the opening and 
#          closing tags, and where each row is associate to a couple of tags.
#          Tags are considered successively following the order of these matrix rows;
#          that is skipping with the first row is performed, then with the remaining
#          lines, skipping witht the second row is performed, and so on.>>
#{include} << Tags to indicate a file (including possible path) by
#             a \samp{character(1)} to include at this point its contents
#             as a text file with the same tags specifications.
#             Including files can be recursive. >>
#{tags} << Character matrix with two columns indicating the opening and closing
#          tags of the lists at different levels. Its row number gives the maximum
#          number of levels. Opening tags must be different.>>
#{sep} << Character sequence used to split the character vectors of every
# line. Notice that \samp{LF} is always considered as a separator.>>
#{rsep} << Indicates if repetitions of \samp{sep} must be considered as 
#          significant or not; when significant \samp{""} values are introduced.
#          If \samp{no_action} then the repetitions will be ignored.>>
#{stag} << Two character strings indicating the tag to define different \samp{sep} for a given
#           [sub]list. These two correspond to \samp{stag[c(1,3)]} of \samp{list2file} function.>>
#VALUE
# a list [of lists [of lists [...] ] ] of character (possibly named) vectors
# or matrices or arrays.
#EXAMPLE
# sink("rbsa.list.txt")
# cat("# comments can be included as well\n")
# cat("<<A>>\n");
# cat("[[a]]/*/v 1*un deux trois\n");
# cat("[[b]]/*/v 1*2*3\n");
# cat("un uno one\n");
# cat("deux dos two\n");
# cat("trois tres three\n");
# cat("<<B>>\n");
# cat("[[a]] un deux trois\n");
# cat("# the following three are interesting\n");
# cat("[[b]] un  uno  one\n");
# cat(" deux dos two\n");
# cat("trois tres three\n");
# cat("<<C>> 1 2 3\n");
# sink();
# rrrfile2list("rbsa.list.txt");
# unlink("rbsa.list.txt");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 08_12_10
#REVISED 14_06_23
#--------------------------------------------
{
# everything, including the check are deported
# to the called functions
#
# from file to character
res <- rrrfile2text(file,
                 path=path,clean=clean,ended=ended,
                 comment=comment,skip=skip,
                 include=include
                );
#
# from character to list
res <- rrrtext2list(res,                      #
                 tags=tags,sep=sep,
                 rsep=rsep,stag=stag
                );
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3translate <- function(words,translation,
                           left=TRUE,eline=TRUE)
#TITLE translate words
#DESCRIPTION 
# returns a \samp{character()} deduced from
# \samp{words} replacing some of the words
# by the translation proposed into \samp{translation}.
# A component of \samp{words} can comprise several
# common words. 
#DETAILS
# When a word has no translation, it is left as it is
# when \samp{left} is \samp{TRUE}, suppressed if not.
# Initial spacing is not respected, ending blanks are
# suppressed and no repeated blanks are separating the
# final words.
#KEYWORDS 
#INPUTS
#{words} << A \samp{character()} of words .>>
#{translation}<< A named \samp{character} vector containing the
# translations through the names:
# \samp{names(translation)[i]} becomes \samp{translation[i]}.>>
#[INPUTS]
#{left} << When \samp{TRUE} non translated words are left.>>
#{eline} << Indicates if resulting empty lines must be left.>>
#VALUE
# The resulting \samp{character} after translation.
#EXAMPLE
# tra <- c(a="A",e="E",i="I",o="O",u="U")
# rrrtext3translate(letters,tra);
# rrrtext3translate(letters,tra,FALSE);
# rrrtext3translate(letters,tra,FALSE,FALSE);
# rrrtext3translate(rrrrbsa0$text3$v,c(The="xxx",the="xxx"));
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 14_06_14
#REVISED 14_06_14
#--------------------------------------------
{
  # checking
  rrrobject9(words,"character");
  rrrobject9(translation,"character");
  rrrobject9(left,"logical",1);
  rrrobject9(eline,"logical",1);
  # standardizing the words to translate
  words <- rrrtext3preparation(words,"bBR");
  # degenerate cases
  if (length(translation)==0) {
    return(words);
  }
  # more checking
  if (length(unique(names(translation)))!=length(translation)) {
    rrrerreur(translation,"Some words have different translations");
  }
  # translating
  for (ww in rrrbf(words)) {
    www <- words[ww];
    www <- rrrtext3preparation(www,"V");
    for (vv in rrrbf(www)) {
      ou <- which(www[vv]==names(translation));
      if (length(ou)==0) {
        if (!left) {
          www[vv] <- "";
        }
      } else {
        www[vv] <- translation[ou];
      }
    }
    words[ww] <- paste(www,collapse=" ");
  }
  if (eline) {
    words <- rrrtext3preparation(words,"R");
  }
  # returning
  words;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3titre <- function(tit,box="un",
                           lbef=0,sbef=5,
                           saft=0,laft=1,
                           charbox=c("+","|","+","-",
                                     "+","|","+","-",
                                     "."),
                           alig=2,
                           caret=TRUE,
                           imp=TRUE,
                       monitor=rrrrbsa0$monitor$v)
#TITLE  prints or prepares a title
#DESCRIPTION
# prints or prepares the character string \samp{tit} (can comprise
# several components) with more or less emphasis.
#DETAILS
#KEYWORDS print
#INPUTS
#{tit}<<the title to print, each component will be printed in a
#       different line.>>
#[INPUTS]
#{box} <<defines the type of box to put around the title. Possible
#        options are:
#         \cr\samp{"no"} for nothing,
#         \cr\samp{"un"} for underlined,
#         \cr\samp{"ov"} for overlined,
#         \cr\samp{"par"} for parentherized,
#         \cr\samp{"cor"} for the four corners,
#         \cr\samp{"unov"} for underlined and overlined,
#         \cr\samp{"box"} for a complete box.>>
#{lbef} <<if \samp{numeric(1)} defines the number of empty
#         lines before the title; when \samp{character}
#         provides the lines to add before the title.>>
#{sbef} <<either \samp{numeric(1)} or \samp{character(1)}
#         indicating the number of spaces or the characters
#         to introduce before the box.>>
#{saft} <<same as \samp{sbef} but after.>>         
#{laft} <<same as \samp{lbef} but after.>>         
#{charbox} <<nine single characters to define the box starting
#            after noon, and the extra 9th being the filling
#            character for multiple line titles.>>
#{alig} <<The aligment to be done by \samp{rrrform3justify}: 
#         1 for left,2 for center and 3 for right.>>
#{caret} <<Indicates if \samp{\\n} must be added at the end
#          of each created line.>>
#{imp} << Printing is performed and nothing is returned.
#                If FALSE, the character string is returned 
#                (including possible new lines)>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# According to \samp{imp}: nothing when printing is performed,
# a character string 
#EXAMPLE
# rrrform3titre("Some Title");
# rrrform3titre(c("The title","can comprise","several lines"),
#            box="box",lbef=4,laft=2);
# rrrform3titre(c("And the box","can be incomplete", "as well!"),
#            box="cor");
# rrrform3titre(c("The title","can comprise","several lines"),
#            box="box",lbef=4,laft=2,
#            charbox=c("*","+","/","=","*","+","/","="," "));
#REFERENCE
#SEE ALSO
#CALLING {rrrform3repeat}
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_09_24
#REVISED 10_07_25
#--------------------------------------------
{
# constants
#            1    2    3     4     5      6     7
boxes <- c("no","un","ov","par","cor","unov","box");
# checking
if (monitor$chk$v) {
    rrrobject9(caret,"logical",1);
    rrrobject9(tit,"character",-1);
    rrrobject9(box,"character",1);
    if (!rrrbelong9(box,boxes,TRUE)) {
        rrrerreur(list(box,boxes),"'box' must be one of 'boxes'");
    }
    rrrobject9(lbef,c("numeric","character"),c(1,Inf));
    rrrobject9(laft,c("numeric","character"),c(1,Inf));
    rrrobject9(sbef,c("numeric","character"),1);
    rrrobject9(lbef,c("numeric","character"),1);
    rrrobject9(charbox,"character",9);
    if (!all(nchar(charbox)==1)) {
        rrrerreur(charbox,"For the moment only unicharacter in 'charbox'");
    }
    rrrobject9(imp,"logical",1);
}
# the null case
if (rrrvoid9(tit)) { return(character(0));}
# preparing
if (is.numeric(lbef)) { lbef <- rep(" ",max(0,round(lbef[1])));}
if (is.numeric(laft)) { laft <- rep(" ",max(0,round(laft[1])));}
if (is.numeric(sbef)) { sbef <- rrrform3repeat(" ",max(0,round(sbef[1])));}
if (is.numeric(saft)) { saft <- rrrform3repeat(" ",max(0,round(saft[1])));}
#
if (box %in% boxes[c(1,2,4)])   { charbox[c(7,8,1)]   <- NA;}
if (box %in% boxes[c(1,3,4)])   { charbox[c(3,4,5)]   <- NA;}
if (box %in% boxes[c(1,2,3,6)]) { charbox[c(1:3,5:7)] <- "";}
if (box %in% boxes[c(5)])       { charbox[c(4,8,2,6)] <- " ";}
#
lmax <- max(sapply(tit,nchar));
tit <- rrrform3justify(tit,lmax,format=alig,carac=charbox[9]);
# producing
res <- character(0);
# first lines
for (ii in rrrbf(lbef)) { res <- c(res,lbef[ii]);}
# overline
if (box %in% boxes[c(3,5,6,7)]) {
    lili <- sbef;
    lili <- paste0(lili,charbox[7]);
    lili <- paste0(lili,paste(rep(charbox[8],lmax),collapse=""));
    lili <- paste0(lili,charbox[1]);
    lili <- paste0(lili,saft);
    res <- c(res,lili);
}
# title lines
for (ii in rrrbf(tit)) {
    lili <- sbef;
    lili <- paste0(lili,charbox[6]);
    lili <- paste0(lili,tit[ii]);
    lili <- paste0(lili,charbox[2]);
    lili <- paste0(lili,saft);
    res <- c(res,lili);
}
# underline
if (box %in% boxes[c(2,5,6,7)]) {
    lili <- sbef;
    lili <- paste0(lili,charbox[5]);
    lili <- paste0(lili,paste(rep(charbox[4],lmax),collapse=""));
    lili <- paste0(lili,charbox[3]);
    lili <- paste0(lili,saft);
    res <- c(res,lili);
}
# last lines
for (ii in rrrbf(laft)) { res <- c(res,laft[ii]);}
# adding carriage returns
if (caret) { for (ii in rrrbf(res)) {
    res[ii] <- paste0(res[ii],"\n");
}}
# returning
if (imp) {
    cat(res,sep="");
    invisible();
} else { return(res);}
# 
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3display <- function(x,pau=FALSE,cat=FALSE,...)
#TITLE  displays with its name any object
#DESCRIPTION
#  displays any object after giving the name of
# the variable containing it. A rrrpause can be introduced
# to give the opportunity to scrutinize the result.  
#DETAILS
#KEYWORDS print
#INPUTS
#{x}<<The object to print.>>
#[INPUTS]
#{pau} << Must a rrrpause be performed after the display?>>
#{cat} << Must the printing be done with 'cat' instead of print?>>
#{\dots} <<possible arguments for the print function.>>
#VALUE
# a print (or cat) is done and \samp{x} is returned
#EXAMPLE
# uu <- "azerty";
# rrrform3display(uu);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 08_09_17
#REVISED 09_10_12
#--------------------------------------------
{
cat("<< Displaying ",deparse(substitute(x))," >>\n");
if (cat) { cat(x,"\n");
} else { print(x,...);}
if (pau) { rrrpause("affichage");}
# returning
x;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext2list <- function(text,
                      tags=rrrrbsa0$tag1$v,
                      sep=rrrrbsa0$sep0$v,
                      rsep="no_action",
                      stag=c("/","/"),
                      comment="#",
                      monitor=rrrrbsa0$monitor$v
                     )
#TITLE  transforms a character into a list (of lists) of characters, and conversely
#ALIAS
#DESCRIPTION
# from a conveniently character vector comprising tags, returns nested lists.
# Most often such a character is obtained by reading a file with
# the function \samp{rrrfile2text}, this explains why it is convenient
# to speak about 'lines' rather about 'components of the character
# vector'.\cr
# All lists are named lists and the tags give the names of their
# component. The maximum number of nested levels of lists is given by the
# number of rows of the matrix \samp{tags}. Its corresponding two columns
# providing the opening and closing sequences of the tags. Final lists 
# contain \samp{character} vectors, each component of them being
# on the same line and/or on the following line (before a new tag).\cr
# All tags must start at the very beginning of a line. All separator
# tags must be used sticked to the list tag.\cr
# Lines starting with a \samp{comment} are first eliminated.
#DETAILS
# It is compulsory to tag each level of the lists, this implies that 
# when starting a new list, a character vma (see \samp{rrrtext2vma}
# for the details: vma means vector or matrix or array) is provided meaning 
# that this is the final level for this branch, or a new sublist 
# is started with the corresponding tag, or a new component of the
# list is given at a level less or equal to the present.\cr
# Separator between different character components is given by the
# \samp{sep} argument or indicated by the separator tag (\samp{stag}), 
# in the last case, it can be different from a leaf list to another.\cr
# Be aware that before reading a line for a character translation,
# all starting and ending spaces are eliminated.\cr
#KEYWORDS IO
#INPUTS
#{text} << character vector to be transformed into a list.>>
#[INPUTS]
#{tags} << Character matrix with two columns indicating the opening and closing
#          tags of the lists at different levels. Its row numbers gives the maximum
#          number of levels. Opening tags must be different.>>
#{sep} << Character sequence used to split the character vectors of every
# line. Notice that \samp{LF} is always considered as a separator.>>
#{rsep} << Indicates if repetitions of \samp{sep} must be considered as 
#          significant or not and which null value to introduce.
#          If \samp{no_action} then the repetitions will be ignored
#          if not \samp{rsep} component(s) will be introduced.>>
#{stag} << Two character strings indicating the tag to define different \samp{sep} for a given
#           [sub]list. These two correspond to \samp{stag[c(1,3)]} of \samp{list2file} function.>>
#{comment} << Set a characters indicating that the line is a comment line when found
#             in the first positions.>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# a list [of lists [of lists [...] ] ] of character (possibly named) vectors
# or matrices or arrays.
#EXAMPLE
# text <- c("<<A>>",
#            "[[a]]/*/v 1*un deux trois",
#            "[[b]]/*/v 1*2*3",
#            "un uno one",
#            "deux dos two",
#            "trois tres three",
#            "<<B>>",
#            "[[a]] un deux trois",
#            "[[b]] un  uno  one",
#            " deux dos two",
#            "trois tres three",
#            "<<C>> 1 2 3");
# rrrtext2list(text);
# rrrlist2text(rrrtext2list(text));
#REFERENCE
#SEE ALSO rrrtext2file rrrfile2list
#CALLING
#COMMENT
#FUTURE Extend the possibility of the list by allowing
#       contents at the intermediate levels (not only
#       at the final levels as it is presently).
#AUTHOR J.-B. Denis
#CREATED 08_12_10
#REVISED 14_08_02
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(text,"character",-1,mensaje="text must be a character");
    #
    rrrobject9(tags,"character",c(2,Inf),mensaje="tags must be a CHARACTER matrix with two columns");
    if (!is.matrix(tags)) {rrrerreur(tags,"tags must be a character MATRIX of two columns");}
    if (ncol(tags)!=2) {rrrerreur(tags,"tags must be a character matrix OF TWO COLUMNS");}
    if (length(unique(tags[,1]))!= nrow(tags)) { rrrerreur(tags,"Opening tags are not unique!");}
    #
    rrrobject9(sep,"character",1,mensaje="sep must indicate one string character");
    if (nchar(sep) == 0) { rrrerreur(sep,"sep must not be an empty string");}
    #
    rrrobject9(rsep,"character",1,mensaje="rsep must be a character(1)");
    #
    rrrobject9(stag,"character",2,mensaje="stag must be a character of length two");
    if (any(nchar(stag) == 0)) { rrrerreur(stag,"stag must not have zero length components");}
    if (length(grep(" ",stag)) > 0) { rrrerreur(stag,"stag must not comprise spaces");}
    #
    rrrobject9(comment,"character",1,mensaje="'comment' must be a character of length one");
}
#
# removing the comment lines
if (!rrrvoid9(comment)) {
  text <- rrrfilter8text(text,comment,remove=TRUE,
                         exact=8,lower=FALSE,
                         monitor=monitor);
}
#
# getting the tagged lines 
tagged <- matrix(NA,0,6);
dimnames(tagged) <- list(NULL,c("num","lev","name","leaf","sep","type"));
#
# tagged will contain all useful informations about the tags.
#    num: line number into text
#    lev: level number of the tag
#    name: name of the tag
#    leaf: "yes" when final branch, "no" if not. The last
#          two columns will be significant only for leaf sub-list.
#    sep: separator string to be used to get the vma character
#    type: indicates the type of the associated sub-sub-...-sub list
#          we are dealing with : "lists" indicates that more
#          levels are present, if not it is a leaf of the structure
#          and \samp{rrrrbsa0$vma$v["c","v","V","w","m",...]} one of the type of vma character
#          is expected. 
nbtag <- 0;
for (i in rrrbf(text)) {
    lev <- 0;
    # looking for a tag
    for (j in rrrbc(nrow(tags))) { if (lev == 0) {
        pot <- rrrform3crop(text[i],tags[j,1],"");
        if (nchar(pot) < nchar(text[i])) { lev <- j;}
    }}
    if (lev > 0) {
        # some tag was perhaps identified looking for the name
        put <- strsplit(pot,tags[lev,2],fixed=TRUE)[[1]];
        if (nchar(put[1]) < nchar(pot)) {
            # indeed it was a tag
            fch <- substr(put[2],1,1);
            if (fch %in% c(""," ",NA)) {
                # default tags and separator
                ssep <- sep;
                typing <- rrrrbsa0$vma$v["v"];
                text[i] <- put[2];
            } else {
                # getting the possible explicit separator
                if (substr(put[2],1,nchar(stag[1]))==stag[1]) {
                    # looking for the separator
                    pit <- rrrform3crop(put[2],stag[1]);
                    pat <- strsplit(pit,stag[2],fixed=TRUE)[[1]];
                    ssep <- pat[1];
                    pyt <- pat[2];
                } else {
                    # default separator
                    ssep <- sep;
                    pyt <- put[2];
                }
                # getting the possible information
                ouou <- strsplit(pyt," ",fixed=TRUE)[[1]];
                if (length(ouou) > 1) {
                    text[i] <- paste(ouou[-1],collapse=" ");
                } else {
                    text[i] <- "";
                }
                pyt <- ouou[1];
                # getting the type of the tag
                tta <- which(pyt==rrrrbsa0$vma$v);
                if (length(tta)==0) {
                    # default type
                    typing <- rrrrbsa0$vma$v["v"];
                } else {
                    # explicite type
                    typing <- names(rrrrbsa0$vma$v)[tta];
                }
            }
        tagged <- rbind(tagged,c(i,lev,put[1],"???",ssep,typing));
        } # end of nchar(put[1]) < nchar(pot)
    } # end of lev > 0
} # end of i in rrrbf(text)
#
#
# in case of no list
if (nrow(tagged) == 0) { return(vector("list",0));}
# detecting the final lists
tagged[,"leaf"] <- "yes";
for (ii in rrrbc(nrow(tagged)-1)) {
    lev1 <- as.numeric(tagged[  ii,"lev"]);
    lev2 <- as.numeric(tagged[ii+1,"lev"]);
    if (lev2-lev1==1) { tagged[ii,"leaf"] <- "no";}
}
#
#   matrix 'tagged' is rrrnow available
#
# checking the sequential progress of the levels
for (ii in rrrbc(nrow(tagged)-1)) {
    lev1 <- as.numeric(tagged[  ii,"lev"]);
    lev2 <- as.numeric(tagged[ii+1,"lev"]);
    if ((lev2-lev1) > 1) {
        rrrerreur(tagged,"The progression of levels is not accepted: all levels must be introduced");
    }
}
# checking that no leaf tags are followed by another tag
for (ii in rrrbc(nrow(tagged)-1)) { if (tagged[ii,"leaf"]=="no") {
    actuel <- as.numeric(tagged[  ii,"num"]);
    suivan <- as.numeric(tagged[ii+1,"num"]);
    if ((suivan-actuel) != 1) {
        rrrerreur(list(tagged,ii),"In line 'ii', a 'no-leaf' tags is not immediately followed by another tag (only terminal levels can have contents)");
    }
}}
#
if (tagged[1,"lev"] != "1") {
    rrrerreur(tagged,"The first level must be 1");
}
# replacing the introduced NA
for (ii in rrrbf(text)) { if (is.na(text[ii])) { text[ii] <- "";}}
#
#
# constructing the character strings associated to each leaf list
#
#
oul1 <- which(tagged[,"leaf"]=="yes");
oul2 <- as.numeric(tagged[oul1,"num"]);
lulu <- character(length(oul1));
oul2 <- c(oul2,length(text)+1);
for (ss in rrrbf(lulu)) {
    leaf <- oul1[ss];
    ssep <- tagged[leaf,"sep"];
    for (sss in oul2[ss]:(oul2[ss+1]-1)) {
        nou <- rrrform3crop(text[sss]," ");
        if (nchar(nou)>0) {
            if (nchar(lulu[ss])>0) {  
                # here multiple lines are concatenated
                lulu[ss] <- paste(lulu[ss],nou,sep=ssep);
            } else {
                lulu[ss] <- nou;
            }
        }
    }
}
#
# building the list of lists...
res <- vector("list",0);
vari <- "res"; niv <- 0; ulul <- 0;
for (ll in rrrbc(nrow(tagged))) {
    # getting the name
    leve <- as.numeric(tagged[ll,"lev"]);
    if (niv<leve) {
        if (niv+1!=leve) { stop("Error(1) found into 'file2test'");}
        vari <- paste0(vari,"[['",tagged[ll,"name"],"']]");
    }
    if (niv>=leve) {
        va <- strsplit(vari,"[['",fixed=TRUE)[[1]];
        vari <- paste(va[1:leve],collapse="[['");
        nom <- tagged[ll,"name"];
        #???if (nom == "<NA>") { nom <- NULL;}
        vari <- paste0(vari,"[['",nom,"']]");
    }
    niv <- leve;
    # filling the value
    if (tagged[ll,"leaf"]=="yes") {
        ulul <- ulul+1;
        jsep <- tagged[ll,"sep"];
        # getting the character vector
        choc <- strsplit(lulu[ulul],jsep,fixed=TRUE)[[1]];
	if (rsep=="no_action") {
	    choc <- choc[choc!=""];
	}
        # building the code to interpret the character vector
        chacha <- rrrtext2vma(choc,tagged[ll,"type"],rrrrbsa0$sep1$v,monitor=monitor);
        eval(parse(text=paste(vari,"<- chacha")));
    } else {
        eval(parse(text=paste(vari,"<- vector('list',0)")));
    }
}
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtexts4text <- function(text,where,
                       addbeg=TRUE,addend=TRUE,rmeli=TRUE)
#TITLE returns a list of texts from a text
#DESCRIPTION from a text, 
# gives back a list of texts according to the delimitations
# proposed into \samp{where}; typically \samp{where} can be deduced from 
# output of \samp{rrrplaces4text6tags}. Empty lines can be removed.\cr
# The idea is from a text comprising several lines, to get pieces of it
# after eliminating the delimitating tags; when the tags have already been
# identified since they are indicated through their positions with
# the argument \samp{where}. A delimitation is given in each row of \samp{where},
# so according to the values of \samp{addbeg} and \samp{addend} the maximum
# length of the returned list is \samp{nrow(where)-1}, \samp{nrow(where)}
# or \samp{nrow(where)+1}.
#DETAILS
# When \samp{nrow(where)==0} an empty list is returned without
# warning.\cr
# The component number \samp{i} of the resulting list, is the portion of \samp{text}
# in between the two positions indicated by rows \samp{i} and \samp{i+1} of
# \samp{where} of course, it can comprise several lines, that is be a
# character vector.\cr
# The two positions defining a delimitation are included (i.e. second and third columns
# of the \samp{where} matrix) in the tag. Nevertheless a split can be obtained at
# position \samp{7} in row \samp{3} by the following row in matrix \samp{where}:
# \samp{c(3,7,6)}.
# When delimitations are not consistent, a fatal error is issued.\cr
#KEYWORDS 
#INPUTS
#{text}<< A named \samp{character} vector containing the text
# (a component, a line).>>
#{where}<< A two or three column matrix (when a two column matrix, the second column
# is first duplicated as third column and substracted with 1 to generate a simple splitting
# at this point). The first column gives the line number;
# the second column, the first character of the position (excluded in the previous text
# component); the third column, the last character of the position (excluded in the
# next text component). Overlapping of the extraction is not admitted so the difference
# between the second and third column can be greater than one.>>
#[INPUTS]
#{addbeg} <<Add (or not) what is before the first position?>>
#{addend} <<Add (or not) what is after the last position?>>
#{rmeli} <<Must empty lines of the found texts (\samp{""}) must be eliminated?>>
#VALUE
# A list of \samp{character}s each corresponding to an extracted text.
#EXAMPLE
# uu <- paste(letters,collapse=""); UU <- toupper(uu);
# ww <- matrix(c(1,2,2,1,12,13,2,5,10),ncol=3,byrow=TRUE);
# rrrtexts4text(c(uu),ww[1,,drop=FALSE]);
# rrrtexts4text(c(uu),c(1,7,6));
# rrrtexts4text(c(uu),ww[1:2,]);
# rrrtexts4text(c(uu,UU),ww);
# rrrtexts4text(c(uu,"",UU),ww);
# rrrtexts4text(c(uu,"",UU),ww,rmeli=FALSE);
# rrrtexts4text(c(uu),ww);
# rrrtexts4text(c(uu),ww[1,,drop=FALSE],addbeg=FALSE);
# rrrtexts4text(c(uu),ww[1,,drop=FALSE],addbeg=FALSE,addend=FALSE);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_07_10
#REVISED 14_06_13
#--------------------------------------------
{
  # transforming
  if (is.data.frame(where)) {
    where <- as.matrix(where);
  }
  if (!is.matrix(where)) { if (length(where)%in%2:3) {
    where <- matrix(where,ncol=length(where));
  }}
  # checking
  if (length(rrrobject9(where,c("matrix","data.frame"))) > 0) {
    rrrerreur(where,"'where' must be a matix!");
  }
  if ((ncol(where) < 2) | (ncol(where)>3)) {
    rrrerreur(where,"'where' must have two or three columns");
  }
  if (nrow(where) < 1) {
    rrrerreur(where,"'where' must have at least one row");
  }
  # completing with a third column
  if (ncol(where)==2) {
    where <- cbind(where,where[,2]-1);
  }
  # further check
  if (any(where[,3] < where[,2]-1)) {
    rrrerreur(where,"Third column cannot be less that the second one!");
  }
  # dimnaming if necessary
  if (is.null(dimnames(where))) {
    dimnames(where) <- list(paste0("r.",rrrbc(nrow(where))),
                            paste0("c.",rrrbc(ncol(where))));
  } else {
    if (is.null(dimnames(where)[[1]])) {
      dimnames(where)[[1]] <- paste0("r.",rrrbc(nrow(where)));
    }
  }
  # possibly adding the ends
  if (addbeg) {
    where <- rbind(c(1,0,0),where);
    dimnames(where)[[1]][1] <- "<antes>";
  }
  if (addend) {
    ll <- length(text);
    lc <- nchar(text[ll]);
    where <- rbind(where,c(ll,lc+1,lc+1));
    dimnames(where)[[1]][nrow(where)] <- "<despues>";
  }
  # eliminating irrelevant rows
  where <- where[where[,1]>=1,,drop=FALSE];
  where <- where[where[,1]<=length(text),,drop=FALSE];
  #where <- where[where[,2]<1,,drop=FALSE];
  #where <- where[where[,2]>nchar(text[where[,1]]),,drop=FALSE];
  #where <- where[where[,3]<1,,drop=FALSE];
  #where <- where[where[,3]>nchar(text[where[,1]]),,drop=FALSE];
  if (nrow(where) < 1) {
    rrrerreur(where,"'where' must have at least one valid row");
  }
  # sorting
  rrr <- rrrtext3ij2n(where[,1:2,drop=FALSE],text);
  where <- where[order(rrr),,drop=FALSE];
  # initializing
  nbt <- nrow(where)-1;
  res <- vector("list",nbt);
  # looking for each text in turn
  for (itt in rrrbc(nbt)) {
    bb <- where[itt   ,1:3];
    be <- where[itt+1,1:3];
    afaire <- TRUE;
    if (bb[1]  > be[1]) { afaire <- FALSE;}
    if (bb[1] == be[1]) {
      if (bb[3]>be[2]) { afaire <- FALSE;}
    }
    if (afaire) {
      if (bb[1] == be[1]) {
        res[[itt]] <- substr(text[bb[1]],bb[3]+1,be[2]-1);
      } else {
        res[[itt]] <- substr(text[bb[1]],bb[3]+1,nchar(text[bb[1]]));
        for (iii in rrrbd(bb[1]+1,be[1]-1)) {
          res[[itt]] <- c(res[[itt]],text[iii]);
        }
        res[[itt]] <- c(res[[itt]],
                    substr(text[be[1]],      1,be[2]-1));
      }
    } else {
      res[[itt]] <- character(0);
    }
    # possibly removing empty lines
    if (rmeli) {
      for (ll in rev(rrrbf(res[[itt]]))) {
        if (nchar(res[[itt]][ll])==0) {
          res[[itt]] <- res[[itt]][-ll];
        }
      }
    }
  }
  # naming
  names(res) <- dimnames(where)[[1]][-(nrow(where))];
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3ij2n <- function(ij,text,flexible=TRUE)
#TITLE gives the absolute positions of (row,column) positions within a /text/
#DESCRIPTION 
# returns numerical values equal to the rank of a position when all 
# lines of the \samp{text} are concatenated without separators.
#DETAILS
# Definition of the positions is flexible (outside positions
# is interpreted as minimum / maximum positions) or not, according to \samp{flexible}.
#KEYWORDS 
#INPUTS
#{ij} << A \samp{numeric(2)} indicating (line,column) of the position.
#          May also be a matrix with two columns.>>
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#[INPUTS]
#{flexible} << When \samp{TRUE} values of \samp{n} outside the acceptable
#              range are accepted and replaced with first or last positions.>>
#VALUE
# The resulting ranks associated to the rows of \samp{ij}.
#EXAMPLE
# aa <- c(paste(letters,collapse=""),paste(LETTERS,collapse=""),paste(0:9,collapse=""));
# rrrtext3ij2n(c(1,3),aa);
# rrrtext3ij2n(c(2,1),aa);
# rrrtext3ij2n(rbind(c(1,3),c(2,1)),aa);
#ALIAS text-indexing
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_17
#REVISED 13_10_29
#--------------------------------------------
{
  # checking and normalizing
  if (is.matrix(ij)) {
    if (ncol(ij) != 2) {
      stop("'ij' must have two columns");
      rrrobject9(ij,"numeric");
    }
  } else {
    rrrobject9(ij,"numeric",len=2);
    ij <- matrix(ij,ncol=2);
  }
  rrrobject9(text,"character");
  # degenerate case
  if (nrow(ij)==0) { return(numeric(0));}
  # possibly adjusting
  if (flexible) {
    i1 <- pmin(length(text),pmax(1,ij[,1]));
    j1 <- pmin(nchar(text[i1]),pmax(1,ij[,2]));
  } else {
    rrrobject9(ij[,1],"numeric",con=c(1,length(text)));
    if ((any(ij[,2] < 1)) | (any(ij[,2]>nchar(text[ij[,1]])))) {
      rrrerreur(cbind(ij,nchar(text[ij[,1]])),"Bad column numbering");
    }
  }
  # computing 
  csu <- c(0,cumsum(nchar(text)));
  # returning
  csu[i1] + j1;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrpause <- function(what="",mensaje=NULL,top=0,ans=NULL) 
#TITLE  pauses the program until an answer is given
#DESCRIPTION
# This function issues a rrrpause with a message allowing 
# to stop the process or to continue it (without or with and 
# according to the answer given by the user).
#DETAILS
# The answer provided by the user is interpreted as
#  a \samp{character(1)}.\cr
# When \samp{top} is \samp{NA} the call is ignored.\cr
# When \samp{top} is \samp{NULL} no stopping is possible
#      but the answer is returned. This answer can be
#      anything when \samp{is.na(ans)} but a void answer
#      is refused. This answer can be anything including
#      nothing when \samp{is.null(ans)}. The answer is
#      returned.\cr
# When \samp{top} is \samp{numeric}, \samp{ans} is not considered.
#       an empty answer allows the process to follow
#       and no void answers stop the process.\cr
# When \samp{top} and \samp{ans} are \samp{character}s, 
#      the only accepted answers are a component of \samp{top} (to 
#      stop the process) or a component of \samp{ans} (to be returned).\cr
# When \samp{top} is a \samp{character} and \samp{is.null(ans) | is.na(ans)}
#      stopping is made when the answer is one of the components of
#      \samp{top}, other answers will be returned.\cr
#KEYWORDS IO
#INPUTS
#[INPUTS]
#{what} << Short message commenting the rrrpause>>
#{mensaje} << A possible longer message to give details. When
#          \samp{NULL} a shortened message is issued.>>
#{top} << Precises the rules to follow for the execution to be stopped
#         (see the \samp{details} section).>>
#{ans} << Precises the accepted answers
#         (see the \samp{details} section).>>
#VALUE
# When the process is not stopped, returns the answer provided
# by the user.
#EXAMPLE
# rrrpause("Time for lunch!",top=NA);
# \dontrun{rrrpause("Time for lunch?")}
# \dontrun{rrrpause("Look at your results before the process continue",top=NULL)}
# \dontrun{rrrpause("Can we continue?",top=c("n","N","no","NO"),ans=c("y","Y","yes","YES"))}
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 07_10_25
#REVISED 14_07_31
#--------------------------------------------
{
  # no action
  if(length(top)==1) { if (is.na(top)) { return(invisible());}}
  # writing the question
  if (!rrrvoid9(mensaje)) {cat(">>> (",mensaje,")\n");}
  if (!rrrvoid9(what))  {cat(">>> --------------> ",what,"\n");}
  # precising the possibilities and getting the answer
  if (is.numeric(top)) {
    # raw rrrpause
    cat(">>> 'Enter' to continue | any key(s) +'Enter' to stop \n");
    quoi <- scan(what="character",nmax=1);
    if (length(quoi) != 0) {
      stop("(...YOU decided to stop the job...)",call.=FALSE);
    } else {
      return(invisible());
    }
  } else {
    if (is.null(top)) {
      # no stopping is possible
      if (is.null(ans)) {
        # any answer is possible
        cat(">>> Give your answer and 'Enter' \n");
        quoi <- scan(what="character",nmax=1);
        return(quoi);
      } else {
        if (is.na(ans[1])) {
          # any non void answer is possible
          cat(">>> Give a non-empty answer and 'Enter' \n");
          quoi <- scan(what="character",nmax=1);
          while (length(quoi) == 0) {
            cat("EMPTY ANSWER is not admitted\n");
            quoi <- scan(what="character",nmax=1);
          }
          return(quoi);
        } else {
          # a restricted answer is required
          cat(">>> Give an answer among the possible and 'Enter' \n");
          cat(">>> Possible are: /",paste(ans,collapse="/"),"/\n",sep="");
          quoi <- scan(what="character",nmax=1);
          while((length(quoi)==0) || (!(quoi %in% as.character(ans)))) {
            cat("POSSIBLE ANSWERS are in the following list:\n");
            cat(">>> /",paste(ans,collapse="/"),"/\n",sep="");
            quoi <- scan(what="character",nmax=1);
          }
          return(quoi);
        }
      }     
    } else {
      # possibly stopping
      if (is.null(ans)) {
        # any answer is possible
        cat(">>> Give an answer and 'Enter' \n");
        cat(">>> To stop: /",paste(top,collapse="/"),"/\n",sep="");
        quoi <- scan(what="character",nmax=1);
        if ((length(quoi)>0) && (quoi %in% top)) {
          stop("(...YOU decided to stop the job...)",call.=FALSE);
        } else {
          return(quoi);
        }
      } else {
        if (is.na(ans[1])) {
          # any non void answer is possible
          cat(">>> Give a non empty answer and 'Enter' \n");
          cat(">>> To stop: /",paste(top,collapse="/"),"/\n",sep="");
          quoi <- scan(what="character",nmax=1);
          while (length(quoi) == 0) {
            cat("EMPTY ANSWER is not admitted\n");
            quoi <- scan(what="character",nmax=1);
          }
          if (quoi %in% top) {
            stop("(...YOU decided to stop the job...)",call.=FALSE);
          }
          return(quoi);
        } else {
          # a restricted answer is required
          topans <- union(top,ans);
          if (length(topans)!=(length(ans)+length(top))) {
            stop("Possible answers and stopping answers are not disjoint!");
          }
          cat(">>> Give an answer among the possible and 'Enter' \n");
          cat(">>> To stop: /",paste(top,collapse="/"),"/ ",sep="");
          cat("or: /",paste(ans,collapse="/"),"/\n",sep="");
          quoi <- scan(what="character",nmax=1);
          while((length(quoi)==0) || (!(quoi %in% as.character(topans)))) {
            cat("POSSIBLE ANSWERS are in the following list:\n");
            cat(">>> To stop: /",paste(top,collapse="/"),"/ ",sep="");
            cat("or: /",paste(ans,collapse="/"),"/\n",sep="");
            quoi <- scan(what="character",nmax=1);
          }
          if (quoi %in% top) {
            stop("(...YOU decided to stop the job...)",call.=FALSE);
          }
          return(quoi);
        }
      }
    }
  }
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrfilter8text <- function(text,pattern,remove=TRUE,
                         exact=FALSE,how="a",lower=FALSE,
                         monitor=rrrrbsa0$monitor$v)
#TITLE filters components of a character.
#DESCRIPTION
# For each component of a character vector, checks with function 
# \samp{rrrbelong9} if it satisfies the conditions given by the other
# arguments, the component is removed (or selected).
#DETAILS
#KEYWORDS utilities
#INPUTS
#{text} << The \samp{character} to be filtered.>>
#{pattern} <<(\samp{character}) the character string(s) to be found
#             in each component.>>
#[INPUTS]
#{remove} << Must the conforming components be removed (or selected)?>>
#{exact} << When exact, one component must
# be strictly identical, if not a subtring is sufficient.>>
#{how} << Indicates what to do when \samp{length(sch)>1}. The choice are 
# \samp{v}: a logical vector gives back each check independently;
# \samp{1}: returns \samp{TRUE} when at least one of the component belongs
# to the series \samp{ch} and \samp{a} when all components must comply to get TRUE.>>
#{lower} << Must the comparisons being done after case lowering?>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE A \samp{character} comprising the non-removed (or selected) components
#      of the initial components.
#EXAMPLE
# rrrfilter8text(rrrrbsa0$text2$v,"little",exact=FALSE,lower=TRUE,remove=FALSE)
# rrrfilter8text(rrrrbsa0$text2$v,"On",exact=8,lower=TRUE,remove=FALSE)
# rrrfilter8text(rrrrbsa0$text2$v,"On",        lower=TRUE,remove=FALSE)
# rrrfilter8text(rrrrbsa0$text2$v,"On",exact=8,lower=TRUE,remove=TRUE)
#REFERENCE
#SEE ALSO rrrbelong9
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 04_08_03
#REVISED 14_08_03
#--------------------------------------------
{
  # first checking
  if (monitor$chk$v) {
    rrrobject9(remove,"logical",1);
    rrrobject9(text,"character");
  }
  # degenerate cases
  if (length(pattern)==0) {
    if (remove) {
      return(text);
    } else {
      return(character(0));
    }
  }
  if (length(text)==0) { return(text);}
  # second checking
  if (monitor$chk$v) {
    # other checking are made within 'rrrbelong9'
  }
  # getting the conforming components
  coco <- vector("logical",length(text));
  for (ii in rrrbf(text)) {
    coco[ii] <- rrrbelong9(pattern,text[ii],
                        exact=exact,how=how,lower=lower,
                        monitor=monitor);
  }
  if (remove) {
    coco <- (!coco);
  }
  # returning
  text[coco];
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3repeat <- function(cha="-",nb=10,imp=FALSE,cr=imp)
#TITLE  prints a repeated given string
#DESCRIPTION
# Component by component,  prints the concatenation
# of the given string(s) repeated \samp{nb} times.
#DETAILS
#KEYWORDS print
#INPUTS
#[INPUTS]
#{cha} << The string to repeat>> 
#{nb} << Number of repetitions>> 
#{imp} << Printing when TRUE or returning (default)>>
#{cr} << Must a line feed be added?>>
#VALUE
# character string or printing according to \samp{imp}
#EXAMPLE
# rrrform3repeat('-+',20,TRUE)
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
for (jbd in rrrbc(nb)) { res <- paste0(res,cha);}
if (cr) { res <- paste0(res,"\n");}
if(!imp) { return(res);}
cat(res);
invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3justify <- function(chaine,
                          nbc=8,
                          format=3,
                          tronc=TRUE,
                          carac=" ",
                       monitor=rrrrbsa0$monitor$v)
#TITLE  formats a character string
#DESCRIPTION
# Formats character string(s).
# The main use of this function is to produce
# aligned columns lists while not printing
# them at the same time.
#DETAILS
# The justification is made component by component,
# no concatenation between them is done.
#KEYWORDS print
#INPUTS
#{chaine}<<the character string to be printed, can be a vector.>>
#[INPUTS]
#{nbc} << Desired number of characters for the result; when
#         \samp{chain} is a vector can be a vector of the same length>>
#{format} << Indicates the type of alignment:\cr
#   0 no aligment (no supplementary character added to reach \samp{nbc})\cr
#   1 to the left side\cr
#   2 centered\cr
#   3 to the right side>>
#{tronc} << If true, no more than
#     \samp{nbc} characters are returned and
# possibly the string is truncated. In that
# case, \samp{$} is introduced to indicate the fact.>>
#{carac} << Character to use for enlarging the string>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE a character string
#EXAMPLE
# rrrform3justify("vers")
# rrrform3justify("versification",5)
# rrrform3justify(letters[1:5],format=2,carac="+");
#REFERENCE
#SEE ALSO rrrtext3preparation
#CALLING
#COMMENT
#FUTURE Improve the behavior when \samp{nchar(carac)>1}.
#AUTHOR J.-B. Denis
#CREATED 1999_05_25
#REVISED   14_08_27
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(chaine,c("character","numeric"),-1);
    rrrobject9(nbc,"integer",c(1,Inf));
}
# the null case
if (length(chaine)==0) { return(character(0));}
# preparing
nbc[nbc < 1] <- 8;
if (length(nbc) < length(chaine)) {
    nnbc <- rep(nbc,length(chaine));
} else {
    nnbc <- nbc;
}
#
itr <- "$"; # truncation indicator
rres <- cchaine <- chaine;
for (rr in rrrbf(rres)) {
    res <- rres[rr];
    nbc <- nnbc[rr];
    chaine <- cchaine[rr];
    # truncation
    if ( (nchar(res) > nbc) & tronc ) {
     if (format <= 1) {
      res <- substring(chaine,1,nbc-1);
      res <- paste0(res,itr);
      }
     else {
      if (format == 2) {
       otg <- (nchar(chaine) - nbc) %/% 2;
       res <- substring(chaine,1+otg);
       res <- substring(res,1,nbc-2);
       res <- paste0(itr,res,itr);
       }
      else {
       res <- substring(chaine,1+nchar(chaine)-nbc+1,
			nchar(chaine));
       res <- paste0(itr,res);
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
rrrlist4text <- function(text,content=c(1,length(text)),ccom="#",
                      lndel=c("{","}"),lvdel=c("<<",">>"))
#TITLE Analyses and prepares a list from a /text/ 
#DESCRIPTION from a text with convenient tags
# returns a named list. Labelling tags must be 
# in first position of the lines, possibly after removing
# some precised character.
#DETAILS
# Useless lines can exist before and after the list, and also
# between the list items.
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# from which to extract the list.>>
#[INPUTS]
#{content} << Indicates which components of \samp{text} to consider.
#             Usually \samp{numeric(2)} to indicate the interval of lines to consider.
#             When \samp{numeric(1)}, only this line. When \samp{0} or \samp{diff(content)<0} 
#             a \samp{list(0)} is returned.>>
#{ccom} <<Character(s) at the beginning of lines to possibly be removed.
#         When \samp{character(0)} no removing is performed.>>
#{lndel} <<Starting and ending delimiters for the name of the item (at the beginning of the line).>>
#{lvdel} <<Starting and ending delimiters for the description of the item.>>
#VALUE
# A named list.
#EXAMPLE
# rrrlist4text(c("{a}<<Il etait une fois>>",
#                "{b}<<   un petit et rouge chaperon>>{c}",
#                "<< pas bon>>",
#                "{d}    <<et le m'echant loup...>>"))
# rrrlist4text(c("{a}","<<","pour voir",">>"));
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_10_14
#REVISED 14_06_10
#--------------------------------------------
{
  # initializing
  res <- vector("list",0);
  # checking
  rrrobject9(content,
             "numeric",len=1:2,con=c(0,Inf),
             mensaje="Bad 'content' argument (1)");
  if (content[1]==0)   { return(res);}
  if (diff(content)<0) { return(res);}
  if (length(content)!=2) { content <- rep(content[1],2);}
  rrrobject9(content,"numeric",len=2,con=c(1,length(text)),
               mensaje="Bad 'content' argument (2)");
  #
  rrrobject9(ccom,
             "character",len=0:1,
             mensaje="Bad 'ccom' argument");
  # restricting the treatment
  text <- text[content[1]:content[2]];
  ## removing starting tags
  if (length(ccom)==1) {
    nbc <- nchar(ccom);
    if (nbc>0) {
      for (tt in rrrbf(text)) {
        tex <- text[tt];
        if (substr(tex,1,nbc)==ccom) {
          text[tt] <- substr(tex,nbc+1,nchar(tex));
        }
      }
    }
  }
  #
  rrrobject9(lndel,
             "character",len=2,
             mensaje="Bad 'lndel' argument");
  if (any(lndel=="")) {
    rrrerreur(lndel,"bad labelling tags");
  }
  rrrobject9(lvdel,
             "character",len=2,
             mensaje="Bad 'lvdel' argument");
  if (any(lvdel=="")) {
    rrrerreur(lndel,"bad valuing tags");
  }
  # extracting the list
  #
  # localizing the list delimiters
  oudon <- rrrplaces4text6tags(text,tags=c(lndel[1],lndel[2],
                                        lvdel[1],lvdel[2]));
  # getting the lines of item names
  ou1 <- which((oudon[,"tag"] == lndel[1]) & (oudon[,"co1"] == 1));
  ou2 <- which(oudon[,"tag"] == lndel[2]);
  oun <- intersect(oudon[ou1,"li1"],oudon[ou2,"li1"]);
  nbi <- length(oun);
  # getting the item names
  inames <- rep("",nbi);
  for (ii in rrrbc(nbi)) {
    tex <- text[oun[ii]];
    inames[ii] <- substr(tex,nchar(lndel[1])+1,as.numeric(regexpr(lndel[2],tex))-1);
  }
  # getting the value for each item name
  ivalues <- rep("",nbi);
  for (ii in rrrbc(nbi)) {
    # interval where must be the values
    i1 <- oun[ii];
    i2 <- c(oun,length(text))[ii+1];
    # first position for opening the value item
    v1 <- which(oudon[,"tag"] == lvdel[1]);
    v1 <- oudon[v1,"li1"];
    v1 <- v1[(v1>=i1)&(v1<=i2)];
    if (length(v1) == 0) {
      ivalues[ii] <- NA;
    } else {
      # first position for closing the value item
      v1 <- v1[1];
      w1 <- which((oudon[,"li1"]==v1)&(oudon[,"tag"]==lvdel[1]));
      if (length(w1)==0) { rrrerreur(v1,"Internal Error 1!");}
      w1 <- oudon[w1[1],"co1"];
      v2 <- which(oudon[,"tag"] == lvdel[2]);
      v2 <- oudon[v2,"li1"];
      v2 <- v2[(v2>=v1)&(v2<=i2)];
      if (length(v2) == 0) {
        ivalues[ii] <- NA;
      } else {
        v2 <- v2[1];
        w2 <- which((oudon[,"li1"]==v2)&(oudon[,"tag"]==lvdel[2]));
        if (length(w2)==0) { rrrerreur(v2,"Internal Error 2!");}
        w2 <- oudon[w2[1],"co2"];
        ivalues[ii] <- list(rrrtext3stext(text,c(v1,w1,v2,w2)));
      }
    }
    # removing the remaining braces
    resu <- ivalues[[ii]];
    if (length(resu) > 0) {
      if (nchar(resu[1]) == nchar(lvdel[1])) {
        resu <- resu[-1];
      } else {
        resu[1] <- substr(resu[1],nchar(lvdel[1])+1,nchar(resu[1]));
      }
      der <- length(resu);
      if (nchar(resu[der]) == nchar(lvdel[2])) {
        resu <- resu[-der];
      } else {
        resu[der] <- substr(resu[der],1,nchar(resu[der])-nchar(lvdel[2]));
      }
    }
    #
    res[[ii]] <- resu;
  }
  names(res) <- inames;
  
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrinterv7belonging <- function(x,int,monitor=rrrrbsa0$monitor$v)
#TITLE  checks if a series of values belong to a series of intervals
#DESCRIPTION
# computes and returns the indicator vector of the positions of
# values with respect to intervals.
#DETAILS
# This function is compatible with real infinite values
#KEYWORDS misc
#INPUTS
#{x} <<Vector of value(s) to be scrutinized.>>
#{int} <<Series of interval(s) to be considered.
# Either a \samp{numeric(2)} or a matrix with two columns.
# Empty intervals (\samp{numeric(0)} are not admitted.>>
#[INPUTS]
#{monitor} <<List providing the monitoring constants, see \code{rrrrbsa0$monitor$v}
#            to know the contents.>>
#VALUE
# A matrix with rows associated to the \code{x} values and
# columns associated to the \code{int} intervals giving
# \code{-2,-1,0,1,2} according to whether \code{x} is less than,
# equal to the lower bound, inside, equal to the upper bound or
# greater than the interval.
#EXAMPLE
# rrrinterv7belonging(1:5,1:2);
# rrrinterv7belonging(1:5,matrix(1:10,ncol=2));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 10_11_17
#REVISED 10_12_13
#--------------------------------------------
{
    # checking
    if (monitor$chk$v) {
        rrrobject9(  x,"numeric",-1,mensaje="  'x' must be numeric");
        rrrobject9(int,"numeric",-1,mensaje="'int' must be numeric");
        if (!is.matrix(int)) {
            if (length(int) != 2) {
                rrrerreur(int,"When 'int' is not a matrix, it must be a numeric(2)");
            }
            if (is.nan(diff(int))) {rrrerreur(int,"This is not an accepted interval");}
            if (diff(int)<0) { rrrerreur(int,"'int' does not define an interval!");}
        } else {
            if (ncol(int)!=2) {
                rrrerreur(int,"When 'int' is a matrix, it must comprise 2 columnes");
            }
            ru <- int[,2] - int[,1];
            if (any(is.nan(ru))) { rrrerreur(int,"Some rows are not accepted as intervals");}
            if (any((ru<0))) {
                rrrerreur(int,"Not all rows of 'int' define an interval");
            }
        }
    }
    # getting a uniform presentation
    if (!is.matrix(int)) { int <- matrix(int,ncol=2);}
    # preparing the result
    nbx <- length(x); nbint <- nrow(int);
    res <- matrix(NA,nbx,nbint);
    dimnames(res) <- list(names(x),dimnames(int)[[1]]);
    # degenerate case
    if (length(res)==0) { return(res);}
    # ancillary functions
    be0 <- function(x,int0) {
      if (is.finite(int0)) {
        ss <- sign(x-int0);
      } else {
        ss <- rep(-sign(int0),length(x));
        ss[x==int0] <- 0;
      }
      ss;
    }
    bel <- function(x,int) {
        be0(x,int[1]) + be0(x,int[2]);
    }
    # computation
    for (ii in rrrbc(nrow(int))) {
        res[,ii] <- bel(x,int[ii,]);
    }
    # returning
    res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrbd <- function(n1,n2)
#TITLE  sequence for insertions
#DESCRIPTION
# This function returns \samp{n1:n2} when \samp{n1<=n2} and
#         \samp{numeric(0)} otherwise.
# Quite useful when some insertion must be done within
# a sequence
#DETAILS
#ALIAS
#KEYWORDS iteration helpful
#INPUTS
#{n1}    <<first element>>
#{n2}    <<second element>>
#[INPUTS]
#VALUE
# \samp{n1:n2} if \samp{n1<n2}
# else \samp{numeric(0)}.
#EXAMPLE
# xx <- 1:5;
# for (ii in 1:6) { print(c(xx[rrrbd(1,ii-1)],10,xx[rrrbd(ii,5)]));}
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
# rrrbc for BouCle (loop in French)
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
rrrtext2vma <- function(cha,what=rrrrbsa0$vma$v["v"],
                     xsep=rrrrbsa0$sep1$v,nat="C",
                       monitor=rrrrbsa0$monitor$v)
#TITLE  transforms a character into a vector (or matrix, or array), and conversely
#DESCRIPTION
# from a \samp{character} vector, returns a vector, or a matrix, or
# an array of _characters_ with possibly names, or dimames. The information
# can be supplied in different ways for each of the three possibilities.
# It is advised to try the proposed examples.
#DETAILS
# The processing is done in character mode but the result can be
# transformed into numerical or logical values with the help of argument \samp{nat}.
# \cr
# In fact \samp{rrrrbsa0$vma$v} coding is used for the argument \samp{what}. 
# This allows to easily modify the coding.
#KEYWORDS IO
#INPUTS
#{cha} << The character to transform.>>
#[INPUTS]
#{what} << Indicates which structure to return: either
# a vector, a matrix or an array.
#
# \cr ~~ \cr
#
#   For vectors, the possibilities are c/C/u/v/U/V (in fact the content of
#     \samp{rrrrbsa0$vma$v["c"]}, \samp{rrrrbsa0$vma$v["C"]},... but for the sake of the
#     simplicity, the names will be used because they are by default identical to 
#     the value; the same will be done for the other types):
#    \cr\samp{c} for a no named character(1); collapsing
#             is done with \samp{rrrrbsa0$sep0$v}.
#    \cr\samp{C} for a no named character() of any length
#             (components are separated with \samp{xsep} which are
#              removed from the result); collapsing
#             is done with \samp{rrrrbsa0$sep0$v}.
#    \cr\samp{v} or \samp{u} for a no named vector;
#    \cr\samp{V} for a named vector with
#          all names before the values; then an even number
#          of components must be provided.
#    \cr\samp{U} for a named vector with
#          names interlaced with the value (name_i, value_i); then an even number
#          of components must be provided.
#
# \cr ~~ \cr
#
#   For matrices, the possibilities are m/n/o/p/M/N/O/P:
#    \cr\samp{m} for a no named matrix given by rows, two adjacent rows
#          being separated with \samp{xsep} sequence, introduced as one of the
#          component of \samp{cha}, then for a 2x3 matrix, the length of \samp{cha}
#          will be 6+2 = 8.
#    \cr\samp{n} for a matrix with only the columns names. The expected sequence is
#          the names of columns, then the values as for \samp{m}; then for a 2x3
#          matrix, the length of \samp{cha} will be 3+1+8=12.
#    \cr\samp{o} for a matrix with only rows named. The expected sequence is
#          name of row, values of rows... Then 2x3 will imply a length of 8+2=10.
#    \cr\samp{p} when names for columns and rows, in a mixed way... Then 2x3 will imply
#          a length of 14.
#    \cr
#     When \samp{M}, \samp{N},
#          \samp{O} or \samp{P},
#          the same but the matrix will be transposed after
#          being read; said in another way, the values are given by columns.
#
# \cr ~~ \cr
#
#   For arrays, the possibilities are a/A/b/B:
#    \cr\samp{a} for a no named array, the dimensions, \samp{xsep}, the values in
#    the classical order (varying the first index the fastest). 2x3 will give
#    a length of 2+1+6=9.
#    \cr\samp{A} for a dimnamed array, the dimensions, \samp{xsep}, the dimnames of each
#    dimension in the right order also separated and finished with \samp{xsep}. 
#    2x3 will gives a length of 2+1+2+1+3+1+6=16.
#    \cr\samp{b} for a named dimensions array, the dimensions, \samp{xsep}, the names for the
#    dimension in the right order not separated and finished with \samp{xsep}. 
#    2x3 will gives a length of 2+1+2+1+6=12.
#    \cr\samp{B} for a named and dimnamed array, the dimensions, \samp{xsep}, the names for the
#    dimension in the right order not separated and finished with \samp{xsep}, then the dimnames separated
#    before the values. 
#    2x3 will gives a length of\cr (2+1)+(2+1)+(2+1+3+1)+(6)=19. >>
#{xsep} << Character sequence used to separate the character vector into blocks
#    giving information about the structure (see the examples).>>
#{nat} << Nature of the returned structure. Can be \samp{C} for character, \samp{N}
#         for numeric or \samp{L} for logical.>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE a vector or a matrix or an array according to the arguments.
#EXAMPLE
#####
## vectors
# rrrtext2vma(letters,"c");
# rrrtext2vma(letters,"C",xsep="e");
# rrrtext2vma(letters);
# rrrtext2vma(letters,"V");
# rrrtext2vma(letters,"u");
# rrrtext2vma(c(LETTERS,letters),rrrrbsa0$vma$v["V"]);
# rrrtext2vma(c("A","a","B","b","C","c"),rrrrbsa0$vma$v["U"]);
#####
## matrices
# rrrtext2vma(c(1:3,"//",4:6),rrrrbsa0$vma$v["m"]);
# rrrtext2vma(c(1:3,"//",4:6),rrrrbsa0$vma$v["M"]);
# rrrtext2vma(c(LETTERS[1:3],"//",1:3,"//",4:6),rrrrbsa0$vma$v["n"]);
# rrrtext2vma(c(LETTERS[1:3],"//",1:3,"//",4:6),"N");
# rrrtext2vma(c("a",1:3,"//","b",4:6),"o");
# rrrtext2vma(c(c(LETTERS[1:3],"//","a",1:3,"//","b",4:6)),rrrrbsa0$vma$v["p"]);
#####
## arrays
# rrrtext2vma(c(2:4,"//",1:24),"a");
# rrrtext2vma(c(2:4,"//","one","two","//",LETTERS[1:3],"//",
#          letters[1:4],"//",1:24),"A");
# rrrtext2vma(c(2:4,"//","one","two","//",LETTERS[1:3],"//",
#          letters[1:4],"//",1:24),"A",nat="N");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_03_28
#REVISED 10_08_13
#--------------------------------------------
{
# flexibility
if (rrrvoid9(cha)) { cha <- character(0);}
# of use
ssep <- which(cha==xsep);
# checking
if (monitor$chk$v) {
    rrrobject9(cha,"character",-1,mensaje="cha must be a character");
    #
    rrrobject9(what,"character",1,
               mensaje="what (character(1)) must indicate the type of desired output");
    if (!(what %in% rrrrbsa0$vma$v)) {
        rrrerreur(what,"'what' not in the list of possible values...");
    }
    #
    rrrobject9(xsep,"character",1,mensaje="must indicate the character string of separation");
    #
    if (what %in% rrrrbsa0$vma$v[c("V","U")]) {
        if ((length(cha) %% 2) != 0) {
            rrrerreur(list(cha,what),"Here the length of 'cha' must be even");
        }
    }
    #
    if (what %in% rrrrbsa0$vma$v[c("n","N","m","M","o","O")]) {
        if (length(ssep) >= 1) {
            nbc <- ssep[1] - 1;
            nbr <- (length(cha)+1) / (nbc+1);
            if ((nbr!=round(nbr)) | (length(cha) != (nbc*nbr+nbr-1))) {
                rrrerreur(list(what,length(cha),nbc,nbr),"Dimensions not consistent");
            }
        }
    }
    #
    if (what %in% rrrrbsa0$vma$v[c("p","P")]) {
        if (length(ssep) >= 1) {
            nbc <- ssep[1] + 1;
            X <- 2 + length(cha);
            if ((X %% nbc) != 0) {
                rrrerreur(list(what,length(cha),nbc),"Dimensions not consistent");
            }
        }
    }
    #
    if (what %in% rrrrbsa0$vma$v[c("a","A","b","B")]) { 
	# an array must be returned
        if (length(ssep)==0) {
            rrrerreur(list(cha,xsep),"For arrays, dims must be provided");
        }
        didi <- cha[1:(ssep[1]-1)];
        didi <- as.numeric(didi);
        if (sum(is.na(didi))>0) {
            rrrerreur(list(cha,xsep),"Dimensions are not numeric !");
        }
        if (any(didi<0)) {
            rrrerreur(didi,"Negative dimensions were provided");
        }
        if (any(didi!=round(didi))) {
            rrrerreur(didi,"Non integer dimensions were provided");
        }
        X <- length(didi) + 1 + prod(didi);
        if (what %in% rrrrbsa0$vma$v[c("A","B")]) { X <- X + sum(didi) + length(didi); }
        if (what %in% rrrrbsa0$vma$v[c("b","B")]) { X <- X + length(didi) + 1; }
	if (length(cha) != X) {
            rrrerreur(list(cha,what,xsep),"Inconsistency for an array");
	}
    }
    #
    rrrobject9(nat,"character",1,mensaje="'nat' must be a character(1)");
}
#
res <- character(0);
#
if (what %in% rrrrbsa0$vma$v["c"]) {
    # a single character must be returned
    # rrrrbsa0$sep0$v is used as collapsor
    res <- paste(cha,collapse=rrrrbsa0$sep0$v);
}
#
if (what %in% rrrrbsa0$vma$v["C"]) {
    # a character must be returned
    # rrrrbsa0$sep0$v is used as separator
    res <- paste(cha,collapse=rrrrbsa0$sep0$v);
    res <- strsplit(res,xsep,fixed=TRUE)[[1]];
    if (length(res)>1) { for (ii in 2:length(res)) {
        res[ii  ] <- rrrform3crop(res[ii  ],rrrrbsa0$sep0$v,"",1);
        res[ii-1] <- rrrform3crop(res[ii-1],"",rrrrbsa0$sep0$v,1);
    }}
}
#
if (what %in% rrrrbsa0$vma$v[c("v","V","u","U")]) {
    # a vector must be returned
    if (what %in% c(rrrrbsa0$vma$v[c("v","u")])) {
        res <- cha;
    } else {
        nb <- floor(length(cha)/2);
        if (what == rrrrbsa0$vma$v["V"]) { nam <- rep(c(TRUE,FALSE),each=nb);
        } else { nam <- rep(c(TRUE,FALSE),nb); }
    res <- cha[!nam];
    names(res) <- cha[nam];
    }
}
#
if (what %in% rrrrbsa0$vma$v[c("m","M","n","N","o","O","p","P")]) {
    # a matrix must be returned
    nbc <- ssep[1] - 1;
    cha <- c(cha,xsep);
    if (what %in% c(rrrrbsa0$vma$v[c("p","P")])) {
        cha <- c(" ",cha);
    }
    if (what %in% rrrrbsa0$vma$v[c("p","P")]) { add <- 1;
    } else { add <- 0; }
    cha <- matrix(cha,ncol=nbc+1+add,byrow=TRUE);
    cha <- cha[,-ncol(cha),drop=FALSE];
    if (what %in% rrrrbsa0$vma$v[c("m","M")]) {
        res <- cha;
    }
    if (what %in% rrrrbsa0$vma$v[c("n","N")]) {
        res <- cha[-1,,drop=FALSE];
        dimnames(res) <- list(NULL,cha[1,,drop=FALSE]);
    }
    if (what %in% rrrrbsa0$vma$v[c("o","O")]) {
        res <- cha[,-1,drop=FALSE];
        dimnames(res) <- list(cha[,1,drop=FALSE],NULL);
    }
    if (what %in% rrrrbsa0$vma$v[c("p","P")]) {
        res <- cha[-1,-1,drop=FALSE];
        dimnames(res) <- list(cha[-1,1,drop=FALSE],cha[1,-1,drop=FALSE]);
    }
    if (what %in% rrrrbsa0$vma$v[c("M","N","O","P")]) {
        res <- t(res);
    }
}
#
if (what %in% rrrrbsa0$vma$v[c("a","A","b","B")]) {
    if (length(ssep) == 0) { rrrerreur(cha,"For array, dimensions must be provided");}
    # an array must be returned
    didi <- cha[1:(ssep[1]-1)];
    didi <- as.numeric(didi);
    nbdi <- length(didi);
    #
    if (what == rrrrbsa0$vma$v["a"]) { vvv <- cha[-(1:ssep[1])];}
    if (what == rrrrbsa0$vma$v["A"]) { vvv <- cha[-(1:ssep[1+nbdi])];}
    if (what == rrrrbsa0$vma$v["b"]) { vvv <- cha[-(1:ssep[2])];}
    if (what == rrrrbsa0$vma$v["B"]) { vvv <- cha[-(1:ssep[2+nbdi])];}
    #
    res <- array(vvv,dim=didi);
    #
    if (what %in% rrrrbsa0$vma$v[c("A","B")]) {
        ndi <- vector("list",0);
        for (jj in rrrbf(didi)) {
            jjj <- jj + (what == rrrrbsa0$vma$v["B"]);
            if (ssep[jjj+1]-ssep[jjj]>1) {
                ndi[[jj]] <- cha[(ssep[jjj]+1):(ssep[jjj+1]-1)];
            }
        }
        dimnames(res) <- ndi;
    }
    if (what %in% rrrrbsa0$vma$v["b"]) {
        ndi <- vector("list",0);
        for (jj in rrrbf(didi)) {
            ndi[[jj]] <- 1:didi[jj];
        }
        dimnames(res) <- ndi;
    }
    if (what %in% rrrrbsa0$vma$v[c("b","B")]) { 
        names(dimnames(res)) <- cha[(ssep[1]+1):(ssep[2]-1)];
    }
}
# transtyping
if (nat %in% c("N","L")) {
    rrr <- as.numeric(res);
    if (nat == "L") { rrr <- as.logical(rrr);}
    attributes(rrr) <- attributes(res);
    res <- rrr;
}
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3n2ij <- function(n,text,flexible=TRUE)
#TITLE gives the (row,column) position from the absolute position within a /text/
#DESCRIPTION 
# returns the (line, column) position associated to an absolute position when the
# lines of the \samp{text} are concatanated.\cr
# In fact, the reverse of \samp{rrrtext3ij2n}.
#DETAILS
#KEYWORDS 
#INPUTS
#{n} << A \samp{numeric(1)} indicating the absolute position.
#          May also be a vector. Values outside the possible range
#          are bounded to the extrema according to \samp{flexible}.  >>
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#[INPUTS]
#{flexible} << When \samp{TRUE} values of \samp{n} outside the accetable
#              range are accepted and replaced with first and last positions.>>
#VALUE
# The resulting positions associated to \samp{n}. When \samp{length(n)==1}
# it is a \samp{numeric(2)} if not a matrix with two columns.
#EXAMPLE
# aa <- c(paste(letters,collapse=""),paste(LETTERS,collapse=""),paste(0:9,collapse=""));
# rrrtext3n2ij(rrrtext3ij2n(c(1,3),aa),aa);
# rrrtext3n2ij(rrrtext3ij2n(c(2,1),aa),aa);
# rrrtext3n2ij(rrrtext3ij2n(rbind(c(1,3),c(2,1)),aa),aa);
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 14_06_12
#REVISED 14_06_13
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  if (flexible) {
    rrrobject9(n,"integer");
    n <- pmax(n,1);
    n <- pmin(n,sum(nchar(text)));
  } else {
    rrrobject9(n,"integer",con=c(1,sum(nchar(text))));
  }
  # initializing
  res <- matrix(NA,length(n),2);
  # degenerate cases
  if (length(n)==0) { return(res);}
  if (length(text)==0) { return(res);}
  # computing the limits
  limi <- cumsum(c(0,nchar(text)));
  # computing the line numbers
  iin <- apply(outer(n,limi,">"),1,sum);
  # computing the column numbers
  jjn <- n - limi[iin];
  #
  if (length(n)==1) {
    res <- c(iin,jjn);
  } else {
    res <- cbind(iin,jjn);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3stext <- function(text,sub=c(l1=1,p1=1,l2=Inf,p2=Inf))
#TITLE returns a subtext
#DESCRIPTION from a text,
# gives back a portion according to the \samp{sub} argument.
# Delimiting positions are included, so \samp{sub=c(1,1,1,1)}
# gives back the first character of \samp{text}.
#DETAILS
# Definition of the subtext is flexible (outside positions
# are interpreted as minimum / maximum positions).
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#[INPUTS]
#{sub} << A \samp{numeric(4)} indicating the portion of the \samp{text} to consider; can be a matrix but only one row is allowed.
# See the function \samp{rrrtext3interval} for details.>>
#VALUE
# The resulting sub-text
#EXAMPLE
# aa <- c(paste(letters,collapse=""),paste(LETTERS,collapse=""),paste(0:9,collapse=""));
# rrrtext3stext(aa);
# rrrtext3stext(aa,c(1,12,1,15));
# rrrtext3stext(aa,c(1,12,2,15));
# rrrtext3stext(aa,c(1,12,3,15));
# # In this last example, the last character of line 2 is kept (flexibility)!
# rrrtext3stext(aa,c(2,68,3,15));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_16
#REVISED 14_06_12
#--------------------------------------------
{
  # checking finished by rrrtext3interval
  if (is.matrix(sub)) {
    if (nrow(sub)>1) {
      rrrerreur(sub,"Only one interval is allowed");
    }
  }
  # getting the interval
  sub <- rrrtext3interval(text,sub);
  # extracting
  ii <- 1;
  if (is.na(sub[ii,1])) {
    res <- character(0);
  } else {
    res <- text[rrrbd(sub[ii,1],sub[ii,3])];
    if (sub[ii,1]==sub[ii,3]) {
      res <- substr(res,sub[ii,2],sub[ii,4]);
    } else {
      res[1] <- substr(res[1],sub[ii,2],nchar(res[1]));
      res[length(res)] <- substr(res[length(res)],1,sub[ii,4]);
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3acceptance <- function(line,posi,
                           b.accepted=NULL,a.accepted=NULL,
                           b.rejected=NULL,a.rejected=NULL,
                           b.line=TRUE,e.line=TRUE)
#TITLE acceptance of word within a line 
#DESCRIPTION 
# Indicates the acceptance of a sequence within a line
# based on some positive and negative criteria.
#DETAILS
# \samp{posi[1]>=posi[2]} case is accepted. Outside
# position are also accepted being valued to the bounds.
#KEYWORDS 
#INPUTS
#{line}<< A \samp{character(1)} within which is the sequence.>>
#{posi}<< The two positions within \samp{line} of the sequence.>>
#[INPUTS]
#{b.accepted} << A \samp{character} providing the surrounding characters 
#                before the words to determine the acceptance of the
#                sequence.>>
#{a.accepted} << A \samp{character} providing the surrounding characters 
#                after the words to determine the acceptance of the
#                sequence.>>
#{b.rejected} << A \samp{character} providing the surrounding characters 
#                before the words to determine the rejection of the
#                sequence.>>
#{a.rejected} << A \samp{character} providing the surrounding characters 
#                after the words to determine the rejection of the
#                sequence.>>
#{b.line} << Must sequences at the beginning of a line be accepted?>>
#{e.line} << Must sequences at the end of a line be accepted?>>
#VALUE
# \samp{TRUE} or \samp{FALSE}.
#EXAMPLE
# rrrtext3acceptance("Et pourquoi pas ?",c(1,3));
# rrrtext3acceptance("Et pourquoi pas ?",c(4,11),b.accepted=" ");
# rrrtext3acceptance("Et pourquoi pas ?",c(4,11),b.accepted=letters);
# rrrtext3acceptance("Et pourquoi pas ?",c(5,11),b.accepted=letters);
# \dontrun{rrrtext3acceptance("Et",1:2,b.accepted="letters",b.rejected=LETTERS)}
#REFERENCE
#SEE ALSO rrrtext3replace
#CALLING
#COMMENT This function was written to be used within \samp{rrrtext3replace}.
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 14_09_03
#REVISED 14_09_03
#--------------------------------------------
{
  # checking
  rrrobject9(line,"character");
  rrrobject9(posi,"integer",2);
  posi <- pmax(pmin(posi,nchar(line)),1);
  rrrobject9(b.accepted,c("null","character"));
  rrrobject9(b.rejected,c("null","character"));
  rrrobject9(a.accepted,c("null","character"));
  rrrobject9(a.rejected,c("null","character"));
  if ((!is.null(b.accepted)) & (!is.null(b.rejected))) {
    rrrerreur(list(b.accepted,b.rejected),"'b.accepted' and 'b.rejected' cannot be active simultaneously");
  }
  if ((!is.null(a.accepted)) & (!is.null(a.rejected))) {
    rrrerreur(list(a.accepted,a.rejected),"'a.accepted' and 'a.rejected' cannot be active simultaneously");
  }
  rrrobject9(b.line,"logical");
  rrrobject9(e.line,"logical");
  # looking for the before acceptability
  befo <- FALSE;
  if (b.line & (posi[1]==1)) {
    befo <- TRUE;
  } else {
    if (!is.null(b.accepted)) {
      for (mot in b.accepted) {
        if (nchar(mot) < posi[1]) {
          if (substr(line,posi[1]-nchar(mot),posi[1]-1)==mot) { befo <- TRUE;}
        }
      }
    } else {
      befo <- TRUE;
      if (!is.null(b.rejected)) {
        for (mot in b.rejected) {
          if (nchar(mot) < posi[1]) {
            if (substr(line,posi[1]-nchar(mot),posi[1]-1)==mot) { befo <- FALSE;}
          }
        }
      }
    }
  }
  # looking for the before acceptability
  afte <- FALSE;
  if (befo) {
    if (e.line & (posi[2]==nchar(line))) {
      afte <- TRUE;
    } else {
      if (!is.null(a.accepted)) {
        for (mot in a.accepted) {
          if (posi[2]+nchar(mot) <= nchar(line)) {
            if (substr(line,posi[2]+1,posi[2]+nchar(mot))==mot) { afte <- TRUE;}
          }
        }
      } else {
        afte <- TRUE;
        if (!is.null(a.rejected)) {
          for (mot in a.rejected) {
            if (posi[2]+nchar(mot) <= nchar(line)) {
              if (substr(line,posi[2]+1,posi[2]+nchar(mot))==mot) { afte <- FALSE;}
            }
          }
        }
      }
    }
  } 
  # returning
  (befo & afte);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext2file <- function(text,file,path="",
                      ap=FALSE,
                       monitor=rrrrbsa0$monitor$v
                     )
#TITLE  writes a character onto a file
#DESCRIPTION
# The reverse operation of \samp{rrrfile2text}. This function is not
# very tricky, it was written for completeness.
#KEYWORDS IO
#INPUTS
#{text} << The \samp{character} to be written into
#          \samp{file}.>>
#{file} << file to be written. According to \samp{ap} when already exist, 
# it will be destroyed and recreated or completed.
# When the \samp{file} is \samp{character(0)}, no file is considered but the
# corresonding \samp{character} is returned.
# >>
#[INPUTS]
#{path} << Directory containing the file.>>
#{ap} <<Must the file be appended if it already exist?>>
#{monitor} <<List of constants indicating the monitoring choices,
#            see the \samp{rrrrbsa0$monitor$v} provided object as an example.>>
#VALUE
# Nothing but a file is created or completed when everything is right -or-
# according to the value of \samp{file}, a character is returned.
#EXAMPLE
# rrrtext2file(letters,"toto.txt");
# unlink("toto.txt");
# rrrtext2file(letters,character(0));
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 08_02_13
#REVISED 14_06_24
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(text,"character",-1,mensaje="'text' must be a character");
    #
    rrrobject9(file,"character",c(0,1),mensaje="file must indicate the name of one file");
    #
    if (!rrrvoid9(path)) {
      path <- rrrdipa(path);
      rrrobject9(path,"character",1,mensaje="path must indicate the name of one directory");
      #
      if (is.na(file.info(path)$size)) { rrrerreur(path,"This directory seems not to exist.");}
    }
    #
    rrrobject9(ap,"logical",1,mensaje="'ap' must be a single logical");
}
#
# The returning case
if (rrrvoid9(file)) { return(text);}
#
# opening the proposed file
if (!is.null(path)) { file <- paste0(path,file);}
sink(file,append=ap);
#
# writting each component
for (ii in rrrbf(text)) {
    cat(text[ii],"\n");
}
#
# closing the file
sink();
#
# returning
invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrtext3interval <- function(text,sub=c(l1=1,p1=1,l2=length(text),p2=Inf))
#TITLE transforms /sub/ into a consistent interval for the /text/
#DESCRIPTION from a text and an interval
# gives back a precised interval taking into account the lengths
# of \samp{text} and possible inconsistencies. When the interval is empty
# \samp{numeric(0)} is returned.
#DETAILS
# Definition of the subtext is flexible (outside positions
# are interpreted as minimum / maximum positions).
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line).>>
#[INPUTS]
#{sub} << A \samp{numeric(4)} indicating the portion of the \samp{text} to consider.
# More precisely (first line, first position within the line, last line and 
# last position in the last line. Can be a matrix with four columns as well.>>
#VALUE
# The resulting interval: \samp{matrix[nsub,4]} where \samp{nsub} is the number
# of intervals. When one interval is not consistent, \samp{NA} are introduced.
#EXAMPLE
# aa <- c(paste(letters,collapse=""),paste(LETTERS,collapse=""),paste(0:9,collapse=""));
# rrrtext3interval(aa);
# rrrtext3interval(aa,c(1,1,12,15));
# rrrtext3interval(aa,c(1,2,12,15));
# rrrtext3interval(aa,c(10,2,1,15));
# rrrtext3interval(aa,c(1,3,12,15));
# rrrtext3interval(aa,c(2,3,68,15));
# rrrtext3interval(aa,matrix(c(1,1,12,15,
#                           1,2,12,15,
#                           10,2,1,15,
#                           1,3,12,15,
#                           2,3,68,15),byrow=TRUE,ncol=4));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_17
#REVISED 13_10_29
#--------------------------------------------
{
  # checking
  rrrobject9(text,"character");
  if (is.matrix(sub)) {
    if (ncol(sub)!=4) {
      stop("A four columns matrix is expected");
    }
  } else {
      rrrobject9(sub,"numeric",len=4);
      sub <- matrix(sub,ncol=4);
  }
  # initialization
  res <- matrix(NA,nrow=nrow(sub),ncol=ncol(sub));
  # getting the lines
  i1 <- pmin(pmax(1,sub[,1]),length(text));
  i2 <- pmax(pmin(length(text),sub[,3]),1);
  # getting the positions
  p1 <- pmin(pmax(1,sub[,2]),nchar(text[i1]));
  p2 <- pmax(pmin(nchar(text[i2]),sub[,4]),1);
  # preparing the result
  res <- cbind(i1,p1,i2,p2);
  res[(i1>i2) | ((i1==i2)&(p1>p2)),] <- NA;
  dimnames(res) <- list(NULL,c("l1","p1","l2","p2"));
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrform3parag <- function(texte,titre=-1,
                            wid=60,fli=NULL,sep=1,
                            jus=1,trunc=TRUE,
                            ed="  ",ef="",
                            monitor=rrrrbsa0$monitor$v,
                            imp=TRUE)
#TITLE  prints or prepares paragraphes from a character vector.
#DESCRIPTION
# prints or prepares a character string \samp{texte}
# as a small formatted text.\cr
#         Each component is supposed to be a
#         paragraph but the first one can be
#         considered as a title.
#DETAILS
#KEYWORDS print
#INPUTS
#{texte}<<The text to print (character vector).>>
#[INPUTS]
#{titre} << When > -2 indicates that the first
#           component is a title (if not a
#           simple paragraph). Then the absolute value
#           of \samp{titre} gives the emphasize to put
#           on the title. Notice that the title is
#           not splitted in several lines as are
#           the other components according to 'wid'
#           value.>>
#{wid} << The maximum width for each line (in characters) without
#        including indentation and frames.>>
#{fli} << When NULL, the first line of
#        each paragraph (in fact the second if there is 
#        a title) is issued as a standard line. If not,
#        \samp{fli[1]} spaces are added before and the considered
#        width is \samp{fli[2]} (not including the added spaces).
#        Also this means that your already wrote \samp{fli[3]}
#        characters on the first line; this last possibility
#        can be used only when there are no title and for the 
#        the first component.
#        For instance, a French paragraph will be issued
#        with \samp{fli=c(5,wid,0)}. The possibility of modifying
#        wid for the first line can be of use when adding
#        the name of an item first.>>
#{sep} << Number of lines separating each paragraph.>>
#{jus} << Type of justification (1=left, 2=centred,
#             3=left).>>
#{trunc} << Must truncation be done when a word
#          is greater than the proposed wid?>>
#{ed} << Framing at the beginning of each line.>>
#{ef} << Framing at the end of each line.>>
#{monitor} << List of monitoring constants, see \samp{rrrrbsa0$monitor$v} to
#             know its structure.>>
#{imp} << Printing is performed and nothing is returned.
#                If FALSE, the character string is returned 
#                (including possible new lines)>>
#VALUE
# either nothing or a character string according to \samp{imp}
#EXAMPLE
# rrrform3parag(c("My Title","My important first sentence.","Last one!"));
# rrrform3parag(rrrrbsa0$text2$v,wid=20);
#REFERENCE
#SEE ALSO
#CALLING {rrrform3repeat}
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 08_08_01
#REVISED 09_09_29
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(texte,"character");
    rrrobject9(wid,"numeric",1,c(1,Inf));
    if (!is.null(fli)) {
      rrrobject9(fli,"numeric",3,c(0,Inf));
    }
    rrrobject9(sep,"integer",1,c(1,Inf));
    rrrobject9(jus,"integer",1,c(1,3));
    rrrobject9(trunc,"logical",1);
    rrrobject9(ed,"character",1);
    rrrobject9(ef,"character",1);
    rrrobject9(imp,"logical",1);
}
sep <- min(max(0,sep),5);
if (titre > -2) {
    # the possible title
    res <- rrrform3title(texte[1],empha=titre,indent=nchar(ed),imp=imp);
    if (imp) { cat(rrrform3repeat("\n",sep+1));
    } else { res <- c(res,rrrform3repeat("\n",sep));}
    texte <- texte[-1];
} else { res <- character(0);}
#
for (i in rrrbf(texte)) {
    # paragraph after paragraph
    mots <- strsplit(texte[i]," ")[[1]];
    nlig <- 0;
    while (length(mots) > 0) {
        # the paragraph is not empty
        nlig <- nlig+1;
        # is it the first line and must it be different?
        spfl <- ((nlig==1)&(!rrrvoid9(fli)));
        if (spfl) {
            if (monitor$chk$v) {
                rrrobject9(fli,"numeric",c(3,Inf));
            }
            wiid <- fli[2];
            lili <- rrrform3repeat(" ",fli[1]);
            trop <- fli[3];
        } else {
            wiid <- wid;
            #lili <- character(0);
            lili <- "";
            trop <- 0;
        }
        if (rrrvoid9(wiid)) { wiid <- 60;}
	while ( ((length(mots) > 0) &&
		 (nchar(paste(lili,collapse=" ")) < (wiid-nchar(mots[1])))
		) || (
		 (length(lili) == 0)
		)
	      ) {
	    lili <- c(lili,mots[1]);
	    mots <- mots[-1];
	}
	#cat("<<",length(mots),">>\n");
        lili <- paste(lili,collapse=" ");
        #cat("{{",nchar(lili),lili,"}}\n");
        lili <- rrrform3justify(lili,wid-trop,jus,trunc);
        if (spfl) { lili <- paste0(lili,ef,"\n");
        } else { lili <- paste0(ed,lili,ef,"\n");}
	if (imp) { cat(lili);
	} else {
	    res <- c(res,lili);
	    lili <- character(0);
	}
    }
    if (imp) { cat(rrrform3repeat("\n",sep));
    } else { res <- c(res,rrrform3repeat("\n",sep));}
}
# returning
if (!imp) { return(res);}
cat(res);
invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrdisplay8k <- function(constants,what="exhibit")
#TITLE  returns information about a /rbsa/ constant
#DESCRIPTION
# Exhibits a \samp{rbsa} constant 
# or returns a \samp{character} vector for the names and
# the definitions of \samp{constants} or a \samp{list} for
# the values.
#DETAILS
# A \samp{rbsa} type constant is a named list, each component being
# a list of two components: \samp{$d} for the definition which must be
# a character (concatenated if not a scalar) and \samp{$v} which
# can take any value.
#KEYWORDS misc helpful
#INPUTS
#{constants} << The \samp{rbsa} constant list (see the details for the description).>>
#[INPUTS]
#{what}    <<a character(1) indicating what must be returned. The acceptable
#           values are \samp{exhibit}, \samp{names}, \samp{definitions} and \samp{values}.
#           For user convenience, the argument is lowercased before checked,
#           and only the first character is taken into consideration.>>
#VALUE
# Information is displayed and nothing is returned for exhibition
# else a named \samp{character} or a named \samp{list}. 
#EXAMPLE
# rrrdisplay8k(rrrrbsa0)
# rrrdisplay8k(rrrrbsa0$monitor$v,"v")
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 14_05_21
#REVISED 14_06_03
#--------------------------------------------
{
  # checking
  rrrobject9(what,"character",1,mensaje="Awaiting for a character(1)");
  what <- tolower(what);
  what <- substr(what,1,1);
  rrrobject9(what,"character",1,c("d","v","n","e"),
          "Must be 'names', 'descriptions', 'values' or 'exhibit'");
  # looking if it is a constant or a list of constants
  nana <- names(constants);
  if (length(unique(c(nana,"v","d")))==2) {
    # must be a constant
    lcon <- list(k=constants);
  } else {
    # must be a list of constant
    lcon <- constants;
  }
  # filling values
  if (what == "v") {
    res <- vector("list",length(lcon));
    names(res) <- names(lcon);
    for (vv in rrrbf(lcon)) {
      res[[vv]] <- lcon[[vv]]$v;
    }
  }
  # filling descriptions
  if (what == "d") {
    res <- vector("character",length(lcon));
    names(res) <- names(lcon);
    for (vv in rrrbf(lcon)) {
      res[vv] <- paste(lcon[[vv]]$d,collapse=" / ");
    }
  }
  # filling names
  if (what == "n") {
    res <- names(lcon);
  }
  # displaying
  if (what == "e") {
    for (cc in rrrbf(lcon)) {
      rrrform3title(names(lcon)[cc],4);
      print(lcon[[cc]]$d);
      print(lcon[[cc]]$v);
    }
    res <- invisible();
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rrrrbsa7list9 <- function(lili,monitor=rrrrbsa0$monitor$v)
#TITLE checks whether a list is rbsa-compatible list
#DESCRIPTION
# To be a rbsa-list, a list must satisfy the following two properties:
# (i) all components and sub-components are named.
# (ii) all components and sub-components are either a list
# or a character vector/matrix/array (i.e. vma components); they are 
# the leafs of the rbsa-list.\cr
# To be handled by \samp{list2file} or \samp{rrrfile2list} functions, a list must
# rbsa-compatible.
#DETAILS
#KEYWORDS IO
#INPUTS
#{lili} << The list structure to be checked.>>
#[INPUTS]
#{monitor} << List of monitoring constants, see \samp{rrrrbsa0$monitor$v} to
#             know its structure.>>
#VALUE
# \samp{TRUE} or \samp{FALSE} according to the results of the checks.
#EXAMPLE
# rrrrbsa7list9(rrrrbsa0$lis1$v);
# rrrrbsa7list9(list(rrrrbsa0$lis1$v));
# rrrrbsa7list9(list(rrrrbsa0$lis1$v,nu=NULL));
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 10_04_19
#REVISED 10_04_20
#--------------------------------------------
{
# checking
if (monitor$chk$v) {
    rrrobject9(lili,"list",-1,mensaje="lili must be a list");
}
# exploring the list
eee <- rrrexplore8list(lili);
# checking the presence on unamed components
#          and the types of the list
comment <- character(0);
for (ii in rrrbc(nrow(eee))) {
    if ((eee[ii,"name"]=="<NA>")|(eee[ii,"name"]=="")) {
        comment <- c(comment,
                     paste("The component",eee[ii,"numbers"],
                           "has got no name:",
                           paste("'",eee[ii,"names"],"'",sep=""))
                    );
    }
    coco <- rrrget8comp7list(lili,eee[ii,,drop=FALSE],monitor=monitor);
    if (length(coco) > 0) { coco <- coco[[1]];}
    if (!(is.list(coco) |
          is.numeric(coco) |
          is.character(coco) | 
          is.matrix(coco) | 
          is.array(coco))) {
        comment <- c(comment,
                     paste("The component (",eee[ii,"numbers"],
                           ") with name: '",eee[ii,"names"],
                           "is not list/vector/matrix/array")
                    );
    }
}
# preparing the result
res <- rrrvoid9(comment);
# returning
res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
