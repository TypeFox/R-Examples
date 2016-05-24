
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
code7objects4text6tags <- function(text,otags=documair0$tags$v,
                                        bra=c("deb","fin"),
                                        addbra=FALSE)
#TITLE returns the code of a series of objects
#DESCRIPTION from  \samp{text}, 
# gives back a list, each component of it are \samp{text}
# containing the code associated to one object.\cr
# The opening and closing tags are taken from \samp{otags[[bra[1]]}
# and \samp{otags[[bra[2]]]}.
#DETAILS
# The name of each object is supposed to be on the line
# after the opening tag (i.e. the second line) and separated with
# blanks or nothing at the beginning and blanks or \samp{=} or
# \samp{<-} at the end.
#KEYWORDS 
#INPUTS
#{text}<< A \samp{character} vector containing the text
# (a component, a line). If every opening tag is not closed with
# a closing tag, a fatal error is issued.>>
#[INPUTS]
#{otags}<< \samp{list()} defining the tags to be used.
# Its structure must be identical to its of \samp{documair0$tags$v}.>>
#{bra} << \samp{character(2)} indicating the name of the two components
#         of \samp{otags} which must be considered as opening and closing
#         brackets.>>
#{addbra}<< Must brackets be added as additional lines?>>
#VALUE
# A named list of the found objects
#EXAMPLE
# code7objects4text6tags(documair0$text4$v,documair0$tgs$v);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_28
#REVISED 13_12_30
#--------------------------------------------
{
  # some checking
  rrrobject9(text,"character",mensaje="'code7objects4text6tags' awaits for a 'character' as first argument");
  rrrobject9(otags,"list",mensaje="'code7objects4text6tags' awaits for a 'list' as second argument");
  rrrobject9(bra,"character",2,names(otags),mensaje="'bra' not included in 'names(otags)'");
  # getting the tags values
  brack <- c(otags[[bra[1]]]$v,otags[[bra[2]]]$v);
  # getting the different pieces
  res <- rrrtext3brackets(text,bra=brack,
                       col1=c(1,1),col2=c(1,1),rembra=!addbra);
  # removing empty lines
  for (cc in rrrbf(res)) {
    res[[cc]] <- rrrtext3preparation(res[[cc]],"R");
  }
  # degenerate case
  if (length(res)==0) { return(res);}
  # naming them
  uu <- rep("",length(res));
  for (cc in rrrbf(res)) {
    uuu <- res[[cc]][1+addbra];
    if (is.na(uuu)) {
      print(uuu);
      stop("No first line found to get the object name!");
    }
    # removing the starting blanks
    uuu <- rrrform3crop(uuu," ","");
    # replacing affectation sign with blank
    uuu <- sub( "="," ",uuu,fixed=TRUE);
    uuu <- sub("<-"," ",uuu,fixed=TRUE);
    # getting the name
    uu[cc] <- strsplit(uuu," ",fixed=TRUE)[[1]][1];
  }
  names(res) <- uu;
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
parse8code <- function(code,
                       otags=documair0$tags$v,
                       tags=names(documair0$tags$v))
#TITLE returns the tag contents of a code associated with an object
#DESCRIPTION from  \samp{code}, 
# gives back a named list with \samp{tags}, each component can be (mostly) \samp{text}
# or a \samp{boolean} or a \samp{list}
# containing the content associated to a tags described in \samp{otags}.\cr
# More, the two components \samp{$name} and \samp{$usage} are deduced from the analysis
# of the first line before the tag \samp{documair0$tags$v[["tit"]]$v}. This must be in first position!
#DETAILS
#KEYWORDS 
#INPUTS
#{code}<< A \samp{character} vector containing the code.
# Tagging must be consistent with the general rule of /documair/.>>
#[INPUTS]
#{otags}<< \samp{list()} defining the tags to be used.
# Its structure must be identical to its of \samp{documair0$tags$v}.>>
#{tags} << \samp{character()} indicating the name of the components
#         of \samp{otags} which must be considered as tags.>>
#VALUE
# A named list of texts associated to each tag except the last one
# which is just a closing tag.
#EXAMPLE
# parse8code(documair0$text3$v,
#            documair0$tgs$v,
#            names(documair0$tgs$v));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_30
#REVISED 14_05_21
#--------------------------------------------
{
  # internal constants
  care <- "|||";
  # some checking
  rrrobject9(code,"character",mensaje="'parse8code' awaits for a 'character' as first argument");
  rrrobject9(otags,"list",mensaje="'parse8code' awaits for a 'list' as second argument");
  rrrobject9(tags,"character",-1,names(otags),mensaje="'tags' not included in 'names(otags)'");
  # getting name and usage components from the first lines
  outi <- rrrtext3places8word(code,otags[["tit"]]$v,column=c(1,1));
  if (dim(outi)[1] != 1) {
    print(code);
    print(outi);
    stop("The 'tit' tags was not found in the code.");
  }
  entete <- code[rrrbc(outi[1,1]-1)];
  # concatanation preserving line breaks
  entete <- paste(entete,collapse=care);
  # removing all blanks
  entete <- gsub(" ","",entete,fixed=TRUE);
  # replacing affectation sign with blank
  if (grep("<-",entete,fixed=TRUE)==0) {
    entete <- sub( "="," ",entete,fixed=TRUE);
  } else {
    entete <- sub("<-"," ",entete,fixed=TRUE);
  }
  # getting the name
  uuu <- strsplit(entete," ",fixed=TRUE)[[1]];
  nom <- uuu[1];
  # getting the usage if any
  #gtest <- grep("function(",uuu[2],fixed=TRUE);
  gtest <- (tolower(substr(uuu[2],1,9))=="function(");
  usa <- NULL;
  if (!is.na(gtest)) {
    uuu <- strsplit(uuu[2],"(",fixed=TRUE)[[1]];
    usa <- paste(nom,"(",
                 paste(uuu[rrrbc(length(uuu)-1)+1],collapse="("),
                 sep="");
    # reintroducing line breaks
    sali <- paste("\n",paste(rep(" ",nchar(nom)+3),collapse=""),sep="");
    usa <- paste(strsplit(usa,care,fixed=TRUE)[[1]],collapse=sali);
  }
  # getting the tags values / contents / preparation
  tagos <- extract8object(otags[tags],"v");
  conte <- as.numeric(extract8object(otags[tags],"c"));
  names(conte) <- names(tagos);
  prepa <- extract8object(otags[tags],"p");
  # splitting the text according to the tags
  tagpl <- rrrplaces4text6tags(code,tagos);
  # checking that no duplicate tags was found
  tagus <- extract8object(otags[tags],"u");
  tagus <- tagus[which(tagus==1)];
  for (tat in names(tagus)) {
    if (paste0("1.",tat) %in% dimnames(tagpl)[[1]]) {
      print(dimnames(tagpl)[[1]]);
      stop("Looks like duplicated tag(s) exists");
    }
  }
  res <- rrrtexts4text(code,tagpl[,c(1,2,4),drop=FALSE],
                      addbeg=FALSE,addend=FALSE);
  # preparing each tag content
  for (tt in names(res)) {
    # identifying the tag to decode
    if (!(tt %in% tags)) {
      # removing some possible numbering
      ttt <- strsplit(tt,".",fixed=TRUE)[[1]][2];
      if (!(ttt %in% tags)) {
        print(tags);
        stop(paste("Programming error:",tt,"or",ttt,"not in tags"));
      }
    } else {
      ttt <- tt;
    }
    # checking and preparation
    ppp <- prepa[ttt];
    if (is.na(ppp)) {
      stop(paste("no 'preparation' indication found! for",ttt));
    }
    ccc <- conte[ttt];
    if (is.na(ccc)) {
      stop(paste("no 'content' indication found! for",ttt));
    }
    if (ccc == 0) {
      # the presence case
      res[[tt]] <- TRUE;
    } else {
      # applying the required formatting
      ccc <- c(1,min(ccc,length(res[[tt]])));
      pte <- rrrtext3preparation(res[[tt]],preparation=ppp,content=ccc);
      res[[tt]] <- pte;
    }
  }
  # adding missing logical
  talo <- setdiff(tags[conte==0],names(res));
  for (ll in talo) {
    llo <- length(res)+1;
    res[[llo]] <- "FALSE";
    names(res)[llo] <- ll;
  }
  # adding name and usage
  res$name <- nom;
  res$usage <- usa;
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
make8rd <- function(parsed,addi=TRUE,alias=NULL)
#TITLE returns the Rd documentation of a parsed object
#DESCRIPTION from  \samp{parsed} (typically a list coming from
# the function \samp{parse8code}) elaborates the equivalent of
# an \samp{Rd} file.\cr
# The only compulsory components are \samp{name},
# \samp{tit} (for title) and \samp{des} (for description).
#DETAILS
# Be aware that the names of \samp{documair0$tags$v} are directly used;
# if you modify them, you must modify this function.
#KEYWORDS 
#INPUTS
#{parsed}<< A \samp{list} coming from the \samp{parse8code}
# or equivalent.>>
#[INPUTS]
#{addi} << Add the additional sections. >>
#{alias} << If present, the alias names for the object.>>
#VALUE
# A \samp{character}.
#EXAMPLE
# par <- parse8code(documair0$text3$v,
#                   documair0$tgs$v,
#                   names(documair0$tgs$v));
# make8rd(par);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_11_04
#REVISED 14_07_03
#--------------------------------------------
{
  # some checking
  obli <- c("name","des","tit");
  rrrobject9(parsed,"list",mensaje="'make8rd' awaits for a 'list' as first argument");
  if (length(unique(obli,names(parsed))) > length(names(parsed))) {
    rrrform3display(parsed);
    rrrerreur(obli,"At least, one of these components is missing in 'parsed'!");
  }
  # initializing
  res <- c("%",paste0("% created by /documair/ at ",rrrnow("dm")),"%");
  # alias
  if (is.null(parsed[["ali"]])) {
    parsed[["ali"]] <- parsed[["name"]];
  }
  # compulsory components
  stan <- c(name="name",alias="ali",title="tit");
  for (ix in rrrbf(stan)) {
    nx <- names(stan)[ix];
    vx <- stan[ix];
    if (is.null(parsed[[vx]])) {
      print(names(parsed));
      stop(paste0("No ",nx,"[",vx,"] found!"));
    }
    mettre <- parsed[[vx]];
    res <- c(res,paste0("\\",nx,"{",mettre,"}"));
    res <- c(res,"%");
  }
  # optional components
  stan <- c(usage="usage",examples="exa");
  for (ix in rrrbf(stan)) {
    nx <- names(stan)[ix];
    vx <- stan[ix];
    if ((!is.null(parsed[[vx]]))&(length(parsed[[vx]])>0)) {
      mettre <- parsed[[vx]];
      if (nx == "usage") {
        for (im in rrrbf(mettre)) {
          imv <- mettre[im];
          sanom <- strsplit(imv,"(",fixed=TRUE)[[1]];
          if (length(sanom)>1) {
            # it is a fonction
            poner <- rrrform3crop(sanom[1]);
            poner <- strsplit(poner,".",fixed=TRUE)[[1]];
            if (length(poner)>1) {
              # there are dot(s) in the function name
              premier <- paste0("{",paste(poner[-length(poner)],collapse="."),"}");
              second <- paste0("{",poner[length(poner)],"}");
              sanom[1] <- paste0("\\method",premier,second);
            }
            mettre[im] <- paste(sanom,collapse="(");
          }
        }
      }
      res <- c(res,paste0("\\",nx,"{"));
      res <- c(res,paste0("  ",mettre));
      res <- c(res,"}");
      res <- c(res,"%");
    }
  }
  # formatted components
  form <- c(details="det",description="des",value="val");
  for (ix in rrrbf(form)) {
    nx <- names(form)[ix];
    vx <- form[ix];
    if (!is.null(parsed[[vx]])) { if (nchar(paste(parsed[[vx]],collapse=""))>0) {
      vz <- rrrparse8text(parsed[[vx]]);
      mettre <- "";
      for (iz in rrrbf(vz)) {
        nz <- substr(names(vz)[iz],1,1);
        if (nz == "P") {
          mettre <- c(mettre,vz[[iz]]);
        } else {
          if (nz == "D") {
            mettre <- c(mettre,"\\describe{");
            for (jz in rrrbf(vz[[iz]])) { 
              mettre <- c(mettre,paste0("\\item{",names(vz[[iz]])[jz],"}"));
              mettre <- c(mettre,vz[[iz]][jz][[1]]);
            }
            mettre <- c(mettre,"}");
          }
          if (nz == "I") {
            mettre <- c(mettre,"\\itemize{");
            for (jz in rrrbf(vz[[iz]])) { 
              mettre <- c(mettre,paste0("\\item"));
              mettre <- c(mettre,vz[[iz]][jz][[1]]);
            }
            mettre <- c(mettre,"}");
          }
          if (nz == "N") {
            mettre <- c(mettre,"\\enumerate{");
            for (jz in rrrbf(vz[[iz]])) { 
              mettre <- c(mettre,paste0("\\item"));
              mettre <- c(mettre,vz[[iz]][jz][[1]]);
            }
            mettre <- c(mettre,"}");
          }
        }
        mettre <- c(mettre,"");
      }
      res <- c(res,paste0("\\",nx,"{"));
      res <- c(res,paste0("  ",mettre));
      res <- c(res,"}");
      res <- c(res,"%");
    }}
  }
  # arguments (identical ones are deleted without warning)
  argus <- c(parsed[["in0"]],parsed[["in1"]]);
  if (length(argus) > 0) {
    # analyzing the arguments
    vz <- rrrlist4text(argus);
    # eliminating the possible duplicates
    nz <- unique(names(vz));
    if (length(nz) < length(vz)) {
      for (nn in nz) {
        zz <- which(nn==names(vz));
        if (length(zz) > 1) {
          vz <- vz[-zz[-1]];
        }
      }
    }
    # writing the code
    res <- c(res,"","\\arguments{");
    for (jz in rrrbf(vz)) {
      lequel <- names(vz)[jz];
      res <- c(res,paste0("    \\item{",lequel,"}{",
                          paste0(vz[[lequel]],collapse=" "),"}"));
    }
    res <- c(res,"}","","%");
  }
  # keywords
  for (kk in rrrbf(parsed$kwd)) {
    res <- c(res,paste0("\\keyword{",parsed$kwd[kk],"}"));
  }
  #
  # additional information
  if (addi) {
    if (!rrrvoid9(parsed$fut)) {
      res <- c(res,"\\section{Future}{");
      res <- c(res,parsed$fut);
      res <- c(res,"}");
    }
    if (!rrrvoid9(parsed$sal)) {
      res <- c(res,"\\seealso{");
      res <- c(res,paste0("\\code{\\link{",parsed$sal,"}}"));
      res <- c(res,"}");
    }
  }
 #
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
