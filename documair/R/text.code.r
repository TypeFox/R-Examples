
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
components9 <- function(noms,
                        seps=as.character(0:9),
                        void=c("","."),
                        rnom=TRUE)
#TITLE returns the components
#DESCRIPTION from a character returns the sorted list
# of elementary components of all components of \samp{noms}.
# Elementary components are given by splitting with digits
# (or any patterns introduced in \samp{seps}).
#DETAILS
# A very raw use is made on the \samp{seps} vector, so ambiguous
# vectors can give unexpected results.
# For instance better to avoid \samp{seps=c("aa","aaa")}!
# "|" is forbidden.
#KEYWORDS 
#INPUTS
#{noms}<< The \samp{character} to be analyzed.>>
#[INPUTS]
#{seps}<< Arrays of the seperators to consider.>>
#{void}<< Series of elementary components to eliminate>>
#{rnom}<< Indicates if complete names must be eliminated
#         from the series of elementary components.>>
#VALUE
# The sorted elementary components in a \samp{character};
# when there is no one digits, returns \samp{character(0)}.
#EXAMPLE
# components9(c("aa12bb","uuu","XX3","YY0ZZ"));
#REFERENCE
#SEE ALSO
#CALLING 
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 14_01_28
#REVISED 14_06_05
#--------------------------------------------
{
  # checking
  rrrobject9(noms,"character");
  if (length(grep("|",noms,fixed=TRUE)) > 0) {
    rrrerreur(noms,"'noms' must not comprise the '|' character");
  }
  rrrobject9(seps,"character");
  rrrobject9(void,"character");
  rrrobject9(rnom,"logical",1);
  # initialization
  res <- character(0);
  # replacing separators with a unique separator
  for (ss in seps) {
    noms <- gsub(ss,"|",noms);
  }
  # getting every elementary component
  for (nn in noms) {
    nom <- strsplit(nn,"|",fixed=TRUE)[[1]];
    if (length(nom)>1) {
      res <- unique(c(res,nom));
    } else {
      if (length(nom)==1) {
        if (nchar(nom) < nchar(nn)) {
          res <- unique(c(res,nom));
        } else {
          if (!rnom) { res <- unique(c(res,nom));}
        }
      }
    }
  }
  # eliminating undesired elementary components
  for (ec in void) {
    if (ec %in% res) {
      res <- res[-which(res==ec)];
    }
  }
  # sorting
  res <- sort(res);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
extract8object <- function(obj,what="v") 
#TITLE returns a vector from a standard object
#DESCRIPTION
# from objects of the style of \samp{documair} package
# returns a named vector of the choosen component
#DETAILS
# The object is supposed to be a list with all its 
# components being a list, every one having as component
# \samp{what}. See the example.\cr
# The type of the vector is deduced from the first component.
#KEYWORDS 
#INPUTS
#{obj} <<The object to be extracted.>>
#[INPUTS]
#{what} <<named vector of the components to extract.>>
#VALUE
# A [possibly named] vector. When \samp{length(obj)==0}
# then \samp{NULL} is returned.
#EXAMPLE
# uu <- list(a=list(u=1,v="A"),b=list(u=2,v="B"));
# extract8object(uu);
# extract8object(uu,"u");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_10_11
#REVISED 13_10_31
#--------------------------------------------
{
  # degenerate case
  if (length(obj)==0) { return(NULL);}
  # some checks
  rrrobject9(obj,"list",
             mensaje="'extract8object' works on 'lists'");
  # looping on the list components
  for (cc in rrrbf(obj)) {
    xx <- obj[[cc]];
    if (is.null(xx[[what]])) {
      print(names(xx));
      cat("The ",cc,"th component of 'obj' doesn't have name ",
          what,"\n",sep="");
      stop();
    } else {
      if (cc==1) {
        res <- xx[[what]];
      } else {
        res <- c(res,xx[[what]]);
      }
    }
  }
  # returning
  names(res) <- names(obj);
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
