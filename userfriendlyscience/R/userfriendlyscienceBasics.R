###########################################################
###########################################################
###
### Collection of very basic functions
###
### File created by Gjalt-Jorn Peters. Questions? You can
### contact me through http://behaviorchange.eu.
###
###########################################################
###########################################################

### Function to remove zero at start of number
noZero <- function (str) {
  return(gsub("0\\.", ".", str));  
}

### Function to format Pearson r
formatR <- function (r, digits) {
  return(noZero(round(r, digits)));
}

### The regular ifelse cannot return objects
ifelseObj <- function(condition, ifTrue, ifFalse) {
  if (condition) {
    return(ifTrue);
  }
  else {
    return(ifFalse);
  }
}

### Basically what Marc Schwartz suggested at Thu Jul 1 19:10:28 CEST 2010
### on the R-help mailing list, see https://stat.ethz.ch/pipermail/r-help/2010-July/244299.html
is.odd <- function(vector) {
  return((vector %% 2) != 0);
}
is.even <- function(vector) {
  return((vector %% 2) == 0);
}

### Case insensitive '%in' variant
`%IN%` <- function(find, table) {
  return(toupper(find) %in% toupper(table));
}

### Paste0 but then immediately displaying on screen
cat0 <- function(..., sep="") {
  return(cat(..., sep=sep));
}

### Check whether elements are true, specifying how 'NA' should be seen
isTrue <- function(x, na = FALSE) {
  naValues <- ifelse(rep(na, length(x)),
                     is.na(x),
                     rep(FALSE, length(x)));
  return(ifelse(is.na(x), naValues, x==TRUE));
}

### Check whether something is a number
is.nr <- function(x) {
  if (!is.null(x)) {
    if (!is.na(x)) {
      if (is.numeric(x)) {
        return(TRUE);
      }
    }
  }
  return(FALSE);
}
