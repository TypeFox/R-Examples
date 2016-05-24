### translate ISO 3166-1 alpha-2 codes to country names
### default result is a character vector of same length, which may contain regular expressions
iso.expand <- function(a,regex=TRUE){
  iso3166 <- get("iso3166")
  AA <- toupper(a)
  if (all(nchar(a)==2)) codes <- iso3166$a2
  else if (all(nchar(a)==3)) codes <- iso3166$a3
  else stop("All codes must be equal length, 2 or 3 characters.")
  nn <- lapply(AA,function(x) iso3166$mapname[which(codes == x)])
  if (regex) {
    nn2 <- lapply(nn,function(x) if (length(x)>1) paste("(^",x,")",sep="",collapse="|") else x)
    unlist(nn2)
  } else unlist(nn)
}

### list all countries that fall under sovereignty of 'sov'
sov.expand <- function(sov,regex=TRUE){
  iso3166 <- get("iso3166")
  sov <- tolower(sov)
  sel <- tolower(iso3166$sovereignty)
  nn <- lapply(sov,function(x) iso3166$mapname[which(sel == x)])
  if (regex) {
    nn2 <- lapply(nn,function(x) if (length(x)>1) paste("(^",x,")",sep="",collapse="|") else x)
    unlist(nn2)
  } else unlist(nn)
}

### the inverse:
### The subtlety lies in dealing with special cases like "China:Hong Kong" which has "HK", not "CN"
### But I think it's OK if "China" returns "CN", not "CN"+"HK"+"MC".
### Similar for Norway/Svalbard, Finland/Aland
### But of course, "China:Hong Kong" must return "HK" etc.

### we does this with a reverse grep
### downsides: the names in x must be complete
iso.alpha <- function(x,n=2){
  iso3166 <- get("iso3166")

## part 1: reverse fit will find all special cases, but names in vector x must be complete
  nam1 <- lapply(seq_along(iso3166$mapname),
                    function(nn) {regex <- paste("(^",iso3166$mapname[nn],")",sep="")
                                  ttt <- grep(regex, x, perl=TRUE, ignore.case=TRUE);
                                  if (length(ttt)>0) cbind(nn,ttt) else NULL})
  fli <- do.call(rbind,nam1)
  sel <- fli[match(seq_along(x),fli[,2]),1]

## part 2:
## try for partial fit. if it gives a single result, use it.
  if (any(is.na(sel))) {
    sel2 <- which(is.na(sel))
    nam2 <- unlist(lapply(x[sel2],
              function(nn) {regx <- paste("(^",nn,")",sep="") ; 
                            ttt <- grep(regx, iso3166$mapname, perl=TRUE, ignore.case=TRUE);
                            if (length(ttt)==1) ttt else NA}))
    sel[sel2] <- nam2
  } 

  if (n==2) iso3166$a2[sel]
  else if (n==3) iso3166$a3[sel]
  else stop("n must be 2 or 3.")
}

