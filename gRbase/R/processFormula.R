processFormula <- function (formula, data, marginal, type = c("Discrete", "Continuous"),
    v.sep = "*", g.sep = "+")
{
    get.var.of.type <- function(type) {
        varNames(data)[varTypes(data) == type]
    }
    used.var <- get.var.of.type(type)


    if (!inherits(formula, "formula")) {
      formula <- list2rhsFormula(formula)
    }

    list.formula <- rhsFormula2list(formula)

    pow <- extract.power(formula)

    if (!is.numeric(pow)){

      ## Check if abbreviations are valid
      ##
      if (any(is.na(pmatch(unlist(list.formula), used.var, duplicates.ok=TRUE))))
        stop("An invalid variable specification has been found\n")

      ## Replace abbreviations with variable names
      ##
      list.formula <-
        lapply(list.formula,
               function(x) {
                 i <- pmatch(x,used.var)
                 used.var[i] })

      formula     <- list2rhsFormula(list.formula)
      str.formula <- paste(deparse(formula[[2]]), collapse = "")
    } else {
      if (!missing(marginal)) {
        used.var <- intersect(marginal, used.var)
      }
      if (pow == -1)
        str.formula <- paste(used.var, collapse = v.sep, sep = "")
      else {
        pow <- min(c(pow, length(used.var)))
        tmp <- selectOrder(used.var, pow)
        str.formula <- paste(unlist(lapply(tmp, paste, collapse = v.sep)),
                      collapse = g.sep, sep = "")
      }
      formula      <- formula(paste("~", str.formula, sep = ""))
      list.formula <- rhsFormula2list(formula)
    }


    num.formula <- lapply(list.formula, function(l) {
        charmatch(l, used.var)
    })

    value <- list(formula = formula, str.formula = str.formula, num.formula = num.formula,
        list.formula = list.formula, gmData = data, varnames = used.var)
    value
}


## Turn a right-hand-sided formula into a list (anything on the left hand side is ignored)
##
## January 2011

rhsFormula2list <- rhsf2list <- function(f){
    if ( is.character( f ) ){
        list( f )
    } else {
        if ( is.numeric( f ) ){
            lapply( list( f ), "as.character" )
        } else {
            if ( is.list( f ) ){
                lapply(f, "as.character")
            } else {
                ## We assume a formula...
                ## FIXME: Was:   ##.xxx. <- f[[2]]
                ## Changed to
                .xxx. <- f[[ length( f ) ]]
                f1 <- unlist(strsplit(paste(deparse(.xxx.), collapse="")," *\\+ *"))
                f2 <- unlist(lapply(f1, strsplit, " *\\* *| *: *| *\\| *"),recursive=FALSE)
                f2
            }
        }
    }
}



## Turn list into right-hand-sided formula
##
## July 2008
list2rhsFormula <- list2rhsf <- function(f){
  if (inherits(f,"formula"))
    return(f)
  as.formula(paste("~",paste(unlist(lapply(f,paste, collapse='*')),collapse="+")),
             .GlobalEnv)
}

allSubsets <- function(x,g.sep="+"){
  if (length(x)==1)
    return(x)
  else {
    val <- x[1]
    for (i in 2:length(x)){
      v <- paste(val,x[i],sep=g.sep)
      val <- c(val,x[i],v)
    }
    val <- strsplit(val,paste("\\",g.sep,sep=""))
    return(val)
  }
}

selectOrder  <- function(x,order=2){
  v <- allSubsets(x)
  value <- v[lapply(v,length)==as.numeric(order)]
  return(value)
}

extract.power<-function(fff){
  mimf  <- paste(as.formula(fff))[2]
  mimf.split <- unlist(strsplit(mimf,""))
  if(length(grep("[a-z]", mimf))>0){
    pow <- mimf
  } else {
    has.hat <- match("^",mimf.split)
    sub <- unlist(strsplit(mimf,"\\^"))
    if (!is.na(has.hat)){
      pow <- ifelse (sub[2]==".", -1, as.numeric(sub[2]))
    }
    else {
      pow <- length(unlist(strsplit(sub,"\\.")))
    }
  }
  return(pow)
}



















## ## SHD, July 2008
## rhsFormula2list <- rhsf2list <- function(f){
##   if (inherits(f,"list")){
##     return(f)
##   } else {
##     if (inherits(f,"character")){
##       return(list(f))
##     }
##   }
## ##   rhs <- paste(f)[2]
## ##   f1 <- strsplit(rhs,"\\+")[[1]]
##   .xxx. <- f[[2]]
##   f1 <- unlist(strsplit(paste(deparse(.xxx.), collapse="")," *\\+ *"))
##   f2 <- unlist(lapply(f1, strsplit, "\\*|:"),rec=FALSE)
##   f2 <- lapply(f2, function(x) gsub(" +","",x))
##   f2
## }


# Notes, functions and examples for generator lists (model formulae
# for hierarchical loglinear models) in R.
# Includes functions dual.rep, add.edge, delete.edge and ..is.graphical.
# I havent tried to optimise the functions in any way, merely to get
# versions which work.
# David Edwards, 12.5.2004.
# adapted to gRbase by Claus Dethlefsen 16.08.04

# The approach could also be used (with some extra work) for
# hierarchical mixed 'mim' models,
# but with 'extended' hierarchical mixed models there is a problem of
# how to handle quadratic terms like X^2.
# Representing such generators as vectors would seem problematic.

# nb: there is an implicit assumption that all models contain at least
# all main effects, eg. the minimal model
# for variable set A,B,C has formula A+B+C.

# ----------------------------------------------------------
# ----------------------------------------------------------
#                         Generators
# ----------------------------------------------------------
# ----------------------------------------------------------

# A generator is implemented as a vector (regarded as a set)
#
#g1 <- c("a", "bc", "x")
#g2 <- c("a", "x")
#g3 <- c("sex", "Age")
#g4 <- 1:5
#
# Thus the following built-in R commands may be used with generators:
#
#    union(g1, g2)
#    intersect(g1, g2)
#    setdiff(g1, g2)
#    setequal(g1, g2)
#    is.element(g1, g2)
#
# setdiff(g1,g2) returns g1 \ g2.
# Comparing two generators, is.element(g1, g2) returns a boolean
# vector of the same length as g1, indicating whether the element of
# g1 is contained in g2. Thus
# g1 <= g2 <=> all(is.element(g1, g2))
# g1 == g2 <=> setequal(g1, g2)
#



# A function to write a generator as a string:
#
..showg <- function(g, v.sep="*") {
  if (length(g)==0)
    s<-'<empty>'
  else {s <- g[1];
        if (length(g)>1)
          for (i in 2:length(g))
            s <- paste(s, g[i], sep=v.sep)
      }
  s
}

# Sometimes we may need the empty set
#g5 <- vector()
#..showg(g5)

# A function to read a generator as a string.
# nb: s is a character (vector), but must have length one.
#
..readg <- function(s, v.sep="*") {
  g <- character(0)
  s <- paste(s, v.sep, sep="") # add a separator
  s <- gsub(" ","",s) # strip spaces
  k1 <- 1
  for (k in ((k1+1):nchar(s))) {
    if (substring(s,k,k) == v.sep) {
      g <- c(g, substring(s,k1,(k-1)))
      k1 <- k+1 }
  }
  g
}

# -------------------------------------------------------------------------------------------------------
#               Generator Lists
# -------------------------------------------------------------------------------------------------------

#f1 <- list(g1, g2, g3)
#f2 <- list(1:10, c(2,3,5), c(3,5,7))


# Function to write a generator list as a formula
showf <- function(f, g.sep="+", v.sep="*") {
   if (length(f)==0) s <- '<empty formula>' else {
   s <- ..showg(f[[1]])
   if (length(f)>1)
     for (i in 2:length(f))
       s <- paste(s, ..showg(f[[i]],v.sep), sep=g.sep)}
   s
}

# Function to read a generator list as a string
# nb: length of string must be one

readf <- function(s, v.sep="*", g.sep="+") {
  gens <- ..readg(s, g.sep)
  l <- list(length(gens))
  for (i in 1:length(gens)) l[[i]] <- ..readg(gens[i], v.sep)
  l
}

#f3 <- readf("A.B.C+B.C.D+C.D.E")
#showf(f3)

# To get the union of all generators i a list l:

..varset <- function(f) unique(unlist(f))

# To find out whether a generator g is contained in (is a subset of
# some element of) a list l:

..in.list <- function(g, l)
  any(unlist(lapply(l, function(xx) all(is.element(g, xx)))))

# To find out whether a generator g contains an element of a list l:

# any(unlist(lapply(l, function(xx) any(is.element(xx, g))))

# A function to find out whether the k.th generator in a list
# is contained in any of the others:

..is.cont <- function(k, l) {
  g <- l[[k]]
  a <- sapply(l, function(xx) all(is.element(g, xx)))
  a[k] <- FALSE
  any(a)
}

# A function to find out whether the k.th generator in a list
# contains any of the others

..contains <- function(k, l) {
  g <- l[[k]]
  a <- unlist(lapply(l, function(xx) all(is.element(xx, g))))
  a[k] <- FALSE
  any(a)
}






# A function to return dual representation for a list.
# See description in Edwards & Havranek, Biometrika (1985), 72, 2, p.341.

# minimal=T: returns dual representation given usual, ie list of
# minimal generators not contained in a generator in the input list.
# minimal=F: returns usual representation given dual, ie list of
# maximal generators not containing a generator in the input list.

## old version gave list()
#dual.rep <- function(glist, S, minimal=TRUE) {
# # S is total varset - often but by no means always given by
# # unique(unlist(g.list))
# list.save <- list()
# if (length(glist)>0) for (v in 1:length(glist)) {
#   m1 <- list.save
#   if (minimal) m2 <- as.list(setdiff(S,glist[[v]])) else m2 <- as.list(glist[[v]])
#   if (v==1) list.save <- m2 else {
#      list.save <- remove.redundant(unlist(lapply(m1, function(g) lapply(m2, union, g)),recursive=FALSE),FALSE)}}
# if (!minimal) list.save <- lapply(list.save, function(g) setdiff(S, g))
# list.save
#}


#m1 <- readf('A.B+A.C')
#showf(m1)
#showf(dual.rep(m1, ..varset(m1)))

#m2 <- readf('B.D+A.D+C.D')
#showf(m2)
#showf(dual.rep(m2, ..varset(m2)))

# The dual of the dual should be the same as the original
#showf(dual.rep(dual.rep(m2, ..varset(m2)), ..varset(m2), F))

# Function to delete 'edge' from a generator list, by (i) converting
# generator list to dual representation,
# (ii) appending the 'edge', (iii) converting back to usual
# representation. 'Edge' is given as vector of length 2:
# it can also have length >2, ie. be a higher-order interaction.

.delete.edge <- function(m, edge) {
  S <- ..varset(m)
  dr <- dual.rep(m, S)
  dr <- c(dr, list(edge))
  removeRedundant(dual.rep(dr, S, FALSE))
}

#m2 <- readf('B.D+A.D+C.D')
#showf(m2)
#showf(delete.edge(m2, c('A','B')))
#showf(delete.edge(m2, c('A','D')))
#m3 <- readf('A.B.C.D')
#showf(delete.edge(m3, c('A','B','C')))

# Function to add 'edge' from a generator list, by (i) converting
# generator list to dual representation,
# (ii) removing the 'edge', (iii) converting back to usual
# representation. 'Edge' is given as vector of length 2:
# it can also have length >2, ie. be a higher-order interaction.

.add.edge <- function(m, edge) {
  S <- ..varset(m)
  dr <- dual.rep(m, S)
  k <- length(dr)
  if (k>0) {for (i in 1:k) if (setequal(dr[[i]], edge)) dr[[i]] <- vector()}
#  if (k>0) {for (i in 1:k) if (setequal(dr[[i]], edge)) dr[[i]] <- NULL}
  dr <- removeRedundant(dr, FALSE)
  dual.rep(dr, S, FALSE)
}

#m2 <- readf('B.D+A.D+C.D')
#showf(m2)
#showf(add.edge(m2, c('A','B')))
#showf(add.edge(m2, c('A','C')))


# From main effect model on 5 vertices
#showf(add.edge(1:5, c(1,2)))
#showf(add.edge(c('a','b','c','d','e'), c('a','b')))

# From saturated model on 5 vertices
#showf(delete.edge(list(1:5), c(1,2)))
#showf(add.edge(list(c('a','b','c','d','e')), c('a','b')))

# Exploiting the fact that a model is graphical iff all its dual
# generators have length 2,
# we get a neat function for graphicalness:

..is.graphical <- function(m) {
   dr <- dual.rep(m, ..varset(m))
   lengths <- lapply(dr, length)
   all(lengths==2)
}

#..is.graphical(readf('A.B+A.C+B.C'))
#..is.graphical(readf('A.B.C'))
#..is.graphical(readf('A.B+B.C'))




########################################################################
###
### Functions removed from here because faster versions exist elsewhere
###
########################################################################


# dual.rep <- function(glist, S, minimal=TRUE) {
#  # S is total varset - often but by no means always given by unique(unlist(g.list))
#  list.save <- list()
#  #if (length(glist)==0) list.save <- list(S)
#  if (length(glist)==1 & is.logical(glist[[1]])) list.save <- list(S)
#  else {
#    for (v in 1:length(glist)) {
#      m1 <- list.save
#    if (minimal) m2 <- as.list(setdiff(S,glist[[v]])) else m2 <- as.list(glist[[v]])
#    if (v==1) list.save <- m2 else {
#       list.save <- removeRedundant(unlist(lapply(m1, function(g)
#                                                   lapply(m2, union,
#                                                          g)),recursive=FALSE),FALSE)}}
#  if (!minimal) list.save <- lapply(list.save, function(g) setdiff(S,
#                                                                   g))}
#  list.save }


## subsetof <- function(g1, g2) all(is.element(g1, g2))
## SHD: Faster version in setopsR.R

# A function to remove redundant generators.  If maximal=T, returns
# the maximal generators, if =F, the minimal generators.  This seems
# difficult to vectorize: any suggestions?  This function is a prime
# candidate for a C routine.

## SHD: Faster version in setopsC.R

# remove.redundant <- function(f, maximal=TRUE) {
#   k <- length(f)
#   new.f <- f
#   if (k>1) {
#     for (i in 1:(k-1)) {
#       g1 <- f[[i]]
#       if (length(g1)>0) {
#         for (j in (i+1):k) {
#           g2 <- f[[j]]
#           if (length(g2)>0) {
#            if (setequal(g1,g2)) f[[j]]<-vector() else {
#              kk <- 0
#              if (subsetof(g1, g2)) {if (maximal) kk <- i else kk <- j}
#              if (subsetof(g2, g1)) {if (maximal) kk <- j else kk <- i}
#              if (kk>0) f[[kk]] <- vector()
#       } # else
#      }} # length(g2)>0; for j in (i+1):k
#     }}  # length(g1)>0; for i in 1:(k-1)
#     f <- f[lapply(f, length)>0]
#   } # if k>1
#   f
# }
















## .processFormula <- function(formula, data, marginal,
##                            type=c("Discrete","Continuous"), v.sep="*", g.sep="+"){

##   get.var.of.type <- function(type){varNames(data)[varTypes(data)==type]}

##   if (!inherits(formula,"formula")){
##     formula <- list2rhsFormula(formula)
##   }


##   used.var <- get.var.of.type(type)
##   pow <- extract.power(formula)
## #  print(pow); print(used.var)

##   if (is.numeric(pow)){
##     if (!missing(marginal)){
##       used.var <- intersect(marginal,used.var)
## #      print(used.var); print(marginal)

##     }
##     if (pow==-1)
##       mimf <- paste(used.var,collapse=v.sep,sep="")
##     else{
##       pow <- min(c(pow, length(used.var)))
##       tmp <- selectOrder(used.var, pow)
##       mimf <- paste(unlist(lapply(tmp, paste, collapse=v.sep)),collapse=g.sep,sep="")

##     }
##   } else {
##     mf    <- as.formula(formula)
##     mimf <-  paste(deparse(mf[[2]]), collapse="")
##   }
##   cat("mf:\n"); print(mf)
##   cat("mimf:\n"); print(mimf)
##   formula <- formula(paste("~",mimf,sep=""))

##   interactions <- strsplit(mimf,paste("\\",g.sep,sep=""))[[1]]
##   interactions <- gsub(" ","",interactions,fixed=TRUE)
##                                         #  interactions <- gsub(g.sep,"",interactions)

##   if (v.sep == "*") v.sep <- "[*|:]"
##   varformula <- strsplit(interactions, v.sep)

##   numformula   <- lapply(varformula, function(l){ match(l,used.var) })

##   value <- list(formula=formula, mimformula=mimf,
##                 numformula=numformula,
##                 varformula=varformula,
##                 gmData=data, varnames=used.var)
##   value
## }


