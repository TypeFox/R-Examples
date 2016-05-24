am.zandrel <-
function (mdf,df, k, l, v, x, y, cohortparts, components,relmat,ctable) {
  # am.zandrel()
  # construct the animal model matrices z[r],rel[r] from a data frame (df)
  #  or construct z[] but extract rel[] from a list (mdf)
  # setup am object as a list
  #
  # must specify k = no fixed effects, l = no traits
  #              v = no of animal random effects ( var and cov components)
  #          cohortparts = vector of df ccol nos for parts of cohort defining common env
  #	   components = vector listing animal components to be fitted
  #
  # determines m = no individuals in df
  #            n = no of individuals with data
  # return all in list object called am

  m <- length(df$Id)  # no of individuals in pedigree
  n <- nrow(y)        # no of individuals with y data and X codes

#  check pedigree in dataframe valid
  if(pedcheck(df) > 0){
   stop("dmm(): pedigree not valid:\n")
  }

  if(is.null(cohortparts)) {
    cohortcode <- NULL
    ncohortcodes <- 1
  }
  else {
    cohortcode <- make.cohort(df,cohortparts)
    cohortcodes <- unique(cohortcode)  # list of codes
    cohortcodes <- cohortcodes[!is.na(cohortcodes)]  # NA's removed
    ncohortcodes <- length(cohortcodes)
#    ncohortcodes <- length(unique(cohortcode)) counts NA !!
    cat("ncohortcodes = ",ncohortcodes,"\n")
  }
# setup Z,R matrices - for e, a, m, ...  effects
  z <- list(i=matrix(0,n,m),m=matrix(0,n,m),c=matrix(0,n,ncohortcodes))
  rel <- list(e=NULL,a=NULL,aa=NULL,d=NULL,ad=NULL,dd=NULL,s=NULL)
  miss <- match(dimnames(df)[[1]], dimnames(y)[[1]])
#
# z$i - incidence of individuals measured in pedigree - (n x m)
#     - always do z$i
  inrow <- 0.0  # inrow counts individuals with no missing data(or X codes)
  for (i in 1:m) {   # i indexes individuals by row in df
    if (!is.na(miss[i])) {
      inrow <- inrow + 1   # inrow indexes individuals by row in x,y,z
      for (j in 1:m) {
        z$i[inrow,j] <- 0.0
      }
      z$i[inrow,i] <- 1.0
    }
  }
# Additive genetic effects - need rel$a matrix
 if(any(!is.na(match(components,ctable$addg)))){
   # rel$a - additive relationship matrix
  if(relmat == "inline") {
   rel$a <- am.arel(df)
  }
  else if(relmat == "withdf"){
   rel$a <- as.matrix(mdf$rel$a)
  }
 }
# Residual or E effect - use z$i for residual Z matrix
# rel$e - residual relationship matrix - usually I
  if(relmat == "inline"){
    rel$e <- diag(m)
  }
  else if(relmat == "withdf") {
    rel$e <- as.matrix(mdf$rel$e)
  }
# Non additive genetic effects
  # dominance
  if(any(!is.na(match(components,ctable$domg)))){
    rel$d <- as.matrix(mdf$rel$d)
  }
  # epistasis-additive
  if(any(!is.na(match(components,ctable$epiaddg)))){
    rel$aa <- as.matrix(mdf$rel$aa)
  }
  # epistasis-dominance
  if(any(!is.na(match(components,ctable$epidomg)))){
    rel$ad <- as.matrix(mdf$rel$ad)
    rel$dd <- as.matrix(mdf$rel$dd)
  }
  # sex-linked
  if(any(!is.na(match(components,ctable$sexlinaddg)))){
     rel$s <- as.matrix(mdf$rel$s)
  }

# Maternal effects
# z$m - incidence of dams of measured individuals in pedigree - ( n x m)
 if(any(!is.na(match(components,ctable$mat)))){
  inrow <- 0.0
# miss <- match(dimnames(df)[[1]], dimnames(y)[[1]])
  for(i in 1:m) {
    if(!is.na(miss[i])) {
      inrow <- inrow + 1   # inrow indexes individuals by row in x,y,z
      for (j in 1:m) {
        z$m[inrow,j] <- 0.0
      }
      if(!is.na(df$DId[i])){
        z$m[inrow,df$DId[i]] <- 1.0
      }
    }
  }
 }

# cat("z$m:\n")
# print(z$m)
# Common cohort env effect ( same cohort)
# z$c - incidence of measured individuals in cohorts (n x ncohorts)
  if(any(!is.na(match(components,ctable$cohort)))){
    cohortcode <- as.factor(cohortcode)
#  aov() method - discarded due to problems with nrow(z$c)
#   ce.aov <- aov(df$Ymat ~ -1 + cohortcode,x=T)
#   z$c <- ce.aov$x
# match method - similar to z$i and z$m
    inrow <- 0
#   miss <- match(dimnames(cohortcode)[[1]], dimnames(y)[[1]])
  # miss gives mth cohortcode[] matches nth y[] value
    for(i in 1:m){
      if(!is.na(miss[i])) {
        inrow <- inrow + 1 # inrow indexes individuals by row in x,y,z
        for(j in 1:ncohortcodes) {
          z$c[inrow,j] <- 0
        }
#       if(!is.na(cohortcode[i])){
          whichcohort <- match(cohortcode[i], levels(cohortcode))
          z$c[inrow,whichcohort] <- 1.0
#       }
      }
    }

    if(nrow(z$i) != nrow(z$c)) {
      cat("nrows in z$i = ",nrow(z$i),"\n")
      cat("nrows in z$c = ",nrow(z$c),"\n")
      stop("these must be equal:\n")
    }
#    cat("z$c:\n")
#    print(z$c)
# rel$c = pseudo relationship matrix for common env - not needed
  }

# Maternal permanent environmental effect (same dam but not same cohort)
# z$m - dams in pedigree - as above
#  rel$m - pseudo relationship matrix for same dam and not same cohort - not needed

# construct the am list object
  am <- list(m=m,n=n,k=k,l=l,v=v,x=x,y=y,z=z,rel=rel,components=components)
 
  return(am)
}
