
`multmoments` <- 
function (moment, current.matrix, current.cell, moment.rep, row_col) 
{
# A recursive function to compute the representation of a multivariate moment
#    using upper - triangular matrices
#
# moment: a vector of integers representing the moment
#     eg, c(3, 2, 4) for a 3-dimensional normal vector
#     corresponding to the moment E[(X1^3)(X2^2)(X3^4)]

# current.matrix: input/output matrix under consideration in recursion
#     this is an upper-triangular integer matrix (l(i, j))
#     organized by row, with   lmom * (lmom + 1) / 2   total elements  
#
# current.cell: input/output cell in current matrix under consideration in recursion

# moment.rep: input/output current set of representations
#     function adds each satisfying matrix to moment.rep

# row_col: 2x(lmom * (lmom + 1) / 2 matrix giving rows and columns for square matrix
#     for each cell so that they don't have to be calculated each time

#  algorithm:
#  loop through cell values from 0 to min(moments[row] - rowsum, moments[col] - colsum)
#    if the this new matrix satisfies moment criterion, 
#       add to moment.rep and return
#    if the current matrix is too great in any dimension, 
#       return
#    if the current matrix is < moment in any dimension, and
#       at most moment(i) for all other indexs i, continue

summomentmatrix <- function (moment.matrix) 
{
#  compute the row/col sums of moment.matrix
#  uses Matrix package
#  construct a square upper diagonal matrix from moment.matrix

#  length of moment based on representation:
lmom <- (sqrt(8 * length(moment.matrix) + 1) - 1) / 2 
#  make a square matrix for use with row and column sum functions
tempmatrix <- matrix(c(rep(0, lmom^2)), nrow=lmom)
tempmatrix[1, ] <- moment.matrix[1:lmom]
endrow <- lmom
if (lmom > 1)
  {
  for (irow in 2:lmom)
     {startrow <- endrow + 1
      endrow <- endrow + (lmom - irow + 1)
      tempmatrix[irow, ] <- cbind(c(rep(0, (irow - 1)), moment.matrix[startrow:endrow])) }
     
  summomentmatrix <- colSums(tempmatrix) + rowSums(tempmatrix) 
  }
if (lmom == 1){summomentmatrix <- 2 * sum(moment.matrix)}     

return(summomentmatrix)}


lmom <- length(moment)
totcells <- lmom * (lmom + 1) / 2   # total cells in a moment representation
 
thisrow <- row_col[1, current.cell]
thiscol <- row_col[2, current.cell]
 
moment.row <- moment[thisrow]
moment.col <- moment[thiscol]
rowcells <- row_col[1, row_col[1, ] == thisrow]
colcells <- row_col[2, row_col[2, ] == thiscol]
 
#  determine sums for row and columm
#  maxvalue is the largest value that can be added to a
#    a row/column sum and still be no more than criterion
#  maxvalue will be the minimum of the moments minus these sums

rowsum <- summomentmatrix(current.matrix)[thisrow]
colsum <- summomentmatrix(current.matrix)[thiscol] 
 
maxvalue <- min(moment[thisrow] - rowsum, moment[thiscol] - colsum)
 
for (ivalue in 0:maxvalue)
   {         # for:loop through all possible values for cell
 
   current.matrix[current.cell] <- ivalue
   current.sum <- summomentmatrix(current.matrix)
             # determine current sum for this matrix

#  if matrix fulfills criterion, add to reps and return
   if (sum(moment == current.sum) == lmom)
      {moment.rep <- rbind(moment.rep, current.matrix)
       return(moment.rep)}

#   if the sum is too large, return because any other sum will also be too large
   if (sum(current.sum > moment) > 0){return(moment.rep)}

#   if at least one term in current sum is smaller than moment, 
#       go down one cell unless this is the last cell

   if ((sum(current.sum == moment) < lmom) & (current.cell != totcells))     
      {cc <- current.cell + 1
# recursive step:
       moment.rep <- multmoments(moment, current.matrix, cc, moment.rep, row_col)}
                        
   }        # end of for 
 
#   note: the three conditions above are exhaustive except for last cell 
#      In that case there is nothing to return
return(moment.rep)
}

`callmultmoments` <- 
function (moment) 
{
#  function to compute the representation of a multivariate moment
# 
#  moment: a vector of integers representing the moment
#          eg, c(3, 2, 4) for a 3-dimensional normal vector
#          corresponding to the moment (E[X1^3)(X2^2)(X3^4)]
#  sum of exponents must be even; otherwise moment is 0

#  returns a list of 3 components:
#    1: $moment  -  the input moment vector 
#    2: $representation  -  a matrix containing the representation 
#       in terms of upper-triangular matrices
#    3: $coefficients  -  the coefficients corresponding to the rows of the representation
#  if sum odd, returns  -1 and prints "Sum of powers is odd. Moment is 0."
#  if any component is negative, returns  -2 and prints "All components of the moment must be non-negative."
#  if any component is not an integer, returns  -3 and prints "All components of the moment must be integers."

subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}

nestedreps <- function (input.vector, inner.rep, outer.rep) 
{
# replicates input.vector, first by inner.rep, then by outer.rep

# input.vector: vector to replicate
# inner.rep:    count of replicates of elements
# outer.rep:    count of replicates of resulting vector

temp.vector <- NULL
for (ivec in 1:length(input.vector))
     {temp.vector <- c(temp.vector, rep(input.vector[ivec], inner.rep))}

total.vector <- rep(temp.vector, outer.rep) 

return(total.vector)}


mrepnames <- function (ndim) 
{
# get colnames for a representation
combs <- expand.grid(1:ndim, 1:ndim)
char.combs <- paste("S(", combs[, 2], ", ", combs[, 1], ")")
for (elem in 1:(ndim^2))
  {char.combs[elem] <- subblank(char.combs[elem])}
m <- sort(matrix((1:ndim^2), nrow=ndim, byrow=TRUE)[ !lower.tri(matrix((1:ndim^2), nrow=ndim, byrow=TRUE))])
return(char.combs[m])
}


lmom <- length(moment)
lrep <- (lmom^2 + lmom) / 2  # length of representation using upper-triangular matrices
moment.rep <- matrix(rep(0, lrep), nrow=1)   # initial null representation to be augmented
if (sum(trunc(moment) == moment) < length(moment)) 
   {print("All components of the moment must be integers.") 
   return( -3)}
if (sum(moment<0) > 0)
  {print("All components of the moment must be non-negative.")
   return( -2)}
if (trunc(sum(moment) / 2) != sum(moment) / 2)
  {print("Sum of powers is odd. Moment is 0.")
   return( -1)}
     
icells <- matrix(rep(0, lmom^2), nrow=lmom)
m <- matrix(rep(1, lmom^2), nrow=lmom) 
limits <- c((lmom:1)%*%(m * !(lower.tri(m))))
# sum of row lengths: lmom, lmom + lmom - 1, ... , for use in row_col
 
row_col <- matrix(rep(0, (2 * lrep)), nrow=2)
#  2x(nm * (nm + 1) / 2 matrix giving rows and columns for each cell
#       so that they don't have to be calculated each time in multmoment

for (icell in (1:lrep))
   {
   row_col[1, icell] <- min((1:lmom)[icell<=limits])
   if (row_col[1, icell] == 1){row_col[2, icell] <- icell}
   if (row_col[1, icell]>1){row_col[2, icell] <- icell - limits[row_col[1, icell] - 1] + 
                                             row_col[1, icell] - 1 }
  }

# initial current.matrix and current.cell
current.matrix <- c(rep(0, lrep))  
current.cell <- 1
 
# call recursive function to determine upper-triangular representations

moment.rep <- multmoments(moment, current.matrix, current.cell, moment.rep, row_col)
if (dim(moment.rep)[1] == 2)    # get rid of initial 0 matrix representation
   { moment.rep <- matrix(moment.rep[2, ], nrow=1) }
if (dim(moment.rep)[1] > 2)
   {moment.rep <- moment.rep[2:dim(moment.rep)[1], ]}
rownames(moment.rep) <- 1:(dim(moment.rep)[1])

##################################################################
# now determine coefficients for upper-triangular representations

l.representation <- moment.rep
lmom <- length(moment)
nrep = dim(l.representation)[1]  
totlength <- lmom^2
rep.coefficients <- c(1:nrep)  # coefficients corresponding to nrep representations

#  multiplier for all terms
overallcoeff <- ((1 / 2)^(sum(moment) / 2)) * prod(factorial(moment)) / factorial((sum(moment) / 2)) 
 
for (irep in 1:(dim(l.representation)[1]))
  {
#  loop through all matrices
  thisrep <- l.representation[irep, ] 

#  determine the coefficient for each term based on switching equivalent terms

#  "base" gives the number of switches that can be made to each element of the l-matrix
#  diagonal elements are not switchable, but are included to allow subtraction below

  base <- c(rep(1, lmom * (lmom + 1) / 2))
  base[1] <- 1   # first diagonal element  -  not switchable

  totreps <- 1   # total number of transpostions
# if there is only one element, it must be the diagonal, so is not switchable  -  skip
if (lmom > 1){
  base[1] <- 1   # first diagonal element
  for (cell in 2:length(base))
    {
     icol = row_col[1, cell]  #  determine if diagonal element
     irow = row_col[2, cell]
     if (irow == icol){base[cell] <- 1}  # diagonal  -  not switchable
     if (icol != irow)
         {totreps <- totreps * (1 + thisrep[cell])
         base[cell] <- 1 + thisrep[cell]} 
    }  #  done with computing base and total transpositions (totreps)
  }
  
mcoeff <- 1  #  sum of multinomial coefficients
if (totreps > 1){ 

#  baserep will represent the lower diagonal (including diagonal)
#  of the augmented matrices
     baserep <- matrix(rep(0, totreps * length(base)), nrow=totreps)
     basegt1 <- base[base>1]
     nbase <- 0
     for (ibase in 1:length(base))
       {if (base[ibase] > 1)
           {nbase = nbase + 1
            if (nbase == 1){baserep[, ibase] <- nestedreps(c(0:(basegt1[nbase] - 1)), 1, totreps / prod(basegt1[1:nbase])) }    
            if (nbase > 1) {baserep[, ibase] <- nestedreps(c(0:(basegt1[nbase] - 1)), prod(basegt1[1:(nbase - 1)]), totreps / prod(basegt1[1:nbase])) }
           }
       }
 
#  now go through each transposition
     if ( !is.na(totreps) & totreps != 1)
    {
     mcoeff <- 0
     for (jrep in (1:totreps)) # check each transposition
          {newrep <- baserep[jrep, ]       #  added lower diagonal elements
           addrep <- sort(newrep, decreasing=TRUE)[1:(lmom * (lmom - 1) / 2)] 
           fulnrep <- c((thisrep - newrep), addrep)
           thiscoeff <- ((length(fulnrep))^sum(fulnrep)) * dmultinom(x=fulnrep, prob=rep(1.0, length(fulnrep)))
           mcoeff <- mcoeff + thiscoeff
#  the multinomial coefficient is obtained from the multinomial distribution
#  multiply by an appropriate power to get rid of probability
        }  
    }
}  
  if (is.na(totreps)){mcoeff <- 1}
if (totreps == 1)
   {mcoeff <- (length(thisrep))^sum(thisrep) * dmultinom(x=thisrep, prob=rep(1.0, length(thisrep)))}

#  determine full coefficient  -  round because all coefficients should be integers
#                                (Note - this statement has not been proved)
rep.coefficients[irep] <- round(overallcoeff * mcoeff)

  cell <- 0
  for (irow in (1:length(moment)))
    {
    for (icol in (irow:length(moment)))
      {cell <- cell + 1
 #   exponent of term
       expo <- l.representation[irep, cell]
 
       }
    }
  }
output <- list(moment, moment.rep, rep.coefficients)
names(output) <- c('moment', 'representation', 'coefficients')
colnames(output$representation) <- mrepnames(length(moment)) 
names(output$coefficients) <- paste("rep", (1:length(output$coefficients)))
class(output) <- "moment" 
return(output)}

`toLatex.moment` <- 
function (object, ...) 
{
#  build latex code for the l-matrix representation of the moment
#  object is the representation of the l-matrices for moment
#  each row is such an l-matrix

# object: list from callmultmoment (class moment)
#     with first component the moment itself, 
#     the second component the set of upper-triangular
#          representations of the moment, 
#     and third component, their correpsonding coefficients

#  note that Latex backslashes are doubled to allow writing to file 
#      with writeLines

subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}
sortmatrix <- function (xmat) 
{
# sort matrix xmat by successive columns
# output is the index vector, not the sorted matrix

ncol <- dim(xmat)[2]
 
sortstring <- "order("
for (icol in 1:ncol)
  {
sortstring <- paste(sortstring, "xmat[, ", icol, "]")
if (icol != ncol){sortstring <- paste(sortstring, ", ")}
  }
sortstring <- paste(sortstring, ", decreasing=TRUE)")
matindex <- eval(parse(text=sortstring))

return(matindex)}


moment.fullrep <- object

maxchars.latex <- 300
# maximum latex characters to print on one line
# this includes formatting characters (subscript notation, etc)

#  extract the components for convenience
moment <- moment.fullrep[[1]]
moment.rep <- moment.fullrep[[2]][sortmatrix(moment.fullrep[[2]]), ]
#       latex representation will be in sorted order
coefficients.rep <- moment.fullrep[[3]][sortmatrix(moment.fullrep[[2]])]


doubquote <- subblank(paste("\\", "\\"))

if ( !is.matrix(moment.rep)){moment.rep <- matrix(moment.rep, nrow=1)}
numrep <- dim(moment.rep)[1]

latex.moment <- rep(" ", numrep + 1)  
 
#   write the left hand side, ie, E[X1 ... Xn] =

latex.moment[1] <- "E["

for (imoment in(1:(length(moment))) )
  {latex.moment[1] <- paste(latex.moment[1], 
         "X_{", imoment, "}^{",  moment[imoment],  "}")  }

latex.moment[1] <- paste(latex.moment[1], "] =", doubquote)
latex.moment[1] <- subblank(latex.moment[1])
 
#  write the right hand side, ie, the set of terms

if (sum(moment==0) == length(moment))
  {latex.moment[2] <- "1"
   return(latex.moment)}

totchars <- 0        # used with totchars.latex
for (irep in (1:numrep))
  {
   mcoeff <- coefficients.rep[irep]
   thisrep <- moment.rep[irep, ] 

  if (mcoeff != 1){latex.moment[irep + 1] <- as.character(mcoeff)}
#              omit coefficient if it is 1
  cell <- 0
  for (irow in (1:length(moment)))
    {
    for (icol in (irow:length(moment)))
      {cell <- cell + 1
#   exponent of term, that is, the "l" value
       exponent <- moment.rep[irep, cell]
 
#   if exponent is 1, omit it as obvious
       if (exponent == 1){latex.moment[irep + 1] <- paste(latex.moment[irep + 1], 
          "\\sigma_{ ",  irow,  ", ",  icol,  "}")}

        if (exponent > 1){latex.moment[irep + 1] <- paste(latex.moment[irep + 1], 
          "\\sigma_{ ",  irow,  ", ",  icol, "}^{",  exponent, "} ")}
       }  # end of for
     }  # end of for
  totchars <- totchars + nchar(latex.moment[irep + 1])
  if (irep < numrep){latex.moment[irep + 1] <- paste(latex.moment[irep + 1], " + ")}
  if (totchars > maxchars.latex) 
      {totchars <- 0
      latex.moment[irep + 1] <- paste(latex.moment[irep + 1], doubquote)}
      latex.moment[irep + 1] <-  subblank(latex.moment[irep + 1])
  }  # end of for

return(latex.moment)}


`simulate.moment` <- 
function(object, nsim, seed=NULL, Mean, Sigma, ...){

# function: method to calculate moment of the multivariate normal distribution
#           using Monte-Carlo integration (Rizzo, 2008)
# object is an object of class moment
# nsim is the number of samples to generate
# seed is the seed for the random number generator
# Mean is the mean of the (X1, ..., Xn)
# Sigma is the variance-covariance of (X1^k1, ..., Xn^kn), dimension nXn


# requires package mvtnorm for function rmvnorm

moment.fullrep <- object
if (is.numeric(seed)){set.seed(seed)}

if (class(moment.fullrep) == "moment"){thismoment <- moment.fullrep$moment}
if (class(moment.fullrep) != "moment")
   {print("moment must be of class 'moment'")
    return(-1)}   
    
ndim <- length(thismoment)                                                                                        
sample <- rmvnorm(n=nsim, mean=Mean, sigma=matrix(Sigma,nrow=length(Mean)))
exponents <- matrix(rep(thismoment, nsim), nrow=nsim, byrow=TRUE)
powers <- sample^exponents
prods <- rep(1, nsim)        #  calculate product of powers of Xs
for (icol in (1:ndim))
   {prods <- prods * powers[, icol]}
moment.value <- mean(prods)
 
return(moment.value)}



`evaluate` <- function(object, sigma) 
   {UseMethod("evaluate",object)}

`evaluate.moment` <- 
function (object, sigma) 
{
#  evaluate the moment using the representation from callmultmoment
#      at the upper-triangular value of sigma
  
#  object from callmultmoment
#      list with first component the moment itself
#      the second component the set of upper-triangular 
#      matrices representing the moment
#      and third component, their corresponding coefficients
#
#  sigma is the upper-triangular matrix of covariance terms
#      at which the moment is to be evaluated
#
#  returns the value of the moment at this sigma

moment <- object[[1]]
moment.rep <- object[[2]]
coefficients.rep <- object[[3]]

#  evaluate the moment by adding the value at each representation
#  this is the product of all sigma[i, j]^l[i, j] 
#     if sigma and l are thought of as square matrices and l is the representation

moment.value <- 0
for (irep in 1:(dim(moment.rep)[1]))
  {moment.value <- moment.value + coefficients.rep[irep] * prod(sigma^moment.rep[irep, ])}

return(as.vector(moment.value))}



`print.moment` <- 
function(x, ...){

# function: method to print a moment of the multivariate normal distribution
# x is an object of class moment

subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}



moment <- x$moment
coef <- as.numeric(x$coefficient)
representation  <- x$representation
express <- "E["
for (imom in 1:length(moment))
  {term <- subblank(paste("X",imom,"^",moment[imom])) 
   express <- paste(express,term)}
express <- paste(express,"]:") 
cat(express, "\n")
print(cbind(coef,representation))
 
invisible(x)}

`make.all.moments` <- 
function (moment,verbose=TRUE) 
{
# make all moments less than or equal to input moment
# put these into the symmoments namespace

subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}

make.moment.name <- function (moment) 
{
# returns character name of moment
# moment is the moment as a vector, eg, c(1,2,3)

allchar <- c(0:9,letters,toupper(letters))
ndim <- length(moment)
moment.name <- subblank(paste("m",allchar[moment[1]+1]))
if (ndim > 1)
  {for (idim in 2:ndim)
     {Ak <- allchar[moment[idim]+1]
     moment.name <- subblank(paste(moment.name,Ak))
     }
  }

return(moment.name)
}

make.moment.vector <- function (moment) 
{
# returns character vector of moment
# moment is the moment as a vector, eg, c(1,2,3)

ndim <- length(moment)
moment.vector <- subblank(paste("c(",moment[1]))
if (ndim > 1)
  {
   moment.vector <- subblank(paste(moment.vector,","))
   for (idim in 2:ndim)
     {Ak <- moment[idim]
     moment.vector <- subblank(paste(moment.vector,Ak))
     if (idim < ndim)
       {moment.vector <- subblank(paste(moment.vector,","))}
     }
  }
moment.vector <- subblank(paste(moment.vector,")"))
return(moment.vector)
}
toSquare <- function(L.ut)
{
n <- (-1 + sqrt(1+8*length(L.ut)))/2
L <- 0*diag(n)
start <- 1
for (irow in 1:n)
  {
  L[irow,] <- c(rep(0,(irow-1)),L.ut[start:(start+n-irow)])
  start <- start+n-irow + 1
  }
return(L)}


create.envir <- TRUE
if (!exists('symmoments'))
  {symmoments <- NULL}
if (exists('symmoments'))
  {
   if (class(symmoments) == 'environment')
      {create.envir <- FALSE}
  }
if (create.envir)
  {
   print('Environment symmoments must exist to receive the moment objects.')
   return('Please create it using   symmoments <- new.env()')
   }

product = prod(moment+1)
cumproduct = cumprod(moment+1)
mdim = length(moment)
thisone <- rep(0,mdim)
for (mcount in 0:(product-1))
{  #1
  remain <- mcount
  if (mdim > 1)
    {remain <- mcount
     for (mindex in mdim:2)
       {
        denom <- cumproduct[mindex-1]
        thisone[mindex] <- trunc(remain/denom)
        remain <- remain - denom*thisone[mindex]
       }
     thisone[1] <- remain
     }
  if (mdim == 1)
     {thisone[1] <- mcount}


  if (sum(thisone)%%2 ==0)
    { #2
     notmoment <- 0
     Tmoment <- make.moment.name(thisone)
     if (exists(eval(parse(text=subblank(paste("'",Tmoment,"'")))),envir=symmoments,inherits=FALSE))
        { #3
         if (class(eval(parse(text=subblank(paste('symmoments$',Tmoment))))) == 'moment')
           {
            notmoment <- 0
            if (verbose)
               {print(paste(subblank(subblank(paste("symmoments$",Tmoment)))," exists"))}
            }  # done with this one
         if (class(eval(parse(text=subblank(paste("symmoments$",Tmoment))))) != 'moment')
           {notmoment <- 1}

        } #3
     if (!exists(eval(parse(text=subblank(paste("'",Tmoment,"'")))),envir=symmoments,inherits=FALSE) | notmoment == 1)  
        { #3
           sortmoment <- sort(thisone)
           if (sum(sortmoment==thisone) == length(thisone))  
              {  #4
               thisvec <- make.moment.vector(thisone)
               if (verbose)
                  {print(paste("Starting ",subblank(paste('symmoments$',Tmoment))))}
               eval(parse(text=subblank(paste("symmoments$",Tmoment," <- callmultmoments(",thisvec,")"))))
              }  #4  canonical, so just make it
           if (sum(sortmoment==thisone) != length(thisone))  #  unsorted
              { #4
               Smoment <- make.moment.name(sortmoment)
               if (exists(eval(parse(text=subblank(paste("'",Smoment,"'")))),envir=symmoments,inherits=FALSE))
                  { #5
                   if (class(eval(parse(text=subblank(paste("symmoments$",Smoment))))) == 'moment')  # canonical moment exists
                     { #6
                      thisvec <- make.moment.vector(thisone)
                      if (verbose)
                         {print(paste("Starting ",subblank(paste("symmoments$",Tmoment))))}
                      eval(parse(text=subblank(paste("symmoments$",Tmoment," <- tounsorted(",thisvec,",symmoments$",Smoment,")"))))
                     } #6
                      if (class(eval(parse(text=subblank(paste("symmoments$",Smoment))))) != 'moment')
                         {notmoment <- 1}
                  } #5

               if (!exists(eval(parse(text=subblank(paste("'",Tmoment,"'")))),envir=symmoments,inherits=FALSE) | notmoment == 1)
                  { #5
                   thisvec <- make.moment.vector(sortmoment)
                   if (verbose)
                      {print(paste("Starting ",subblank(paste("symmoments$",Smoment))," to create",subblank(paste("symmoments$",Tmoment))))}
                   eval(parse(text=subblank(paste("symmoments$",Smoment," <- callmultmoments(",thisvec,")"))))
                   thisvec <- make.moment.vector(thisone)
                   if (verbose)
                      {print(paste("Starting ",subblank(paste('symmoments$',Tmoment))))}
                   eval(parse(text=subblank(paste("symmoments$",Tmoment," <- tounsorted(",thisvec,",symmoments$",Smoment,")"))))
                  } #5
              } #4
        } #3
    } #2
} #1
return(NULL)
}

`toLatex_noncentral` <-
function (moment,envir='symmoments') 
{
# Compute the Latex representation of a noncentral moment 


exists.envir <- FALSE
if (exists(envir))
  {
   if (class(eval(parse(text=envir))) == 'environment')
      {exists.envir <- TRUE}
  }
if (!exists.envir)
  {return(paste('There is no environment named ', envir)) }


subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}

make.moment.name <- function (moment) 
{
# returns character name of moment
# moment is the moment as a vector, eg, c(1,2,3)

allchar <- c(0:9,letters,toupper(letters))
ndim <- length(moment)
moment.name <- subblank(paste("m",allchar[moment[1]+1]))
if (ndim > 1)
  {for (idim in 2:ndim)
     {Ak <- allchar[moment[idim]+1]
     moment.name <- subblank(paste(moment.name,Ak))
     }
  }

return(moment.name)
}
muproduct.Latex <- function (exponent) 
{
# create Latex expression for a mu-product with this exponent
ndim <- length(exponent)
term <- " "
if (exponent[1] > 0)
  {
   if (exponent[1] == 1)
     {term <- "\\mu_{1}"}     
   if (exponent[1] > 1)  
     {term <- paste("\\mu_{1}^{",exponent[1],"}")}
   }
if (ndim > 1)
  {
   for (idim in 2:ndim)
     {
      if (exponent[idim] == 1)
         {term <- paste(term,"\\mu_{",idim,"}")}
      if (exponent[idim] > 1)
         {term <- paste(term,"\\mu_{",idim,"}^{",exponent[idim],"}")}
     }
  }
if (sum(exponent) == 0)
   {term <- " "}
term <- subblank(term)

return(term)
}

doubquote <- subblank(paste("\\", "\\"))
fullmoment <- "E["
for (imoment in(1:(length(moment))) )
  {fullmoment <- paste(fullmoment, 
         "X_{", imoment, "}^{",  moment[imoment],  "}")  }
fullmoment <- paste(fullmoment, "\\mid \\mu, \\Sigma] =", doubquote)
fullmoment <- subblank(fullmoment)

product = prod(moment+1)
cumproduct = cumprod(moment+1)
mdim = length(moment)
thisone <- rep(0,mdim)

# special case: X**0
if (mdim == 1 & sum(moment) == 0)
  fullmoment <- paste(fullmoment," 1 \\\\")

if (mdim > 1 | sum(moment) > 0)
{ 

for (mcount in 0:(product-1))
{
  if (mdim > 1)
    {remain <- mcount
     for (mindex in mdim:2)
       {
        denom <- cumproduct[mindex-1]
        thisone[mindex] <- trunc(remain/denom)
        remain <- remain - denom*thisone[mindex]
       }
     thisone[1] <- remain
     }
  if (mdim == 1)
     {thisone[1] <- mcount}

  if (sum(thisone)%%2 ==0)
    {
     sortmoment <- sort(thisone)
     Tmoment <- make.moment.name(sortmoment)
     if (eval(parse(text=subblank(paste("!exists('",Tmoment,"',envir=",envir,",inherits=FALSE)")))))   
       {return(paste(Tmoment,' does not exist in environment ',envir))}

     Tmoment <- subblank(paste(envir,'$',Tmoment))

     eval(parse(text=paste("thismoment <- ",Tmoment)))       

     if (sum(thismoment$moment == thisone) < length(thismoment$moment))
        {thismoment <- tounsorted(thisone,thismoment)   
        }

      muproduct <- muproduct.Latex(moment-thismoment$moment)
 
      combinproduct <- prod(choose(moment,thismoment$moment))
      if (combinproduct == 1)
         {Acombinproduct <- " "}
      if (combinproduct > 1)
         {Acombinproduct <- subblank(paste(combinproduct))}
      thisLatex <- toLatex(thismoment)[-1] # get rid of initial expression
 
      if (length(thisLatex) > 1)

         {thisLatex[1] <- paste("(",thisLatex[1])
          thisLatex[length(thisLatex)] <- subblank(paste(thisLatex[length(thisLatex)],")"))
          thisLatex <- c(thisLatex," \\\\")
          }

      if (length(thisLatex) == 1)
         {#print(paste(thisLatex,thisLatex=="1"))
          if (subblank(thisLatex) != "1")
            {#print(paste("length",length(thisLatex),thisLatex))
             thisLatex[1] <- paste("(",thisLatex[1])
             thisLatex[length(thisLatex)] <- subblank(paste(thisLatex[length(thisLatex)],")"))
             thisLatex <- c(thisLatex," \\\\")
            }

          else {thisLatex <- c("\\\\")}
          }  
      if (mcount == 0)
         {thisLatex <- c(paste(Acombinproduct,"\\;",muproduct),thisLatex)}
      if (mcount > 0)
         {thisLatex <- c(paste("+",Acombinproduct,"\\;",muproduct),thisLatex)}
 
    fullmoment <-c(fullmoment,thisLatex)
   }
}
}
return(fullmoment)
}

`evaluate_noncentral` <- 
function (moment,mu,sigma,envir='symmoments') 
{
# Evaluate noncentral moment with mean mu and covariance matrix sigma
# moment and mu are vectors
# sigma is a vector of the upper diagonal
# envir is a character variable containing the name of the environment containing the required central moments

exists.envir <- FALSE
if (exists(envir))
  {
   if (class(eval(parse(text=envir))) == 'environment')
      {exists.envir <- TRUE}
  }
if (!exists.envir)
  {return(paste('There is no environment named ', envir)) }

subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}

make.moment.name <- function (moment) 
{
# returns character name of moment
# moment is the moment as a vector, eg, c(1,2,3)

allchar <- c(0:9,letters,toupper(letters))
ndim <- length(moment)
moment.name <- subblank(paste("m",allchar[moment[1]+1]))
if (ndim > 1)
  {for (idim in 2:ndim)
     {Ak <- allchar[moment[idim]+1]
     moment.name <- subblank(paste(moment.name,Ak))
     }
  }

return(moment.name)
}

product = prod(moment+1)
cumproduct = cumprod(moment+1)
mdim = length(moment)
thisone <- rep(0,mdim)
value <- 0

for (mcount in 0:(product-1))
{
  if (mdim > 1)
    {remain <- mcount
     for (mindex in mdim:2)
       {
        denom <- cumproduct[mindex-1]
        thisone[mindex] <- trunc(remain/denom)
        remain <- remain - denom*thisone[mindex]
       }
     thisone[1] <- remain
     }
  if (mdim == 1)
     {thisone[1] <- mcount}

if (sum(thisone)%%2 ==0)
    {
     sortmoment <- sort(thisone)
     Tmoment <- make.moment.name(sortmoment)
     if (eval(parse(text=subblank(paste("!exists('",Tmoment,"',envir=",envir,",inherits=FALSE)")))))   
       {return(paste(Tmoment,' does not exist in environment ',envir))}
     Tmoment <- subblank(paste(envir,'$',Tmoment))
     eval(parse(text=paste("thismoment <- ",Tmoment)))       
     if (mdim > 1)
       {
        if (sum(thismoment$moment == thisone) < length(thismoment$moment))
         {thismoment <- tounsorted(thisone,thismoment)} 
        }
      muproduct <- prod(mu^(moment-thismoment$moment))
      combinproduct <- prod(choose(moment,thismoment$moment))
      thisvalue <- combinproduct*evaluate.moment(thismoment,sigma)*muproduct
      value <- value + thisvalue
   }
}
return(value)
}
`convert.mpoly` <- 
function (poly) 
{
# convert between a mpoly object and a list giving the corresponding moments and coefficients

if (class(poly)=="mpoly") 
 {
   mpoly.list <- unclass(poly)
   matrix.size <- length(mpoly.list)

   poly.size <- prod(matrix.size)
   coeff <- rep(0,poly.size)
   variables <- "coef"
   vcount <- 0
   for (mcount in 1:poly.size)
     {   
      coeff[mcount] <- unlist(mpoly.list[[mcount]])['coef']
      variables <- c(variables,names(unlist(mpoly.list[mcount])))
     } 
   variables <- unique(variables)
   if (length(variables) == 1 & variables[1] == 'coef')
     {
      powers <-  matrix(0,ncol=1,nrow=1)
      return(list(coeff=coeff,powers=powers))
     }

   if (length(variables) > 1) 
      {
      variables <- variables[variables !='coef']
      ndim <- length(variables)
      powers <- matrix(rep(0,ndim*poly.size),ncol=ndim)

      for (mcount in 1:poly.size)
        {   
         thisone <- rep(0,ndim)
         thisexp <- unlist(unlist(mpoly.list[mcount]))
         thesevars <- c(names(unlist(mpoly.list[mcount])))
         nvars <- length(thesevars) - 1
         thesevars <- thesevars[thesevars != 'coef']
         nvars <- length(thesevars)
         if (nvars > 0)
            {for (idim in 1:ndim)
              {
               for (ivar in 1:nvars)
                  {
                   if (thesevars[ivar] == variables[idim])
                      {
                      thisone[idim] <- thisexp[ivar]
                      }
                  }
               }
            }
         powers[mcount,] <- thisone
        } 
     return(list(coeff=coeff,powers=powers))
     }
   }


if (class(poly)!="mpoly")  # assume conversion from matrix (ie, list) to mpoly
 { 
  n.powers <- dim(poly$powers)[1]
  if (!is.null(n.powers) | 0==0)
    {
    mpoly.list <- vector("list",length=n.powers)
    moment.size <- dim(poly$powers)[2]
  
    variables <- gsub(" ","",paste("X",1:moment.size))
     for (ipower in 1:n.powers)
      { 
       toeval <- " "
       thisterm <- poly$powers[ipower,]
       if (sum(thisterm) == 0)
         {mpoly.list[1] <- list(c(coef=poly$coeff[ipower]))} 

       if (sum(thisterm) > 0)
          {  
           thesevars <- variables[thisterm > 0]
           thesevalues <- thisterm[thisterm > 0] 
           toeval <- paste("c(",thesevars[1],"=",thesevalues[1])

           if (length(thesevars) > 1)
             { 
             for (ivar in 2:length(thesevars))
               {toeval <- paste(toeval,",",thesevars[ivar],"=",thesevalues[ivar])} 
             } 
             toeval    <- paste(toeval,", coef = ",poly$coeff[ipower],")")
             eval(parse(text=paste("mpoly.list[[ipower]] <- ",toeval)))  
          } 
      } 
      output <- mpoly.list
    return(output)
    }
  } 
  
}

toSquare <- function (L.ut) 
{
    n <- (-1 + sqrt(1 + 8 * length(L.ut)))/2
    L <- 0 * diag(n)
    start <- 1
    for (irow in 1:n) {
        L[irow, ] <- c(rep(0, (irow - 1)), L.ut[start:(start + 
            n - irow)])
        start <- start + n - irow + 1
    }
    return(L)
}


`convert.multipol` <- function (poly) 
{
# convert between a multipol object and a list giving the corresponding moments and coefficients

subblank <- function (inputstring) 
{
# remove blanks from a character string
    outputstring <- sub(" ", "", inputstring)
    if (outputstring == inputstring) {
        return(outputstring)
    }
    if (outputstring != inputstring) {
        outputstring <- subblank(outputstring)
    }
    return(outputstring)
}
if (class(poly)=="multipol") 
 {
   multipol.array <- as.array(poly)

   array.size <- dim(multipol.array)
   poly.size <- prod(array.size)
   moment.size <- length(dim(multipol.array)[dim(multipol.array)>1]) # problem of multipol array component x**0
   if (is.null(moment.size) | moment.size==0)
      {moment.size <- 1}

   coeff <- rep(0,poly.size)
   powers <- matrix(nrow=poly.size,ncol=moment.size)
   cumproduct = cumprod(dim(multipol.array))
   thisone <- rep(1,moment.size)

   nonzero <- 0
   for (mcount in 1:poly.size)
     {  
      thisone <- rep(1,moment.size)
      remain <- mcount
      if (moment.size == 1)
         {thisone <- c(remain)}
      if (moment.size > 1)
         {
         for (mindex in moment.size:2)
           {  
            thisone[mindex] <- trunc((remain-1)/cumproduct[mindex-1])
            remain <- remain - cumproduct[mindex-1]*thisone[mindex]
            thisone[mindex] <- thisone[mindex] + 1
           }  
         }
      thisone[1] <- remain

      powers[mcount,] <- thisone - rep(1,moment.size)  ###  account for indexing in multipol
      thistemp = paste(thisone[1])
      if (moment.size > 1)
        {
        for (mindex in 2:moment.size)
          {thistemp = paste(thistemp,",",thisone[mindex])}
        }
        thistemp = subblank(thistemp)
        coefftemp <- NULL
        eval(parse(text=paste("coefftemp <- multipol.array[",thistemp,"]")))
        if (coefftemp != 0)
           {
            nonzero <- nonzero + 1
            eval(parse(text=paste("coeff[nonzero] <- multipol.array[",thistemp,"]")))
            powers[nonzero,] <- thisone - rep(1,moment.size)  ###  account for indexing in multipol
           }
      } 
  coeff <- coeff[1:nonzero]
  powers <- as.matrix(powers[1:nonzero,],ncol=dim(powers)[2])
  output <- list(coeff=coeff,powers=powers)
  return(output)
 }


if (class(poly)!="multipol")  # assume conversion from matrix (ie, list) to multipol
 {
  n.powers <- dim(poly$powers)[1]
  moment.size <- dim(poly$powers)[2]
  if (dim(poly$powers)[1] == 1)
     {dimmult <- c(1,poly$powers + rep(1,moment.size))} 
  if (dim(poly$powers)[1]> 1)
     {dimmult <- apply(poly$powers, 2, max) + rep(1,moment.size)}
  n.elements <- prod(dimmult)
  multipol.array <- array(rep(0,n.elements),dim=dimmult)
  output <- as.multipol(multipol.array)
  for (ipower in 1:n.powers)
    {
      thistemp = paste(1+poly$powers[ipower,1])
      if (moment.size > 1)
        {
        for (mindex in 2:moment.size)
          {thistemp = paste(thistemp,",",(1+poly$powers[ipower,mindex]))}
        }
        thistemp = subblank(thistemp)
        eval(parse(text=paste("multipol.array[",thistemp,"] <- poly$coeff[",ipower,"]")))
   }
    output <- as.multipol(multipol.array) 

  return(output)
 }


}

`evaluate_expected.polynomial` <- 
function (poly,mu,sigma,envir='symmoments') 
{
# compute expected value of a multidimensional polynomial
# poly is either a multipol objects or
# a list with components "powers", which lists the moments involved,
# and "coeff, their coefficients
# mu is the mean as a vector
# sigma is the variance-covariance matrix as the vector of the upper diagonal part, including diagonal
# envir is a character variable containing the name of the environment containing the required central moments

exists.envir <- FALSE
if (exists(envir))
  {
   if (class(eval(parse(text=envir))) == 'environment')
      {exists.envir <- TRUE}
  }
if (!exists.envir)
  {return(paste('There is no environment named ', envir)) }


temp.poly <- poly
if (class(poly) == "multipol")
  {temp.poly <- convert.multipol(poly)}
if (class(poly) == "mpoly")
  {temp.poly <- convert.mpoly(poly)}


ndim <- dim(temp.poly$powers)[2]
npowers <- dim(temp.poly$powers)[1]
powers <- temp.poly$powers
coeff <- temp.poly$coeff

value <- 0
for (imom in 1:npowers)  
   {
    if (coeff[imom] != 0)
       {
        thisvalue <- coeff[imom]*evaluate_noncentral(powers[imom,],mu,sigma,envir)
        value <- value + thisvalue
       }
   }

return(value)
}



`toMoment` <- 
function (inputobject, tip.label = NULL) 
{
    if (class(inputobject) != "matching" & !is.matrix(inputobject)) {
        if (class(inputobject) != "L-Newick") {
            temp <- inputobject
            Newick.out <- inputobject
        }
        if (class(inputobject) == "L-Newick") {
            temp <- inputobject$Newick
            Newick.out <- inputobject$Newick
        }
        temp <- gsub(";", " ", temp, fixed = TRUE)
        temp <- gsub("(", " ", temp, fixed = TRUE)
        temp <- gsub(")", " ", temp, fixed = TRUE)
        temp <- gsub(",", " ", temp, fixed = TRUE)
        temp <- gsub("  ", " ", temp, fixed = TRUE)
        tips <- strsplit(temp, " +")[[1]]
        tips <- tips[tips != ""]
        n <- length(tips)
        neworder <- 1:n
        if (!is.null(tip.label)) {
            for (itips in 1:n) {
                neworder[itips] <- grep(tip.label[itips], tips)
            }
        }
        tips <- tips[neworder]
        np <- n
        temp.tips <- rep(" ", n)
        temp.Newick <- gsub(";", " ", inputobject, fixed = TRUE)
        for (itips in 1:n) {
            temp.tips[itips] <- gsub(" ", "", paste(";", as.character(itips), 
                ";"), fixed = TRUE)
            temp.Newick <- gsub(tips[itips], temp.tips[itips], 
                temp.Newick, fixed = TRUE)
        }
        n.temp.tips <- n
        temp.tip.names <- temp.tips
        for (itips in 1:n) {
            temp.tips[itips] <- gsub(" ", "", temp.tips[itips], 
                fixed = TRUE)
        }
        couples.n <- 1
        while (couples.n > 0) {
            couples.n <- 0
            for (cspec in 1:np) {
                for (rspec in 1:np) {
                  couple <- gsub(" ", "", paste("(", temp.tips[cspec], 
                    ",", temp.tips[rspec], ")"), fixed = TRUE)
                  if (length(grep(couple, temp.Newick)) > 0) {
                    couples.n <- couples.n + 1
                    n.temp.tips <- n.temp.tips + 1
                    temp.tip.names <- c(temp.tip.names, couple)
                    temp.tips <- c(temp.tips, gsub(" ", "", paste(";", 
                      as.character(n.temp.tips), ";"), fixed = TRUE))
                    temp.Newick <- gsub(couple, gsub(" ", "", 
                      paste(";", as.character(n.temp.tips), ";"), 
                      fixed = TRUE), temp.Newick, fixed = TRUE)
                  }
                }
            }
            np <- np + couples.n
        }
        temp.tip.names.n <- length(temp.tip.names)
        L <- 0 * diag(2 * (n - 1))
        for (icouple in (n + 1):temp.tip.names.n) {
            temp <- gsub("(", " ", temp.tip.names[icouple], fixed = TRUE)
            temp <- gsub(")", " ", temp, fixed = TRUE)
            temp <- gsub(",", " ", temp, fixed = TRUE)
            temp <- gsub(";", " ", temp, fixed = TRUE)
            temp <- gsub("  ", " ", temp, fixed = TRUE)
            ind <- strsplit(temp, " +")[[1]]
            ind <- ind[ind != ""]
            indn <- as.integer(ind)
            indn <- ind[order(indn)]
            irow <- as.integer(indn[1])
            icol <- as.integer(indn[2])
            L[irow, icol] <- 1
        }
        L.ut <- NULL
        for (irow in 1:(2 * (n - 1))) {
            L.ut <- c(L.ut, L[irow, irow:(2 * (n - 1))])
        }
        SMNM <- list(L = L, L.ut = L.ut, Newick = Newick.out, 
            tip.label = tips, tip.label.n = n)
        class(SMNM) <- "L-matrix"
        return(SMNM)
    }
    if (class(inputobject) == "matching" | is.matrix(inputobject)) {
        temp.tip.label <- NULL
        if (class(inputobject) == "matching") {
            temp.matching <- inputobject$matching[, c(1, 2)]
            if (!is.null(inputobject$tip.label)) {
                temp.tip.label <- inputobject$tip.label
            }
        }
        if (is.matrix(inputobject)) {
            temp.matching <- inputobject[, c(1, 2)]
        }
        nL <- max(temp.matching)
        n <- 1 + nL/2
        L <- 0 * diag(nL)
        if (is.matrix(temp.matching)) {
            L[temp.matching] <- 1
        }
        if (!is.matrix(temp.matching)) {
            L[matrix(temp.matching, nrow = 1)] <- 1
        }
        L.ut <- NULL
        for (irow in 1:(2 * (n - 1))) {
            L.ut <- c(L.ut, L[irow, irow:(2 * (n - 1))])
        }
        if (!is.null(temp.tip.label)) {
            Newick <- toNewick(L, type = "square", tip.label = temp.tip.label)
        }
        if (is.null(temp.tip.label)) {
            Newick <- toNewick(L, type = "square")
        }
        if (is.null(tip.label)) {
            alphab <- gsub(" ", "", paste(rep(letters, times = rep(676, 
                26)), paste(rep(letters, times = rep(26, 26)), 
                letters)))
            if (n <= 26) {
                tips <- letters[1:n]
            }
            if (n > 26) {
                tips <- alphab[1:n]
            }
        }
        if (!is.null(tip.label)) {
            tips <- tip.label
        }
        if (is.null(tip.label) & class(inputobject) == "matching") {
            if (!is.null(inputobject$tip.label)) {
                tips <- inputobject$tip.label
            }
        }
        SMNM <- list(L = L, L.ut = L.ut, Newick = Newick$Newick, 
            tip.label = tips, tip.label.n = n)
        class(SMNM) <- "L-matrix"
        return(SMNM)
    }
}
`toNewick` <- 
function (L, type = NULL, tip.label = NULL) 
{
toSquare <- function (L.ut) 
{
    n <- (-1 + sqrt(1 + 8 * length(L.ut)))/2
    L <- 0 * diag(n)
    start <- 1
    for (irow in 1:n) {
        L[irow, ] <- c(rep(0, (irow - 1)), L.ut[start:(start + 
            n - irow)])
        start <- start + n - irow + 1
    }
    return(L)
}

    if (class(L) == "L-matrix") {
        L.matrix <- L$L
    }
    if (!is.null(type)) {
        if (type == "square") {
            L.matrix <- L
        }
        if (type == "ut") {
            L.matrix <- toSquare(L)
        }
    }
    nL <- dim(L.matrix)[1]
    n <- 1 + nL/2
    alphab <- gsub(" ", "", paste(rep(letters, times = rep(676, 
        26)), paste(rep(letters, times = rep(26, 26)), letters)))
    if (!is.null(tip.label)) {
        tips <- tip.label
    }
    if (is.null(tip.label) & class(L) == "L-matrix") {
        tips <- L$tip.label
    }
    if (is.null(tip.label) & class(L) != "L-matrix") {
        if (n <= 26) {
            tips <- letters[1:n]
        }
        if (n > 26) {
            tips <- alphab[1:n]
        }
    }
    nodes <- tips
    for (icol in 2:nL) {
        irow <- (1:icol)[L.matrix[1:icol, icol] == 1]
        if (length(irow) == 1) {
            ccol <- nodes[icol]
            crow <- nodes[irow]
            if (icol <= n) {
                couple <- gsub(" ", "", paste("(", crow, ",", 
                  ccol, ")"), fixed = TRUE)
            }
            if (icol > n) {
                couple <- gsub(" ", "", paste("(", ccol, ",", 
                  crow, ")"), fixed = TRUE)
            }
            nodes <- c(nodes, couple)
        }
    }
    nodes.n <- length(nodes)
    temp.Newick <- gsub(" ", "", paste(nodes[nodes.n], ";"))
    nodes <- gsub("!", "", nodes, fixed = TRUE)
    Newick <- list(Newick = temp.Newick, tip.label = tips, tip.label.n = n, 
        L = L.matrix, node.label = nodes, nodes.label.n = nodes.n)
    class(Newick) <- "L-Newick"
    return(Newick)
}
`toMatching` <- 
function (L, type = NULL, tip.label = NULL) 
{
toSquare <- function (L.ut) 
{
    n <- (-1 + sqrt(1 + 8 * length(L.ut)))/2
    L <- 0 * diag(n)
    start <- 1
    for (irow in 1:n) {
        L[irow, ] <- c(rep(0, (irow - 1)), L.ut[start:(start + 
            n - irow)])
        start <- start + n - irow + 1
    }
    return(L)
}

    if (class(L) == "L-matrix") {
        L.matrix <- L$L
        if (!is.null(L$tip.label) & is.null(tip.label)) {
            tip.label <- L$tip.label
        }
    }
    if (!is.null(type)) {
        if (type == "square") {
            L.matrix <- L
        }
        if (type == "ut") {
            L.matrix <- toSquare(L)
        }
    }
    nL <- dim(L.matrix)[1]
    n <- 1 + nL/2
    temp.matching <- matrix(rep(0, 3 * (n - 1)), ncol = 3)
    if (!is.null(tip.label)) {
        temp.tip.label <- tip.label
    }
    if (is.null(tip.label)) {
        alphab <- gsub(" ", "", paste(rep(letters, times = rep(676, 
            26)), paste(rep(letters, times = rep(26, 26)), letters)))
        if (!is.null(tip.label)) {
            temp.tip.label <- tip.label
        }
        if (is.null(tip.label)) {
            if (n <= 26) {
                temp.tip.label <- letters[1:n]
            }
            if (n > 26) {
                temp.tip.label <- alphab[1:n]
            }
        }
    }
    temp.node.label <- temp.tip.label
    nrows <- 0
    for (icol in 2:nL) {
        irow <- (1:icol)[L.matrix[1:icol, icol] == 1]
        if (length(irow) == 1) {
            ccol <- temp.node.label[icol]
            crow <- temp.node.label[irow]
            if (icol <= n) {
                couple <- gsub(" ", "", paste("(", crow, ",", 
                  ccol, ")"), fixed = TRUE)
            }
            if (icol > n) {
                couple <- gsub(" ", "", paste("(", ccol, ",", 
                  crow, ")"), fixed = TRUE)
            }
            nrows <- nrows + 1
            temp.node.label <- c(temp.node.label, couple)
            temp.matching[nrows, ] <- c(irow, icol, length(temp.node.label))
        }
    }
    temp.node.label <- gsub("!", "", temp.node.label, fixed = TRUE)
    match.obj <- list(matching = temp.matching, tip.label = temp.tip.label, 
        node.label = temp.node.label)
    class(match.obj) <- "matching"
    return(match.obj)
}

`tounsorted` <- 
function (moment,sorted.moment) 
{
# converts a sorted moment to a specified unsorted moment
# eg,   m(2,3,5) ->  m(5,2,3)
# each row is sorted separately
# the rows may be ordered differently

# moment: noncanonical moment to obtain 
#         moment is in vector form, eg, c(3,1,2)

# sorted.moment: canonical moment of class "moment"
#         this moment is monotone in its powers;


toSquare <- function (L.ut) 
{
    n <- (-1 + sqrt(1 + 8 * length(L.ut)))/2
    L <- 0 * diag(n)
    start <- 1
    for (irow in 1:n) {
        L[irow, ] <- c(rep(0, (irow - 1)), L.ut[start:(start + 
            n - irow)])
        start <- start + n - irow + 1
    }
    return(L)
}



nestedreps <- function (input.vector, inner.rep, outer.rep) 
{
# replicates input.vector, first by inner.rep, then by outer.rep

# input.vector: vector to replicate
# inner.rep:    count of replicates of elements
# outer.rep:    count of replicates of resulting vector

temp.vector <- NULL
for (ivec in 1:length(input.vector))
     {temp.vector <- c(temp.vector, rep(input.vector[ivec], inner.rep))}

total.vector <- rep(temp.vector, outer.rep) 

return(total.vector)}


output.moment <- sorted.moment
output.moment$moment <- moment
lmom <- length(moment)
r <- order(moment)
lrep <- lmom*(lmom+1)/2
representation <- sorted.moment$representation

nrep <- dim(representation)[1]
overallcoeff <- ((1/2)^(sum(moment)/2))*prod(fact(moment))/fact(sum(moment)/2)

m <- matrix(rep(1, lmom^2), nrow=lmom) 

limits <- c((lmom:1)%*%(m * !(lower.tri(m))))
# sum of row lengths: lmom, lmom + lmom - 1, ... , for use in row_col
row_col <- matrix(rep(0, (2 * lrep)), nrow=2)
for (icell in (1:lrep))
   {
   row_col[1, icell] <- min((1:lmom)[icell<=limits])
   if (row_col[1, icell] == 1){row_col[2, icell] <- icell}
   if (row_col[1, icell]>1){row_col[2, icell] <- icell - limits[row_col[1, icell] - 1] + 
                                             row_col[1, icell] - 1 }
  }
#  2x(nm * (nm + 1) / 2 matrix giving rows and columns for each cell
#    compute here so that they don't have to be calculated each time 


for (irep in 1:nrep)
  {

   utri <- toSquare(representation[irep,])
   noncanonical.matrix <- 0*utri

   for (irow in 1:lmom)
      {
       jrow <- r[irow]
       for (icol in irow:lmom)
          {
           jcol <- r[icol]
          if (jrow <= jcol)
             {noncanonical.matrix[jrow,jcol] <- utri[irow,icol]}
          if (jrow > jcol)
               {noncanonical.matrix[jcol,jrow] <- utri[irow,icol]}
          }
      }
   thisrep <- t(noncanonical.matrix)[t(!lower.tri(noncanonical.matrix))]


  output.moment$representation[irep,] <- thisrep


#  determine the coefficient for each term based on switching equivalent terms
#          this is taken from callmultmoments

#  "base" gives the number of switches that can be made to each element of the l-matrix
#  diagonal elements are not switchable, but are included to allow subtraction below

  base <- c(rep(1, lmom * (lmom + 1) / 2))
  base[1] <- 1   # first diagonal element  -  not switchable

  totreps <- 1   # total number of transpostions
# if there is only one element, it must be the diagonal, so is not switchable  -  skip
if (lmom > 1){
  base[1] <- 1   # first diagonal element
  for (cell in 2:length(base))
    {
     icol = row_col[1, cell]  #  determine if diagonal element
     irow = row_col[2, cell]
     if (irow == icol){base[cell] <- 1}  # diagonal  -  not switchable
     if (icol != irow)
         {totreps <- totreps * (1 + thisrep[cell])
         base[cell] <- 1 + thisrep[cell]} 
    }  #  done with computing base and total transpositions (totreps)
  }
 
mcoeff <- 1  #  sum of multinomial coefficients
if (totreps > 1){ 

#  baserep will represent the lower diagonal (including diagonal)
#  of the augmented matrices
     baserep <- matrix(rep(0, totreps * length(base)), nrow=totreps)
     basegt1 <- base[base>1]
     nbase <- 0
     for (ibase in 1:length(base))
       {if (base[ibase] > 1)
           {nbase = nbase + 1
            if (nbase == 1){baserep[, ibase] <- nestedreps(c(0:(basegt1[nbase] - 1)), 1, totreps / prod(basegt1[1:nbase])) }    
            if (nbase > 1) {baserep[, ibase] <- nestedreps(c(0:(basegt1[nbase] - 1)), prod(basegt1[1:(nbase - 1)]), totreps / prod(basegt1[1:nbase])) }
           }
       }
 
#  now go through each transposition
     if ( !is.na(totreps) & totreps != 1)
    {
     mcoeff <- 0
     for (jrep in (1:totreps)) # check each transposition
          {newrep <- baserep[jrep, ]       #  added lower diagonal elements
           addrep <- sort(newrep, decreasing=TRUE)[1:(lmom * (lmom - 1) / 2)] 
           fulnrep <- c((thisrep - newrep), addrep)
           thiscoeff <- ((length(fulnrep))^sum(fulnrep)) * dmultinom(x=fulnrep, prob=rep(1.0, length(fulnrep)))
           mcoeff <- mcoeff + thiscoeff
#  the multinomial coefficient is obtained from the multinomial distribution
#  multiply by an appropriate power to get rid of probability
        }  
    }
}  
  if (is.na(totreps)){mcoeff <- 1}
if (totreps == 1)
   {mcoeff <- (length(thisrep))^sum(thisrep) * dmultinom(x=thisrep, prob=rep(1.0, length(thisrep)))}

#  determine full coefficient  -  round because all coefficients should be integers
#                                (Note - this statement has not been proved)

     output.moment$coefficients[irep] <- round(overallcoeff * mcoeff)

  }  # end of representations
   

return(output.moment)}  

`integrate.polynomial` <- 
function (poly,mu,sigma,lower=NULL,upper=NULL) 
{
# integrate polynomial moment against MVN

# poly: either a multipol objects or 
#       a multipol defined by a list with moment powers and coefficients
# mu: mean of multivariate normal as vector
# sigma: variance-covariance matrix of multivariate normal
# lower, upper: vectors giving limits of integration
#    if one is NULL, then make it the mean +/- 6 * SD

if (is.null(lower)) 
  {lower <- mu - 6*sqrt(diag(sigma))}
if (is.null(upper))
   {upper <- mu + 6*sqrt(diag(sigma))}

thispoly <- poly
if (class(poly) == "multipol")
  {thispoly <- convert.multipol(poly)}
if (class(poly) == "mpoly")
  {thispoly <- convert.mpoly(poly)}

ndim <- dim(thispoly$powers)[2]
npowers <- dim(thispoly$powers)[1]
powers <- thispoly$powers
coeff <- thispoly$coeff
value <- 0
   
    f <- function(x)
     {
     y <- x[1]^powers[imom,1]
     for (idim in 2:ndim)
       {
        y <- y*x[idim]^powers[imom,idim]
       }
     y <- y*dmvnorm(x,mean=mu,sigma=sigma, log=FALSE)
     return(y)
    }

for (imom in 1:npowers)   
   {thisvalue <- adaptIntegrate(f,lower,upper)$integral
    value <- value + coeff[imom]*thisvalue
   }
return(value)}


