semimetric.BAGIDIS =function(DATA1,
                             DATA2=DATA1,
                             p = 2, 
                             wk=NULL, 
                             Param=0.5,
                             wdw= min(ncol(DATA1),30),  
                             Evol =FALSE,
                             Overlap = wdw-1,
                             method = c('TS','BD')
                             ){ 
#===============================================================================
# FUNCTION : semimetric.BAGIDIS
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial
# Sciences, UCLouvain. 
# Last update: 05/05/2012 
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# This function computes the Bagidis semidistance between curves.
# If several curves are provided, it returns a matrix of semidistances.
# The function is an alias for either semimetric.BAGIDIS.TS or semimetric.BAGIDIS.BD, depending on the value of the parameter 'method'.
# 
#-------------------------------------------------------------------------------
# USAGE: 
#semimetric.BAGIDIS(  DATA1,
#                     DATA2=DATA1,
#                     p = 2, 
#                     wk=NULL, 
#                     Param=0.5,
#                     wdw= min(ncol(DATA1),30),  
#                     Evol =FALSE,
#                     Overlap = wdw-1,
#                     method = 'TS'
#                     )
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# DATA1 and DATA2:  matrices containing the series to be compared row by row.  Each row of DATA1 has its semidistance being computed with every row of DATA2. If only DATA1 is provided, then DATA2 is taken to be equal to DATA1. We must have ncol(DATA1)= ncol(DATA2).    
# p: the kind of norm to be used for computing the partial diastance in the B-D plane.  Must be numeric or Inf.
# wk: a vector of weights of length ncol(DATA1)-1. If not provided, wk = log(N+1-(1:N))/log(N+1) with N= ncol(DATA1)-1 , as this is the default of function BAGIDIS.dist that is required for the computation.
# Param : the balance parameter between the differences along the breakpoint axis and along the detail axis. Param must be in [0;1]. Param= 1 means that only breakpoints differences are taken into account. Param=0 means that only details differences are taken into account.
# wdw: In case distances are measured between "long" series, it could be advantageous to make use of a windowed semimetric. wdw encode the length of the window in which the semimetric will be computed between the subseries. By default there is no windowing if the length of the series ( =ncol(DATA1) ) is smaller than or equal to 30, and a windows length of 30 otherwise.
# Evol : Logical. In case a windowing is applied, should the matrices of local (windowed) dissimilarities be returned. Default is FALSE.
# Overlap : In case a windowing of the series is applied, Overlap determines how the subseries overlap each other. By default, a one-step-sliding distance is computed.
# method: either 'TS' or 'BD' : the method for computing the matrix of semi-distances in case of multiple series. Results are identical. 'TS' recompute the BUUHWE transform for each pairwise comparison, 'BD' computes all signatures beforehand and store them before computing the distances. 'TS' requires more time, 'BD' requires more storage.  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# VALUE:
# dissimilarity.matrix :  Matrix of semidistances between the nrow(DATA1) series of DATA1 and the nrow(DATA2) series of DATA2. Dimensions: nrow(DATA1) x nrow(DATA2) . 
# dissimilarity.evol :  Array of local matrices of semidistances between the windowed series of DATA1 and DATA2. Dimensions: nrow(DATA1) x nrow(DATA2) x Nb_subseries. Nb_subseries is determined by the three quantities nrow(DATA1), wdw and Overlap.   
#-------------------------------------------------------------------------------
# WARNING/NOTE :
# - Using semimetric.BAGIDIS requires functions BUUHWE and BAGIDIS.dist, semimetric.BAGIDIS.BD, semimetric.BAGIDIS.TS, Breakpoints, Details, DefineSubseries  .
# - Computation time is affected by the number of rows in DATA1 and DATA2. If nrow(DATA1)> nrow(DATA2), it increases the number of operations to be computed. If nrow(DATA2)>nrow(DATA1), it increases the memory usage.
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS,a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences. 
#-------------------------------------------------------------------------------
# SEE ALSO:
# BUUHWE , BAGIDIS.dist ,  semimetric.BAGIDIS.BD, semimetric.BAGIDIS.TS ...
#-------------------------------------------------------------------------------
# EXAMPLES:
# x= 1:10
# y=2:11
# A=rbind(x,y)
# semimetric.BAGIDIS(A)
  #          [,1]     [,2]
  # [1,] 0.000000 3.578209
  # [2,] 3.578209 0.000000
# B= rbind(x,x,y)
# semimetric.BAGIDIS(A,B)
  #          [,1]     [,2]     [,3]
  # [1,] 0.000000 0.000000 3.578209
  # [2,] 3.578209 3.578209 0.000000
# x= 1:30
# y= 1:30
# A= rbind(x,y)
# B= rbind(x,x, y)
# semimetric.BAGIDIS(A,B, wdw =15, Evol =TRUE, Overlap =0)
  # $dissimilarity.matrix 
  #      [,1] [,2] [,3]
  # [1,]    0    0    0
  # [2,]    0    0    0
  # 
  # $dissimilarity.evol
  # , , 1
  #
  #      [,1] [,2] [,3]
  # [1,]    0    0    0
  # [2,]    0    0    0
  #
  # , , 2
  #
  #      [,1] [,2] [,3]
  # [1,]    0    0    0
  # [2,]    0    0    0
#===============================================================================
#===============================================================================





#-------------------------------------------------------------------------------
# 1. Tranforming vectors into matrices if needed
#-------------------------------------------------------------------------------

if (is.vector(DATA1)){
DATA1 = matrix(DATA1, nrow=1,byrow=TRUE)}
if  (is.vector(DATA2) ){ 
DATA2 =matrix(DATA2,nrow=1,byrow=TRUE) }


#-------------------------------------------------------------------------------
# 2. Checking for incorrect values of the arguments
#  - if all series have the same length :   ncol(DATA1) = ncol(DATA2)
#  - if wdw is in [2;ncol(DATA1)]
#  - if Overlap is in [0; wdw-1]
#  - if Param is included in [0;1]
#  - if length(p) =1 
#  - if wk=NULL or length(wk) =wdw -1 
#  - if Evol is boolean
#-------------------------------------------------------------------------------

if (ncol(DATA1) != ncol(DATA2)){
stop('ERROR: Lengths of the series in DATA1 and DATA2 should be equals')}

if (wdw > ncol(DATA1)){
wdw= min(ncol(DATA1),30)
warning('WARNING: The length of the window (wdw) cannot be higher than the length of the series. wdw is set to its default value.')
}

if (wdw < 2){
stop('ERROR: The length of the window (wdw) cannot be smaller than 2.')
}

if (Overlap > wdw-1 ){
stop('ERROR: The Overlap parameter cannot be higher than the wdw-1.')
}

if (Overlap < 0){
stop('ERROR: the Overlap parameter must be non negative.')
}

if (Param < 0 | Param > 1){
stop('ERROR: the balance parameter Param must be included in [0;1].')}


if (length(p) !=1){
stop('ERROR: p must be a single numeric value, not a vector.')}

if (!is.null(wk) & length(wk)!= (wdw-1) ){
stop('ERROR: wk must be a numeric vector of length wdw-1, or must be set to NULL (so that its default value is applied).')}

if (!is.logical(Evol)){
warning('WARNING: Evol must be boolean. It is set to its default value (FALSE).')} 


 if(! method %in% c('BD','TS')) {
 method ='TS'
 warning("WARNING: method must be 'BD' or 'TS'. It is set to its default value (TS).")
 }



if (method=='BD'){

if(identical(DATA1,DATA2)){idem =TRUE}  else {idem=FALSE}
 
if (!is.matrix(DATA1)) {DATA1=matrix(DATA1,nrow=1,ncol =length(DATA1))}
if (!idem) {DATA2=as.matrix(DATA2,nrow=1)}



Subseries1=NULL
for (i in 1:nrow(DATA1)){
Subseries1 = rbind(Subseries1,DefineSubseries(DATA1[i,], wdw, Overlap)) }


if (!idem) {
Subseries2 =NULL
for (i in 1:nrow(DATA2)){
Subseries2 = rbind(Subseries2,DefineSubseries(DATA2[i,], wdw, Overlap)) }} 

Dataset1.BUUHWE = apply(Subseries1, 1, BUUHWE)
if (!idem) {Dataset2.BUUHWE = apply(Subseries2, 1, BUUHWE)} 

Breakpoints1 = Breakpoints(Dataset1.BUUHWE)  # col by col
Details1 = Details(Dataset1.BUUHWE)
if (!idem) {Breakpoints2 = Breakpoints(Dataset2.BUUHWE) 
Details2 = Details(Dataset2.BUUHWE)} 

if(idem){
Breakpoints2 = Breakpoints1
Details2 = Details1
}


NbSubseries = floor(nrow(Subseries1)/nrow(DATA1))

out = semimetric.BAGIDIS.BD(Details1, 
                            Breakpoints1,
                            Details2, 
                            Breakpoints2,
                            NbSubseries,
                             p, 
                             wk, 
                             Param,                             
                             Evol
                             )


} else if (method =='TS'){ # EO If BD

out = semimetric.BAGIDIS.TS(DATA1, 
                            DATA2,
                             p, 
                             wk, 
                             Param,
                             wdw,  
                             Evol,
                             Overlap
                             )

} # EO If TS


return(out)

                             
}#EOF