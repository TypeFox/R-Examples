semimetric.BAGIDIS.TS= function(DATA1,
                             DATA2=DATA1,
                             p = 2, 
                             wk=NULL, 
                             Param=0.5,
                             wdw= min(ncol(DATA1),30),  
                             Evol =FALSE,
                             Overlap = wdw-1 
                             ){  
#===============================================================================
# FUNCTION : semimetric.BAGIDIS
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 16/08/2010 17:31:15
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# This function computes the Bagidis semidistance between series of measurements.
# If several series are provided, the function returns a matrix of semidistances.
# A user-friendly interface for this function is given by semimetric.BAGIDIS with the option method='BD'.
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
#                     Overlap = wdw-1
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
#-------------------------------------------------------------------------------
# DETAILS:
# The \textsc{Bagidis} semimetric has been introduced by \citet{BagidisPaper1} so as to measure differences between curves that are characterized by some sharp features. It is a functional data-driven and wavelet-based method that is highly adaptive to the curves being considered. Its main originality is its ability to consider simultaneously horizontal and vertical variations of patterns occurring in the curves. The idea is as follows: As a first step, we expand each series in the Unbalanced Haar wavelet basis that is best suited to it, as obtained by the BUUHWT algorithme described in \citet{Piotr}.  Those unbalanced Haar wavelet bases are orthonormal bases that are made up of one constant vector, and a set of Haar-like (i.e. up-and-down shaped) orthonormal wavelets whose discontinuity point (hereafter the breakpoint) between the positive and negative parts is not necessarily located at the middle of its support. The BUUHWT algorithm allows for organising the basis vectors according to a principle of hierarchy of the patterns that forms the curve. 

# Let's denote such an expansion of a series $x^{(i)}$ as follows:
#$$ \mathbf{x}^{(i)}= \sum_{k=0}^{N-1} d_k^{(i)} \boldsymbol{\psi}_k^{(i)}, $$
#where the coefficients $d_k^{(i)}$ are the projections of  $ \mathbf{x}^{(i)}$ on the corresponding basis vectors $\boldsymbol{\psi}_k^{(i)}$ (i.e. the \textit{detail} coefficients) and  where the set of vectors $\{ \boldsymbol{\psi}_k^{(i)}\}_{k=0 \ldots N-1}$ is the Unbalanced Haar wavelet basis that is best suited to the series $\mathbf{x}^{(i)}$, as obtained using the BUUHWT algorithm.
#Let us also denote $b_k^{(i)}$, the breakpoint $b$ of the wavelet $\boldsymbol{\psi}_k^{(i)}$, $k= 1 \ldots N-1$. 
#\citet{Piotr} has shown that the ordered set of breakpoints $\{ b_k^{(i)}\}_{k=1 \ldots N-1}$ determines the basis $\{ \boldsymbol{\psi}_k^{(i)}\}_{k=0 \ldots N-1}$ uniquely.	As a consequence, the set of points $ \{ y_k^{(i)}\}_{k=1}^{N-1}$ = ($b_k^{(i)}, d_k^{(i)}$)$_{k=1 \ldots N-1}$ determines the shape of the series $\mathbf{x}^{(i)}$ uniquely (i.e., it determines the series, except for a change of the mean level of the series, that is encoded by the additional coefficient $d_0^{(i)}$). 

#Given that, the \textsc{Bagidis} semimetric between $\mathbf{x}^{(1)}$ and $\mathbf{x}^{(2)}$ is defined as follows: 
#\be \label{expr_semidist}
#d_{p\lambda}(\mathbf{x}^{(1)},\mathbf{x}^{(2)}) = \sum_{k=1}^{N-1} w_k \left\| \mathbf{y}_k^{(1)} -  \mathbf{y}_k^{(2)} \right\| _{p\lambda}, \quad \mathrm{with} \quad p= 1,2,\dots, \infty
#\ee
#with 
#\be \label{expr_semidistBD}
#\left\| \mathbf{y}_k^{(1)} -  \mathbf{y}_k^{(2)} \right\| _{p\lambda} = \left(\lambda  \left| b_k^{(1)} - b_k^{(2)} \right|^p + (1-\lambda) \left|  d_k^{(1)} - d_k^{(2)}\right|^p \right)^{1/p}.
#\ee
#and 
#$$w_k = \frac{ \log(N+1-k)}{\log(N+1)}, \quad k=1 \ldots N-1. $$

#In \ref{expr_semidistBD}, the parameter $\lambda$ actually corresponds to a scaling parameter in the planes $\{ b_k^{(i)}, d_k^{(i)}\}_{k=1 \ldots N-1}$, and hence also in the original units of the problem. This dissimilarity has been proven to be a semimetric in \citet{BagidisPaper1}.
#
##\paragraph{Avoiding Feature Confusion by Using a Sliding Distance.}
##Increasing the length of the series we want to compare often implies also increasing the number of features that appear in the series. If we aim at comparing the fine structure of the series, this could be problematic as feature confusion could occur in the BUUHWT expansions of the different series. 
#In that case, it appears thus useful to make use of a localized version of our semi-distance, when dealing with long series.
#The principle is as follows. We consider a sliding window of length $\Delta$.  For each localization of that window, we obtain the BUUHWT expansion of the so-defined subseries and we compute the semi-distance of expression (\ref{expr_semidist}) or expression (\ref{def_semidist_lambda}). 
#One can then take the mean of those local semi-distances so as to obtain a global measure of dissimilarity.
#
#\paragraph{Localizing Differences between Long Series.}
#When dealing with long series, one can also be interested to point the abscissas where the series do significantly differ.  Using a sliding semi-distance  allows us to consider the curve of semi-distances between two series, so that we can automatically identify the abscissas where two series are close to each other (as seen using the distance of equation (\ref{expr_semidist}) or equation (\ref{def_semidist_lambda})) and the ones where there are more distant. This could be a helpful diagnostic tool for comparing observational curves with a target curve for instance. 
#
#For sure, the choice of the length $\Delta$ of the windows should be problem-dependent: $\Delta$ defines the range in which a given feature could possibly move along the abscissa axis while remaining identified as a unique feature from one series to another. Experience shows little sensitivity to small variations of that parameter, so that only an order of magnitude of $\Delta$ should actually be provided. 
#-------------------------------------------------------------------------------
# VALUE:
# dissimilarity.matrix :  Matrix of semidistances between the nrow(DATA1) series of DATA1 and the nrow(DATA2) series of DATA2. Dimensions: nrow(DATA1) x nrow(DATA2) . 
# dissimilarity.evol :  Array of local matrices of semidistances between the windowed series of DATA1 and DATA2. Dimensions: nrow(DATA1) x nrow(DATA2) x Nb_subseries. Nb_subseries is determined by the three quantities nrow(DATA1), wdw and Overlap.   
#-------------------------------------------------------------------------------
# WARNING/NOTE :
# - Using semimetric.BAGIDIS requires functions BUUHWE and BAGIDIS.dist .
# - Computation time is affected by the number of rows in DATA1 and DATA2. If nrow(DATA1)> nrow(DATA2), it increases the number of operations to be computed. If nrow(DATA2)>nrow(DATA1), it increases the memory usage.
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS,a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences. 
#-------------------------------------------------------------------------------
# SEE ALSO:
# BUUHWE , BAGIDIS.dist , ...
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
#
#if (is.vector(DATA1)){
#DATA1 = matrix(DATA1, nrow=1,byrow=TRUE)}
#if  (is.vector(DATA2) ){ 
#DATA2 =matrix(DATA2,nrow=1,byrow=TRUE) }
#
##-------------------------------------------------------------------------------
## 2. Checking for incorrect values of the arguments
##  - if all series have the same length :   ncol(DATA1) = ncol(DATA2)
##  - if wdw is in [2;ncol(DATA1)]
##  - if Overlap is in [0; wdw-1]
##  - if Param is included in [0;1]
##  - if length(p) =1 
##  - if wk=NULL or length(wk) =wdw -1 
##  - if Evol is boolean
##-------------------------------------------------------------------------------
#
#if (ncol(DATA1) != ncol(DATA2)){
#stop('ERROR: Lengths of the series in DATA1 and DATA2 should be equals')}
#
#if (wdw > ncol(DATA1)){
#wdw= min(ncol(DATA1),30)
#warning('WARNING: The length of the window (wdw) cannot be higher than the length of the series. wdw is set to its default value.')
#}
#
#if (wdw < 2){
#stop('ERROR: The length of the window (wdw) cannot be smaller than 2.')
#}
#
#if (Overlap > wdw-1 ){
#stop('ERROR: The Overlap parameter cannot be higher than the wdw-1.')
#}
#
#if (Overlap < 0){
#stop('ERROR: the Overlap parameter must be non negative.')
#}
#
#if (Param < 0 | Param > 1){
#stop('ERROR: the balance parameter Param must be included in [0;1].')}
#
#
#if (length(p) !=1){
#stop('ERROR: p must be a single numeric value, not a vector.')}
#
#if (!is.null(wk) & length(wk)!= (wdw-1) ){
#stop('ERROR: wk must be a numeric vector of length wdw-1, or must be set to NULL (so that its default value is applied).')}
#
#if (!is.logical(Evol)){
#warning('WARNING: Evol must be boolean. It is set to its default value (FALSE).')} 
 
#-------------------------------------------------------------------------------
# 3. Determining how many subseries will have to be selected from each series.
# Note: If  wdw = ncol(DATA1), NbSubseries =1.  
#-------------------------------------------------------------------------------

Decalage = wdw-Overlap       # Default: =1.

NbSubseries= (ncol(DATA1) - wdw +1) # Nb of subseries in case of Decalage = 1.
NbSubseries_overlap= floor((ncol(DATA1) - wdw)/(wdw-Overlap))+1 #ceiling(NbSubseries/Decalage ) # Nb of subseries for which the semidistance will be calculated. If Overlap is set to its default value, NbSubseries_overlap = NbSubseries.

if (NbSubseries_overlap !=  (ncol(DATA1) - wdw)/(wdw-Overlap) +1) {
warning('WARNING: the length of the series (ncol(DATA1)) is not a multiple of the Overlap parameter. The series are thus truncated to the closer multiple of Overlap. ')
}


#-------------------------------------------------------------------------------
# 4. Initializing the dissimilarity matrices 
#-------------------------------------------------------------------------------

dissimilarity.matrix = matrix(nrow=nrow(DATA1),ncol=nrow(DATA2))

if (Evol) { 
dissimilarity.evol = array(dim=c(nrow(DATA1),nrow(DATA2),NbSubseries_overlap) )}

#-------------------------------------------------------------------------------
# 5. Transformation of the series of DATA2 into subseries and BUUHWE transform of those subseries.
#-------------------------------------------------------------------------------

## Defining a matrix of subseries whose BUUHWE transform has to be computed:
## Dimensions of this matrix wdw*(NbSubseries_overlap*nrow(DATA2))
## Subseries to be transformed are encoded column by column.
        
        DefineSubseries =function(Y,NbSubseries,wdw,Decalage){                          
        ### Transforming  the series j in subseries
        Y.Windowed_tmp = matrix(NA,nrow=wdw, ncol=NbSubseries)
        for (k in 1:NbSubseries){Y.Windowed_tmp[,k]= Y[k:(k+wdw-1)]}        
        ### Select subseries depending on the Ovelap parameter
        Y.Windowed = Y.Windowed_tmp[,seq(1,ncol(Y.Windowed_tmp),by =Decalage)]
        ## Subseries are returned as columns of the matrix Y.Windowed.
        return(Y.Windowed) }

        # Apply function DefineSubseries to each row of DATA2. 
        
DATA2_Subseries = matrix(apply(DATA2,1,DefineSubseries, NbSubseries,wdw,Decalage), nrow =wdw, ncol = NbSubseries_overlap*nrow(DATA2))
        
## BUUHWE transform of each column of  DATA2_Subseries, and encoding as an array of list.
## This array has nrow(DATA2) columns (each columns corresponds to subseries of one single series in DATA2) and NbSubseries_overlap rows.  

BUUHWE.DATA2 = array(apply(DATA2_Subseries,2,BUUHWE), dim=c(NbSubseries_overlap,nrow(DATA2)))       



#-------------------------------------------------------------------------------
# 6. Computation of the dissimilarities
#-------------------------------------------------------------------------------


### Definition of a function that computes the dissimilarities for one row of DATA1 and one row of DATA2
    
    ComputeDissimVec = function(BUUHWE.DATA2.rowj, BUUHWE.DATA1.rowi, NbSubseries_overlap, p, wk, Param){ 
        Dissim = rep(NA, NbSubseries_overlap)
        for(k in 1:NbSubseries_overlap){Dissim[k]= BAGIDIS.dist(BUUHWE.DATA1.rowi[[k]],
                                           BUUHWE.DATA2.rowj[[k]],
                                           p=p,
                                           wk =wk,
                                           Param=Param)} # eo FOR ( index in 1:NbSubseries_overlap) 
          ## --> The vector Dissim contains the sequence of partial 
          ##     dissimilarities between series i and j. 
    return(Dissim)    
    } 
    

### Computation of the dissimilarities
        
if (identical(DATA1,DATA2)){ #----------------------------------if DATA1 = DATA2

for (i in 1:nrow(DATA1)){  #### We work row by row in DATA1 
                                 ### for each row of DATA1,
                                 ### we select successively BUUHWE transforms 
                                 ### of each Row of DATA2 and compute                
                                 ### the local semidistance.
                                 ### That is, we select  successively 
                                 ### each column of BUUHWE.DATA2 

BUUHWE.DATA1.rowi =  BUUHWE.DATA2[,i]

DissimRowsData2toDATA1rowi = apply(array(BUUHWE.DATA2[,i:ncol(BUUHWE.DATA2)],dim =c(NbSubseries_overlap,nrow(DATA2)-i+1)), 2, ComputeDissimVec,BUUHWE.DATA1.rowi, NbSubseries_overlap, p, wk, Param) 
 
dissimilarity.matrix[i,i:ncol(BUUHWE.DATA2)]=  apply(matrix(DissimRowsData2toDATA1rowi, ncol =nrow(DATA2)-i+1) ,2 ,sum)        
dissimilarity.matrix[i:ncol(BUUHWE.DATA2),i]=  apply(matrix(DissimRowsData2toDATA1rowi, ncol =nrow(DATA2)-i+1) ,2 ,sum)        

if (Evol) { 
dissimilarity.evol[i,i:ncol(BUUHWE.DATA2),] = t(DissimRowsData2toDATA1rowi)       
dissimilarity.evol[i:ncol(BUUHWE.DATA2),i,] = t(DissimRowsData2toDATA1rowi) }          
     
} # eo FOR 

} else { #----------------------------------------------------------------- ELSE

for (i in 1:nrow(DATA1)){  #### We work row by row in DATA1 
                                 ### for each row of DATA1,
                                 ### we select successively BUUHWE transforms 
                                 ### of each Row of DATA2 and compute                
                                 ### the local semidistance.
                                 ### That is, we select  successively 
                                 ### each column of BUUHWE.DATA2 

### Transformation of series i in subseries
### Subseries are returned as colums of matrix X.Windowed. This matrix has dimension wdw*NbSubseries_overlap.
    X= as.numeric(DATA1[i,])
    X.Windowed = DefineSubseries(X,NbSubseries,wdw,Decalage)
       
###  BUUHWE  transform of subseries on row i
    
    BUUHWE.DATA1.rowi  = apply(matrix(X.Windowed,ncol=NbSubseries_overlap)  , MARGIN=2,FUN= BUUHWE)    
           ## list of lists


### Computes dissimilarities between row i of DATA1 and every rows of DATA2 

DissimRowsData2toDATA1rowi = apply(BUUHWE.DATA2, 2, ComputeDissimVec,BUUHWE.DATA1.rowi, NbSubseries_overlap, p, wk, Param) 
       ## Column j contains the NbSubseriesOverlap local dissimilarities of row j of DATA2 and row i of DATA1.  Dimensions: NbSubseries_overlap * ncol(DATA2)

### Encoding dissimilarities of row i of DATA1 to rows of DATA2 in dissimilarity.matrix and dissimilarity.evol

dissimilarity.matrix[i,]=  apply(matrix(DissimRowsData2toDATA1rowi, ncol =nrow(DATA2)),2 ,sum)        

if (Evol) { 
dissimilarity.evol[i,,] = t(DissimRowsData2toDATA1rowi) }       
   
} # eo FOR (row in DATA1)    


} # eo IF (DATA1 ?= DATA2)

#-------------------------------------------------------------------------------
#  7. Return values of the function
#-------------------------------------------------------------------------------

if (!Evol){
return(dissimilarity.matrix)}
else{
return(list('dissimilarity.matrix' =dissimilarity.matrix,
              'dissimilarity.evol' =dissimilarity.evol))}

}#EOF

#===============================================================================




