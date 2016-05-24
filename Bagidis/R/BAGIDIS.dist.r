   
BAGIDIS.dist= function(BUUHWE.out.1,
                       BUUHWE.out.2, 
                       p = 2, 
                       wk=NULL, 
                       Param=0.5){
#===============================================================================
# FUNCTION : BAGIDIS.dist
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 8/15/2010 10:59:17 AM
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# Computation of the BAGIDIS semidistance between two series, encoded as argument of the function by their BUUHWE expansion. 
#-------------------------------------------------------------------------------
# USAGE: 
# BAGIDIS.dist( BUUHWE.out.1,
#              BUUHWE.out.2, 
#              p = 2, 
#              wk=NULL, 
#              Param=0.5)
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# BUUHWE.out.1 : BUUHWE expansion of a series, as obtained from function BUUHWE.
# BUUHWE.out.2 : BUUHWE expansion of a second series, as obtained from function BUUHWE.
# p: the kind of norm to be used for computing the partial distance in the B-D plane. Must be numeric or Inf.
# wk: a vector of weights having length of the series -1. If not provided, wk = log(N+1-(1:N))/log(N+1) with N= ncol(DATA1)-1 (this is the default of function BAGIDIS.dist).
# Param : the balance parameter between the differences along the breakpoint axis and along the detail axis. Param must be in [0;1]. Param= 1 means that only breakpoints differences are taken into account. Param=0 means that only details differences are taken into account.
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
#-------------------------------------------------------------------------------
# VALUE:
# a numeric being the BAGIDIS semidistance between the two series. 
#-------------------------------------------------------------------------------
# WARNING/NOTE :
# There is no check for the validity of the parameters in this function as it will be mainly used as a call within function semimetric.BAGIDIS, that already performs that check. This allows for avoiding identical checks to be performs repetetively when BAGIDIS.dist is called form semimetric.dist. Conditions for the parameters are:
#  - at least BUUHWE.out.1 and BUUHWE.out.2 must be provided as arguments
#  - if both series have the same length (i.e length(BUUHWE.out.1$detail) = length(BUUHWE.out.2$detail) ) 
#  - if Param is included in [0;1]
#  - if length(p) =1. p must be numeric or Inf. 
#  - if wk=NULL or length(wk) = length(BUUHWE.out.1$detail) -1 
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS,a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences. 
#-------------------------------------------------------------------------------
# SEE ALSO:  
# semimetric.BAGIDIS, BUUHWE, ...
#-------------------------------------------------------------------------------
# EXAMPLES:
#> x= 1:10
#> y=10:1
#> BAGIDIS.dist(BUUHWE(x),BUUHWE(y))
#[1] 30.01798
#> BAGIDIS.dist(BUUHWE(x),BUUHWE(y), p=Inf, wk= 1:9)
#[1] 337.3122
#===============================================================================
#===============================================================================

#-------------------------------------------------------------------------------
# 1. Select details and breakpoints in both series
#-------------------------------------------------------------------------------

if (p!= Inf){
det.S1  = ((1-Param)^(1/p)) * BUUHWE.out.1$detail[-1]
break.S1= (Param^(1/p)) * BUUHWE.out.1$split.abs[-1]

det.S2  =  ((1-Param)^(1/p)) * BUUHWE.out.2$detail[-1]
break.S2= (Param^(1/p)) * BUUHWE.out.2$split.abs[-1] }
else
{    # if p=Inf
det.S1  = (1-Param) * BUUHWE.out.1$detail[-1]
break.S1= Param * BUUHWE.out.1$split.abs[-1]

det.S2  =  (1-Param) * BUUHWE.out.2$detail[-1]      
break.S2= Param * BUUHWE.out.2$split.abs[-1] }

#-------------------------------------------------------------------------------
# 2. Define weights (if unspecified in the external call of the function )
#-------------------------------------------------------------------------------

if (is.null(wk)){
N =length(det.S1)                                                       
wk = log(N+1-(1:N))/log(N+1) }

#-------------------------------------------------------------------------------
# 3. Computation of partial dissimilarities (in norm p)
#-------------------------------------------------------------------------------

if (p!= Inf){
dk =  ( (abs(break.S1-break.S2))^p + (abs(det.S1-det.S2))^p  )^(1/p)
} else {   #if p= Inf
dk = max( abs(break.S1-break.S2), abs(det.S1-det.S2) )    }
#print(dk)
#-------------------------------------------------------------------------------
# 4. Computation of the dissimilarity (semidistance)
#-------------------------------------------------------------------------------

dissimilarity =sum(wk*dk)

#-------------------------------------------------------------------------------
# 5. return the value of the semidistance
#-------------------------------------------------------------------------------

return(dissimilarity) } # EOF
