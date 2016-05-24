BD.plot = function(x, 
                   y= NULL,
                   BUUHWE.out.x =BUUHWE(x), 
                   BUUHWE.out.y =BUUHWE(y),
                   French =FALSE, 
                   col = c('black','red')){  

#===============================================================================
# FUNCTION : BD.plot
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 8/15/2010 10:59:17 AM
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# function for graphical representation of the Bottom Up Unbalanced Haar Wavelet Expansion of a series in the Breakpoints-Details plane.
#-------------------------------------------------------------------------------
# USAGE: 
# BD.plot(x, y= NULL, BUUHWE.out.x =BUUHWE(x), BUUHWE.out.y =BUUHWE(y), French =FALSE, col = c('black','red'))
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# x :  a numeric vector (a series) whose BUUHWE expansion has to be computed and represented and the B-D plane 
# y:  an optional second numeric vector (a series) whose BUUHWE expansion has to be computed and superimposed to the one of x in the B-D plane.  length(y) must be equal to length(x).
# BUUHWE.out.x : output of \code{BUUHWE(x)}, included as an optional argument for saving computation time within the function if it has already been computed and saved outside the function; otherwise BUUHWE transform of x is compute at the call of function \code{BD.plot}.
# BUUHWE.out.y: output of \code{BUUHWE(y)}, included as an optional argument for saving computation time within the function if it has already been computed and saved outside the function; otherwise BUUHWE transform of y is compute at the call of function \code{BD.plot}.
# French: logical. Should labels be written in french?(default=english)
# col: vector of size one or two, indicating the colors for representing series x and - if needed -series y in the B-D plane. Defaulf is black for series x and red for series y. 
#-------------------------------------------------------------------------------
# DETAILS:
#
# Unbalanced Haar wavelet bases are orthonormal bases that are made up of one constant vector and a set of Haar-like (i.e. \textit{up-and-down}-shaped ) orthonormal wavelets whose discontinuity point between positive and negative parts is not necessarily located at the middle of its support. Using the notation of \citet{Piotr}, the general mathematical expression of one of those Haar-like wavelets is given by
# \begin{equation}\label{HaarVec}
# \boldsymbol{\phi}_{e,b,s}(t) = (\frac{1}{b-s+1}- \frac{1}{e-s+1})^{1/2} \, \Ind_{s \le t \le b} 
#                + (\frac{1}{e-b}- \frac{1}{e-s+1})^{1/2} \, \Ind_{b+1 \le t \le e},  \nonumber
# \end{equation}
# where $t =1 \ldots N$ is a discrete index along the abscissa axis, and where $s < b < e$ are values along this axis that determine the particular shape of one wavelet ($s$, $b$ and $e$ stand for \textit{start}, \textit{breakpoint} and \textit{end} respectively). Each wavelet $\boldsymbol{\phi}_{e,b,s}(t)$ is thus associated with a level change from one observation (or group of observations) to the consecutive one, and the projection of the series $\mathbf{x}(t)$ on the wavelet $\boldsymbol{\phi}_{e,b,s}(t)$ encodes the importance of the related  level change in the series. Different choices of $N-1$ sets of values $s$, $b$ and $e$ leading to orthonormal wavelets define a whole family of bases $\{ \boldsymbol{\phi}_k\}_{k=0 \ldots N-1}$.
# 
#\paragraph{The Basis Pursuit Algorithm and the Property of Hierarchy.} In 2007, \citet{Piotr} proposed an algorithm for building the unbalanced Haar wavelet basis $\{ \boldsymbol{\phi}_k\}_{k=0 \ldots N-1}$ that is best suited to a given series, according to the principle of hierarchy - namely, the vectors of this basis and their associated coefficients are ordered following the importance of the level change they encode  for describing the global shape of the series. He called it the \textit{bottom-up unbalanced Haar wavelet transform} (here-after BUUHWT). The resulting expansion is organized in a hierarchical way and avoids the dyadic restriction that is typical for classical wavelets. The family of unbalanced Haar wavelets is thus really adaptive to the shape of the series.
#{Representing the Series in the \textit{b-d} Plane.} Let us denote the ``optimal'' Unbalanced Haar wavelet expansion of a series $ \mathbf{x}^{(i)}$ as follows:
#$$ \mathbf{x}^{(i)}= \sum_{k=0}^{N-1} d_k^{(i)} \boldsymbol{\psi}_k^{(i)}, $$
#where the coefficients $d_k^{(i)}$ are the projections of  $ \mathbf{x}^{(i)}$ on the corresponding basis vectors $\boldsymbol{\psi}_k^{(i)}$ (i.e. the \textit{detail} coefficients) and  where the set of vectors $\{ \boldsymbol{\psi}_k^{(i)}\}_{k=0 \ldots N-1}$ is the Unbalanced Haar wavelet basis that is best suited to the series $\mathbf{x}^{(i)}$, as obtained using the BUUHWT algorithm.
#Let us also denote $b_k^{(i)}$, the breakpoint $b$ of the wavelet $\boldsymbol{\psi}_k^{(i)}$, $k= 1 \ldots N-1$, as defined in equation (\ref{HaarVec}). An interesting property of the basis $\{ \boldsymbol{\psi}_k^{(i)}\}_{k=0 \ldots N-1}$, that has been proved by  \citet{Piotr}, is the following: \\
#\noindent
#\textbf{Property:} The ordered set of breakpoints $\{ b_k^{(i)}\}_{k=1 \ldots N-1}$ determines the basis $\{ \boldsymbol{\psi}_k^{(i)}\}_{k=0 \ldots N-1}$ uniquely.	

#\noindent  Consequently, the set of pairs ($b_k^{(i)}, d_k^{(i)}$)$_{k=1 \ldots N-1}$ determines the shape of the series $\mathbf{x}^{(i)}$ uniquely (i.e., it determines the series, except for a change of the mean level of the series, that is encoded by the additional coefficient $d_0^{(i)}$). 
#This leads us to the possibility of representing any series $\mathbf{x}$ in the \textit{b-d} plane, i.e. the plane formed by the breakpoints and the details coefficients. 

#
#-------------------------------------------------------------------------------
# VALUE:
# A graphical representation of the expansion of a series in its unbalanced Haar wavelet basis. The  series is plotted in the plane that is defined by the values of its breakpoints and its detail coefficients. Points are numbered according to their rank in the hierarchy.
# 
#-------------------------------------------------------------------------------
# WARNING/NOTE :
#  This function requires the use of function \code{\link{BUUHWE}}.
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS, a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences.  http://www.stat.ucl.ac.be/ISpub/dp/2010/DP1030.pdf
# Fryzlewicz P., 2007,  "Unbalanced Haar Technique for Non Parametric Function Estimation," Journal of the American Statistical Association, 102, 1318-1327.
#-------------------------------------------------------------------------------
# SEE ALSO:  
# \code{\link{BUUHWE}}, \code{\link{BUUHWE.plot}}, ...
#-------------------------------------------------------------------------------
# EXAMPLES:
# x= c(1,7,3,0,-2,6,4,0,2)
# BD.plot(x)
#===============================================================================
#===============================================================================





#==============================================================
# LAYOUT
#==============================================================

def.par <- par(no.readonly = TRUE)
                              # save default, for resetting...
  
# definition du layout
mylayout =rbind( 1 , 2, 2, 2 )
nf <- layout(mylayout,  TRUE)
    #layout.show(nf)  #visualisation du layout
    
# Definition des labels (selon la langue choisie)
if (!French){
ylabSeries ='Series'
xlabSeries ='Abscissa'
Myxlab = 'Breakpoints'
Myylab = 'Details'}
else{
ylabSeries ='Serie'
xlabSeries ='Abscisses'                     
Myxlab = 'Points de discontinuite'
Myylab = 'Details'  }

    

#==============================================================
# 1. Representation de la serie
#==============================================================
# Parametres graphiques
par(mar = c(bottom= 5, left =5, top =1, right=5) )

if (!is.null(y)) {  # si 2 series
  ylimits= c(min(x,y),max(x,y))}
else  {ylimits= c(min(x),max(x)) }
    
# Representation de la serie
plot(1:length(x), x,
    type='b', pch =20, cex=1.5, col =col[1],
    ylim= ylimits,
    xlab =xlabSeries,ylab = ylabSeries, cex.lab =1.45)

if (!is.null(y)) {  # si 2 series
    lines(1:length(x),y,type='b',pch=20,cex=1.5,col=col[2])}

#==============================================================
# 2. Representation dans le plan BD
#==============================================================

# Parametres graphiques
par(mar = c(bottom= 5, left =5, top =1, right=5) )

# Selection des details et breakpoints, et definition des labels
details =BUUHWE.out.x$detail[-1]
breakpoints =BUUHWE.out.x$split.abs[-1]
id = 1:length(details)

# Representation
plot(breakpoints,details, pch =20, 
     ylim = c(min(details)-1,max(details)+1),
     xlab = Myxlab, ylab = Myylab, cex.lab=1.45, col =col[1])
text(breakpoints,details,labels =id,pos=3,offset =0.6,
     cex=1.2,font=2, col =col[1])

for( i in 1: (length(details)-1) ){
segments( breakpoints[i],details[i], 
          breakpoints[i+1],details[i+1],
          col=col[1]) }
          
if (!is.null(y)) {  #si 2 series
  details =BUUHWE.out.y$detail[-1]
  breakpoints =BUUHWE.out.y$split.abs[-1]

  points(breakpoints,details, pch =20, col =col[2])
  text(breakpoints,details,labels =id,pos=3,offset =0.6,
       cex=1.2,font=2, col =col[2])
  
  for( i in 1: (length(details)-1) ){
  segments( breakpoints[i],details[i], 
            breakpoints[i+1],details[i+1],
            col=col[2]) }

}          
  
#==============================================================
# BACK TO CLASSICAL LAYOUT
#==============================================================
  
  par(def.par)#- reset to default         
  
  } # EOF 