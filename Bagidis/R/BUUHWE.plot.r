BUUHWE.plot =function(BUUHWE.out,row.max=length(BUUHWE.out$series),French=FALSE, Color=TRUE){

#===============================================================================
# FUNCTION : BUUHWE.plot
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 8/15/2010 10:59:17 AM
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# function for graphical representation of the Bottom Up Unbalanced Haar Wavelet Expansion of a series.
#-------------------------------------------------------------------------------
# USAGE: 
# BUUHWE.plot(BUUHWE.out,row.max=length(BUUHWE.out$series),French=FALSE, Color=TRUE)
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# BUUHWE.out : output of BUUHWE(.).
# row.max: number of basis vectors to be represented. By defauld, all vectors  are shown.
# French: logical. should labels be written in french?(default=english)
# Color: logical. should the representation be in color (alternative: black and white)? (default = color) 
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
#
#-------------------------------------------------------------------------------
# VALUE:
# A graphical representation of the expansion of a series in its unbalanced Haar wavelet basis. 
# The series is represented on the top of the graph, with the time axis at the bottom of the main figure. This main figure is the representation of the Unbalanced Haar basis, row by row, with rank k = 0 (constant vector ) on the top. Vectors are expressed as a function of time (horizontaly) and the height of that function is given by the "blue to red" color scale. Vertically, the detail coefficients are represented rank by rank at the same level as the corresponding basis vectors. This is the expansion of the series in its associated Unbalanced Haar basis, starting from rang k= 0 on the top of the gaph - wich is represented missing for graphical purpose and continuing downwards with rang k increasing. If you sum those detail coeffcients times their associated basis vectors, you get back the initial series. 
#-------------------------------------------------------------------------------
# WARNING/NOTE :
# The use of BUUHWE.plot requires function Blue2Red.Col.r  
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS, a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences. 
# Fryzlewicz P., 2007,  "Unbalanced Haar Technique for Non Parametric Function Estimation", Journal of the American Statistical Association, 102, 1318-1327.
#-------------------------------------------------------------------------------
# SEE ALSO:  
# BUUHWE, BD.plot, ...
#-------------------------------------------------------------------------------
# EXAMPLES:

#===============================================================================
#===============================================================================

#-------------------------------------------------------------------------------
# The series
#-------------------------------------------------------------------------------

x = BUUHWE.out$series

# Defining labels (depending on the language chosen)
if(French){
MyxlabBase = "Abscisses"
MyylabBase = "Vecteurs de base"
MyylabTS =   'Serie'
MyLegLab = 'Code couleurs'
} else{
MyxlabBase = "Abscissa"
MyylabBase = "Basis vectors"
MyylabTS =   'Series'
MyLegLab = 'Color key' }

if (Color){ # representation in color or in black and white?   ------IF COLOR --

#-------------------------------------------------------------------------------
# LAYOUT
#-------------------------------------------------------------------------------
def.par <- par(no.readonly = TRUE) # save default, for resetting...

#definition of the layout
mylayout =rbind( c(3,3,3,5),
               c(1,1,1,2),
               c(1,1,1,2),
               c(1,1,1,2),          
               c(4,4,4,4))
nf <- layout(mylayout,  TRUE)
layout.show(nf)  #visualisation of the layout
#par(bg ='grey')
#===============================================================================
#### FIGURE 1: WAVELET BASIS ####
#===============================================================================
par(mar = c(bottom= 5, left =5, top =0, right=0))

nvec = length(x)
nbrow = min(nvec,row.max)
nz = 99       ### has to be even so as to have 0 = white
mymatrix = BUUHWE.out$basis[,seq(nbrow,1, by=-1)]
#xnames =   paste("t", seq(1,nrow(matrix), by=1), sep='')
myX = 1:nvec
myY = 0:(nbrow-1)
#ynames =   paste("psi", seq(nrow-1,0, by=-1))
ynames =   paste("expression(paste(psi,", seq(nbrow-1,1, by=-1),")),")
ynames =   c(ynames,paste("expression(paste(psi,", 0,"))"))
ynames = c("c(",ynames,")")
image(z =mymatrix,
      x =myX,
      y =myY,
      zlim =  c(-max(abs(mymatrix)),max(abs(mymatrix))),
      axes =FALSE,
      xlab = MyxlabBase,
      ylab = MyylabBase,
      col =Blue2Red.Col(nz+1),
      bty ='l',
      cex.lab=1.45
      )
axis(side =2, at= c(-1,myY,nbrow) , tick =TRUE, labels = c(NA,eval(parse(text =ynames)),NA) , las =1,cex.axis=1.25) 
axis(side =1, at= c(0,myX,max(myX)+1) , tick =TRUE, labels = c(NA,myX,NA), las =1,cex.axis=1.25)      

#===============================================================================
#### FIGURE 2: SERIES of DETAILS ####
#===============================================================================
 par(mar = c(bottom= 5, left =0, top =0, right=1) )
 plot(c(NA,BUUHWE.out$detail[-1])[1:nbrow],nbrow:1,type ='o',pch = 19, cex =2,
      xlim = c(-max(abs(BUUHWE.out$detail)),max(abs(BUUHWE.out$detail))), yaxt ='n', bty = ']',fg ='black', col.axis ='blue', col.lab ='blue', xlab =' Coefficients', col='blue',bty ='u',cex.axis=1.25,cex.lab=1.45)
 points(0,nbrow,pch ='x',cex =2)    
 abline(v=0,col='red')
 
# lines( c(NA,uh.out$correctdetail[-1]) ,nvec:1,type ='o', cex =2,pch = 19,col='green')
#===============================================================================
#### FIGURE 3: TIME SERIES ####
#===============================================================================
par(mar = c(bottom= 0, left =5, top =1, right=0) )
plot(1:nvec,x, type='o',pch=19,cex =2, xaxt = 'n',bty ='o',fg ='black', col.axis ='blue', col.lab ='blue', ylab =MyylabTS,col ='blue',cex.axis=1.25,cex.lab=1.45) 


#===============================================================================
#### FIGURE 4: LEGEND  OF THE WAVELET BASIS REPRESENTATION ####
#===============================================================================

par(mar = c(bottom= 5, left =13, top =1, right=9) )
mc =  -max(abs(mymatrix))
md =  max(abs(mymatrix))
ma = 1
mb = nz+1
mylegX = (ma:mb*(md-mc)/(mb-ma)) -(ma*(md-mc)/(mb-ma) ) + mc 

image(x=mylegX,y=1,z=as.matrix(1:(nz+1)),col=Blue2Red.Col(nz+1),yaxt='n',ylab='',xlab =MyLegLab,font.lab=2,cex.lab=1.3, cex.axis=1.2)


#===============================================================================
#### FIGURE 5: CUSTOMIZE MAIN GRAPH LAYOUT ####
#===============================================================================
par(mar = c(bottom= 0, left =0, top =1, right=1) )
plot(1,1, bty=']',type='n',xaxt = 'n',yaxt = 'n',xlim=c(0,4),ylim=c(0,4))
#
box(which = "inner", lty = "solid")

#===============================================================================
# BACK TO CLASSICAL LAYOUT
#===============================================================================

par(def.par)#- reset to default   

  } #eo If "Color" 
else   {  #------------------------------------ IF BLACK AND WHITE -------------

#===============================================================================
# LAYOUT
#===============================================================================
def.par <- par(no.readonly = TRUE) # save default, for resetting...

#definition of the layout
mylayout =rbind( c(3,3,3,3,4),
               c(1,1,1,1,2),
               c(1,1,1,1,2),
               c(1,1,1,1,2),
               c(1,1,1,1,2))
nf <- layout(mylayout,  TRUE)
#layout.show(nf)  #visualisation du layout
#par(bg ='grey')

#===============================================================================
#### FIGURE 1:WAVELET BASIS ####
#===============================================================================

# graphical parameters
par(mar = c(bottom= 5, left =5, top =0, right=0))

# number of basis vectors to be represented
nvec = length(x)
nbrow = min(nvec,row.max)

# preliminary for the definition and label of the axis
myX = 1:nvec
myY = 0:(nbrow-1)

ynames =   paste("expression(paste(psi,", seq(nbrow-1,1, by=-1),")),")
ynames =   c(ynames,paste("expression(paste(psi,", 0,"))"))
ynames = c("c(",ynames,")")

# representation of the wavelet basis

  ## Representation of the first vector psi_0
    ### draw the line (with a level change'step' in between the two concerned abscissa --> shift -0.5).
  
  abscisse =   (1:(length(BUUHWE.out$basis[,1])+1))-0.5
  ordonnee =  BUUHWE.out$basis[,1]   +
              +rep(ncol(BUUHWE.out$basis)-1,length(BUUHWE.out$basis[,1])) 
  ordonnee = c(ordonnee, ordonnee[length(ordonnee)])
  
  plot(abscisse, 
       ordonnee,
       type ='s', lwd=2 ,
       ylim =c(-0.5,ncol(BUUHWE.out$basis)-0.5), 
       xlim =c(0.5,length(BUUHWE.out$basis[,1])+0.5) ,
       xlab=MyxlabBase, ylab =MyylabBase,  cex.lab=1.45,
       xaxt = 'n',yaxt = 'n')
  
    ### Representation of the point of evaluation of the wavelet 
 
  points((1:length(BUUHWE.out$basis[,1])), 
         BUUHWE.out$basis[,1]  +
         +rep(ncol(BUUHWE.out$basis)-1,length(BUUHWE.out$basis[,1])),
         pch=19)
    
    ### Horizontal axis "0" of wavelet psi_0   (+ #vertical axis)
  
  abline(h = ncol(BUUHWE.out$basis)-1, lty ='dotted')
  #abline(v = ncol(BUUHWE.out$basis), lty ='dotted')
  
  ## idem for representing vectors psi_1 to psi_nbrow.  
    
  for (i in 2:ncol(BUUHWE.out$basis)){
  
  abline(h = ncol(BUUHWE.out$basis)-i, lty ='dotted') 
  #abline(v = ncol(BUUHWE.out$basis)-(i-1), lty ='dotted') 
  
  abscisse =   (1:(length(BUUHWE.out$basis[,i])+1))-0.5
  ordonnee =  BUUHWE.out$basis[,i] +
              +rep(ncol(BUUHWE.out$basis)-i ,length(BUUHWE.out$basis[,i])) 
  ordonnee = c(ordonnee, ordonnee[length(ordonnee)])
  
  lines (abscisse,ordonnee ,type ='s',lwd=2)
  
  points((1:length(BUUHWE.out$basis[,i])), 
          BUUHWE.out$basis[,i]+rep(ncol(BUUHWE.out$basis)-i,
          length(BUUHWE.out$basis[,i])), 
          pch=19,)
  
  } # eoFor
                                
# representation and label of the axe 
  
axis(side =2, at= c(-1,myY,nbrow) , tick =TRUE,
     las =1,cex.axis=1.25, las =1,cex.axis=1.25, 
     labels = c(NA,eval(parse(text =ynames)),NA) ) 
     
axis(side =1, at= c(0,myX,max(myX)+1) , tick =TRUE, 
      las =1,cex.axis=1.25,
      labels = c(NA,myX,NA))      

#===============================================================================
#### FIGURE 2: SERIES DETAILS ####
#===============================================================================
 
 # graphical parameters
 par(mar = c(bottom= 5, left =0, top =0, right=1) )
 
 # Representation of the series
 plot(c(NA,BUUHWE.out$detail[-1])[1:nbrow],(nbrow-1):0,
        type ='o',pch = 19, cex =1.5, 
        ylim = c(-0.5,ncol(BUUHWE.out$basis)-0.5),
        xlim = c(-max(abs(BUUHWE.out$detail)),max(abs(BUUHWE.out$detail))),
        yaxt ='n', bty = ']',fg ='black', col.axis ='black', col.lab ='black',          xlab =' Coefficients', 
        col='black',bty ='u',cex.axis=1.25,cex.lab=1.45)

 points(0,nbrow-1,pch ='x',cex =2)    

 abline(v=0,col='black')
 
# lines( c(NA,uh.out$correctdetail[-1]) ,nvec:1,type ='o', cex =2,pch = 19,col='green')


#=======================================================================
#### FIGURE 3: TIME SERIES ####
#=======================================================================

# graphical parameters
par(mar = c(bottom= 0, left =5, top =1, right=0) )

#representation of the series
plot(1:nvec,x, 
    xlim =c(0.5,length(BUUHWE.out$basis[,1])+0.5) ,    
    type='o',pch=19,cex =1.5, xaxt = 'n',bty ='o',fg ='black', 
    col.axis ='black', col.lab ='black', 
    ylab =MyylabTS, col ='black',
    cex.axis=1.25,cex.lab=1.45) 


#===============================================================================
#### FIGURE 5: CUSTOMIZE MAIN GRAPH LAYOUT ####
#===============================================================================
par(mar = c(bottom= 0, left =0, top =1, right=1) )
plot(1,1, bty=']',type='n',xaxt = 'n',yaxt = 'n',xlim=c(0,4),ylim=c(0,4))
#
box(which = "inner", lty = "solid")

#===============================================================================
# BACK TO CLASSICAL LAYOUT
#===============================================================================

par(def.par)#- reset to default   
  
  } # eo Else if


} # Eof
  
