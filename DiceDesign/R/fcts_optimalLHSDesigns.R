
############################################################################
# Elementary functions for:
# discrepSA_LHS,discrepESE_LHS, maximinSA_LHS and maximinESE_LHS
############################################################################

# EP  = Elementary Permutation

#####discrepancyTERM_L2_STAR#####
#---------------------------------------------------------------------------|
#args :  i,j : two points x_i and x_j from the design m                     |      
#        m   : the design                                                   |
#depends :  both "term_oneDD" and "term_twoDD"                              |
#---------------------------------------------------------------------------|

discrepancyTERM_L2_STAR <- function(i,j,m)
{
  n<-nrow(m)
  d<-ncol(m)
  if(i!=j)
  { t<-c()
    for (l in 1:d)
    {t<-c(t,1-max(m[i,l],m[j,l]))}
    t<-(prod(t))/(n^2)  
  }
  else 
  {t<-term_oneDD(i,m)/(n^2)-((2^(1-d))/n)*term_twoDD(i,m)}
  return(t)
}

#####term_twoDD#####

term_twoDD=function(i,m)
  
{
  t2<-1-m[i,]^2
  t2<-prod(t2)
  return(t2)
}

#####term_oneDD#####

term_oneDD=function(i,m)
  
{
  t1<-1-m[i,]
  t1<-prod(t1)
  return(t1)
}


#####FUNCTION PERFORMING ELEMENTARY PERMUTATION (EP) IN LHD##
#####USED IN SA ALGORITHMS#####

#####lhs_EP#####
#---------------------------------------------------------------------------|
#args :  m = the design                                                     |
#out  :  l = list including design after EP, ligns and columns defining EP  |          
#---------------------------------------------------------------------------|

lhs_EP<-function(m) 
{  
  G<-m
  if (!is.matrix(G)) G<-m[[1]]
  d<-ncol(G)
  n<-nrow(G)
  ligns<-trunc(runif(2,0,n))+1      
  column<-trunc(runif(1,0,d))+1
  x<-G[ligns[1],column]
  G[ligns[1],column]<-G[ligns[2],column]
  G[ligns[2],column]<-x 
  l<-list(G,ligns[1],ligns[2],column)   
  return(l)
} 

#####phiP_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the power in phiP criterion                                |
#out     l     : a list with new design and its mindist value               |
#---------------------------------------------------------------------------|

phiP_EP_ESE<-function(m,k,p) 
{  
  d<-ncol(m)
  n<-nrow(m)
  G<-m
  ligns<-trunc(runif(2,0,n))+1      
  x<-G[ligns[1],k]
  G[ligns[1],k]<-G[ligns[2],k]
  G[ligns[2],k]<-x 
  l<-list(phiP(G,p),G)   
  return(l)
} 


#####INTERMEDIATE FUNCTIONS#####

alpha_oneDD<-function(i1,i2,k,m)
{
  
  alpha<-(1-m[i2,k])/(1-m[i1,k])
  alpha
}

beta_oneDD<-function(i1,i2,k,m)
{
  
  beta<-(1-m[i2,k]^2)/(1-m[i1,k]^2)
  beta
}


gamma_oneDD<-function(i1,i2,k,j,m)
{
  
  gamma<-(1-max(m[i2,k],m[j,k]))/(1-max(m[i1,k],m[j,k]))
  gamma
}

#####DESIGN'S L2 STAR DISCREPANCY VALUE AFTER AN EP#####

#####discrepancyL2_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the L2_star_discrepancy of m                 |
#out     l : list with dL2_2 (the square of the L2_star_discrepancy         |
#             of the EP design and the matrix corresponding to the EP design|                    
#depends :  term_oneDD, term_twoDD, alpha_oneDD, beta_oneDD,                |
#           gamma_oneDD, discrepancyTERM_L2_STAR                            |
#---------------------------------------------------------------------------|

discrepancyL2_EP_ESE<-function(m,k,p)
{
  G<-m
  i<-trunc(runif(2,0,nrow(m)))+1   # On genère un autre LHS à partir de M, en le perturbant de façon minimale
  x<-G[i[1],k]
  G[i[1],k]<-G[i[2],k]
  G[i[2],k]<-x 
  i1<-i[1]
  i2<-i[2]
  
  
  n<-nrow(m)
  d<-ncol(m)
  dL2_2<-p+(term_oneDD(i1,m)*alpha_oneDD(i1,i2,k,m))/(n^2)-beta_oneDD(i1,i2,k,m)*term_twoDD(i1,m)*((2^(1-d))/n)-discrepancyTERM_L2_STAR(i1,i1,m)
  dL2_2<-dL2_2+(term_oneDD(i2,m)/(alpha_oneDD(i1,i2,k,m)*(n^2)))-((term_twoDD(i2,m)*((2^(1-d))/n))/(beta_oneDD(i1,i2,k,m)))-discrepancyTERM_L2_STAR(i2,i2,m)
  y<-rep(0,n-2)
  w<-0
  for (j in 1:n)
  {if (j!=i1 & j!=i2)
  {
    w<-w+1
    y[w]<-gamma_oneDD(i1,i2,k,j,m)*discrepancyTERM_L2_STAR(i1,j,m)-discrepancyTERM_L2_STAR(i1,j,m)+(discrepancyTERM_L2_STAR(i2,j,m)/(gamma_oneDD(i1,i2,k,j,m)))-discrepancyTERM_L2_STAR(i2,j,m)
  }
   
  }
  
  y<-sum(y)
  dL2_2<-dL2_2+2*y
  l<-list(dL2_2,G)
  return(l)
}

#####DESIGN'S L2 DISCREPANCY VALUE AFTER AN EP#####


#####discrepancyL2_EP#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        i1,i2 : the two ligns of the lhs elementary permutation            |              
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the L2_star_discrepancy of m                 |
#out     dL2_2 : the square of the L2_star_discrepancy of the EP design     |
#depends :  term_oneDD, term_twoDD, alpha_oneDD, beta_oneDD,                |
#           gamma_oneDD, discrepancyTERM_L2_STAR                            |
#---------------------------------------------------------------------------|

discrepancyL2_EP<-function(m,i1,i2,k,p)
{
  n<-nrow(m)
  d<-ncol(m)
  dL2_2<-p+(term_oneDD(i1,m)*alpha_oneDD(i1,i2,k,m))/(n^2)-beta_oneDD(i1,i2,k,m)*term_twoDD(i1,m)*((2^(1-d))/n)-discrepancyTERM_L2_STAR(i1,i1,m)
  dL2_2<-dL2_2+(term_oneDD(i2,m)/(alpha_oneDD(i1,i2,k,m)*(n^2)))-((term_twoDD(i2,m)*((2^(1-d))/n))/(beta_oneDD(i1,i2,k,m)))-discrepancyTERM_L2_STAR(i2,i2,m)
  y<-rep(0,n-2)
  w<-0
  for (j in 1:n)
  {if (j!=i1 & j!=i2)
  {
    w<-w+1
    y[w]<-gamma_oneDD(i1,i2,k,j,m)*discrepancyTERM_L2_STAR(i1,j,m)-discrepancyTERM_L2_STAR(i1,j,m)+(discrepancyTERM_L2_STAR(i2,j,m)/(gamma_oneDD(i1,i2,k,j,m)))-discrepancyTERM_L2_STAR(i2,j,m)
  }
   
  }
  
  y<-sum(y)
  dL2_2<-dL2_2+2*y
  
  return(dL2_2)
}

#####gammaWDD#####

gammaWDD=function(i1,i2,k,j,m)
{
  gamma<-((3/2)-(abs(m[i2,k]-m[j,k])*(1-abs(m[i2,k]-m[j,k]))))/((3/2)-(abs(m[i1,k]-m[j,k])*(1-abs(m[i1,k]-m[j,k]))))
  gamma
}

#####ccWDD#####

ccWDD=function(i,j,m)
{
  n<-nrow(m)
  
  if(i!=j)
  {c<-(3/2)-((abs(m[i,]-m[j,]))*(1-abs(m[i,]-m[j,])))
   c<-(prod(c))/(n^2)}
  c
}


#####discrepancyW2_EP#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        i1,i2 : the two ligns of the lhs elementary permutation            |              
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the W2_star_discrepancy of m                 |
#out     dW2_2 : the square of the W2_star_discrepancy of the new design    |
#depends       : ccWDD, gammaWDD                                            |
#---------------------------------------------------------------------------|

discrepancyW2_EP=function(m,i1,i2,k,p)        
  #m=le LHS initial, i1 et i2=les deux lignes de la permutation, k la colonne,  DW=discrépance Wrap-around de m au carré
{
  
  n<-nrow(m)
  dW2_2<-rep(0,n-2)
  w<-0
  for (j in 1:n)
  {if (j!=i1 & j!=i2)
    
  {w<-w+1
   dW2_2[w]<-gammaWDD(i1,i2,k,j,m)*ccWDD(i1,j,m)-ccWDD(i1,j,m)+ccWDD(i2,j,m)/(gammaWDD(i1,i2,k,j,m))-ccWDD(i2,j,m)}
   
  }
  
  dW2_2<-sum(dW2_2)
  dW2_2<-p+2*dW2_2
  return(dW2_2)
  
}

#####discrepancyW2_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the W2_star_discrepancy of m                 |
#out     l     : a list with new design and it W2 discrepancy value         |
#depends       : ccWDD, gammaWDD                                            |
#---------------------------------------------------------------------------|

discrepancyW2_EP_ESE=function(m,k,p)    
{
  n<-nrow(m)
  G<-m
  i<-trunc(runif(2,0,n))+1  
  i1=i[1]
  i2=i[2]
  x<-G[i1,k]
  G[i1,k]<-G[i2,k]
  G[i2,k]<-x 
  
  
  
  dW2<-rep(0,n-2)
  w<-0
  for (j in 1:n)
  {if (j!=i1 & j!=i2)
    
  {w<-w+1
   dW2[w]<-gammaWDD(i1,i2,k,j,m)*ccWDD(i1,j,m)-ccWDD(i1,j,m)+ccWDD(i2,j,m)/(gammaWDD(i1,i2,k,j,m))-ccWDD(i2,j,m)}
   
  }
  
  dW2<-sum(dW2)
  dW2<-sqrt(p+2*dW2)
  l<-list(dW2,G)
  return(l)
  
}



#####alphaDD#####

alphaDD<-function(i1,i2,k,m)
{
  
  alpha<-(1+abs(m[i2,k]))/(1+abs(m[i1,k]))
  alpha
}

#####betaDD#####

betaDD<-function(i1,i2,k,m)
{
  
  beta<-(2-abs(m[i2,k]))/(2-abs(m[i1,k]))
  beta
}

#####gammaDD#####

gammaDD<-function(i1,i2,k,j,m)
{
  
  gamma<-(2+abs(m[i2,k])+abs(m[j,k])-abs(m[i2,k]-m[j,k]))/(2+abs(m[i1,k])+abs(m[j,k])-abs(m[i1,k]-m[j,k]))
  gamma
}

#####gDD#####

gDD<-function(i,m)
{
  g<-1+abs(m[i,])
  g<-prod(g)
  g
}

#####hDD#####

hDD<-function(i,m)
{
  h<-1+(abs(m[i,])/2)-0.5*(m[i,]^2)
  h<-prod(h)
  h
}

#####ccDD#####

ccDD<-function(i,j,m)
{
  n<-nrow(m)
  
  if(i!=j)
  {c<-0.5*(2+abs(m[i,])+abs(m[j,])-abs(m[i,]-m[j,]))
   c<-(prod(c))/(n^2)}
  else {
    c<-gDD(i,m)/(n^2)-2*hDD(i,m)/n}
  c
}


#####discrepancyC2_EP#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        i1,i2 : the two ligns of the lhs elementary permutation            |              
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the C2_star_discrepancy of m                 |
#out     dC2_2 : the square of the C2_star_discrepancy of the new design    |
#depends       : alphaDD, betaDD, gammaDD, gDD, hDD, ccDD                   |
#---------------------------------------------------------------------------|

discrepancyC2_EP=function(m,i1,i2,k,p)        #m=le LHS initial,     i1 et i2=les deux lignes de la permutation, k la colonne,  DW=discrépance Wrap-around de m au carré
{
  
  n<-nrow(m)
  m<-m-0.5
  dC2_2<-p+(gDD(i1,m)*alphaDD(i1,i2,k,m))/(n^2)-2*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)*hDD(i1,m)/n
  dC2_2<-dC2_2-ccDD(i1,i1,m)-ccDD(i2,i2,m) 
  dC2_2<-dC2_2+(gDD(i2,m)/((n^2)*alphaDD(i1,i2,k,m)))-(2*hDD(i2,m)/(n*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)))
  
  
  x<-rep(0,n-2)
  w<-0
  for (j in 1:n)
  {if (j!=i1 & j!=i2)
  {
    w<-w+1
    x[w]<-gammaDD(i1,i2,k,j,m)*ccDD(i1,j,m)-ccDD(i1,j,m)+ccDD(i2,j,m)/(gammaDD(i1,i2,k,j,m))-ccDD(i2,j,m)}
  } 
  
  x<-sum(x)
  dC2_2<-dC2_2+2*x
  return(dC2_2)
  
  
}


#####discrepancyC2_EP_ESE#####
#---------------------------------------------------------------------------|
#args :  m     : the design before the EP                                   |
#        k     : the column of lhs elementary permutation                   |
#        p     : the square of the C2_star_discrepancy of m                 |
#out     l     : a list with new design and it C2 discrepancy value         |
#depends       : alphaDD, betaDD, gammaDD, gDD, hDD, ccDD                   |
#---------------------------------------------------------------------------|

discrepancyC2_EP_ESE=function(m,k,p)        #m=le LHS initial,     i1 et i2=les deux lignes de la permutation, k la colonne,  DW=discrépance Wrap-around de m au carré
{
  n<-nrow(m)
  G<-m
  i<-trunc(runif(2,0,n))+1 
  i1=i[1]
  i2=i[2] 
  x<-G[i1,k]
  G[i1,k]<-G[i2,k]
  G[i2,k]<-x 
  
  m<-m-0.5
  dC2<-p+(gDD(i1,m)*alphaDD(i1,i2,k,m))/(n^2)-2*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)*hDD(i1,m)/n
  dC2<-dC2-ccDD(i1,i1,m)-ccDD(i2,i2,m) 
  dC2<-dC2+(gDD(i2,m)/((n^2)*alphaDD(i1,i2,k,m)))-(2*hDD(i2,m)/(n*alphaDD(i1,i2,k,m)*betaDD(i1,i2,k,m)))
  
  
  x<-rep(0,n-2)
  w<-0
  for (j in 1:n)
  {if (j!=i1 & j!=i2)
  {
    w<-w+1
    x[w]<-gammaDD(i1,i2,k,j,m)*ccDD(i1,j,m)-ccDD(i1,j,m)+ccDD(i2,j,m)/(gammaDD(i1,i2,k,j,m))-ccDD(i2,j,m)}
  } 
  
  x<-sum(x)
  dC2<-sqrt(dC2+2*x)
  l=list(dC2,G)
  return(l)  
}
