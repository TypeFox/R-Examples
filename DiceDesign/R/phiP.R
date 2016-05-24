
#####phiP#####

# Compute the phiP criterion 
# (Lp norm of the sum of the inverses of the design inter-point distances)
# Reference: Pronzato, L. and Muller, W.,2012, Design of computer experiments: 
#              space filling and beyond, Statistics and Computing, 22:681-701.
# A higher phiP corresponds to a more regular scaterring of design points

#---------------------------------------------------------------------------|
#args :  design     : the design                                            |
#        p          : the "p" in the Lp norm which is taken (default=50)    |
#out  : the phiP criterion                                                  |
#---------------------------------------------------------------------------|

phiP<-function(design,p=50)   
{
  D<-dist(design)
  D<-D^(-p)
  fi_p<-(sum(D))^(1/p)
  return(fi_p)  
}
