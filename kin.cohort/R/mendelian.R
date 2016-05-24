`mendelian` <-
function(p){
# function for computing mendelian conditional probabilities
#
# p = allele frequency
# tp = 2 by 2 by 2 array of conditional probabilities
#      tp(i,j,k) gives the conditional probability of genotype j
#       given the proband has genotype i for k-th type of relationship
#       i = 1(non-carrier), 2(carrier)
#       j = 1(non-carriers), 2(carriers)
#       k = 1(parent | offspring), 2(sibling),  3 (grand-parent | parent-cibs) 
q <- 1-p
tp <- array(data=0, dim=c(3, 3, 3))

# T   p    q     0           O    p2  2pq  q2       I  1   0   0
#     p/2  1/2   q/2              p2  2pq  q2          0   1   0
#     0    p     q                p2  2pq  q2          0   0   1
#

# parent or offspring
# T

tp[1,1,1] = q;       #  aa|aa
tp[1,2,1] = p;       #  Aa|aa 
tp[1,3,1] = 0;       #  AA|aa 

tp[2,1,1] = q/2      #  aa|Aa
tp[2,2,1] = 1/2      #  AA|Aa
tp[2,3,1] = p/2      #  Aa|Aa

tp[3,1,1] = 0;       #  aa|AA
tp[3,2,1] = q;       #  Aa|AA
tp[3,3,1] = p;       #  AA|AA

# sibling
# S=I/4 + T/2 +O/4

tp[1,1,2] = (1/4)*(1+q)^2;           # aa|aa
tp[1,2,2] = (1/2)*p*(1+q);           # Aa|aa
tp[1,3,2] = p*p/4;                   # AA|aa 

tp[2,1,2] = (1/4)*q*(1+q);           # aa|Aa
tp[2,2,2] = (1/2)*(1+p*q);           # Aa|Aa
tp[2,3,2] = (1/4)*p*(1+p)            # AA|Aa

tp[3,1,2] = q*q/4                    # aa|AA 
tp[3,2,2] = (1/2)*q*(1+p);           # Aa|AA  
tp[3,3,2] = (1/4)*(p+1)^2;           # AA|AA 

# grand-parents | parent-sibs
# TT=T/2 + O/2

tp[1,1,3] = (1/2)*q*(q+1)
tp[1,2,3] =  p/2 + p*q
tp[1,3,3] = (1/2)*p^2

tp[2,1,3] = q/4 + q*q/2
tp[2,2,3] = p*q + 1/4
tp[2,3,3] = p/4 + p*p/2

tp[3,1,3] = q*q/2
tp[3,2,3] = p*q + q/2
tp[3,3,3] = p*p/2 + p/2

# cousins   T^3 = T/4 + O*3/4
#

tp
}

