#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified


########################################################
##################### fisherextest #####################
########################################################
#FISHEREXTEST - Fisher's Exact Probability Test
#   FISHEREXTEST performs the Fisher exact probability test for a table of 
#   frequency data cross-classified according to two categorical variables, 
#   each of which has two levels or subcategories (2x2). It is a non-parametric
#   statistical test used to determine if there are nonrandom associations 
#   between the two categorical variables. Fisher's exact test is used to 
#   calculate an exact P-value with small number of expected frequencies, for
#   which the Chi-square test is not appropriate (in case the total number of
#   observations is less than 20 or there is a cell-value less than 5). This
#   test is based upon the hypergeometric probability. The test was proposed in
#   the 1934 edition of the famous Ronald Aylmer Fisher's book 'Statistical 
#   Methods for Research Workers'.
#
#   So, according to the next 2x2 table design,
#
#                    Var.1
#                --------------
#                  a       b      r1=a+b
#         Var.2
#                  c       d      r2=c+d
#                --------------
#                c1=a+c  c2=b+d  n=c1+c2
#
#   The Fisher's exact test it is conditioned to all the a+b values (all the 
#   possible cell combinations that would still result in the marginal frequencies
#   as highlighted) such that,
#
#          H(a+b) = {X: X  H and a+b = r1}.
#
#   Then,
#
#          P(a,b|a+b,n) = P(X|X  H(r1)) = b(a;c1,p)*b(b;c2,p)/b(a+b;n,p).
#
#   This binomials relationship reduces each 2x2 table to the exact hypergeometric 
#   distribution to compute the P-value. Now considering only the binomial 
#   coefficients,
#
#          P(X|X H(r1)) = C(c1,a)*C(c2,b)/C(n,a+b)
#
#   Thereby Fisher's exact P-values are readily evaluated as,
#
#          P = SUM(P(X|X H(r1))).        
#
#   IMPORTANT: Due that Matlab could not work for factorials greater than 170. We
#   use the function sum([log(x+1)......]), in order to avoid further calculation 
#   problems.
#
#   Syntax: function Fisherextest(a,b,c,d) 
#      
#   Inputs:
#         a,b,c,d - observed frequency cells
#
#   Output:
#         A table with three p-values [decide to use Left, Right or 2-Tail before
#         collecting (or looking at) the data]: 
#           - Left tail: Use this when the alternative to independence is that there is 
#             negative association between the variables. That is, the observations
#             tend to lie in lower left and upper right. 
#           - Right tail: Use this when the alternative to independence is that there is
#             positive association between the variables. That is, the observations 
#             tend to lie in upper left and lower right. 
#           - 2-Tail: Use this when there is no prior alternative. 
#
#   Example: From the example given on the Ina Parks S. Howell's internet homepage 
#            (http://www.fiu.edu/~howellip/Fisher.pdf). Suppose Crane and Egret are
#            two very small collages, the results of the beginning physics course at
#            each of the two schools are given in the follow table.
#
#                                   Physics
#                              Pass         Fail
#                            ---------------------
#                      Crane     8           14
#            Collage                 
#                      Egret     1            3
#                            ---------------------
#                                       
#   Calling on Matlab the function: 
#             Fisherextest(8,14,1,3)
#
#   Answer is:
#
#   Table of the Fisher's exact test results for:
#   a = 8, b = 14, c = 1, d = 3
#   ----------------------------------------------
#                     P-value  
#   ----------------------------------------------
#     Left tail       Right tail         2-tail   
#   ----------------------------------------------
#   0.8412876883     0.5685618729     1.0000000000
#   ----------------------------------------------
#
#  Created by A. Trujillo-Ortiz, R. Hernandez-Walls and A. Castro-Perez
#             Facultad de Ciencias Marinas
#             Universidad Autonoma de Baja California
#             Apdo. Postal 453
#             Ensenada, Baja California
#             Mexico.
#             atrujo@uabc.mx
#             And the special collaboration of the post-graduate students of the 2004:2
#             Multivariate Statistics Course: Laura Rodriguez-Cardozo, Norma Alicia Ramos-Delgado,
#             and Rene Garcia-Sanchez. 
#  Copyright (C) September 27, 2004.
#
#  To cite this file, this would be an appropriate format:
#  Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez, L. Rodriguez-Cardozo 
#    N.A. Ramos-Delgado and R. Garcia-Sanchez. (2004). Fisherextest:Fisher's Exact
#    Probability Test. A MATLAB file. [WWW document]. URL http://www.mathworks.com/
#    matlabcentral/fileexchange/loadFile.do?objectId=5957
#
#  References:
# 
#  Agresti, A. (1992), A Survey of Exact Inference for Contegency Tables. 
#           Statistical Science,7:131-153.
#  Fisher, R.A. (1934), Statistical Methods for Research Workers. Chapter 12. 
#           5th Ed., Scotland:Oliver & Boyd.
#  Howell, I.P.S. (Internet homepage), http://www.fiu.edu/~howellip/Fisher.pdf
#  Zar, J.H. (1999), Biostatistical Analysis (2nd ed.). NJ: Prentice-Hall,
#           Englewood Cliffs. p. 543-555. 
#

fisherextest <- function(a,b,c,d){


tt <- rbind(cbind(a,b),cbind(c,d))
tt <- as.matrix(tt)
t  <- rbind(cbind(a,b),cbind(c,d))
t  <- as.matrix(t)

pp <- t[,1]/rowSums(t)

 if(pp[1]<pp[2]){
   if(a==0){
    t<-t
   }else{
    t<-t[c(2,1), ]
   }
 }

a <-t[1,1]
b <-t[1,2]
c <-t[2,1]
d <-t[2,2]
r1 <- a+b
r2 <- c+d
c1 <- a+c
c2 <- b+d
n <- c1+c2


ap <- 0: min(cbind(r1,c1))
bp <- r1-ap
cp <- c1-ap
dp <- r2-cp
len <- length(ap)

###############################
# Hier wegen NANs in R ########
###############################

if(r1<=0){r1<-1}
if(r2<=0){r2<-1}
if(c1<=0){c1<-1}
if(c2<=0){c2<-1}
if(n<=0){n<-1}

###############################

factor1<-sum(c(log(1:r1),log(1:r2),log(1:c1),log(1:c2),-log(1:n)))

Pap <- matrix(0,1,len)

######### Hier wegen NANs in R ####

for (x in 1:len){
if(ap[x]<=0){ap[x]<-1}
if(bp[x]<=0){bp[x]<-1}
if(cp[x]<=0){cp[x]<-1}
if(dp[x]<=0){dp[x]<-1}
}
###########################################################

for (k in 1:len){
 factor2 <- sum(c(log(1:ap[k]),log(1:bp[k]),log(1:cp[k]),log(1:dp[k])))
 Pap[k] <- factor1-factor2
}


Pap <- exp(Pap)
P1 <- sum(Pap[which(ap<=a)])
if(P1>1){P1<-round(P1)}else{P1<-P1}

P2 <- sum(Pap[which(ap>=a)])
 if(P2>1){P2 <- round(P2)}else{P2<-P2}

factor2 <- sum(c(log(1:a),log(1:b),log(1:c),log(1:d)))
Pa <- exp(factor1-factor2)
P3<- sum(Pap[which(Pap<=Pa)])
if(P3>1){P3 <- round(P3)}else{P3<-P3}


a<- tt[1,1]
b<- tt[1,2]
c<- tt[2,1]
d<- tt[2,2]

p<-P3

return(p)
}
