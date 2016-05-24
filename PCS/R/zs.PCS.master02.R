##### S.PCS.MASTER collects all of the functions necessary for the PCS code
##### Created: 3/26/09
############################################################################
##### Table of Contents
#1. PdofCSt.cyc2: master exact function for d-best
#2. PdofCSt.T1or2: exact d-best for case T=1 or T=2
#3. PofCSGt:       exact G-best
#4. PofCSt:        exact non-d or G best, with t>=1
#5. G.H.Quad:      performs Gauss-Hermite quadrature in PofCSt
#6. PdofCSGt.bootstrap5: bootstrap for d-best and G-best
#7. PdofCSLt.bootstrap5: bootstrap for L-best
#8. PdCSGt.bootstrap.NP2: non-parametric bootstrap for d-best and G-best
#9. tindep: obtain t-statistics inside of PdCSGt.bootstrap.NP2
#10. PCS.exact: wrapper for PdofCSt.cyc2 and PofCSGt
#11. PCS.boot.par: wrapper for PdofCSGt.bootstrap5 and PdofCSLt.bootstrap5
#12. PCS.boot.np: wrapper for PdCSGt.bootstrap.NP2

##### 1 ####################################################################
PdofCSt.cyc2 <- function(theta, T, d, m=20, tol=1e-8) {
#PDOFCST.CYC2 estimates Pd(CSt) by intelligently cycling through
#and grabbing only the probability contributions above the tolerance
#parameters.  This version is debugged and generalized to accomodate
#all cases of T and d.  As such, it uses two subroutines, 
#PofCSt and PdofCSt.T1or2
#PDOFCST.CYC2 is identical to .CYC except the tolerance input parameters have been consolidated, 
#  for ease of insertion into the wrapper [PCS.exact]; and the print lines have been commented out
#theta: vector of statistics (standardized to have common mean 1)
#T: number of populations desired
#d: tolerance within the Tth population researcher is willing to accept populations as same
#m: Gauss-Hermite quadrature parameter, degree of polynomial fitting curve = 2*m
#tol is a single tolerance parameter, which is given to control the following four tolerances:
#  tol.1: tolerance parameter for individual components of "for (g in G:0)" loop
#  tol.2: tolerance parameter for probability chunks of "for (g in G:0)" loop
#  tol.3: tolerance parameter for probability in a new cycle of the "for (g in G:0)" being too small (inside repeat loop); key tolerance parameter
#  tol.4: tolerance parameter to all "PofCSt" calls
#new.out is a 1x3 vector consisting of the following parts:
#  count: the number of the final combination reached in searching for probability components > tol
#  p:     the cumulative probability = PdCSt = the answer
#  PgCSt.sum: the final probability added to the sum  
#Example: PdofCSt.cyc2(theta, T=2, d=0, m=20)


#Initialization for PdofCSt part of code:
#require(statmod)       #Library for gauss.quad.prob
tol.1=tol; tol.2=tol; tol.3=tol; tol.4=tol   #Set tolerance parameters (PdofCSt.cyc allows user control of each parameter)
y = gauss.quad.prob(m, dist="normal")$nodes #Load the constants for GH quadrature
k = length(theta)
cut1 = theta[k-T+1]-d-0.00000001 #Note: the 0.00000001 is a fudge factor to circumvent R's bug of theta[123]=1>=1=cut returning FALSE
cut2 = theta[k-T+1]+d
t1 = sum(ifelse(cut2<theta, 1, 0))                   #Count of fixed values
t3 = sum(ifelse((cut1<=theta)&(theta<=cut2), 1, 0))  #Count of indifferent values.
t2 = T-t1
t3 = t3-t2
G = t3    #G is the old notation, but I'm retaining it for simplicity
GG=G  #G begins as t3(d), but is updated according to cyclic permutation needs
      #GG is a constant that remains original G(d) throughout
if ((k-GG-T) == 0) {lower = NULL
} else lower = theta[1:(k-GG-T)]   #Vector of fixed lower theta's; "If" accounts for the case where all pops. are chosen


#Initialization for cyclic part of code:
i=t2-2; j=0; 
count=0
end.cyc.prob=0
mins=t2:3  #Vector of minimum values for each row of combinations
top.rows=(t2+G):(G+3)  #G+3 = (t2+G) - (t2-2) + 1, Top minus debit plus one
p=0
out=NULL
inds = 1:(GG+t2)  #Indices used in creating Abar

#Preliminary 'if' statement: don't run body if t2=1 or 2 because of the algorithm PdofCSt.cyc; use special subroutine
if (t2==1 || t2==2) { 
  out = PdofCSt.T1or2(theta, T, d, m, tol.4) #Subroutine for case T=1 or 2
  new.out = c(choose(t2+G,t2), out, out)     #Put results into common form
  return(new.out)                            #Return output
} else {  #All cases where t2>2

#Preliminary 'if' statement: don't run body if G = 0; P(CSt) = p, d=0 will be the correct result
if (G == 0) {
  PdCSt = PofCSt(theta,T,m,tol.4)
  new.out=c(1,PdCSt[1],PdCSt[1])
  return(new.out)
} else {  #All cases where G>0

#Master Repeat Loop: Controls all cases t2>2 & G>0
  repeat {
    #Loop working with the specific combination cycles
    PgCSt = 0                               #Initialize probability value
    for (g in G:0) { 
      #Initializations
      res = matrix(NA, nrow=t2, ncol=(g+1))  #Generate empty matrix
      res[1:(t2-2),]=top.rows                #Fill matrix with top rows
      res[t2-1,]=g+2                         #Fill row T-1
      res[t2,]=(g+1):1                       #Fill bottom row
      count = count+res[t2,1]                #Combinatoric index

      #Construct Matrices A and Abar
      A = res      #Set A equal to "res", a portion of true A
      Abar = matrix(NA, nrow=GG, ncol=(g+1))   #Initialize Abar -> needs to be NA's
      for (ii in 1:(g+1)) {                    #Loop that obtains Abar
        bool = (inds %in% A[,ii])            #by selecting elements
        Abar[,ii] = inds[!bool]  #not in A
      } #end for ii
      A = t(A);     Abar = t(Abar)          #Transpose to change form
      A = A+k-GG-T; Abar = Abar+k-GG-T      #Add constant to make correct indices

      #Obtain PgCSt's by creating theta permutations and calling PofCSt
      R = dim(A)[1] #Initialize number of index sets; This is g+1
      PrCSt = 0     #Initialize probability value
      for (r in 1:R) {
        if (t1==0) {theta.top = NULL  #Critical "if" statement, accounting for special case of t1=0       
        } else theta.top = theta[(k-t1+1):k]
        theta.star = c(lower, theta[Abar[r,]], theta[A[r,]], theta.top) #Concatenate the 4 pieces of vector together; note 2nd vector reads indices bottom to top
        PrCSt[r] = PofCSt(theta.star,T,m,tol.4)[1]
        if (PrCSt[r]<tol.1) break #Break when individual components become too small
      } #end for r

      PgCSt[G-g+1] = sum(PrCSt)          #Probability chunk, carefully indexed to be in order and positive
      if (PgCSt[G-g+1] < tol.2) break    #Break when chunks become too small
    } #end for g

    PgCSt.sum = sum(PgCSt)           #Sum of all chunks through a particular cycle
    p = p + PgCSt.sum                #Running sum of all probability
    new.out = c(count, p, PgCSt.sum) #Updated vector of output
#    print(new.out)
    #out = rbind(out, new.out)            #Running summary of all output

    ii = i  #Hold the index
    #Inner Repeat Loop to update top t2-2 rows
    repeat {          
      if ((top.rows[i]-1+j)<mins[i]) { #When moving index hits bottom...
        i=i-1   #change index
        #if (i==0) break
        j=1     #loop flag
        top.rows[i:(t2-2)]=(top.rows[i]-1):(top.rows[i]-t2+i+1) #The next index up, and all those below
      } else { #end if
          if (j==0) {top.rows[i]=top.rows[i]-1; break}
          break } #end else
    } #end repeat
  
    #If statement that breaks if error tolerance is reached
    if (ii>i) {   #When a new sub-cycle begins...
#      print("##### FLAG #####")
#      print(c(ii, i, end.cyc.prob, p))
      if (p-end.cyc.prob < tol.3) {break #Either break    
      } else end.cyc.prob=p              #Or update cycle probability
    } #end if statement

#    print(i)
    #If statement that breaks if the last combination is reached, else resets row index
    if (i==1 && top.rows[i]==mins[i]) break
    else i=t2-2 #Reset i
    j=0        #Reset j, loop flag
    G=top.rows[i]-3 #Update G index

  } #end repeat

  #Calculate the final probability component
  #Note: Since this code is an approximation, it is in many senses best to eliminate
  #this portion of the code because for large k it will be irrelevant.  However, it is
  #here included in order that PdofCSt.cyc can be made to obtain identical results as
  #PdofCSt.T1or2 for all cases, otherwise there is the slight inconsistency.  It will
  #not appreciably increase run time, because it is a single PCS calculation.
  res[1,1]=res[1,1]-1 #Final component, see cyclic.nchoosek code
  A = res      #Set A equal to "res", a portion of true A
  Abar = matrix(NA, nrow=GG, ncol=(g+1))   #Initialize Abar -> needs to be NA's
  for (ii in 1:(g+1)) {                    #Loop that obtains Abar
      bool = (inds %in% A[,ii])            #by selecting elements
      Abar[,ii] = inds[!bool]  #not in A
  } #end for ii
  A = t(A);     Abar = t(Abar)          #Transpose to change form
  A = A+k-GG-T; Abar = Abar+k-GG-T      #Add constant to make correct indices
  R = dim(A)[1] #Initialize number of index sets; This is g+1
  PrCSt = 0     #Initialize probability value
  if (t1==0) {theta.top = NULL  #Critical "if" statement, accounting for special case of t1=0       
  } else theta.top = theta[(k-t1+1):k]
  theta.star = c(lower, theta[Abar[r,]], theta[A[r,]], theta.top) #Concatenate the 4 pieces of vector together; note 2nd vector reads indices bottom to top
  last.component = PofCSt(theta.star,T,m,tol.4)[1]
  p = p + last.component  #Running sum of all probability
  new.out = c(count+1, p, last.component) #Updated vector of output
#  print(new.out)
  return(new.out)

} #end if G==0
} #end if t2==1 || t2==2
} #end function


##### 2 ########################################################
PdofCSt.T1or2 <- function(theta, T, d, m=20, tol=1e-8) {
#Subroutine of PdofCSt.cyc for the special case T=1 or T=2 when G>0.
#This provides the 'exact' solution, up to the quadrature. (This is
#as opposed to T>2 cases, where cases are excluded because of small
#probability contributions.)
#Calculate Pd(CSt) using the subroutine P(CSt)
#Note: This version is identical to my original, except that #require(utils) with "combn"
#					replaces outdated    #require(vsn) with "nchoosek"
# theta is the vector of means
# T is the number of populations desired for selection
# d is the distance tolerance parameter, the innovative part
#E.g. PdofCSt.T1or2(theta, T=2, d=0.3, m=20, tol=1e-8)

#require(utils)                 #Library for combn
#require(statmod)	       #Library for gauss.quad.prob
y = gauss.quad.prob(m, dist="normal")$nodes #Load the constants for GH quadrature
k = length(theta)
cut1 = theta[k-T+1]-d-0.00000001  #Note: the 0.00000001 is a fudge factor to circumvent R's bug of theta[123]=1>=1=cut returning FALSE
cut2 = theta[k-T+1]+d
t1 = sum(ifelse(cut2<theta, 1, 0))                   #Count of fixed values
t3 = sum(ifelse((cut1<=theta)&(theta<=cut2), 1, 0))  #Count of indifferent values.
t2 = T-t1
t3 = t3-t2
#Gd = sum(ifelse(theta>=cut, 1, 0)) - T #Count of values >= cut.  Note: >= is chosen so that d=0 reduces to the usual P(CS1) case
#print(c(cut1,cut2,t1,t2,t3))

if (t3 > 0) {                   #Statement should not run if t3 = 0; P(CSt) = p, d=0 will be the correct result
    A = combn((t3+t2),t2)    #Obtain relative upper indices
    A = t(A)                    #transpose to change form
    A = A+k-(t1+t2+t3)          #Add constant to make correct indices
    Abar = combn((t3+t2),t3) #Obtain relative middle indices
    Abar = t(Abar)              #transpose to change form
    Abar = Abar+k-(t1+t2+t3)    #Add constant to make correct indices

    if ((k-(t1+t2+t3)) == 0) {lower = NULL
    } else lower = theta[1:(k-(t1+t2+t3))]   #Vector of fixed lower theta's; "If" accounts for the case where all pops. are chosen
    R = dim(A)[1]		#Initialize number of index sets
    PdCSt = 0			#Initialize probability value
    for (r in 1:R) {            #The summation in the PCS formula!
      if (t1==0) {theta.top = NULL  #Critical "if" statement, accounting for special case of t1=0       
      } else theta.top = theta[(k-t1+1):k]
      theta.star = c(lower, theta[Abar[R-r+1,]], theta[A[r,]], theta.top) #Concatenate the 4 pieces of vector together; note 2nd vector reads indices bottom to top
      PdCSt[r] = PofCSt(theta.star,T,m,tol)[1]
#      print(c(theta.star[(k-T+1):k], sum(theta.star[(k-T+1):k]), round(PdCSt[r], 4)))
    } #end for r
} else {PdCSt = PofCSt(theta,T,m,tol); PdCSt=PdCSt[2:length(PdCSt)]} #end if

return(sum(PdCSt))
} #end function


##### 3 ########################################################
PofCSGt <- function(theta, T, Gd, m=20, tol=1e-8) {
#POFCSGT calculates the probability that we correctly select the top
#population t elements inside of the top t+G sample elements
#The code is an exact version.  It was developed after PdofCSt.cyc,
#when I saw that this concept may be more useful than the PdofCSt
#concept.  It is expected that the bootstrap which I already created
#will match these results and therefore there is no need to create
#a theoretical version which will obtain results for large k
#Judicious use is made of the core subroutine P(CSt)
# theta is the vector of means
# T is the number of populations desired for selection
# Gd is the number of additional populations admitted
#E.g. PofCSGt(theta, T=2, Gd=2, m=20, tol=1e-6)

#require(utils)                 #Library for combn
#require(statmod)	       #Library for gauss.quad.prob
y = gauss.quad.prob(m, dist="normal")$nodes #Load the constants for GH quadrature
k = length(theta)

if (Gd > 0) {               #Statement should not run if Gd = 0; P(CSt) = p, d=0 will be the correct result
    A = combn((k-T),Gd)  #Obtain relative upper indices
    A = t(A)                #transpose to change form
    R = dim(A)[1]		#Initialize number of index sets

    if ((k-Gd-T) == 0) {Abar = NULL
    } else {
      #Construct Abar matrix
      Abar = matrix(NA, nrow=R, ncol=(k-T-Gd))   #Initialize Abar -> needs to be NA's
      inds = 1:(k-T)                       #Indices, move out of loop in revision
      for (ii in 1:R) {                    #Loop that obtains Abar
        bool = (inds %in% A[ii,])            #by selecting elements
        Abar[ii,] = inds[!bool]  #not in A
      } #end for ii
    } #end if/else

    A.top = rep((k-T+1):k,times=R)
    A.top = matrix(A.top,R,length((k-T+1):k),byrow=TRUE)
    A = cbind(A,A.top)


    PdCSt = 0			#Initialize probability value
    for (r in 1:R) {
      theta.star = c(theta[Abar[R-r+1,]], theta[A[R-r+1,]]) #Note we take them from bottom first
      PdCSt[r] = PofCSt(theta.star,(T+Gd),m,tol)[1]
    } #end for r
} else {PdCSt = PofCSt(theta,T,m,tol); PdCSt=PdCSt[2:length(PdCSt)]} #end if

sum(PdCSt)
} #end function
##############################################################

##### 4 #########################################################
PofCSt <- function(theta, T, m, tol=1e-7) {
# The function calculates P(CSt) using equation (2.4) of Gupta 1998
# y are the constants for the Gaussian-Hermite quadrature, to be applied afterwards
# theta is the vector of means
# T is the number of populations desired to be chosen
#Note: For future error dialoguing I'll need to set the case T=k
#giving PCSt=1.0, and all wrong calls of T to error.  Presently
#the machine error for this is, "if (Qti[k-T+1-i)] <tol) break:
#                                missing value where TRUE/FALSE needed"
# E.g. PofCSt(theta, T, m)
#################################################################
#require(statmod)       #Library for gauss.quad.prob

y = gauss.quad.prob(m, dist="normal")$nodes #Load the constants for GH quadrature
k = length(theta)
len.y = length(y)
#tol = 1e-10
Qti = NA
phi.y = NA

for (i in (k-T):1) {           #Loop over the i's in Qti
 
      y.theta = y + theta[i]   #Vector of constants required for run, allowing simultaneous computation at each quadrature point

      #Begin first CDF product
      if ((i-1) >= 1) {        #if statement required for case i=0
        h = seq(1, (i-1))
        len.h = length(h)
        y.theta.h = rep(y.theta, each=len.h)
        norm.mat = pnorm(y.theta.h-theta[h])
	norm.mat = matrix(norm.mat,len.y, len.h, byrow=TRUE) #Rows are quadrature points, columns are pnorm values for indices h 
        G1 = apply(norm.mat, 1, prod)
      } else {G1 = 1}          #end if

      #Begin second CDF product
      if ((i+1)<=(k-T)) {      #if statement required for case i=k-T
        l = seq((i+1), (k-T))
        len.l = length(l)
        y.theta.l = rep(y.theta, each=len.l)
        norm.mat = pnorm(y.theta.l-theta[l])
	norm.mat = matrix(norm.mat, len.y, len.l, byrow=TRUE) #Rows are quadrature points, columns are pnorm values for indices l
        G2 = apply(norm.mat, 1, prod)
      } else {G2 = 1}

      #Begin third CDF product
      j = seq((k-T+1), k)
      len.j = length(j)
      y.theta.j = rep(y.theta, each=len.j)
      norm.mat = 1-pnorm(y.theta.j-theta[j])
      norm.mat = matrix(norm.mat, len.y, len.j, byrow=TRUE) #Rows are quadrature points, columns are pnorm values for indices j
      G3 = apply(norm.mat, 1, prod)

      G123 = cbind(G1,G2,G3)
      phi.y = apply(G123, 1, prod)     #Vector of the complete product values all quadrature points
    Qti[(k-T+1-i)] = G.H.Quad(phi.y, m)  #Perform numerical quadrature for case i
    if (Qti[(k-T+1-i)] < tol) break      #Run-time reduction based on postulate of monotonicity
    #print(c(round(i), round(Qti[k-T+1-i]),2))
} #end for i
PCSt = sum(Qti)
c(PCSt, Qti)
} #end


##### 5 ########################################################
G.H.Quad <- function(x, m) {
# Performs Gauss-Hermite Quadrature for general m
#require(statmod)
weight = gauss.quad.prob(m, dist="normal")$weights
res = sum(x*weight)
res
} #end
#############################################################


##### From sim_Article2_functions.txt
##### 6 ##############################################################
PdofCSGt.bootstrap5 <- function(theta, T, D, G, B, SDE, dist=c("normal","t"), df=14, trunc=6, est.names=c("O")) {
#PDOFCSGT.BOOTSTRAP5 is revised from PdofCSGt.bootsrap3 (skipped 4 to be even with PofCSLt).
#It includes assorted minor adjustments to prepare it for PCS package publication.  One significant change
#is that the user supply of n and var was changed to SDE and df, in order to be more general.  
#For use with PdofCSGt.simulation3.
#theta = vector of true means, need not be in order
#T = vector of T values desired in table
#D = vector of D values desired in table 
#G = vector of G values desired in table 
#B = number of bootstrap samples
#SDE = standard error of the statistics theta, row-wise
#      e.g. if theta is a vector of means, then the standard error of one of the means is sqrt(var/n)
#dist = The distributional assumption used for estimating PCS
#       "normal" uses a N(theta, var) dist (common or uncommon var)
#       "tcommon" uses a t-dist shifted by theta, with a common df by the Satterthwaite approximation
#       "tuncommon" uses a t-dist shifted by theta, with uncommon df estimated by MOM
#df = the common degrees of freedom for one of t-statistics in the vector theta;
#	the parameter is used only if dist="t"
#trunc = number of standard deviations below minimum selected population to disregard
#        in the estimation of PCS; it is a truncation parameter to increase run time
#        6 seems to give pretty good accuracy; below that reduces accuracy
#est.names = Vector of estimators desired.  O is Olkin, M is McCulloch & Dechter
#    B is Bofinger, and E is edwards (MIN).  They must be entered in the order
#    given, i.e. if you want "E", then the others must be included, in this
#    present version of the code.
#Example: theta = c(4.2, 2.5, 2.3, 1.7, 1.5, 1.0); 
#         PdofCSGt.bootstrap5(theta, c(1,2,3), c(0,2), c(0,2), 10, 1)
####################################################################
beg = Sys.time()
 
#Initilize parameters
T = sort(T); D = sort(D); G = sort(G)
T.len = length(T); D.len=length(D); G.len=length(G)
name.len = length(est.names)
resD = array(NA, c(D.len,T.len,name.len), list(D,T,est.names)) #Indicator matrix for D by T by 4 estimators
resG = array(NA, c(G.len,T.len,name.len), list(G,T,est.names))
num.corrD = array(0, c(D.len,T.len,name.len), list(D,T,est.names)) #Count matrix for D by T by 4 estimators
num.corrG = array(0, c(G.len,T.len,name.len), list(G,T,est.names))

#First part of truncation step (significantly improve run time)
bound = T[T.len]      #Smallest T which will ever be selected
const = 1/SDE	      #Constant modifier of the mean, user supplied, needed when statistics theta are not standardized
theta = sort(theta)
K = length(theta)     #True length of theta

##Obtain different estimator terms
#g.mean = mean(theta)  
#SS = sum((theta-g.mean)^2)
#a.MD = max(0, (1 - ((K-3)*varE/n)/SS));  #M, McCulloch & Dechter (1985)
#p = pchisq(SS*n/varE, (K-1));
#Z.B = qchisq(p/2, (K-1));
#V = max((SS/Z.B - varE/n), 0.01);
#a.B = const*sqrt(1/(1/V + n/varE)) 	#B, Bofinger (1990)
#a.MIN = max(0, (1 - ((K-3)*(K-1)*(varE/n)/(4*K*(theta[K]-g.mean)^2)))); #Edwards (1992)

#Bootstrap loop
for (b in 1:B) {
  #Generate deviates according to distribution "dist"
  if        (dist[1]=="normal")    {xbar=rnorm(K,0)
  } else if (dist[1]=="t")   	    xbar=rt(K,df)

  #Loop through calculation for each estimator, O, M, B, & E
  for (est in 1:name.len) {
    if        (est.names[est]=="O") {Theta=const*theta       #O, Olkin, Sobel, & Tong (1982)
    } #else if (est.names[est]=="M") {Theta=const*a.MD*theta  #M, McCulloch & Dechter (1985)
#    } else if (est.names[est]=="B") {Theta=const*a.B*theta   #B, Bofinger (1990)
#    } else if (est.names[est]=="E") {Theta=const*a.MIN*theta #E, Edwards (1992)
#    } #end if
    x=xbar+Theta  #Obtain deviates according to estimator

    #Second part of truncation step
    x = x[x>(x[K-bound+1]-trunc)] #Cut off parameters below interest
    k = length(x)       #Find shorter length for truncated vector
    Th.i = rank(Theta[(K-k+1):K])    #Vector of 'true' ranks of truncated theta

    X.i = order(x)
    for (j in 1:T.len)   {
 
      #d-parameter for-loop
      for (i in 1:D.len) {         
        t = T[j]; d = D[i]
        cut1 = Theta[K-t+1]-d      #lower bound
        cut2 = Theta[K-t+1]+d      #upper bound
        t1 = sum(ifelse(cut2<Theta, 1, 0))                   #Count of fixed values
        t3 = sum(ifelse((cut1<=Theta)&(Theta<=cut2), 1, 0))  #Count of indifferent values.
        t2 = t-t1
        t3 = t3-t2
        if ((t1+t2+t3)>k) {(resD[i:D.len,j,est]=1)&break #Condition accounts for when truncation is "large"
        } else {        
          #1st part of definition, if all A1 are in s, then remove A1 from s, else next
          x.i  = X.i[(k-t1-t2+1):k]  #Sub-vector of sample ranks, all
          if (t1==0) {th.i=NULL
            } else {
            th.i = Th.i[(k-t1+1):k]  #Sub-vector of true ranks, fixed set A1 (Note: "max" fixes the case t1=0)
            dCSt = match(th.i, x.i)    #Indices of matches (Note: I couldn't get "nomatch=NULL" to work, and so do it this strange way)
            dCSt = dCSt[!is.na(dCSt)]  #Remove the NA's
            if (length(dCSt) == t1) {x.i=x.i[-dCSt]  
            } else (resD[i,j,est]=0)&next
          } #end if t1==0
          #2nd part of definition, if all s~A1 are in A2, then selection is correct
          th.i = Th.i[(k-t1-t2-t3+1):(k-t1)] #Sub-vector of true ranks, set A2
          dCSt = x.i %in% th.i       #Logical vector indicating whether selections were correct
          if (sum(dCSt) == t2) {(resD[i:D.len,j,est]=1)&break  #Flag if all correct are in selected
          } else resD[i,j,est]=0
        } #end if (t1+t2+t3)>k
      } #end for i

      #G-parameter for-loop
      for (i in 1:G.len) {  
        t = T[j]; g = G[i]  #Fix variables for easy reading
        if ((t+g)>k) {(resG[i:G.len,j,est]=1)&break}  #Jump to next column when bottom row is reached, or all populations selected
        th.i = Th.i[(k-t+1):k]   #Sub-vector of true ranks
        x.i  = X.i[(k-t-g+1):k]  #Sub-vector of sample ranks
        GCSt = th.i %in% x.i     #Logical vector indicating whether selections were correct
        if (sum(GCSt) == t) {(resG[i:G.len,j,est]=1)&break  #Flag if all correct are in selected
        } else resG[i,j,est]=0
      } #end for i

    } #end for j

    num.corrD[,,est] = num.corrD[,,est] + resD[,,est] #Running sums
    num.corrG[,,est] = num.corrG[,,est] + resG[,,est]

  } #end for est

} #end for B

est.PdCSt = num.corrD/B
est.PCSGt = num.corrG/B

end = Sys.time()
print(end-beg)

out = list(est.PdCSt,est.PCSGt)
names(out)=c("d", "G")
return(out)
} #end function


################################################################
##### From s.Article3_selection02.txt
##### 7 ###########################################################
PofCSLt.bootstrap5 <- function(theta, T, L, B, SDE, dist=c("normal", "t"), df=14, trunc=6, est.names=c("O")) {
#POFCSLT.BOOTSTRAP5 is built off of PdofCSGt.bootstrap3, but it employs the FDR-type
#selection approach, which is Mahamunulu goal I (p. 1080) with c=L, s=t.  It is
#used by package function "PCS.boot" 
#PofCSLt.bootstrap5 reflects parallel minor updating with PdofCSGt.bootstrap5, for PCS package publication
#theta = vector of true means, need not be in order
#T = vector of T values desired in table, where T = total number selected
#L = vector of minimum number of populations selected correctly
#B = number of bootstrap samples
#SDE = standard error of the statistics theta, row-wise
#      e.g. if theta is a vector of means, then the standard error of one of the means is sqrt(var/n)
#dist = The distributional assumption used for estimating PCS
#       "normal" uses a N(theta, var) dist (common or uncommon var)
#       "tcommon" uses a t-dist shifted by theta, with a common df by the Satterthwaite approximation
#       "tuncommon" uses a t-dist shifted by theta, with uncommon df estimated by MOM
#df = the common degrees of freedom for one of t-statistics in the vector theta;
#	the parameter is used only if dist="t"
#trunc = number of standard deviations below minimum selected population to disregard
#        in the estimation of PCS; it is a truncation parameter to increase run time
#        6 seems to give pretty good accuracy; below that reduces accuracy
#est.names = Vector of estimators desired.  O is Olkin, M is McCulloch & Dechter
#    B is Bofinger, and E is edwards (MIN).  They must be entered in the order
#    given, i.e. if you want "E", then the others must be included, in this
#    present version of the code.
#Example: theta = c(4.2, 2.5, 2.3, 1.7, 1.5, 1.0); 
#         PofCSGLt.bootstrap5(theta, c(1,2,3), c(0,2), 10)
####################################################################
beg = Sys.time()
 
#Initilize parameters
T = sort(T); L = sort(L)
T.len = length(T); L.len=length(L)
name.len = length(est.names)
resL = array(NA, c(L.len,T.len,name.len), list(L,T,est.names)) #Indicator matrix for L by T by 4 estimators
num.corrL = array(0, c(L.len,T.len,name.len), list(L,T,est.names)) #Count matrix for L by T by 4 estimators

#First part of truncation step (significantly improve run time)
bound = T[T.len]      #Smallest T which will ever be selected
const = 1/SDE	      #Constant modifier of the mean, user supplied, needed when statistics theta are not standardized
theta = sort(theta)
K = length(theta)     #True length of theta

##Obtain different estimator terms
#g.mean = mean(theta)  
#SS = sum((theta-g.mean)^2)
#a.MD = max(0, (1 - ((K-3)*varE/n)/SS));  #M, McCulloch & Dechter (1985)
#p = pchisq(SS*n/varE, (K-1));
#Z.B = qchisq(p/2, (K-1));
#V = max((SS/Z.B - varE/n), 0.01);
#a.B = const*sqrt(1/(1/V + n/varE)) #B, Bofinger (1990)
#a.MIN = max(0, (1 - ((K-3)*(K-1)*(varE/n)/(4*K*(theta[K]-g.mean)^2)))); #Edwards (1992)

#Bootstrap loop
for (b in 1:B) {
  #Generate deviates according to distribution "dist"
  if        (dist[1]=="normal")    {xbar=rnorm(K,0)
  } else if (dist[1]=="t")   	 xbar=rt(K,df)

  #Loop through calculation for each estimator, O, M, B, & E
  for (est in 1:name.len) {
    if        (est.names[est]=="O") {Theta=const*theta       #O, Olkin, Sobel, & Tong (1982)
    } #else if (est.names[est]=="M") {Theta=const*a.MD*theta  #M, McCulloch & Dechter (1985)
#    } else if (est.names[est]=="B") {Theta=const*a.B*theta   #B, Bofinger (1990)
#    } else if (est.names[est]=="E") {Theta=const*a.MIN*theta #E, Edwards (1992)
#    } #end if
    x=xbar+Theta  #Obtain deviates according to estimator

    #Second part of truncation step
    x = x[x>(x[K-bound+1]-trunc)] #Cut off parameters below interest
    k = length(x)       #Find shorter length for truncated vector
    Th.i = 1:k          #Vector of 'true' ranks of truncated theta

    X.i = order(x)
    for (i in 1:L.len) {    
 
      #L-parameter for-loop
        for (j in 1:T.len)   {
        t = T[j]; g = L[i]  #Fix variables for easy reading
        if (g>t) {(resL[i:L.len,j,est]=NA)&next}  #Jump to next column when bottom row is reached, or all populations selected
        th.i = Th.i[(k-t+1):k]   #Sub-vector of true ranks
        x.i  = X.i[(k-t+1):k]  #Sub-vector of sample ranks
        LCSt = th.i %in% x.i     #Logical vector indicating whether selections were correct
        if (sum(LCSt) >= g) {(resL[i,j:T.len,est]=1)&break  #Flag if all correct are in selected
        } else resL[i,j,est]=0
      } #end for i

    } #end for j

    num.corrL[,,est] = num.corrL[,,est] + resL[,,est] #Running Sums

  } #end for est

} #end for B

est.PCSLt = num.corrL/B

end = Sys.time()
print(end-beg)

return(est.PCSLt)
} #end function

##### 8 #############################################################
PdCSGt.bootstrap.NP2 <- function(X1, X2, T, D, G, N, trunc=6) {
#PDCSGT.BOOTSTRAP.NP2 is my fully non-parametric bootstrap for a control-
#treatment experiment; finds d-best and G-best tables; differs from
#NP1 code only in that the truncation parameter has been successfully added. 
#X1, X2:  the original data
#T:       the indices of interest
#D:       vector of d-parameters
#G:       vector of G-parameters
#N:       bootstrap sample size
#trunc:   truncation parameter, as in the parametric code
#theta: the true, or supposed true, mean vector, ordered small to large
#Example: PdCSGt.bootstrap.NP(X1, X2, T, D, G, 10)
##############################################################
beg = Sys.time()

#Obtain the means to estimate the true ordering
K = dim(X1)[1]; r1 = dim(X1)[2]; r2 = dim(X2)[2]
theta =  abs(tindep(X1, X2)[,1])   #True t-statistics
ord.i = order(theta)                #Indices of ordered t-stats
theta = theta[ord.i]                #Order theta
X1 = X1[ord.i, ]; X2 = X2[ord.i, ]; #Order data
#Th.i = rank(theta)                  #Vector of true ranks

#Initialize variables
T = sort(T); D = sort(D); G = sort(G)
T.len = length(T); G.len=length(G); D.len=length(D)
resD = matrix(NA, D.len, T.len)
resG = matrix(NA, G.len, T.len)
num.corrD = matrix(0, D.len, T.len)
num.corrG = matrix(0, G.len, T.len)
colnames(num.corrD)=T; rownames(num.corrD)=D
colnames(num.corrG)=T; rownames(num.corrG)=G

#First part of truncation step (significantly improve run time)
bound = T[T.len]      #Largest T which will ever be selected

#Bootstrap
for (n in 1:N) {

  #Obtain non-parametric realization
  X.C = apply(X1, 1, sample, size=r1, replace=TRUE)
  X.C = t(X.C)
  X.T1 = apply(X2, 1, sample, size=r2, replace=TRUE)
  X.T1 = t(X.T1)
  sample = cbind(X.C, X.T1)
  tstat = abs(tindep(X.C, X.T1)[,1]) #Sample t-statistics

    #Second part of truncation step
    tstat = tstat[tstat>(tstat[K-bound+1]-trunc)] #Cut off parameters below interest
    k = length(tstat)       #Find shorter length for truncated vector
    Th.i = rank(theta[(K-k+1):K])    #Vector of 'true' ranks of truncated theta

  X.i = order(tstat)    #Sample t-stat ranks, the realization (deviate)

  #Construct d-best and G-best table based on deviate
  for (j in 1:T.len)   {

    for (i in 1:D.len) {         #d-parameter for-loop
      t = T[j]; d = D[i]
      cut1 = theta[K-t+1]-d      #lower bound
      cut2 = theta[K-t+1]+d      #upper bound
      t1 = sum(ifelse(cut2<theta, 1, 0))                   #Count of fixed values
      t3 = sum(ifelse((cut1<=theta)&(theta<=cut2), 1, 0))  #Count of indifferent values.
      t2 = t-t1
      t3 = t3-t2
      if ((t1+t2+t3)>k) {(resD[i:D.len,j]=1)&break}  #Jump to next column when row of certainty is reached
      #1st part of definition, if all A1 are in s, then remove A1 from s, else next
      x.i  = X.i[(k-t1-t2+1):k]  #Sub-vector of sample ranks, all
      if (t1==0) {th.i=NULL
        } else {
        th.i = Th.i[(k-t1+1):k]  #Sub-vector of true ranks, fixed set A1 (Note: "max" fixes the case t1=0)
        dCSt = match(th.i, x.i)    #Indices of matches (Note: I couldn't get "nomatch=NULL" to work, and so do it this strange way)
        dCSt = dCSt[!is.na(dCSt)]  #Remove the NA's
        if (length(dCSt) == t1) {x.i=x.i[-dCSt]  
        } else (resD[i,j]=0)&next
      } #end if t1==0
      #2nd part of definition, if all s~A1 are in A2, then selection is correct
      th.i = Th.i[(k-t1-t2-t3+1):(k-t1)] #Sub-vector of true ranks, set A2
      dCSt = x.i %in% th.i       #Logical vector indicating whether selections were correct
      if (sum(dCSt) == t2) {(resD[i:D.len,j]=1)&break  #Flag if all correct are in selected & move to next row
      } else resD[i,j]=0
    } #end for i

    for (i in 1:G.len) {  #G-parameter loop
      t = T[j]; g = G[i]  #Fix variables for easy reading
      if ((t+g)>k) {(resG[i:G.len,j]=1)&break}  #Select remaining and jump to next column when row of certainty is reached
      th.i = Th.i[(k-t+1):k]   #Sub-vector of true ranks
      x.i  = X.i[(k-t-g+1):k]  #Sub-vector of sample ranks
      GCSt = th.i %in% x.i     #Logical vector indicating whether selections were correct
      if (sum(GCSt) == t) {(resG[i:G.len,j]=1)&break  #Flag if all correct are in selected & move to next row
      } else resG[i,j]=0
    } #end for i
  } #end for j

  num.corrD = num.corrD + resD #Running sums
  num.corrG = num.corrG + resG
  #print(cbind(theta,tstat,Th.i,X.i))
} #end for N

est.PdCSt = num.corrD/N
est.PCSGt = num.corrG/N
out = list(est.PdCSt,est.PCSGt)
names(out) = c("d","G")

end = Sys.time()
print(end-beg)

return(out)
} #end function

##### 9 ##############################################################################################
##### Make function "tindep"
tindep <- function(X, Y, flag=0) {
#tindep performs a Welch independent t-test on two matrices of data (unequal unknown variances)
#	X is the first kxnx matrix of n1 samples from k populations
#	Y is the second kxny matrix of n2 samples
#	flag is defaulted to not return the multiple comparison pvalues.
#	  if it is changed, then the library "multtest" is invoked and
#	  various step up and step down procedures run
#	output is a matrix of T-values and Pvalues (labeled)
#	Examples: ans=tindep(X,Y, flag=1)
#		  ans=tindep(golub[,1:27], golub[,28:38], flag=1)
##############################
k = dim(X)[1];   		 #number of populations, should be same for both sets
nx = dim(X)[2]; ny = dim(Y)[2]   #number in sample, may vary between sets
Xbar = apply(X, 1, mean) 	 #computes the mean of each row of Dif
Xse = apply(X, 1, var)/nx   	 #computes the sample SE of each row of X
Ybar = apply(Y, 1, mean) 	 #computes the mean of each row of Y
Yse = apply(Y, 1, var)/ny   	 #computes the sample SE of each row of Y
T = (Xbar-Ybar)/sqrt(Xse+Yse)	 #computes t-statistic
df = (Xse + Yse)^2 / (Xse^2/(nx-1) + Yse^2/(ny-1))
Pvalue = 2*(1-pt(abs(T), df= df))  #computes p-value for t-statstic

### Various correction factors for multiple testing, using MULTTEST library
if (flag != 0) {
	#require(multtest)  #load the multtest library
	procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")  #this sequence is taken directly from multtest documentation, p. 5
	res <- mt.rawp2adjp(Pvalue, procs)
	adjp <- res$adjp[order(res$index), ]
	#round(adjp[1:10, ], 2)  #display line
	return(cbind(T,adjp))
	} #end if
else return(cbind(T,Pvalue))
} #end function

#######################################################################################################
##### 10 ###############################################################################################
PCS.exact <- function(theta, t=1, g=NULL, d=NULL, m=20, tol=1e-8) {

### Error dialoguing
#t
t = as.integer(t)
if (length(t)>1) stop("t must be a scalar")
if (!is.integer(t) | t<1) stop("t must be a positive integer")

#g & d
if (is.null(g) & is.null(d)) stop("A value is required for g or d")
if (!is.null(g) & !is.null(d)) stop("Only g or d may be inputted at one time, not both")
if (!is.null(g)) {
	if (length(g)>1) stop("g must be a scalar")
	g = as.integer(g)
	if (g < 0) stop("g must be a non-negative integer")
	}
if (!is.null(d)) {
	if (length(d)>1) stop("d must be a scalar")
	if (d < 0) stop("d must be a non-negative real number")
	}

#m
if (m<=0) stop("m (number of nodes in Gaussian quadrature) must be positive")

#tol
if (tol>0.01) warning("tol (tolerance parameter) is large and may admit gross errors")

###Function call
theta = sort(theta)          #Functions assume theta is ordered
if (is.null(d)) {
	if (length(theta) <= t+g) {out=1; warning("t+g >= length(theta), which implies PCS=1, necessarily"); return(out)}
	out = PofCSGt(theta=theta, T=t, Gd=g, m=m, tol=tol)
}  else {
	if (length(theta) < 2) stop("theta must be a vector of length >= 2")
	out = PdofCSt.cyc2(theta=theta, T=t, d=d, m=m, tol=tol)[2]
}

if (out>1) out=1.00 #Fix results near one that go over due to round-off error, e.g. 1.000003
return(out)
} 


##### 11 ###################################################################################
PCS.boot.par <- function(theta, T=1:1, G=NULL, D=NULL, L=NULL, B=100, SDE=1, dist=c("normal","t"), df=14, trunc=6) {

### Parameter errors/warnings
#theta
if (length(theta) < 2) stop("theta must be a vector of length >= 2")

#T
T = as.integer(T)
if (min(T)<1) stop("The elements of T must be positive integers")
if (max(T)>length(theta)) stop("The elements of T must be less than or equal to the length of theta")

#B
if (B<1) stop("B (bootstrap size) must be greater than one")

#SDE
if (SDE<=0) stop("SDE (standard error) must be positive")

#dist
dist=dist[1]
if (dist != "normal" & dist != "t") stop("dist must be either 'normal' or 't'")

#df
if (df<=0) stop("df (degrees of freedom) must be positive")

#trunc
if (trunc<=0) stop("trunc (truncation parameter) must be positive")
if (trunc<=4) warning("trunc<=4 and may therefore produce inaccuracies in estimating PCS")

### The seven possible input commands:
#Error, no G, D, & L
if (is.null(G) & is.null(D) & is.null(L)) stop("A vector is required for G or D or L")

#D only
if (is.null(G) & !is.null(D) & is.null(L)) {
	if (min(D)<0) stop("The elements of D must be non-negative")
	out = PdofCSGt.bootstrap5(theta, T, D, G=0, B, SDE, dist, df, trunc)
	out = out$d[,,1]; 
	if(length(D)==1) out=matrix(out,nrow=1,ncol=length(T),dimnames=list(D,T))  }

#G only
if (!is.null(G) & is.null(D) & is.null(L)) {
	G = as.integer(G); if (!is.integer(G) | min(G)<0) stop("The elements of G must be non-negative integers")
	out = PdofCSGt.bootstrap5(theta, T, D=0, G, B, SDE, dist, df, trunc)
	out = out$G[,,1]
	if(length(G)==1) out=matrix(out,nrow=1,ncol=length(T),dimnames=list(G,T))  }

#L only
if (is.null(G) & is.null(D) & !is.null(L)) {
	L = as.integer(L); if (!is.integer(L) | min(L)<0) stop("The elements of L must be non-negative integers")
	if (max(L) > max(T)) stop("The maximum element of L must be <= the maximum element of T")
	out = PofCSLt.bootstrap5(theta, T, L, B, SDE, dist, df, trunc)
	out = out[,,1]
	if(length(L)==1) out=matrix(out,nrow=1,ncol=length(T),dimnames=list(L,T))  }

#D & G
if (!is.null(G) & !is.null(D) & is.null(L)) {
	if (min(D)<0) stop("The elements of D must be non-negative")
	G = as.integer(G); if (!is.integer(G) | min(G)<0) stop("The elements of G must be non-negative integers")
	out = PdofCSGt.bootstrap5(theta, T, D, G, B, SDE, dist, df, trunc)
	outD = out$d[,,1]; if(length(D)==1) outD=matrix(outD,nrow=1,ncol=length(T),dimnames=list(D,T))
	outG = out$G[,,1]; if(length(G)==1) outG=matrix(outG,nrow=1,ncol=length(T),dimnames=list(G,T))
	out = list(outG, outD); names(out)=c("G", "D")  }

#G & L
if (!is.null(G) & is.null(D) & !is.null(L)) {
	G = as.integer(G); if (!is.integer(G) | min(G)<0) stop("The elements of G must be non-negative integers")
	L = as.integer(L); if (!is.integer(L) | min(L)<0) stop("The elements of L must be non-negative integers")
	if (max(L) > max(T)) stop("The maximum element of L must be <= the maximum element of T")
	outG = PdofCSGt.bootstrap5(theta, T, D=0, G, B, SDE, dist, df, trunc)
	outL = PofCSLt.bootstrap5(theta, T, L, B, SDE, dist, df, trunc)
	outG = outG$G[,,1]; if(length(G)==1) outG=matrix(outG,nrow=1,ncol=length(T),dimnames=list(G,T))
	outL = outL[,,1];   if(length(L)==1) outL=matrix(outL,nrow=1,ncol=length(T),dimnames=list(L,T))
	out = list(outG, outL); names(out)=c("G", "L")  }

#D & L
if (is.null(G) & !is.null(D) & !is.null(L)) {
	if (min(D)<0) stop("The elements of D must be non-negative")
	L = as.integer(L); if (!is.integer(L) | min(L)<0) stop("The elements of L must be non-negative integers")
	if (max(L) > max(T)) stop("The maximum element of L must be <= the maximum element of T")
	outD = PdofCSGt.bootstrap5(theta, T, D, G=0, B, SDE, dist, df, trunc)
	outL = PofCSLt.bootstrap5(theta, T, L, B, SDE, dist, df, trunc)
	outD = outD$d[,,1]; if(length(D)==1) outD=matrix(outD,nrow=1,ncol=length(T),dimnames=list(D,T))
	outL = outL[,,1];  if(length(L)==1) outL=matrix(outL,nrow=1,ncol=length(T),dimnames=list(L,T))
	out = list(outD, outL); names(out)=c("D", "L")  }

#G, D, & L
if (!is.null(G) & !is.null(D) & !is.null(L)) {
	if (min(D)<0) stop("The elements of D must be non-negative")
	G = as.integer(G); if (!is.integer(G) | min(G)<0) stop("The elements of G must be non-negative integers")
	L = as.integer(L); if (!is.integer(L) | min(L)<0) stop("The elements of L must be non-negative integers")
	if (max(L) > max(T)) stop("The maximum element of L must be <= the maximum element of T")
	out = PdofCSGt.bootstrap5(theta, T, D, G, B, SDE, dist, df, trunc)
	outL = PofCSLt.bootstrap5(theta, T, L, B, SDE, dist, df, trunc)
	outD = out$d[,,1]; if(length(D)==1) outD=matrix(outD,nrow=1,ncol=length(T),dimnames=list(D,T))
	outG = out$G[,,1]; if(length(G)==1) outG=matrix(outG,nrow=1,ncol=length(T),dimnames=list(G,T))
	outL = outL[,,1];  if(length(L)==1) outL=matrix(outL,nrow=1,ncol=length(T),dimnames=list(L,T))
	out = list(outG, outD, outL); names(out)=c("G", "D", "L")  }

return(out)
} 


##### 12 ##########################################################################################
PCS.boot.np <- function(X1, X2, T=1, G=1, D=NULL, B, trunc=6) {
### Error dialoguing
#X1, X2
if (!is.matrix(X1) | !is.matrix(X2)) stop("X1 and X2 must be matrices")
if (nrow(X1) != nrow(X2)) stop("The number of rows in X1 and X2 must be the same")
if (sum(is.na(X1)) | sum(is.na(X2)) >= 1) stop("The entries of X1 and X2 cannot be 'NA'")

#T
T = as.integer(T)
if (min(T)<1) stop("The elements of T must be positive integers")

#N
N=B
if (N<1) stop("N (bootstrap size) must be greater than one")

#trunc
if (trunc<=0) stop("truncation parameter (trunc) must be positive")
if (trunc<=4) warning("trunc<=4 and may therefore produce inaccuracies in estimating PCS")

### Code proper
#Error, no G, D
if (is.null(G) & is.null(D)) stop("A vector is required for G or D")

#D only
if (is.null(G) & !is.null(D)) {
	if (min(D)<0) stop("The elements of D must be non-negative")
	out = PdCSGt.bootstrap.NP2(X1, X2, T, D, G=0, N, trunc=6)
	out = out$d  }

#G only
if (!is.null(G) & is.null(D)) {
	G = as.integer(G); if (!is.integer(G) | min(G)<0) stop("The elements of G must be non-negative integers")
	out = PdCSGt.bootstrap.NP2(X1, X2, T, D=0, G, N, trunc=6)
	out = out$G  }

#D & G
if (!is.null(G) & !is.null(D)) {
	if (min(D)<0) stop("The elements of D must be non-negative")
	G = as.integer(G); if (!is.integer(G) | min(G)<0) stop("The elements of G must be non-negative integers")
	out = PdCSGt.bootstrap.NP2(X1, X2, T, D, G, N, trunc=6)
	out.d = out$d; out.G = out$G; out = list(out.G,out.d); names(out) = c("G","d")  }

return(out)
} #end function
