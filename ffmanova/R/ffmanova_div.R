### $Id: ffmanova_div.R 53 2007-04-20 12:05:00Z bhm $


matlabColon <- function(from, to) { if(from > to) numeric(0) else from:to }
# Author:: Bjoern-Helge Mevik


norm = function(X){
   svd(X, nu = 0, nv = 0)$d[1]
}


fixModelMatrix = function(mOld) {
#print(mOld)
mOld[mOld>0]=1 # Original m has "2" when the variable should be coded via
               # dummy variables
mNew = mOld
varNamesOld = attr(mOld,'dimnames')[[1]]
#varNamesNew = varNamesOld
nVar = length(varNamesOld)
nPower = rep(0,nVar)
index = 1:nVar
for(i in 1:nVar){
   lab = varNamesOld[i]
   lab = unlist(strsplit(lab, "I(", fixed = TRUE))
   if(length(lab)==2){
         lab = lab[2]
         lab = unlist(strsplit(lab, ")", fixed = TRUE))
         labNew = lab
         if(length(lab)==1){
               lab = unlist(strsplit(lab, "^", fixed = TRUE))
               if(length(lab)==1)
                  lab = c(lab,"1")
               if(length(lab)==2){
                     if(lab[2] == sprintf('%d',(as.integer(lab[2])))){ # string-integer-test
                        for(j in 1:nVar){
                           if(lab[1] == varNamesOld[j]){
                                 nPower[i] = as.integer(lab[2])
                                 index[i] = j
                                 #varNamesNew[i] = labNew # Denne er ikke i bruk
                             }
                          }
                       }
               }
          }
    }
 }
for(i in 1:nVar){
       if(index[i]!=i){
          a = rep(0,nVar)
          a[index[i]] = nPower[i]
          a = matrix(a,nVar,1)
          b = mOld[i,,drop=FALSE]
          mNew = mNew + a %*%b
        }
 }
 mNew[index==(1:nVar),,drop=FALSE]
}# end  fixModelMatrix



myrank = function(X,tol_ = 1e-9){ # Ny funksjon ikke matlab
if(dim(X)[2]==0)                  # Kode hentet fra myorth
   return(X)
S = svd(X,nv=0,nu=0)
S = S$d
r = length(S)
tol = max(dim(X)) * S[1] * tol_
if(!tol)
   return(0)
while(S[r] < tol)
   r=r-1
r
}# end myrank
