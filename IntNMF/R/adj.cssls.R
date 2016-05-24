adj.cssls <-
function(CtC, CtA, Pset=NULL){
    # This function was obtained from NMF package (Gaujoux R., BMC Bioinformatics 2010,11:367) at
# https://github.com/renozao/NMF/blob/master/R/algorithms-snmf.R
# and was adjusted for some minor bug fixes.
# The origibal matlab code was proposed by M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
# Solve the set of equations CtA = CtC*K for the variables in set Pset
# using the fast combinatorial approach
K = matrix(0, nrow(CtA), ncol(CtA));
if ( is.null(Pset) || length(Pset)==0 || all(Pset) ){
#K = solve(CtC) %*% CtA;      
 K = ginv(CtC) %*% CtA;  
 #K = pseudoinverse(CtC) %*% CtA;
#K = pinv(CtC)*CtA;
}else{
lVar = nrow(Pset); pRHS = ncol(Pset);
codedPset = as.numeric(2.^(seq(lVar-1,0,-1)) %*% Pset);
sortedPset = sort(codedPset)
sortedEset = order(codedPset)
breaks = diff(sortedPset);
breakIdx = c(0, which(breaks > 0 ), pRHS);
for( k in seq(1,length(breakIdx)-1) ){
cols2solve = sortedEset[ seq(breakIdx[k]+1, breakIdx[k+1])];
vars = Pset[,sortedEset[breakIdx[k]+1]];
#K[vars,cols2solve] = solve(CtC[vars,vars, drop=FALSE]) %*% CtA[vars,cols2solve, drop=FALSE];   
K[vars,cols2solve] = ginv(CtC[vars,vars]) %*% CtA[vars,cols2solve];       
#K[vars,cols2solve] = pseudoinverse(CtC[vars,vars]) %*% CtA[vars,cols2solve];
}
}
# return K
K
}
