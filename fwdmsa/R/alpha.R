alpha=function(X){
# Cronbach's alpha
 J=ncol(X);
 varX=var(X)
 if(sum(varX) > 0 & length(varX) > 1) J/(J-1)* sum(varX - diag(diag(varX)))/sum(varX) else NA
}
