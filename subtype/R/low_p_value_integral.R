low_p_value_integral <-
function(pt1) {
#measures enrichment in low p-values (<lambda) for a given distribution
lambda = 0.2;
ratio = length(which(pt1 < lambda)) / (length(pt1) * lambda);
areauc = 1-mean(pt1); 
return(list=c(ratio=round(ratio, digits=3), areauc=areauc));
}
