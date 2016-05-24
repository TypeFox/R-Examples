clust.kappa <-
function(clus1,clus2,p) {
   
   com=seq(1,p,by=1)
   n11=length(intersect(clus1,clus2))
   s1com=setdiff(com,clus1)
   s2com=setdiff(com,clus2)
   n12=length(intersect(clus1,s2com))
   n21=length(intersect(s1com,clus2))
   n22=length(intersect(s1com,s2com))

   pr_a=(n11+n22)/p

   pr_e=(n11+n12)*(n11+n21)/p^2+(n12+n22)*(n21+n22)/p^2


   if(length(clus1)==0||length(clus2)==0||length(clus1)+length(clus2)==2*p) ka=-1 else
   ka=(pr_a-pr_e)/(1-pr_e)


   
   return(ka)
}
