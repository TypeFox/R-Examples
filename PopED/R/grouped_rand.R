grouped_rand <- function(G,xtopt,dxt,ff1,axt){
# Calculate a new random search xt in a grouped xt manner way.
xt = xtopt
for(k in 1:max(max(G))){
    random_num = randn(1) 
    xt = xt+(dxt/ff1*random_num[,,drop=T]*(axt>0))*(G==k) 
 }
return( xt) 
}
# 
# function xt=grouped_rand(G,xtopt,dxt,ff1,axt)
# % Calculate a new random search xt in a grouped xt manner way.
# size1 = size(xtopt,1);
# size2 = size(xtopt,2);
# xt = zeros(size1,size2);
# for k=1:max(max(G))
# random_num = randn; 
# tmp = ones(size1,size2)*k;
# inters = (G==tmp);
# xt = (xtopt+dxt./ff1.*random_num.*(axt>0)).*(inters==1);
# 
# %     for i=1:size1
# %        for j=1:size2
# %           if (inters(i,j)==1)
#   %              xt(i,j) = xtopt(i,j)+dxt(i,j)/ff1*random_num*(axt(i,j)>0);
# %           end
# %        end
# %     end
# end
# end