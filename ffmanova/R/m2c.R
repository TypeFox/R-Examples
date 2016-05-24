### $Id: m2c.R 45 2006-08-15 13:11:29Z bhm $
# %=============== m2c.m ====================
# % C = m2c(M,df)
# %    partition M into a cell array C
# %
# %   Input:
# %           M(*,*)  -  ordinary unpartitioned matrix
# %           df(1,*) -  number of columns in the cells of C
# %                      default: [1 1 1 ....]
# %   Output:
# %           C{1,*}  -  M partitioned into a cell array
# %
# %   See also: c2m, c2df
# %
# function C = m2c(M,df)
# if nargin < 2
#     df = ones(1,size(M,2));
# end
# C=cell(1,length(df));
# k=0;
# for i=1:length(df)
#     C{i} = M(:,(k+1):(k+df(i)));
#     k=k+df(i);
# end
#####################################################
m2c = function(M,df=rep(1,dim(M)[2])){
CC = vector("list",length(df))
k=0
for(i in 1:length(df)){
      CC[[i]] = M[,matlabColon(k+1,k+df[i]),drop = FALSE];
      k=k+df[i];
   } # end
CC
}# end m2c
