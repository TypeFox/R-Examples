### $Id: c2df.R 45 2006-08-15 13:11:29Z bhm $
# %=============== c2df.m ====================
# % df = c2df(C)
# %   Input: C{1,*}  -  a matrix partitioned into a cell array
# %   Output: df(1,*) -  number of columns in the cells of C
# %
# %   See also: c2m, m2c 
# %
# function df = c2df(C)
# df = size(C{1},2);
# for i=2:length(C)
#    df = [df size(C{i},2)];
# end
###################################################
c2df = function(CC){ # Difference from matlab: returns m only
   df = dim(CC[[1]])[2];
   if(length(CC)>1) 
      for(i in 2:length(CC))
         df = c(df,dim(CC[[i]])[2])      
   df
}# end c2df
