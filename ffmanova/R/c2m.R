### $Id: c2m.R 45 2006-08-15 13:11:29Z bhm $
# %=============== c2m.m ====================
# % [M,df] = c2m(C)
# %   Input:
# %           C{1,*}  -  a matrix partitioned into a cell array
# %   Output:
# %           M(*,*)  -  C as ordinary unpartitioned matrix
# %           df(1,*) -  number of columns in the cells of C
# %
# %   See also: c, c2df
# %
# function [M,df] = c2m(C)
# M = C{1};
# df = size(C{1},2);
# for i=2:length(C)
#    M = [M C{i}];
#    df = [df size(C{i},2)];
# end
##################################################
## Rewritten by bhm
c2m <- function(CC) { # Difference from matlab: returns m only
   do.call("cbind", CC)
}# end c2m
