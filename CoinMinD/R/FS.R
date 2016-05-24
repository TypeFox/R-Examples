FS <-
function(inpmat,alpha)
{
k = length(inpmat)
s = sum(inpmat)
zval = abs(qnorm(1-(alpha/2)))
pi = inpmat/s
FS.LL = pi - (zval/(2*sqrt(s))) 
FS.UL = pi + (zval/(2*sqrt(s))) 
LLA=0
ULA=0
for (r in 1:length(inpmat))
{
if ( FS.LL [r]< 0) LLA[r] = 0 else LLA[r]=FS.LL[r]
if (FS.UL[r] > 1) ULA[r] = 1 else ULA[r]=FS.UL[r]
}
diA=ULA-LLA##FIND LENGTH OF CIs
VOL=round(prod(diA),8)##PRODUCT OF LENGTH OF CIs
cat('Original Intervals\n')
cat('Lower Limit\n')
print(FS.LL)
cat('Upper Limit\n')
print(FS.UL)
cat('Adjusted Intervals\n')
cat('Lower Limit\n')
print(LLA)
cat('Upper Limit\n')
print(ULA)
cat('Volume\n')
print(VOL)
}
