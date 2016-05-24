WALD <-
function(inpmat,alpha)
{
k = length(inpmat)
s = sum(inpmat)
chi = qchisq(1-alpha, df=1)
pi = inpmat/s
WALD.LL = pi - (sqrt(chi*(pi)*(1-pi)/s)) 
WALD.UL = pi + (sqrt(chi*(pi)*(1-pi)/s)) 
LLA=0
ULA=0
for (r in 1:length(inpmat))
{
if ( WALD.LL [r]< 0) LLA[r] = 0 else LLA[r]=WALD.LL[r]
if (WALD.UL[r] > 1) ULA[r] = 1 else ULA[r]=WALD.UL[r]
}
diA=ULA-LLA##FIND LENGTH OF CIs
VOL=round(prod(diA),8)##PRODUCT OF LENGTH OF CIs
cat('Original Intervals\n')
cat('Lower Limit\n')
print(WALD.LL)
cat('Upper Limit\n')
print(WALD.UL)
cat('Adjusted Intervals\n')
cat('Lower Limit\n')
print(LLA)
cat('Upper Limit\n')
print(ULA)
cat('Volume\n')
print(VOL)
}
