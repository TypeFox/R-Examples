GM <-
function(inpmat,alpha)
{
k = length(inpmat)
s = sum(inpmat)
chi = qchisq(1-(alpha/k), df=1)
pi = inpmat/s
GM.UL = (chi + 2*inpmat + sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s))
GM.LL = (chi + 2*inpmat - sqrt(chi*chi + 4*inpmat*chi*(1 - pi)))/(2*(chi+s))
LLA=0
ULA=0
for (r in 1:length(inpmat))
{
if (GM.LL [r]< 0) LLA[r] = 0 else LLA[r]=GM.LL[r]
if (GM.UL[r] > 1) ULA[r] = 1 else ULA[r]=GM.UL[r]
}
diA=ULA-LLA##FIND LENGTH OF CIs
VOL=round(prod(diA),8)##PRODUCT OF LENGTH OF CIs
cat('Original Intervals\n')
cat('Lower Limit\n')
print(GM.LL)
cat('Upper Limit\n')
print(GM.UL)
cat('Adjusted Intervals\n')
cat('Lower Limit\n')
print(LLA)
cat('Upper Limit\n')
print(ULA)
cat('Volume\n')
print(VOL)
}
