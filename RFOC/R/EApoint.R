`EApoint` <-
function()
{
 #   to plot points on an equale area stero net: 
deg2rad = pi/180
rad2deg = 180/pi

x = 0
y = 0
r = sqrt(x^2+y^2)
outx = vector()
outy = vector()
phis = vector()
angs = vector()

k=1
while(r<=1)
{
p = locator(n=1, type='p')
if(length(p$x)<1) break
r = sqrt(p$x^2+p$y^2)
if(r>1) break
iang = rad2deg*2*asin(r/sqrt(2))
phiang = rad2deg*( pi/2 - atan2(p$y,p$x))
print(paste(" ", p$x, p$y, "iANG=", iang, "PHI=",phiang))
outx[k] = p$x
outy[k] = p$y
phis[k] = phiang
angs[k] = iang
k = k+1
}

invisible(list(phi=phis, iang=angs, x=outx, y=outy) ) 

}

