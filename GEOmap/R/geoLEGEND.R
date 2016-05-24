`geoLEGEND` <-
function(names, shades, zx, zy, nx, ny, side=1, cex=0.5)
{

u = par("usr")
pin = par("pin")
pai = par("mai")

dx  = zx*(u[2]-u[1])/pin[1]
dy = zy* (u[4]-u[3])/pin[2]

N = length(names)

if(side==1){
upperleft.x = u[1]
upperleft.y = u[3]
}

if(side==2){
upperleft.x =u[1] - pai[2]* (u[2]-u[1])/pin[1] + 0.1*pai[2]
upperleft.y = u[4]
}

if(side==3){
upperleft.x =u[1] 
upperleft.y = u[4] + pai[3]* (u[4]-u[3])/pin[2] - 0.1*pai[3]
}

if(side==4){
upperleft.x = u[2]
upperleft.y = u[4]
}


for(i in 1:N)
{

ii = i
ixy = RPMG::itoxyz(i, nx, ny,1)

##print(paste(sep=" ", i, ixy$ix, ixy$iy))

x1 = upperleft.x+(ixy$ix-1)*dx
x2 = x1+dx
y1 = upperleft.y -(ixy$iy+1)*dy
y2 = y1+dy

rect(x1,y1,x2,y2, col=shades[ii], xpd=TRUE)
text(mean(c(x1,x2)), mean(c(y1,y2)), labels=names[ii], xpd=TRUE, cex=cex)

}

}
