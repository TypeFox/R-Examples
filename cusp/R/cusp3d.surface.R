cusp3d.surface <- function(
	alpha = c(-5,5),
	beta = c(-3,3),
	y = 41,
	xlim = range(alpha), ylim=range(beta),
	zlim = c(-5,4),
	xlab=expression(alpha), ylab=expression(beta), zlab="equilibrium states",
	main = NULL, sub = NULL,
	phi=20, theta=160, r=sqrt(3), d=1,
	scale = TRUE, expand = 1, 
	hue=240, chroma=35,	surf.alpha=0.75, gamma=1.5,	bcol=NA, lcol="gray", 
	ltheta=90, lphi=70, box=TRUE, axes=FALSE, nticks = 5, 
	ticktype = "simple", floor.lines=TRUE, ...)
{
	ltheta = ltheta * 2 * pi/360
	lphi = lphi * 2 * pi/360
	light = c(cos(ltheta) * sin(lphi), sin(ltheta) * sin(lphi), cos(lphi))
    extprod = function(a, b) c(a[2] * b[3] - a[3] * b[2], 
            a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
	pmat = persp( alpha, beta, matrix(min(zlim),2,2), 
		xlim=xlim, ylim=ylim, zlim=zlim,
		theta=theta,phi=phi,r=r,d=d, 
		xlab=xlab, ylab=ylab, zlab=zlab, axes=axes, box=box, 
		nticks=nticks, ticktype=ticktype,
		main=main, sub=sub)
	# make y-range
	if(length(y)==1){
		root1 = polyroot(c(alpha[1],beta[2],0,-1))
		root1 = Re(root1[abs(Im(root1))<1e-10])
		root2 = polyroot(c(alpha[2],beta[2],0,-1))
		root2 = Re(root2[abs(Im(root2))<1e-10])
		yr = seq(root1, root2, len=y)
	}
	else yr = y
	if(min(yr)<0 & max(yr)>0) {
		yr = sort(c(0,yr)) # always include Y=0 if in range
	}
    plotFloorLines <- function(){
		for(Y in yr){
			# 3D iso-contour coordinates projected onto ground floor projected onto 2D canvas
			y = c(Y,Y)
			b = c(-5,5)
			b[1] = min(beta[2], max(beta[1],(alpha[2]+Y^3)/Y))  # stay inside the box!
			b[2] = max(beta[1], min(beta[2],(alpha[1]+Y^3)/Y))
			a = -b*y+y^3
			# plot floor lines
			lines(trans3d(-a,b,rep(min(zlim),2),pmat), col=lcol, ...)
		}
    }
    plotSurface <- function(){
		Y = min(yr)
		yo = range(Y,Y)
		bo = c(-5,5)
		bo[1] = min(beta[2], max(beta[1],(alpha[2]+Y^3)/Y)) # stay inside the box!
		bo[2] = max(beta[1], min(beta[2],(alpha[1]+Y^3)/Y))
		ao = -bo*yo+yo^3
		# plot surface elements
		for(Y in yr){
			# 3D iso-contour coordinates projected onto 2D canvas
			y = c(Y,Y)
			b = c(-5,5)
			b[1] = min(beta[2], max(beta[1],(alpha[2]+Y^3)/Y))
			b[2] = max(beta[1], min(beta[2],(alpha[1]+Y^3)/Y))
			a = -b*y+y^3
			# determine surface element color
			v1 = c(a[1],b[1],y[1]) - c(ao[2],bo[2],yo[2])
			v2 = c(a[2],b[2],y[2]) - c(ao[1],bo[1],yo[1])
			normal = extprod(v1, v2)
			normal = normal/sqrt(sum(normal^2))
			shade = sum((normal * light)) * 0.5 + 0.5
			shade = if (is.na(shade)) 1 else shade^gamma
			col = hcl(hue, chroma, shade * 100, surf.alpha)
			# treat case Y==0 differently
			if(Y==0) {
				polygon(trans3d(-c(ao,(a)),c(bo,(b)),c(yo,rev(y)),pmat),,col=colo,border=if(!is.na(bcol)&&bcol=="surface") colo else bcol,...)
			}
			else {
				polygon(trans3d(-c(ao,rev(a)),c(bo,rev(b)),c(yo,rev(y)),pmat),,col=col,border=if(!is.na(bcol)&&bcol=="surface") col else bcol,...)
			}
			yo=y; bo=b; ao=a; colo=col
		}
	}
	if(phi>0){
		if(floor.lines)
			plotFloorLines();
		plotSurface();
	}
	else{
		plotSurface();
		if(floor.lines) {
			plotFloorLines();
		}
	}
	invisible(pmat);
}


