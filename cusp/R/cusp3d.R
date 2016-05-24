`cusp3d` <-
function(y, 
alpha = if(!missing(y) && is.list(y)) y$lin[,'alpha'], 
beta=if(!missing(y) && is.list(y)) y$lin[,'beta'],
w = 0.03,
theta = 170,
phi = 35,
B = 4,
Y = 3,
Yfloor = -15,
np = 180, 
n.surface = 30,
surface.plot = TRUE,
surf.alpha = 0.75,
surf.gamma = 1.5,
surf.chroma = 35,
surf.hue = 240,
surf.ltheta = 0,
surf.lphi = 45,
...){
    A = Y^3 - B*Y
    # doorsnedes op beta = B
    v1 = cbind(a=-(function(x,b=B) x^3-b*x)(seq(-Y,Y,len=np)),b=B,y=seq(-Y,Y,len=np))
    v2 = cbind(a=-(function(x,b=-B) x^3-b*x)(seq(-Y,Y,len=np)),b=-B,y=seq(-Y,Y,len=np))
    # doorsnedes op alpha = A = Y^3 - B*Y en alpha = -A
    b1 = min(v1[abs(v1[,1])<=A,][,3]); b2 = max(v1[abs(v1[,1])<=A,][,3])
    e1 = min(v2[abs(v2[,1])<=A,][,3]); e2 = max(v2[abs(v2[,1])<=A,][,3])
    v3 = cbind(a=A,b=(function(y,a=-A)y^2-a/y)(seq(b1,e1,len=np)),y=seq(b1,e1,len=np))
    v4 = cbind(a=-A,b=(function(y,a=A)y^2-a/y)(seq(b2,e2,len=np)),y=seq(b2,e2,len=np))
    # overhellingslijn (die bifurcatie-set bepaald)
    bf = cusp.bifset(seq(0, B,len=np))
    v5 = cbind(a=-bf[,2],b=bf[,1], y=Vectorize(cusp.extrema)(bf[,2],bf[,1])[3,])
    v6 = cbind(a=-bf[,3],b=bf[,1], y=Vectorize(cusp.extrema)(bf[,3],bf[,1])[1,])
    # bifurcatie set op vloer
    v7 = cbind(a=bf[,2], b=bf[,1], y=Yfloor)
    v7 = rbind(cbind(a=bf[,3],b=bf[,1],y=Yfloor),v7[rev(1:NROW(v7)),])
    # vierkant op vloer om bifurcatie-set heen
    v8 = cbind(a=c(A,A,-A,-A),b=c(B, -B,-B,B),y=Yfloor)
    # 'stutten' van vloer naar cusp-oppervlak
    v9 = cbind(a=A,b=B,y=c(Yfloor,b1))
    v10 = cbind(a=A,b=-B,y=c(Yfloor,e1))
    v11 = cbind(a=-A,b=B,y=c(Yfloor,b2))
    v12 = cbind(a=-A,b=-B,y=c(Yfloor,e2))
    
    helper3d = (function(v,p){v=v[abs(v[,1])<=A,];trans3d(v[,1],v[,2],v[,3],p);})
    asRect3d <- function(i,j, a, b, Z=Vectorize(function(a,b)cusp.extrema(a,b)[1])){
        # returns 3d coordinates of a rectangle on a surface specified by Z given coordinate vectors a and b
        # at the i,j-th coordinate
        x = -a[I <- c(i+1,i,i,i+1)]
        y =  b[J <- c(j,j,j+1,j+1)]
        z =  if(is.function(Z)) {Z(-x,y)} else {Z[cbind(I,J)]}
        list(x=x,y=y,z=z)
    }
    plotSurfElem <- function(x, y, z, pmat, lphi=surf.lphi,ltheta=surf.ltheta, 
        hue=surf.hue, chroma=surf.chroma, alpha=surf.alpha, gamma=surf.gamma, bcol=col, ...){
    	# plots a 3d rectangular surface element on the current device using a transform matrix returned by persp
    	# and calculates the fill color taking shading into account using a light source shining from
    	# direction (ltheta, lphi). gamma can be used to exagerate or diminish shading, alpha is transparency channel,
    	# hue and chroma specify the color (see ?hcl)
    	V =  cbind(x,y,z)
    	v1 = V[1,] - V[2,]
    	v2 = V[2,] - V[3,]
    	extprod = function(a,b) c(a[2] *b[3] - a[3] *b[2] , a[3] *b[1] - a[1] *b[3] , a[1] *b[2] - a[2] *b[1] )
    	normal = extprod(v1,v2);normal = normal / sqrt(sum(normal^2))
    	ltheta = ltheta * 2 * pi / 360
    	lphi = lphi * 2 * pi / 360
    	light = c(cos(ltheta)*sin(lphi), sin(ltheta)*sin(lphi),cos(lphi)); 
    	shade =sum(normal*light)*0.5 + 0.5; 
    	shade = if(is.na(shade)) {1} else {shade^gamma}
    	#on.exit(print(shade))
    	col= hcl( hue, chroma, shade*100, alpha)
    	polygon(trans3d(x,y,z,pmat), , col=col, border=bcol, ...)
    }
    vce = Vectorize(cusp.extrema)

    plotLowerSurface <- function(){
    	for(I in i)
    	for(J in j){
    		if(!is.nan(aup[J]) & a[I+1]>aup[J]){
    			if(a[I]<aup[J]){
    				tmpa = rev(c(a[I],aup[J],aup[J+1],a[I]));
    				tmpb = rev(c(b1[J],b1[J],b1[J+1],b1[J+1])); 
    				tmpz = vce(tmpa,tmpb)[1,]
    				.rect <- list(x= -tmpa, y= tmpb,z= tmpz)
    				plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    			}
    			next
    		}
    		.rect = asRect3d(I, J, a=a, b=b1, Z=z1)
    		plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    	}
    }
    plotMiddleSurface <- function(){
    	for(I in i)
    	for(J in j){
    		if(is.nan(aup[J]) | a[I+1]>aup[J] | a[I]<aul[J]) {
    			if(!is.nan(aul[J]) & a[I+1]>=aul[J] & a[I]<aul[J]){
    				tmpa = (c(aul[J],min(a[I+1],aup[J]),min(a[I+1],aup[J+1]),aul[J+1]))
    				tmpb = (c(b1[J],b1[J],b1[J+1],b1[J+1]))
    				tmpz = vce(tmpa,tmpb)[2,]
    				.rect <- list(x= -tmpa, y= tmpb,z= tmpz)
    				plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    			}
    			if(!is.nan(aup[J]) & a[I] <= aup[J] & a[I+1]>aup[J]){
    				tmpa = rev(c(max(a[I],aul[J]), aup[J], aup[J+1], max(a[I],aul[J+1])))
    				tmpb = rev(c(b1[J],b1[J],b1[J+1],b1[J+1]))
    				tmpz = vce(tmpa,tmpb)[2,]; #print(rbind(tmpa,tmpb)); print(vce(tmpa,tmpb)); print("")
    				.rect <- list(x= -tmpa, y= tmpb,z= tmpz)
    				plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    			}
    			next
    		}
    		.rect = asRect3d(I, J, a=a, b=b1, Z=z2)
    		plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    	}
    }
    plotUpperSurface <- function(){
    	for(I in i)
    	for(J in j){
    		if(is.nan(aup[J]) | a[I]<aul[J]) {
    			if(!is.nan(aul[J]) & a[I+1]>=aul[J]){
    				tmpa = rev(c(aul[J],a[I+1],a[I+1],aul[J+1]))
    				tmpb = rev(c(b1[J],b1[J],b1[J+1],b1[J+1]))
    				tmpz = vce(tmpa,tmpb)[3,];
    				.rect <- list(x= -tmpa, y= tmpb,z= tmpz)
    				plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    			}
    			next
    		}
    		.rect = asRect3d(I, J, a=a, b=b1, Z=z3)
    		plotSurfElem(.rect$x, .rect$y, .rect$z, pmat=r, ...)
    	}
    }
    drawPoint = function(alpha, beta, w, y, col="silver", ...){
        .a = alpha+w*crcl$a;
        .b = beta +w*crcl$b;
        .bf= cusp.bifset(.b)
        .mna = .bf[,'alpha.l']
        .mxa = .bf[,'alpha.u']; 
		rts = vce(alpha,beta)
		upps = y>=0 # y lies on upper sheet (if in bifurcation set)?
        ibs = .a^2/4 - .b^3/27 > 0
        .na = if(upps) {
        		ifelse(ibs & !is.nan(.mna) & .a < .mna, .mna, .a) 
        	} 
        	else {
        		ifelse(ibs & !is.nan(.mxa) & .a > .mxa, .mxa, .a) 
        	}
        .ny = vce(.na, .b)[if(upps) 3 else 1, ]
        crcl.xyz = cbind(-.na, .b, .ny)
        
        polygon(helper3d(crcl.xyz, r), col=col, ...)
    }
    drawPoints = Vectorize(drawPoint)

    # data point shape coordinates
    .x = seq(0,2*pi,len=max(300*w,15))
    crcl = list(a = A*sin(.x), b = B*cos(.x))

	a = sort(c(seq(-A,A,len=n.surface-1),0))
	b1 = sort(c(seq(-B,B,len=n.surface-1),0))
	z1 = outer(a,b1,function(a,b)vce(a,b)[1,])
	z2 = outer(a,b1,function(a,b)vce(a,b)[2,])
	z3 = outer(a,b1,function(a,b)vce(a,b)[3,])
	i = 2:length(a)-1
	j = 2:length(b1)-1
	aup <- cusp.bifset(b1)[,'alpha.u']
	aul <- cusp.bifset(b1)[,'alpha.l']

    # setup new plot window
    r = persp(matrix(-Yfloor,2,2),zlim=c(Yfloor,Y),xlim=c(-A,A),ylim=c(-B,B),box=FALSE,xlab='alpha',
        ylab='beta',zlab='y',border='transparent',d=1,theta=theta,phi=phi, ...)
    
    # first the floor and frame
    polygon(helper3d(v7,r),col='gray')
    polygon(helper3d(v8,r),lty=3)
    for(v in paste('v',9:12,sep='')) lines(helper3d(get(v), r), lty=3)
    # the cusp surface
    if(surface.plot){
        plotLowerSurface()
    }
    
    ibs = alpha^2/4 - beta^3/27 > 0 # inside bifurcation set
    .bf <- cusp.bifset(beta); 
    .bf.au <- .bf[,3]
    .bf.al <- .bf[,1]
    .ce <- t(vce(alpha, beta))
    .y <- if(is.list(y)) {y$y} else {y}
    idc =  alpha < ifelse(is.nan(.bf.au), 0, .bf.au) & abs(.y-.ce[,1])<abs(.y-.ce[,3]) #.y < .ce[,2]
    idc =  .y < 0
    # data points lower surface
    if(!missing(y)) {
       drawPoints(alpha=alpha[idc], beta=beta[idc], w=w, y=(if(is.list(y)) y$y else y)[idc], col='gray')
    }
    # middle and upper surface
    if(surface.plot){
        plotMiddleSurface()
    }
    for(v in paste('v',5:6,sep='')) lines(helper3d(get(v), r))
    if(surface.plot){
        plotUpperSurface()
    }
    for(v in paste('v',1:4,sep='')) lines(helper3d(get(v), r))
    # data points upper surface
    if(!missing(y)) {
       drawPoints(alpha=alpha[!idc], beta=beta[!idc], w=w, y=(if(is.list(y)) y$y else y)[!idc], col='gray')
    }

    invisible(r)
}

