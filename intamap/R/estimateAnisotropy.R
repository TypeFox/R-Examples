#######################################################################
#Function	:EstimateAnisotropy function 								
#Description	:Receives a Intamap Object or  calls estimateAnisotropySc 	
#				function and returns a new KrigingObject which 			
#				includes anisotropy parameters estimation.													
#																		
#Input		:object		& The Intamap object Object or Spatial points data 
#										& frame			
#					:anisPar	& Only used if object is Spatial and includes the
#										& dependant variable
#Output		:object		& (i) The modified intamap object in case 
#										& of Intamap object input with a list element anisPar 
#										& anisPar$(list(ratio=R,direction=theta.deg))		
#										& (ii) Only the anisPar list element in case of Spatial
#										& input.								
#######################################################################
estimateAnisotropy<-function(object, depVar,formulaString){
pl=FALSE
  if (is(object,"Spatial")) {
    observations = object
  	if (missing(formulaString)) {
	    if (!missing(depVar)) {
	      formulaString = as.formula(paste(depVar,"~1"))
      }else{
				 formulaString = as.formula("value~1")
			}
    }else{
			formulaString=as.formula(formulaString)
		}
	}else {
    observations = object$observations
    formulaString = object$formulaString
	}
		
	depVar = as.character(formulaString[[2]])
	if (!(formulaString[[3]] == 1)) {
  	m = autofitVariogram(formulaString,observations)$var_model
    g <- gstat(NULL, "value", formulaString, observations, model = m)
     blue0 <- predict(g, newdata = observations, BLUE = TRUE)
     residual <- observations[[depVar]] - blue0$value.pred
  } else residual = observations[[depVar]]
	
	# params = object$params
	xy<-as.matrix(coordinates(observations))
#

if (FALSE) {  # This is removed because of the akima licensing issue
  anisPar<-try(estimateAnisotropySc(xy[,1],xy[,2],residual,method="linear",pl=pl),TRUE)
	if(inherits(anisPar,"try-error")){
		warning("Slope tensors estimation error. Double size anisotropy estimation grid used instead.")
		anisPar<-try(estimateAnisotropySc(xy[,1],xy[,2],residual,len=2*length(xy),method="linear",pl=pl),TRUE)
	}
	if(inherits(anisPar,"try-error") && length(xy)<=400 ){
			warning("Slope tensors estimation error. Anisotropy interpolation method was switched to  biharmonics spline.")
		anisPar<-try(estimateAnisotropySc(xy[,1],xy[,2],residual,len=2*length(xy),method="v4",pl=pl),TRUE)
	}
	if(inherits(anisPar,"try-error") && length(xy)>400 ){
			warning("Slope tensors estimation error. Switch anisotropy interpolation to biharmonics spline method.")
		anisPar<-try(estimateAnisotropySc(xy[,1],xy[,2],residual,method="v4",pl=pl),TRUE)
	}
	if(inherits(anisPar,"try-error")){
			warning("Slope tensors estimation error. Override anistropy estimation.")
		anisPar<-list(R=1,theta.deg=0,Q=cbind(0,0,0),doRotation=FALSE)
	}
}
  anisPar<-try(estimateAnisotropySc(xy[,1],xy[,2],residual,len = min(2*length(xy), 400),method="v4",pl=pl),TRUE)
	if(inherits(anisPar,"try-error")){
			warning("Slope tensors estimation error. Override anistropy estimation.")
		anisPar<-list(R=1,theta.deg=0,Q=cbind(0,0,0),doRotation=FALSE)
	}



  if (is(object,"Spatial")) {
    anisPar
  } else {
    object$anisPar = anisPar
    object
  }
}

#######################################################################
#Function			: rotateAnisotropicData function 								
#Description	:	Receives an Intamap Object containing
#								an element observations containing the data coordinates 
#								and a anisPar element containing  the anisotropy parameters 
#								that will be used for the data rotation. 		 	
#Input		:object		& The original Intamap  or Spatial Object.				
#					:anisPar  & The anisotropy Parameters to be used for transformation 
#										& of the data.
#
#Output		:object		& with the transformed coordinates (intamap input) or
#										& the transformed coordinates(spatial input)	
#######################################################################
rotateAnisotropicData<-function(object,anisPar){
  #if (inherits(object,"Spatial")) {
  #EJP: changed inherits into is
  if (is(object,"Spatial")) {
    locations = object
  } else {
    if (missing(anisPar)) anisPar = object$anisPar
    locations = object$observations
  }
  if (missing(anisPar) || is.null(anisPar)) stop("Argument anisPar missing or equal to NULL")
	
	xy<-as.matrix(coordinates(locations))
	dataNames=colnames(xy)
	x<-as.matrix(xy[,1])
	y<-as.matrix(xy[,2])
	
	theta=(anisPar$direction)*pi/180
	R=anisPar$ratio

	#rotate data
	x_new=(x)*cos(theta)+y*sin(theta)
	y_new=R*(-x*sin(theta)+y*cos(theta))

	coords=as.data.frame(cbind(x_new,y_new))
	colnames(coords)<-dataNames


  coordinates(coords) = as.formula(paste("~",dataNames[1],"+",dataNames[2]))

  if ("data" %in% names(getSlots(class(locations)))) {
	  
    coords = SpatialPointsDataFrame(coords, data = locations@data)
  }
	if (!is(locations,"Spatial")) {
    object$observations = coords
    object
  } else 
    coords
}

#old version that supports only the Intamap Object.
rotateAnisotropicDataOld<-function(object){
	if ("formulaString" %in% names(object)) {
    formulaString = object$formulaString 
  } else formulaString = as.formula("value ~ 1")	
	depVar=as.character(formulaString[[2]])
	
	xy<-as.matrix(coordinates(object$observations))
	dataNames=colnames(xy)
	x<-as.matrix(xy[,1])
	y<-as.matrix(xy[,2])
	
	theta=(object$anisPar$direction)*pi/180
	R=object$anisPar$ratio

	#rotate data
	x_new=(x)*cos(theta)+y*sin(theta)
	y_new=R*(-x*sin(theta)+y*cos(theta))

	coords=cbind(x_new,y_new)
	colnames(coords)<-dataNames
	object$observations = SpatialPointsDataFrame(coords,data=object$observations@data)
	
	return(object)
}

######################### Internal functions ##########################
#######################################################################
#Function	estimateAnisotropySc										
#													
#Description : 	This is the main fucntion used in order to calculate 
#		anisotropy parameters from scattered data
#		
#		Given a Cartesian coordinates system  of axes x and y and an ellipsoid	
#		of correlation isolevels has principal directions  a and b. The anisotropy 
#		rotation angle(theta.deg) express the  angle between principal direction a and axes x.
#			
#												
#Inputs		&x			&The x-axis coordinate of field			
# 			&y			&The y-axis coordinate of field			
#			&r			&The field value in (xi,yi) point		
#			&len		&The length of the interpolated field		
#			&method		&The method that we will use for interpolation 
#			&			&"cubic", "linear", "v4"		
#			&min.x		&minimum x value 	
#			&max.x		&maximum x value		
#			&min.y		&minimum y value		
#			&max.y		&maximum x value		
#			&deb		&toggles debugging mode	Only Used Interactively, outside intamap 	
#			&pl			&toggles plot						Only Used Interactively, outside intamap 
#			&br			&borders coordinates		Only Used Interactively, outside intamap 
#Outputs	&R			&Anisotropy ratio		
#			&theta.deg	&Anisotropy rotation angle in degrees 
#																	
#Packages : akima package 									 
#													
#######################################################################
estimateAnisotropySc<-function(x, y, r, len=length(x), method="linear", min.x=min(x), max.x=max(x), min.y=min(y), max.y=max(y),deb=FALSE,pl=FALSE,br){
    eps=.Machine$double.eps        
	plot.borders=TRUE
	if(missing(br))	{plot.borders=FALSE}
  	if(len<50) return(list(ratio=1, direction=0,Q=c(eps,0,0),doRotation=FALSE))       

	#mesh creation

        #selection of interpolation method
#        if (method=="cubic")   {
#	ri<-interp(x,y,r,xn,yn,linear=FALSE,extrap=FALSE,duplicate="mean")
#				 ri$z<-t(ri$z)		}	#akima package
#        if (method=="linear")  {
#ri<-interp(x,y,r,xn,yn,linear=TRUE,extrap=FALSE,duplicate="mean")	
#				 ri$z<-t(ri$z)	}				#akima package
#        if  (method=="v4")     {  
  	step=min(c(abs((max.x-min.x)/sqrt(len)),abs((max.y-min.y)/sqrt(len))))
   	xn<-seq(min.x,max.x,by=step)
	  yn<-seq(min.y,max.y,by=step)
  if (length(x) < 1500) {
	  mesh<-meshgrid(xn,yn)
    ri<-biharmonicSplineAnisotropy(x,y,r,mesh$x,mesh$y)     
  } else if (FALSE) {
    dat = data.frame(x = x, y = y, z = r)
#    require(MBA)
    ri <- mba.surf(dat, length(xn), length(yn))$xyz.est
    rri <- mba.surf(dat, max(length(xn),length(yn)), max(length(xn),length(yn)))$xyz.est
    mesh = meshgrid(ri$x, ri$y)
  } else  warning("There is currently no method implemented for anisotropy detection in large data sets")
    		
	#calculate anisotropy parameters over regular grid 
	res=estimateAnisotropyGrid(mesh$x,mesh$y,ri$z)
	
		
	if((pl==TRUE & plot.borders==TRUE)){
	        image(mesh$x[1,],mesh$y[,1],t(ri$z),col=rainbow(20), #color.palette=rainbow,
	        plot.title=title(main=method),xlim=range(br[,1]),ylim=range(br[,2]))
		points(br,pch=".",xlim=(br[,1]),ylim=range(br[,2]))
	}else if ((pl==T & plot.borders==FALSE)){
		 image(mesh$x[1,],mesh$y[,1],t(ri$z),col=rainbow(20), #color.palette=rainbow,
	        plot.title=title(main=method),xlim=range(x),ylim=range(y))
	}
	
	return(list(ratio=res$R, direction=res$theta.deg,Q=res$Q,doRotation=res$dump$doRotation))
}


#######################################################################                              
# function [R , theta_deg ]= estimateAnisotropyGrid(xi ,yi , ri)  		
#     theta_deg , R :variables to b calculated              			
#     xi ,yi    : Coordinates on Grid                        			
#         ri    : Random field at (xi,yi)                    			
#######################################################################

estimateAnisotropyGrid<-function(xi,yi,ri){
	dump=NULL
	#Gradient function
	mgradient<-function(k,stx,sty){

		gg<-function(m,st){
			m<-as.matrix(m)
				n=dim(m)[1]                     #rows
				p=dim(m)[2]                     #columns

				g<-matrix(0,n,p)
				g[1,]<-(m[2,]- m[1,])/st        
				g[n,]<-(m[n,]-m[n-1,])/st       
				g[2:(n-1),]=(m[3:n,]- m[1:(n-2),])/(2*st)                        
				return (g)
		}

		y=gg(k,sty)                       
		x=t(gg(t(k),stx))                 
			return(list(x=x,y=y))             
	}
#################################


  stepx<-xi[1,2]-xi[1,1]
  stepy<-yi[2,1]-yi[1,1]
  divZ<-mgradient(ri,stepx,stepy)

  divX2<-as.vector(divZ$x*divZ$x)
  divY2<-as.vector(divZ$y*divZ$y)
  divXY<-as.vector(divZ$x*divZ$y)

  #get rid of Na
  Q11<-mean(divX2[is.na(divX2)==FALSE])
  Q22<-mean(divY2[is.na(divY2)==FALSE])
  Q12<-mean(divXY[is.na(divXY)==FALSE])
	if(Q11<10^-31 ||	is.na(Q11)||is.na(Q12)||is.na(Q22)||(Q11*Q22)<(Q12^2)){
		dump=NULL
		dump$doRotation=FALSE
		e=simpleError("slope tensor error")
		stop(e)	
	 #return (list(R=1,theta.deg=0,Q=cbind(Q11,Q22,Q12),dump=dump))
	}else{

  Zdiag<-Q22/Q11
  Zoff<-Q12/Q11
  theta<-0.5*atan(2*Zoff/(1-Zdiag))

  # R=ifelse(Zdiag-1<10*.Machine$double.eps,1,sqrt(1+((1-Zdiag)/(Zdiag-(1+Zdiag)*(cos(theta))^2))))
	
R=sqrt(1+((1-Zdiag)/(Zdiag-(1+Zdiag)*(cos(theta))^2)))

	#Test for isotropy 
 	#dump<-jpdf(seq(0,3,by=0.01),seq(-90,90,by=1),R,theta*180/pi,length(xi))
	dump=jpdf(length(xi),R)				
	
	if(is.na(R)){
		R=1 	
		theta=0
		dump$doRotation=FALSE
	}

   #Make sure that all results come in the same interval
	if (R<1){
		R=1/R
			if ((theta+pi/2>-pi/2) & (theta+pi/2)<(pi/2)){
				theta=theta+pi/2
			}else{
				theta=theta-pi/2
			}
	}
	theta.deg<-theta*180/pi
	}
  return (list(R=R,theta.deg=theta.deg,Q=cbind(Q11,Q22,Q12),dump=dump))

}

# Description: function that tests if the predicted anisotropy ratio
#							 is inside the 95% confidence interval of the jpdf function
#							 for statistical isotropy. Replaced jpdf_old for speed reasons
jpdf=function(N,R){

	doRotation=TRUE
	r=6	#comes from 95%interval
	
	floor=sqrt((N-2*sqrt((N-r)*r))/(N-2*r))
	ceil =sqrt((N+2*sqrt((N-r)*r))/(N-2*r))

	if(R>floor && R <ceil) doRotation=FALSE 
	
	return(list(doRotation=doRotation))
}


########################################################
#Function   : jpdf_old function 
#
#Description: calculates the joint pdf function for R
#       and theta over a given region 
#Arguments  : 
#   R   : a mesh type matrix of the R region    
#   theta   : a mesh type matrix of the theta region
#   Rest    : the given ratio value
#   thetaEst: the given theta value
#
#Longer Description:
#   The joint pdf gives the probability of R AND theta
#   to lie within the infinitesimal region dR*dtheta. It
#   is normalized to 1 if integrated to R belonging to [0, inf) and
#   theta belonging to (-45, 45) degrees. At the moment it is
#   not dependent to a specific covariance function model since it
#   is a first approximation. But this is useful since it accounts
#   for the worst-case scenario. In other words, the confidence
#   region it provides is conservative. A version which will
#   account for a specific correlation model is under development.
#
#   A few words about the confidence interval (region). The idea is to
#   define the parametric equation of the curve in the (R,theta)-plane
#   which "encloses" 95% of the probability. The equation has the form
#   c(R,theta)=0 so the points (R, theta) for which the equation holds,
#   define the confidence region boundary. Any point (R', theta') for
#   which c(r',theta')>0 lies within the confidence region. The percentile
#   (here 95%) determined by the parameter r and is the inverse of
#   the chi-square distribution with 2 degrees of freedom
#   and the percentile requested (e.g. 0.95).  
#######################################################
jpdf_old<-function(Rf,thetaf,Rest,thetaEst,N){
    
    
    erfc <- function(x){ 2 * pnorm(x * sqrt(2), lower.tail = FALSE)}
        
    RRest=Rest
    tthetaEst=thetaEst
    
    #so far the jpdf is not numerically stable for large values of N due to
    #the exponential terms ~exp(N). But this is not a problem for the
    #evaluation of the confidence interval. So I moved the... 
    #if(N>1200){N=1200}
    #...where it is necessary.

    Rest=1
    thetaEst=0
    
    
    #create mesh grids
    mesh<-meshgrid(Rf,thetaf)
    R<-mesh$x
    theta<-mesh$y

    #threshold value to decide if the given set is in 95% interval of Jpdf
    threshold=5.99416



    #first convert degrees to Radians
    theta=theta*pi/180
    thetaEst=thetaEst*pi/180

    
    #Matrix elements of H -- constants 
    s1=2*(cos(thetaEst)^2 +Rest^2 * sin(thetaEst)^2 )
    s2=2*(sin(thetaEst)^2 +Rest^2 * cos(thetaEst)^2 )
    s12=2*sin(thetaEst)* cos(thetaEst)*(1-Rest^2)


    #Determinant of H -- constant
    d=s1*s2-s12^2
    
    qd = (R^2 +tan(theta)^2) /(1 + R^2 *tan(theta)^2)
    qo = tan(theta) * (1-R^2)/ (1+R^2 *tan(theta)^2)

    #Jacobian for the transformation (qd,qo) -> (R,theta)
    J1= ((2*R*abs(R^2 -1)*(1/cos(theta))^6))
    J2=(1+R^2 *tan(theta)^2)^3
    J=J1/J2

    a = N/(2*d^2) * (qd^2 * s1^2  + 2 * qd * s12^2 + s2^2 - 4 * qo * s12 * (qd * s1 + s2) - 2 * qo^2 * (d - 2 * s1*s2)) 
    b = -(N/d) * (qd * s1 - 2 * qo * s12 + s2);
    k = (N/d)^(3/2) / (4*sqrt(2)* pi^(3/2));

    #contour plot
    lhs = 2 * (R^4 * s1^2 + 2 * R^2 * s12^2 + s2^2) * threshold + tan(theta)^2 * (-2 * N * (2 * s12^2 + s1 * s2 + R^4 * (2 * s12^2 + s1 * s2) + R^2 * (s1^2 - 4 * s1 * s2 + s2^2)) + 4 * (2 * s12^2 + R^2 * (-2 * s12^2 + (s1 - s2)^2) + s1 * s2 + R^4 * (2 * s12^2 + s1 * s2))* threshold + (-N * (4 * R^2 * s12^2 + (s1 - R^2 * s2)^2) + 2 * (s1^2 + 2 * R^2 * s12^2 + R^4 * s2^2) * threshold) * tan(theta)^2)
    rhs = N *(R^4 * s1^2 + s2^2 + R^2 * (4 * s12^2 - 2 * s1 * s2)) + 4 * (-1 + R^2) * s12 * (N - 2 * threshold) * tan(theta) * (R^2 * s1 + s2 + (s1 + R^2 * s2) * tan(theta)^2)
    cc= lhs - rhs
    
    #(see comment in the beginning of the function)
    if(N>1200){N=1200}
    
    #jpdf
    p = (pi/180)*k*J* (exp(-N/2) * (-4 * sqrt(a)* b + (4*a + b^2)* exp(b^2/(8*a))*sqrt(2*pi)*erfc(b/(2*sqrt(2)*sqrt(a)))))/(8*a^(5/2))
    

    #The parametric equation of the contour interval is cc=0. If cc>=0 we are in the
    #interior of the confidence interval
    booleanC<-(cc>=0)

    maxV<-max(mesh$x[which(booleanC)])
    minV<-min(mesh$x[which(booleanC)])
    print(paste("jpdf 95% R maximum ratio for statistical isotropy: ",maxV))
    doRotation=(RRest>maxV || RRest<minV)
#   doRotation=(max(mesh$x[which(booleanC)])<RRest)
        
    
    
    
#   contour(as.matrix(Rf),as.matrix(thetaf),t(p))
#   contour(Rf,thetaf,t(p))

    return(list(c=cc,p=p,bc=booleanC,doRotation=doRotation,R=maxV))




}





#######################################################################
#biharmonics Spline interpolation			
#   scattered data coordinates -> x, y			
#   data value @(x,y)          -> z			
#   grid                       ->xi, yi			
#######################################################################
biharmonicSplineAnisotropy<-function(x,y,z,xi,yi,extrap=FALSE){
	#sortrows function
	SortMat <- function(Mat, Sort)
	{
 	       m <- do.call("order", as.data.frame(Mat[, Sort]))
	        Mat[m, ]
	}

        #sorting input
        Mat<-cbind(x,y,z)
        Mat<-SortMat(Mat,c(2,1))
        x<-Mat[,1]
        y<-Mat[,2]
        z<-Mat[,3]

        #Move to Complex field

        xy<-as.matrix(x+y*1i)
        d<-matrix(xy,length(xy),length(xy))

        #Calculate distances
        d<-abs(d-t(d))

        #remove zeros so we can use log
	ind=which(d==0)
	if(length(ind)>0){
        d[ind]=1
	}
        d=d^2*(log(d)-1)

        #return replace zeros to diag elements
        d[seq(1,length(d),by=(dim(d)[1]+1))]=0
	
	ind=which(is.na(d)==TRUE)
	if(length(ind)>0){
        	d[ind]=1
	}
        
	#calculate weights
        wweights=solve(d)%*%z
        rm(d)
        m=dim(xi)[1]
        n=dim(xi)[2]
        zi<-matrix(0,m,n)
        xy=t(xy)


        for(i in 1:m){
         for(j in 1:n){
             d = abs(xi[i,j]+(1i*yi[i,j])-xy)
             mask=which(d==0)

             if (length(mask)>0) {d[mask]=1}
             g<-(d)^2*(log(d)-1)
             if(length(mask)>0) { g[mask]=0}
              zi[i,j]=g%*% wweights
         }
        }

if(extrap==FALSE){
	cvh=chull(x,y)
	cvh=c(cvh,cvh[1])
	cvhPolygon=SpatialPolygons(list(Polygons(list(Polygon(data.frame(x=x[cvh],y=y[cvh]))),"r1")))
	interpolationGrid=expand.grid(x=xi[1,],y=yi[,1])
	coordinates(interpolationGrid)=c("x","y")
	
	index=which(is.na(over(interpolationGrid, cvhPolygon)))
	temp=as.vector(t(zi))
	temp[index]=NA
	zi=t(matrix(temp,n,m)	)
	
}

  return(list(x=xi,y=yi,z=zi))

}


	#Meshgrid function
	meshgrid <- function(a,b) {
	  list(
	       x=outer(b*0,a,FUN="+"),
	       y=outer(b,a*0,FUN="+")
	       )
	}


########################End of internal functions######################



