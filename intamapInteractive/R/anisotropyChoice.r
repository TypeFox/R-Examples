#####################################################################################
#
#Function	: AnisotropyChoice(object)
#
#Description	: This method combines segmentation of a large dataset 
#		  and anistopy parameters estimation for scattered data.
#		  
#Arguments 	: An intamap type object
#
#Details	: The function AnisotropyChoice function employs the
#	doSegmentation function to  automatically separate the original 
#	dataset into clusters based on the sampling density and the spatial
#	locations of the data. The results of the segmentation procedure and 
#	the anisotropy analysis per cluster are returned in  a matrix of
#	dimension [cl]x5, where [cl] is the number of clusters . Each row of
#	the matrix contains the cluster index, the anisotropy ratio, the
#	anisotropy direction, the number of cluster points and the area
#	inside the convex hull of the cluster. In addition, a single set of
#	anisotropy parameters is returned in the element \code{anisPar}.
#	These parameters are calculated using weighted averages of the
#	covariance Hessian matrix estimates in each cluster. The weights are
#	based on the area enclosed by the convex hull of each cluster.  
##########################################################################################

anisotropyChoice<-function(object){
	
	params<-object$params
	params$doSegmentation=TRUE
	if(params$doSegmentation==FALSE){
	#This is the case that no Segmentation is done 

	
	object$params$AnisPar<-estimateAnisotropy(object)
	return (object)


	}else if (params$doSegmentation==TRUE){
	#returns list with a single pair of anisotropy parameters for all data(singleCl)
	#and a matrix (cln,5) with the values for 
	
	
	if ("formulaString" %in% names(object)) formulaString = object$formulaString else formulaString = as.formula("value ~ 1")	
	depVar=as.character(formulaString[[2]])

	xy<-as.matrix(coordinates(object$observations))
	z<-object$observations[[depVar]]

		
	#first run segmentation
	xyz<-as.matrix(cbind(xy,z))
	if(dim(xy)[1]<200){
		object$clusters=list(
											index=matrix(1,dim(xy)[1],1),
											clusterNumber=1,
											rmdist=NULL
										)
				rmInd=NULL
				clNumber=1
				index=matrix(1,dim(xy)[1],1)
				ind=index
	}else{


		
	
	rmdist=TRUE
	#Step 0 remove distant points
	if(rmdist==TRUE){
		rmInd<-rmDist(xyz)	
		if(length(rmInd)>0){
			xyz_d<-xyz[-rmInd,]	
		}else{ 
		xyz_d=xyz
		}
	}


	# do Segmentation 
	segmentResult<-segmentData(ddd=xyz_d,pl=FALSE,dev=FALSE)
	
	#This index is 
	ind<-segmentResult$index
	clNumber<-segmentResult$clusterNumber
	#rmDist<-segmentResult$clusters$rmDist

	}


	xyz_new<-xyz
	if (length(rmInd)>0){
		xyz_new<-xyz[-rmInd,]
	}

	Qs<-matrix(0,clNumber,3)
	area<-matrix(0,clNumber,1)

	results=matrix(0,clNumber,5)



	#create an index  size of original data xyz 0 value means that the point 
	#is not assigned to a cluster and should not be included in anisotropy estimation.
	index=matrix(0,length(xyz[,1]),1)
	if(length(rmInd)>0){
	index[-rmInd]=ind	
	}else{
	index=ind
	}
	
	

	object$clusters$index=index
	object$clusters$clusterNumber=clNumber
	

	#find anisotropy for each cluster and calculate Q's matrices 
	#Qs are slope tensors - mean covariances matrices 

	for(i in 1:clNumber){
		#select points that belong to the cluster i
		temp<-xyz_new[which(i==ind),]
				

		tempPar<- intamap:::estimateAnisotropySc(temp[,1],temp[,2],temp[,3])
		
		
		
		Qs[i,]<-tempPar$Q
		
		
		#calculate area of each cluster
		###########################################
		cdat<-chull(temp[,1],temp[,2])

		cdat<-SpatialPoints(data.frame(x=temp[cdat,1],y=temp[cdat,2]))
#		proj4string(cdat)=CRS("+init=epsg:3035")
#		cdat<-spTransform(cdat,CRS("+init=epsg:4326"))
		convHull<-cbind(coordinates(cdat))
		convHull<-rbind(convHull,convHull[1,]) #close Polygon
		
#		conv_Poly=Polygon(as.matrix(convHull))
#		ps=appendPolys(NULL,convHull,1,1,FALSE)
#
#		attr(ps,"projection")="LL"
#		psUTM = convUL(ps)
#		polygonArea=calcArea(psUTM,rollup=1)
#		area[i]=polygonArea$area 

		area[i]=area(convHull)	
		###########################################
	

		#store all results together
		results[i,]=c(i,tempPar$ratio,tempPar$direction,length(which(ind==i)),area[i])		

		
			
	}	

		
		#here we use Q matrices to create a weighted mean for Qs 
		#with respect to the area of each clusters chull
		Qs<-apply(kronecker(matrix(1,1,3),area)*Qs,2,sum)/sum(area)
		
		#the weighted Q matrices (slope tensors)
		Q11=Qs[1]
		Q22=Qs[2]
		Q12=Qs[3]
		
		if(Q11==0){
			return (list(R=1,theta.deg=0))
		}else{


		#calculate the weighted correlation lengths 
  		Zdiag<-Q22/Q11
  		Zoff<-Q12/Q11
  		theta<-0.5*atan(2*Zoff/(1-Zdiag))

   		R<-sqrt(1+((1-Zdiag)/(Zdiag-(1+Zdiag)*(cos(theta))^2)))
	
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

	
		singleCl=list(ratio=R,direction=theta.deg)

		object$anisPar=singleCl
		object$anisPar$clusters=results



		#must return object
		#return(list(clParams=resutls,singleCl=singleCl))
		return(object)



	}else{

#		printf("Wrong  Choise for anisotropy")
		return(-1)

	}




}
area <- function(x, y)
{
        if(missing(y)) {
                if(is.matrix(x) && ncol(x) == 2) {
                        y <- x[, 2]
                        x <- x[, 1]
                }
                else if(!is.null(x$x) && !is.null(x$y)) {
                        y <- x$y
                        x <- x$x
                }
        }
        x <- c(x, x[1])
        y <- c(y, y[1])
        i <- 2:length(x)
        return(0.5 * sum(x[i] * y[i - 1] - x[i - 1] * y[i]))
}

#
rmDist<-function(data){
		std.x<-sd(data[,2])/2
		std.y<-sd(data[,1])/2

		dd<-sqrt(std.x^2+std.y^2)

	
		#find which points are too distant 
		x.ind=which(data[,2]>mean(data[,2])+4*sd(data[,2]) |data[,2]< mean(data[,2])-4*sd(data[,2])  )
		y.ind=which(data[,1]>mean(data[,1])+4*sd(data[,1])|data[,1]< mean(data[,1])-4*sd(data[,1]))
		ind=union(x.ind , y.ind)
	
		#check if these points don't have any point near  
		if(length(ind)>2){
			dat<-as.matrix(data[ind,])
						
			xy=as.matrix(dat[,2]+1i*dat[,1])
			xy<-matrix(xy,length(xy),length(xy))
	
			distances=abs(xy-t(xy))
			index<-which(distances < dd & distances >0,arr.ind=TRUE) 
			if(length(index)>0){
				index<-matrix(index,,2)
				index=union(index[,1],index[,2])
				ind=ind[-index]
			}

		}

		return(ind)
}


