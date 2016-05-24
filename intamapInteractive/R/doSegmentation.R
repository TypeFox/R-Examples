doSegmentation<-function(object){
#	ptm<-proc.time()
	
	print("Clustering - Segmentation")
	
  xy<-coordinates(object$observations)
	dd<-cbind(xy[,1:2],object$observations$value)

	if(dim(xy)[1]<200){
		object$clusters=list(
											index=matrix(1,dim(xy)[1],1),
											clusterNumber=1,
											rmdist=NULL
										)
		return(object)
	}

	rmdist=TRUE
	#Step 0 remove distant points
	if(rmdist==TRUE){	
			ind<-rmDist(dd)	
			if(length(ind)>0)	dd<-dd[-ind,]	
	}

	object$clusters<-segmentData(ddd=dd,pl=FALSE,dev=FALSE,br=object$predictionLocations)

	index<-object$clusters$index 
	

	#create the index that corresponds to the original dataset.
	final_index<-matrix(0,length(dd[,2]),1)
	
	if(length(ind)>0){final_index[-ind]=index
		}else{final_index=index}

	object$clusters$index<-final_index


	object$clusters$rmdist<-ind


	#return single a value for all Europe.

#	print("Finished in : ")
#	print(proc.time()-ptm)
	return(object)
}
################################################################################
#Description : This is the main segmentation function if used outside 
#	 Intamap. Calls sequentially the steps described in doSegmentation
#	 help file. 
#
#Input:
#	ddd: A nx2 or nx3  matrix with the observations. First column should 
#		 the horizontal and the second one the vertical coordinate.
#	pl : boolean that toggles plotting. If TRUE several plots during the  
#		 segmentation procedure are shown. (Only used interactively)
#	dev: boolean that toggles saving the plots. (Only used interactively)
#	
#	soft: parameter that controls the weight of sampling density and 
#		distance from neighbouring clusters.
#	br: mx2 matrix with boorders coordinates in case of plotting. This 
#		this parameter is actually deactivated in line 79. 
#
#Output: A list with the following elements
#	index: a nx1 array with the indices.
#	clusterNumber: The total number of clusters detected
#
#################################################################################
segmentData<-function(ddd,pl=FALSE,dev=FALSE,soft=0.2,br){
	plot.borders=TRUE
	
#	if(pl==TRUE) {	
#		par(ask=FALSE)
#		}else {dev=FALSE}
	if(missing(br)==TRUE){	
		plot.borders=FALSE
	#	br<-borders()
	#	br<-long2INSP(br)
	} 

	#Only for bench 
  	ptm<-proc.time()
    	tim=0         		#toggle print times
    	cl_min_size=25




		
#############################################
#Step 1 density matrix
#############################################
	dat<-list(den=calculateDensity(ddd[,1],ddd[,2],1,4,pl=FALSE,dev),dat=ddd)
		
    #############################
#    	if(tim==1){
#     	 print("Step1")
#     	 tmp<-(proc.time() - ptm)
#     	 print(tmp)
#     	 }
    ############################


	Ng<-dat$den$Ng
	N<-dat$den$N
	stepx<-dat$den$stepx
	stepy<-dat$den$stepy
	xx<-dat$dat[,1]
	yy<-dat$dat[,2]
	d.mtx<-dat$den$density_L
	
	x<-dim(d.mtx)[1]
	y<-dim(d.mtx)[2]
##############################################
#Step 2 Find edges
##############################################
     edg<-detectEdges(d.mtx)
	
#	if(pl==TRUE){	
#	image(edg,col=terrain.colors(12))	
#		if(dev==TRUE){
#			dev.copy2eps(file="edges.eps")
#			}
#	}     	

    ############################
#    if(tim==1){
#      print("Step2")
#      print(proc.time() - tmp -ptm)
#      tmp<-proc.time() - tmp -ptm
#      }
    ############################

#Step 3    Separate the closed perimeters and label them
# inputs:  edg        :a (0,1) matrix indicating which points are edges
#          Ng         :Ng number
# output:  cl$ll      :a list with cluster indentified
#          cl$number  :the number of cluster indentfied
	cl<-findPerimeters(edg,Ng,pl=FALSE,dev)
#if no or just one cluster was identified index all clusters as 1.
	if(cl$number==0 | cl$number==1){
		return(list(index=matrix(1,length(ddd[,1]),1),clusterNumber=1))
	}

    ############################
#    if(tim==1){
#      print("Step3")
#      print(proc.time() - tmp -ptm)
#     tmp<-proc.time() - tmp  -ptm
#      }
    ############################

########################################################
#STEP 4: find internal to the perimeters points and reject clusters that are too small
########################################################

    st4<-fillPerimeters(xx,yy,Ng,N,cl,cl_min_size, dat$den$mav_dens_xy,stepx,stepy,pl=FALSE,dev)      # return(list(clnum=clnum,clust_dens=clust_dens,Mx=Mx,My=My,id=id))

      ############################
#    if(tim==1){
#      print("Step4")
#      print(proc.time() - tmp -ptm)
#      tmp<-proc.time() - tmp -ptm
#      }
    ############################


######################################################
#STEP 5 :: Label points that have not yet been labeled
######################################################

    out<-assignRestToClusters(ddd[,1], ddd[,2], st4$id, st4$clnum, st4$Mx, st4$My,dat$den$mav_dens_xy, st4$clust_dens,soft)
    cc<-cbind(xx,yy)

#    if(tim==1){
#      print("Step5")
#      print(proc.time() - tmp -ptm)
#      tmp<-proc.time() - tmp   -ptm
#      }

	if(tim==1){
    	print("Total")
    	print(proc.time() - ptm)
	}
	
	
	#plot the result
	if(pl==TRUE){
		if (plot.borders){
			plot(br,pch=".")
		
			for (temp in 1:st4$clnum){
				Cluster=ddd[which(out==temp),1:2]
				points(Cluster[,1:2],col=temp,pch="*")
				mtext(line=(temp-1)%%4,col=temp+1,paste("cluster",temp),adj=round(temp/8))
				#mtext(". =borders ")
			}
		}else{

				plot(ddd[,1:2],pch="*")

				for (temp in 1:st4$clnum){
				Cluster=ddd[which(out==temp),1:2]
				points(Cluster[,1:2],col=temp,pch="*")

				mtext(line=(temp-1)%%4,col=temp,paste("cluster",temp),adj=round(temp/8))
				#mtext(". =borders ")
				}


		}
	}
	
	
	#return(list(st4=st4,out=out,cc=cc,dat=dat))
	return(list(index=out,clusterNumber=st4$clnum))
	#return(out)
}

#########################INTERNAL FUNCTIONS#############
###################################################
#STEP 1
calculateDensity<-function(x,y,gd,L,pl,dev){
	

	xy<-cbind(x,y)
	N<-dim(xy)[1]
	d<-dim(xy)[2]

	#calculating the grid step
	stepx<-diff(range(x)) / (round(sqrt(gd*N))-1)*1
	stepy<-diff(range(y)) / (round(sqrt(gd*N))-1)*1

	
	#creating the grid
	xg<-seq(min(x),max(x),by=stepx)
	yg<-seq(min(y),max(y),by=stepy)

	xg<-matrix(xg,length(xg),length(xg))
	yg<-matrix(yg,length(yg),length(yg))

	

	#expand grid by 2 with zero filling on each dimension in order to have space
	#for use image filters
	xg<-rbind(matrix(0,2,dim(xg)[2]),xg,matrix(0,2,dim(xg)[2]))
	xg<-t(cbind(matrix(0,dim(xg)[1],2),xg,matrix(0,dim(xg)[1],2)))

	yg<-rbind(matrix(0,2,dim(yg)[2]),yg,matrix(0,2,dim(yg)[2]))
  	yg<-cbind(matrix(0,dim(yg)[1],2),yg,matrix(0,dim(yg)[1],2))


	#initiate variables
	Ng=dim(xg)[1]
	denst<-matrix(0,Ng,Ng)
	inn=vector("list",N)
	#the radius of the square is going to be equal to the grid size, cetred
	#on each node. Therefore, we are looking for points that are +-stepx/2
	#and +-stepy/2 apart.
	#density_xy=zeros(N, 1);

	for( i in 1:Ng){
		for (j in 1:Ng){
			ind<-indfinal<-NULL
			ind<-which(x<xg[i,j]+stepx/2 & x>xg[i,j]-stepx/2)
			if(length(ind)==0){		denst[i,j]=0
				}else{
				indfinal<-which(y[ind]<yg[i,j]+stepy/2 & y[ind]> yg[i,j]-stepy/2)
				denst[i,j]=length(indfinal)


        			}

		}
	}

	#attribute a density to the initial points
	density_xy<-matrix(0,1,N)
	dens_xy<-NULL
 	index<-vector("list",N)

	for (st in 1:N){
		ind=which(x[st]-x>0 & abs(x[st]-x)<stepx/2)
		indfinal=which(y[st]-y[ind]>=0 & abs(y[st]-y[ind])<stepy/2)
		dens_xy[st]=length(indfinal)+1
		inn[[st]]<-c(indfinal)
	}
	density_L=t(denst)
	mav_dens_xy=dens_xy


	#moving average filter, (smoothing)
	density_L<-applyFilter(density_L,"average")

	
	#Plot the Density
	if(pl==TRUE){
	filled.contour(density_L,color.palette=terrain.colors)
		if(dev==TRUE){
			dev.copy2eps(file="smoothed")
		}	
	}
	return(list(mav_dens_xy=mav_dens_xy,density_L=density_L,Ng=Ng,stepx=stepx,stepy=stepy,N=N,inn=inn))
}





###############################################################################
#STEP 2
#This function applies filters on images
###############################################################################
applyFilter<-function(dat,input,TH){


######INTERNAL FUNCTION ###################
#Create Laplacian of gaussian filter
######################################
	 create.log<-function(){
	   eps=2.2204e-016
	    #first we calculate gaussian
	    g=(-1:1)
	    x=outer(g*0,g,"+")
	    y=outer(-g,g*0,"+")
	    sigma=0.5

	    arg = -(x*x + y*y)/(2*(sigma^2))
	    hh=exp(arg)
	    hh[(hh<eps*max(hh))]=0

	    sumh=sum(hh)
	      if(sumh!=0){
		hh=hh/sumh
	      }
	      #now we calculate laplacian
	      h1=hh*(x*x+y*y - 2*(sigma^2))/(sigma^2)^2
	      h2=h1-sum(h1)/prod(5,5)

	      #sum to zero

	      h2=h2-sum(h2)/length(h2)
	      return(h2)
	}


 #mask selection
	mask<-switch(input,
               "laplace" =matrix(c(0,-1,0,-1,4,-1,0,-1,0),3,3),
               "average"=(matrix(c(1,1,1,1,1,1,1,1,1),3,3)/9),
               "average22"=(matrix(c(-1,-1,-1,-1,9,-1,-1,-1,-1),3,3)),
               # "rows"=matrix(c(1,sqrt(2),1,0,0,0,-1,-sqrt(2),-1),3,3),
               # "columns"=matrix( c( -1 ,0 ,1 ,-sqrt(2) ,0 ,sqrt(2), -1, 0, 1),3,3)
               "rows"=matrix(c(1,1,1,0,0,0,-1,-1,-1),3,3),
               "columns"=matrix( c( -1 ,0 ,1 ,-1 ,0 ,1, -1, 0, 1),3,3),
               "average5"=matrix(1/25,5,5),
               "log"=create.log()

                )

	#Threshold initialization
	#TH<-1
	x<-dim(dat)[1]
	y<-dim(dat)[2]

	ll<-(dim(mask)[2]-1)/2
	zerosx<-matrix(0,ll,y)
	zerosy<-matrix(0,(x+2*ll),ll)
	dat<-rbind(zerosx,dat,zerosx)
	dat<-cbind(zerosy,dat,zerosy)

	rm(zerosx,zerosy)
	out<-dat

	for (i in (1+ll):(x+ll)){
	    for(j in (1+ll):(y+ll)){

		#out[i,j]<-{if(sum(mask * dat[(i-ll):(i+ll),(j-ll):(j+ll)])>TH) 1 else 0 }
		out[i,j]<-sum(mask*dat[(i-ll):(i+ll),(j-ll):(j+ll)])
	    }
	 }
	    return(out[(1+ll):(x+ll),(1+ll):(y+ll)])
	}


detectEdges<-function(dat){
	
	thresh=0
	m<-dim(dat)[1]
	n<-dim(dat)[2]
	
	rr<-2:(m-1)
	cc<-2:(n-1)
	edg<-matrix(0,m,n)
 
 	dat<-applyFilter(dat,"average",0)

	out<-applyFilter(dat,"log",0)
  
	ed=which((out[rr,cc]<0 & out[rr,cc+1]>0) & abs(out[rr,cc]-out[rr,cc+1])>thresh,arr.ind=TRUE)
	edg[ed+1]=1

	ed=which((out[rr,(cc-1)]>0 & out[rr,cc]<0) & abs(out[rr,(cc-1)]-out[rr,cc])>thresh,arr.ind=TRUE)
	edg[ed+1]=1

	ed=which((out[rr,cc]<0 & out[rr+1,cc]>0) & abs(out[rr,cc]-out[rr+1,cc])>thresh,arr.ind=TRUE)
	edg[ed+1]=1

	ed=which((out[(rr-1),cc]>0 & out[rr,cc]<0) & abs(out[(rr-1),cc]-out[rr,cc])>thresh,arr.ind=TRUE)
	edg[ed+1]=1

	ed=which(out[rr,cc]==0,arr.ind=TRUE)
	
	if(length(ed)>0){
		zeros=ed
		zz=which(out[zeros-c(-1,0)]<0 & out[zeros+c(2,1)]>0 & abs (out[zeros-c(-1,0)]-out[zeros +c(2,1)])>2*thresh ,arr.ind=TRUE)
		edg[zz+c(1,1)]=1

		zz=which(out[zeros-c(-1,0)]>0 & out[zeros+c(2,1)]<0 & abs (out[zeros-c(-1,0)]-out[zeros +c(2,1)])>2*thresh ,arr.ind=TRUE)
		edg[zz+c(1,1)]=1

		zz=which(out[zeros-c(0,-1)]<0 & out[zeros+c(1,2)]>0 & abs (out[zeros-c(0,-1)]-out[zeros +c(1,2)])>2*thresh ,arr.ind=TRUE)
		edg[zz+c(1,1)]=1

		zz=which(out[zeros-c(0,-1)]>0 & out[zeros+c(1,2)]<0 & abs (out[zeros-c(0,-1)]-out[zeros +c(1,2)])>2*thresh ,arr.ind=TRUE)
		edg[zz+c(1,1)]=1
	     }

       return(edg)
}


################################################################################
#Step 3
################################################################################
findPerimeters<-function(edg,Ng,pl=FALSE,dev=FALSE){

########################INTERNAL FUNCTIONS TO findPerimeters############
find_neighb<-function(cx,cy,jj,Ng){
        nx<- (((cx[jj]-cx)==0 )| ((cx[jj]-cx)==1) | (cx[jj] -cx==-1))
        ny<- (((cy[jj]-cy)==0 )| ((cy[jj]-cy)==1) | (cy[jj]-cy==-1))

        same_point<-((cx[jj]-cx)==0 & (cy[jj]-cy)==0)

        edgex=(( (cx[jj]==1)|(cx[jj]==2) )|((cx[jj]==Ng) | (cx[jj]==(Ng-1) ) ))
        edgey=(( (cy[jj]==1)|(cy[jj]==2) )|((cy[jj]==Ng) | (cy[jj]==(Ng-1) ) ))

	# if we are at an edge point we will continue searching even with a
	#single neighbour
         np_n<-xor((nx & ny),same_point)
         np_n<-which(np_n==TRUE)
         edge<-(edgex | edgey )
         return(list(np_n=np_n,edge=edge))

}


#we will search and remove dublicate points
rm.duplicate<-function(input){
	 L=length(input)
	 for(i in 1:L){
	   if((input[i]!=-1)){
	      input[(which((input[(i+1):L]-input[i])==0)+i)]=-1
	   }
	 }
	  return(input[which(input!=-1)])
}

#function that removes  old points 
rm.old<-function(input){
	 L=length(input)
	 for(i in 1:L){
	   if((input[i]!=-1)){
	      input[(which((input[(i+1):L]-input[i])==0)+i)]=-1
	   }
	 }
	  return(input[which(input!=-1)])
}

#plot.point<-function(dat,jj){
# plot(dat)
# points(dat[jj,1],dat[jj,2],col=2)
#}
###################END OF INTERNAL FUNCTIONS########################


#	if(missing(dev)){dev=FALSE}
#	if(missing(pl)){pl=FALSE}

	indices<-which(edg==1,arr.ind=TRUE)
	len<-dim(indices)[1]

	#here are the coords of edges
	cx<-indices[,1]
	cy<-indices[,2]

	cl<-matrix(0,length(cx),1)


	continuee=1     #flag indicating that the search for the next region is NOT over
	cl_counter=0    #temporary counter of the number of clusters


    while(continuee==1){

        jj=1                        		#start at the first point
        cl_counter=cl_counter+1    		#increase the number of clusters: the number of clusters at this point equals the number
        cl_ind=NULL                		#the indices of the perimeter points that belong to the cluster (the region perimeter in other words)
        fn=find_neighb(cx, cy, jj, Ng)          #np_n is the number of neighbours of the point under investigation, edgee is a flag indicating
						#whether the point under investigation is an edge point (in an open curve which normally would not occur, however just in case)
        

	if(pl==TRUE){
        	#plot.point(cbind(cx,cy),jj)                 #plots the edge points			#PLOTTING
		dat=cbind(cx,cy)
		plot(dat)
		points(dat[jj,1],dat[jj,2],col=2)
		points(cx[fn$np_n],cy[fn$np_n],col=3)							#PLOTTING
	       	}
	
	 if (length(fn$np_n)>1){                    #then we don't have a single point
            #find all points of the cluster
                cl_ind=c(cl_ind ,t(fn$np_n))
                cl_ind_n=cl_ind                    		#temporary variable
                cluster_end=0                       		#flag indicating that the cluster perimeter is not yet finalized
                size_cl=length(cl_ind)            		#cluster size
                tmp=0                               		#old cluster length
                while (cluster_end==0){
                      size_cl=length(cl_ind)
                          for (k in (1+tmp):length(cl_ind_n)){
                              nfn=find_neighb(cx, cy, cl_ind_n[k], Ng)      #find all neighbouring points     returns->(np_n, edgee)
                              cl_ind_n=c(cl_ind_n, t(nfn$np_n))            	#add new neighbouring points to the cluster perimeter
                              cl_ind_n=rm.duplicate(cl_ind_n)
       	                      if(pl==TRUE){
				points(cx[cl_ind_n],cy[cl_ind_n],col=3)     #visualize data selection		#PLOTTING
                          	}
			}

                    #%We might at this point have introduced a point twice in the perimeter because it has
                    #%multiple neighboring points. So we need to remove duplicate points
                    cl_ind_n=rm.duplicate(cl_ind_n)

                    #keep the old length in order to search only in new points
                    tmp<-length(cl_ind)

                    if ((length(cl_ind)-length(cl_ind_n))==0){
                             cluster_end=1                          #when there are no more neighbours, the perimeter is finalised
                     }
                     cl_ind=cl_ind_n
                }


                #create new variable cluster
                assign(paste('cluster',cl_counter,sep=""),cbind(cx[cl_ind], cy[cl_ind]))

                #remove points already assigned to a cluster
                cx=cx[-cl_ind]
                cy=cy[-cl_ind]
                include_edge_points=1
                np=0
                 if (length(cx)==0){
                      continuee=0;
                 }

         }else{
         cluster_end=1;
                cl_ind=jj;


                #assign points to a cluster variable
                assign(paste('cluster',cl_counter,sep=""),cbind(cx[cl_ind], cy[cl_ind]))

                #remove points already assigned to a cluster
                cx=cx[-cl_ind]
                cy=cy[-cl_ind]
                np=0;
                if (length(cx)==0){
                      continuee=0;
                 }
        }
    }#end of while

    #plot(indices)
    ll=vector("list",cl_counter)
    for (i in 1:cl_counter){
      cl<-get(paste("cluster",i,sep=""))
      if(pl==TRUE){
	points(cl,col=(i+1))						#PLOTTING
	}
      ll[[i]]=cl;
    }
	
#	if(dev==TRUE){dev.copy2eps(file="Edges.eps")	}
    return(list(ll=ll,number=cl_counter))
}


###################################################
#Step 4
###################################################
fillPerimeters<-function(x,y, Ng, N, clust , cl_min_size, mav_dens_xy, stepx, stepy,pl=FALSE,dev=FALSE){

###################Internal Functions to fill Perimeters##########
###################################################
#The cluster perimeter is defined from the points xx and yy.
#We want to identyify wether the point (x,y) is inside this perimeter.
###################################################
isinside_cluster<-function(xx,yy,x,y,stepx,stepy){

	if( ( (x>(max(xx)+stepx/2)) | (x<(min(xx)-stepx/2)) | (y> (max(yy)+stepy/2)) | (y<(min(yy)-stepy/2) ))==TRUE){
		not_inside=1
		}else{
			ind=which(y<yy+stepy/2 & y>yy-stepy/2)
			rx=xx[ind]
			ry=yy[ind]
      	if(length(ind)==0){
      		not_inside=1
      		}else{
      		   if(length(which(x>rx))==0  | length(which(x<rx))==0){
      		    not_inside=1
      		    }else{
      			   not_inside=0
      		    }
          }
        }
    return(not_inside)
}
######################END IF INTERNAL FUNCTIONS############


	if(pl==TRUE){
		plot(x,y,pch="*")							#PLOTTING
	}
	
	ax= (max(x)-min(x)) / (Ng-1)
	ay= (max(y)-min(y)) / (Ng-1)

	b.x=min(x)-ax
	b.y=min(y)-ay
	
	not_inside=1
	id=matrix(0,N,1)

	for(b in 1:N){
		for(fff in 1:clust$number){
			xxx=clust$ll[[fff]][,1]
			yyy=clust$ll[[fff]][,2]
	
			#Decide if a point is inside a perimeter or not.
			not_inside=isinside_cluster(ax*xxx+b.x,ay*yyy+b.y,x[b],y[b],stepx,stepy)
			if(not_inside==0){
				id[b]=fff
				not_inside=1
				break
			}
		}
	}

	if(pl==TRUE){
		for(b in 1:clust$number){
			xxx=clust$ll[[b]][,1]
			yyy=clust$ll[[b]][,2]
			
			points(x[id==b],y[id==b],col=b+3,pch="*")					#PLOTTING
			points(ax*xxx+b.x,ay*yyy+b.y,pch="+",col=3)			#PLOTTING
		}
	}



	#reject cluster with too few points
	clnum=0
	clust_dens=NULL
	Mx=NULL
	My=NULL
	
	for(ff in 1:clust$number){
		if(length(x[which(id==ff)])<cl_min_size){
			id[which(id==ff)]=0
		}else{
			clnum=clnum+1
			temp=which(id==ff)
			id[temp]=clnum

			xxx=clust$ll[[ff]][,1]
			yyy=clust$ll[[ff]][,2]

		      assign(paste("final_cluster",clnum,sep=""),cbind(ax*xxx+b.x,ay*yyy+b.y))
	
			#We evaluate the cluster centres as the centroids of the
			#clusters
			Mx[clnum]=sum(x[temp]*mav_dens_xy[temp])/sum(mav_dens_xy[temp]);
			My[clnum]=sum(y[temp]*mav_dens_xy[temp])/sum(mav_dens_xy[temp]);

			clust_dens[clnum]=mean(mav_dens_xy[temp])

			if(is.na(Mx[clnum])|is.na(My[clnum])){
				Mx=Mx[-clnum]
				My=My[-clnum]
				clust_dens=clust_dens[-clnum]
				id[which(id==ff)]=0
				#	rm(paste("final_cluster",clnum,sep=""))
				clnum=clnum-1	
			}
		}
	}
	
	if(dev==TRUE){dev.copy2eps(file="fillPerimeters.eps")	}
	return(list(clnum=clnum,clust_dens=clust_dens,Mx=Mx,My=My,id=id))
}

assignRestToClusters<-function(x, y, id, clnum, Mx, My, mav_dens_xy, clust_dens, soft){
	
	#INTERNAL function calculates cost
	fcost<-function(x,y,Mx,My,soft,mav,clust,maxden,maxdist){
		cost= (sqrt((x-Mx)^2+(y-My)^2))/maxdist +soft*(sqrt(abs(mav-clust))/maxdens)
		return(cost)
	}

  	N<-length(x)
  	index<-matrix(0,N,1)
	

	#in order to avoid dividing by zero
	if(length(which(clust_dens!=0)==0)){	
		clust_dens=clust_dens+1	
	}

	for(i in 1:N){
   		if(id[i]>0){
       			index[i]=id[i]
   		}else{

		#choose the 2 closest
		ind2<-which(rank(sqrt((x[i]-Mx)^2+(y[i]-My)^2))<2.6)
			
			if(length(ind2)!=0){
				maxdist=max(sqrt((x[i]-Mx[ind2])^2+(y[i]-My[ind2])^2))+0.000000001
				maxdens=max(sqrt(abs(mav_dens_xy[i]-clust_dens[ind2])))+0.000000001
		      
				#calculate clusters' costs and assign point to the minimum one
				cost=fcost(x[i],y[i],Mx[ind2],My[ind2],soft,mav_dens_xy[i],clust_dens[ind2],maxdens,maxdist)
				index[i]=ind2[which.min(cost)]
				}
					else{
						index[i]=which.min(rank(sqrt((x[i]-Mx)^2+(y[i]-My)^2)))
					}
			}
  	}
  return(index)

}#

#borders<-function(){
#	library(maptools)
#	eur<-read.shape('countryLim.shp')	
#	brd<-NULL
#	iter=length(eur$Shapes)
#	
#	for(i in 1:iter){
#	brd<-rbind(brd,eur$Shapes[[i]]$verts)
#	}
#	return(brd)
#}
	
#################################################################################
#function	:rmDist								#
#Description	:Removes distant points from the original dataset		#
#		 in order to increase edge detection accuracy			#
#		 the removed points will not be assigned in step5		#
#################################################################################
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


