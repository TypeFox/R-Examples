PHPP<-function(reaction_number=c(28,36),fba_object,PCS="D glucose",flux_range=NULL,ret_OBJ_mat=FALSE,surf_col="red",divs=5,dimension="3",animate=FALSE,objective=NULL){

candidates<-grep(PCS,fba_object$reaction_list,ignore.case=TRUE)
message("\nProbable Carbon Sources:") 
print(fba_object$reaction_list[candidates])
C_S<-as.numeric(readline("\nSelect the Carbon source from the list above by serial number"))

if(length(which(reaction_number==candidates[C_S]))==0)
	{
	message("Alternate carbon source PhPP")
	fba_object<-CHANGE_RXN_BOUNDS(reaction_number=candidates[C_S],fba_object,lb=0,ub=0)
	}else(message("Primary carbon source PhPP"))
if(length(flux_range)==0)
{
a1<-round(Flux_Ranger(reaction_number[1],fba_object=fba_object,divs=divs)$ramp,3)
b1<-round(Flux_Ranger(reaction_number[2],fba_object=fba_object,divs=divs)$ramp,3)
}else{a1<-round(Flux_Ranger(reaction_number[1],fba_object=fba_object,divs=divs,art_limit_range=flux_range)$ramp,3)
b1<-round(Flux_Ranger(reaction_number[2],fba_object=fba_object,divs=divs,art_limit_range=flux_range)$ramp,3)
}

message("Phenotypic PhasePlane Analysis Info-:")
message("X-axis")
message(fba_object$reaction_list[reaction_number[1]]) 
print(a1)
message("Y-axis")
message(fba_object$reaction_list[reaction_number[2]])
print(b1)
OBJ_MAT=matrix(0,divs,divs)

for(i in 1:length(a1))
	{
	for(j in 1:length(b1))
		{
		temp_mod<-CHANGE_RXN_BOUNDS(reaction_number=reaction_number[1],fba_object,lb=a1[i],ub=a1[i])
		temp_mod<-CHANGE_RXN_BOUNDS(reaction_number=reaction_number[2],temp_mod,lb=b1[j],ub=b1[j])	
		temp_sol<-FBA_solve(temp_mod)		
		if(length(objective)>0)		
		{OBJ_MAT[i,j]<-temp_sol$fluxes[objective]}else{OBJ_MAT[i,j]<-temp_sol$objective}
		}
	}

if(dimension=="2")
	{
		#print("we need access")
		rownames(OBJ_MAT)=a1
		colnames(OBJ_MAT)=b1
		OBJ_MAT<-round(OBJ_MAT,3)
		print(levelplot(OBJ_MAT,contour=T,labels=T,xlab=fba_object$reaction_list[reaction_number[1]],ylab=fba_object$reaction_list[reaction_number[2]]))
	
		if(ret_OBJ_mat==TRUE){
		rownames(OBJ_MAT)=a1
		colnames(OBJ_MAT)=b1
		return(OBJ_MAT)
		}

	}

if(dimension=="3")
	{
	print("3D")
	if(length(objective)>0)
	{
	persp3d(x=a1,y=b1,z=OBJ_MAT,xlab=fba_object$reaction_list[reaction_number[1]],ylab=fba_object$reaction_list[reaction_number[2]],
			zlab=fba_object$reaction_list[objective],color=surf_col, alpha=0.95, back="lines",smooth=FALSE,main="PHPP")} else{persp3d(x=a1,y=b1,z=OBJ_MAT,xlab=fba_object$reaction_list[reaction_number[1]],ylab=fba_object$reaction_list[reaction_number[2]],
			zlab=fba_object$reaction_list[which(fba_object$obj==1)],color=surf_col, alpha=0.95, back="lines",smooth=FALSE,main="PHPP")}
		if(animate==TRUE)
		{
		play3d(spin3d(axis=c(1,0,0), rpm=10), duration=6)
		play3d(spin3d(axis=c(0,1,0), rpm=10), duration=6)
		play3d(spin3d(axis=c(0,0,1), rpm=10), duration=6)
		}

	if(ret_OBJ_mat==TRUE){
		rownames(OBJ_MAT)=a1
		colnames(OBJ_MAT)=b1
		return(OBJ_MAT)
		}


	}
}
