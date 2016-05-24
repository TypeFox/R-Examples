freqmat <-
function(r1,r2,r=NULL,oneside=NULL,twoside=NULL,cross = NULL,homo=NULL,...){	
	if( (is.null(twoside))&(is.null(oneside)) ){
		twoside = TRUE
		oneside = FALSE
	}
	if( !(is.null(twoside)) & (is.null(oneside)) ) oneside = !twoside 
	if( (is.null(twoside)) & !(is.null(oneside)) ) twoside = !oneside
	if( (twoside==TRUE & oneside==TRUE) || (twoside==FALSE & oneside==FALSE))stop("'twoside' and 'onesied' should be different logic value!")	
	if(is.null(cross))cross = FALSE	
	if(is.null(homo))homo = TRUE
	if(!cross) r12 = 0
	if(cross) r12 = r1*r2
	if(is.null(r))r = r1 + r2 - 2*r12
	if(r > (r1+r2) ) stop("r should no bigger than r1+r2")
	if(homo & twoside){
		A = matrix(0,4,3)
      		A[1,1] = A[4,1] = 0.5*(1-r)
		A[2,1] = A[3,1] = 0.5*r	
		A[1,2] = A[4,3] = 0.5*(1-r1-r2+r12)
		A[2,2] = A[3,3] = 0.5*(r2-r12)
		A[3,2] = A[2,3] = 0.5*(r1-r12)
		A[4,2] = A[1,3] = 0.5*r12
		colnames(A) = c("Freq","TT","tt")
		rownames(A) = c("AABB","AAbb","aaBB","aabb")
		B = A[,2:3]/A[,1]
		id = which(A[,1]==0)
		if(length(id)>=1)B[id,]=0
		colnames(B) = c("TT","tt")
		rownames(B) = c("AABB","AAbb","aaBB","aabb")
	}
	if(!homo & twoside){
		A = matrix(0,9,4)
		A[1,1] = A[9,1] = 0.25*(1-r)*(1-r)
		A[2,1] = A[8,1] = A[4,1] = A[6,1] = 0.5*r*(1-r)
		A[3,1] = A[7,1] = 0.25*r*r
		A[5,1] = 0.5*((1-r)*(1-r)+r*r) 		
		A[1,2] = A[9,4] = 0.25*(1-r1-r2+r12)*(1-r1-r2+r12)
		A[1,3] = A[9,3] = 0.5*(1-r1-r2+r12)*r12
		A[1,4] = A[9,2] = 0.25*r12*r12
		A[2,2] = A[8,4] = 0.5*(r2-r12)*(1-r1-r2+r12)
		A[2,3] = A[8,3] = 0.5*((r1-r12)*(1-r1-r2+r12)+(r2-r12)*r12)
		A[2,4] = A[8,2] = 0.5*(r1-r12)*r12
		A[3,2] = A[7,4] = 0.25*(r2-r12)*(r2-r12)
		A[3,3] = A[7,3] = 0.5*(r1-r12)*(r2-r12)
		A[3,4] = A[7,2] = 0.25*(r1-r12)*(r1-r12)
		A[4,2] = A[6,4] = 0.5*(r1-r12)*(1-r1-r2+r12)
		A[4,3] = A[6,3] = 0.5*((r2-r12)*(1-r1-r2+r12)+(r2-r12)*r12)
		A[4,4] = A[6,2] = 0.5*(r2-r12)*r12
		A[5,2] = A[5,4] = 0.5*(r1-r12)*(r2-r12)+0.5*(1-r1-r2+r12)*r12
		A[5,3] = 0.5*((1-r1-r2+r12)*(1-r1-r2+r12)+(r1-r12)*(r1-r12)+(r2-r12)*(r2-r12)+r12*r12)
		colnames(A) = c("Freq","TT","Tt","tt")
		rownames(A) = c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb")		
		B = A[,2:4]/A[,1]
		id = which(A[,1]==0)
		if(length(id)>=1)B[id,]=0 
		colnames(B) = c("TT","Tt","tt")
		rownames(B) = c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb")
	}
	if(homo & oneside){
		A = matrix(0,4,3)
      		A[,1] = 0.5	
		A[1,2] = A[2,3] = 0.5*(1-r1)
		A[2,2] = A[1,3] = 0.5*r1
		A[3,2] = A[4,3] = 0.5*(1-r2)
		A[4,2] = A[3,3] = 0.5*r2
		colnames(A) = c("Freq","TT","tt")
		rownames(A) = c("AA__","aa__","__BB","__bb")
		B = A[,2:3]/A[,1]
		colnames(B) = c("TT","tt")
		rownames(B) = c("AA__","aa__","__BB","__bb")
	}
	if(!homo & oneside){
		A = matrix(0,6,4)
      		A[1,1] = A[3,1] = A[4,1] = A[6,1] = 0.25
		A[2,1] = A[5,1] = 0.5
		A[1,2] = A[3,4] = 0.25*(1-r1)*(1-r1)
		A[1,3] = A[3,3] = A[2,2] = A[2,4] =0.5*(1-r1)*r1
		A[1,4] = A[3,2] = 0.25*r1*r1
		A[2,3] = 0.5*((1-r1)*(1-r1)+r1*r1)
		A[4,2] = A[6,4] = 0.5*(1-r2)*(1-r2)
		A[4,3] = A[6,3] = A[5,2] = A[5,4] =0.5*(1-r2)*r2
		A[5,3] = 0.5*((1-r2)*(1-r2)+r2*r2)
		colnames(A) = c("Freq","TT","Tt","tt")
		rownames(A) = c("AA__","Aa__","aa__","__BB","__Bb","__bb")		
		B = A[,2:4]/A[,1]
		colnames(B) = c("TT","Tt","tt")
		rownames(B) = c("AA__","Aa__","aa__","__BB","__Bb","__bb")
	}
	Res = list(Freq=A,Condional_Freq=B)
	return(Res)
}
