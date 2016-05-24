#	Q        Rotation Quaternions				Q - [q1,q2,q3,q4] (Nx4)
#	EV       Euler Vector and rotation angle (degrees)	EV - [m1,m2,m3,MU] (Nx4)
#	DCM      Orthogonal DCM Rotation Matrix			DCM - 3x3xN
#	EA    Euler angles (12 possible sets) (degrees)		EA - [psi,theta,phi] (Nx3)

# DCM2EA	DCM2EV	DCM2Q
# EA2DCM	EA2EV	EA2Q	EA2EA
# EV2DCM	EV2EA	EV2Q
# Q2DCM		Q2EA	Q2EV	Q2GL

EA2Q<-function(EA, EulerOrder='zyx', ichk=FALSE, ignoreAllChk=FALSE)
{# EA - [psi,theta,phi] yaw, pitch, roll to EV - [m1,m2,m3,MU]
# EA in radians
# ichk = FALSE disables near-singularity warnings.
# Identify singularities (2nd Euler angle out of range)
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if (!is.matrix(EA)) EA <- matrix(EA,ncol=3,byrow=FALSE)
        theta = EA[, 2] # N×1
        if (!ignoreAllChk){
        if (substr(EulerOrder,1,1) != substr(EulerOrder,3,3)) {# Type 1 rotation about three distinct axes
            # Type 1 rotation (rotations about three distinct axes)
            if (any(abs(theta) >= pi/2 )) { stop('Second input Euler angle(s) outside -90 to 90 degree range')
            } else if (ichk && any(abs(theta)>88 * (pi/180))) warning('Warning: Second input Euler angle(s) near a singularity (-90 or 90 degrees).')
        } else {
            # Type 2 rotation (1st and 3rd rotation about same axis)
            if (any((theta<=0) | (theta >= pi))) { stop('Second input Euler angle(s) outside 0 to 180 degree range')
            } else if (ichk && (any(theta < 2 * (pi/180)) | (theta>178 * (pi/180)))) warning('Warning: Second input Euler angle(s) near a singularity (0 or 180 degrees).')
        }}
        # Half angles in radians
        HALF = EA / 2# * (pi/360) # N×3
        Hpsi   = matrix(HALF[,1],ncol=1) # N×1
        Htheta = matrix(HALF[,2],ncol=1) # N×1
        Hphi   = matrix(HALF[,3],ncol=1) # N×1
        # Pre-calculate cosines and sines of the half-angles for conversion.
        c1=cos(Hpsi); c2=cos(Htheta); c3=cos(Hphi)
        s1=sin(Hpsi); s2=sin(Htheta); s3=sin(Hphi)
        c13 =cos(Hpsi+Hphi);  s13 =sin(Hpsi+Hphi)
        c1_3=cos(Hpsi-Hphi);  s1_3=sin(Hpsi-Hphi)
        c3_1=cos(Hphi-Hpsi);  s3_1=sin(Hphi-Hpsi)
        if (EulerOrder=='xyx') Q=cbind(c2*c13, c2*s13,s2*c1_3, s2*s1_3) else
        if (EulerOrder=='yzy') Q=cbind(c2*c13, s2*s1_3,c2*s13, s2*c1_3) else
        if (EulerOrder=='zxz') Q=cbind(c2*c13, s2*c1_3,s2*s1_3, c2*s13) else
        if (EulerOrder=='xzx') Q=cbind(c2*c13, c2*s13,s2*s3_1, s2*c3_1) else
        if (EulerOrder=='yxy') Q=cbind(c2*c13, s2*c3_1,c2*s13,  s2*s3_1) else
        if (EulerOrder=='zyz') Q=cbind(c2*c13, s2*s3_1,s2*c3_1, c2*s13) else
        if (EulerOrder=='xyz') Q=cbind(c1*c2*c3-s1*s2*s3, s1*c2*c3+c1*s2*s3,c1*s2*c3-s1*c2*s3, c1*c2*s3+s1*s2*c3) else
        if (EulerOrder=='yzx') Q=cbind(c1*c2*c3-s1*s2*s3, c1*c2*s3+s1*s2*c3,s1*c2*c3+c1*s2*s3, c1*s2*c3-s1*c2*s3) else
        if (EulerOrder=='zxy') Q=cbind(c1*c2*c3-s1*s2*s3, c1*s2*c3-s1*c2*s3,c1*c2*s3+s1*s2*c3, s1*c2*c3+c1*s2*s3) else
        if (EulerOrder=='xzy') Q=cbind(c1*c2*c3+s1*s2*s3, s1*c2*c3-c1*s2*s3,c1*c2*s3-s1*s2*c3, c1*s2*c3+s1*c2*s3) else
        if (EulerOrder=='yxz') Q=cbind(c1*c2*c3+s1*s2*s3, c1*s2*c3+s1*c2*s3,s1*c2*c3-c1*s2*s3, c1*c2*s3-s1*s2*c3) else
        if (EulerOrder=='zyx') Q=cbind(c1*c2*c3+s1*s2*s3, c1*c2*s3-s1*s2*c3,c1*s2*c3+s1*c2*s3, s1*c2*c3-c1*s2*s3) else
        if (!ignoreAllChk) stop('Invalid input Euler angle order')
Q
}

EV2Q <- function(EV,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EV - [m1,m2,m3,MU] to Q - [q1,q2,q3,q4]
        # Euler vector (EV) and angle MU in radians
        if(is.null(dim(EV))) EV<-matrix(EV,ncol=4,byrow=FALSE)
        EVtmp = matrix(EV[,1:3],ncol=3,byrow=FALSE) # N×3
        halfMU = matrix(EV[,4] / 2,ncol=1) #* (pi/360) (N×1) MU/2 in radians
        # Check that input m's constitute unit vector
        delta = sqrt(matrix(apply(EVtmp^2,1,sum),ncol=1)) - 1 # N×1
        if (!ignoreAllChk) if (any(abs(delta) > tol)) stop('(At least one of the) input Euler vector(s) is not a unit vector')            
        # Quaternion
        SIN = sin(halfMU) # (N×1)
        Q = cbind(EVtmp[,2]*SIN, EVtmp[,3]*SIN, cos(halfMU), EVtmp[,1]*SIN)
Q
}

DCM2Q <- function(DCM,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# DCM - 3x3xN to Q - [q1,q2,q3,q4] 
        # NOTE: Orthogonal matrixes may have determinant -1 or 1
        #       DCMs are special orthogonal matrices, with determinant 1
        improper  = FALSE
        DCM_not_1 = FALSE
        if(any(is.null(dim(DCM)))) stop('DCM must be a matrix or array.')
        if(is.na(dim(DCM)[3]))  N =1 else N = dim(DCM)[3]
        if (N == 1){
            # Computing deviation from orthogonality
            delta = DCM %*% t(DCM) - diag(3) # DCM*DCM' - I
            delta = matrix(delta,ncol=1) # 9×1 <-- 3×3
            # Checking determinant of DCM
            DET = det(DCM)
            if (DET<0) improper=TRUE
            if (ichk && (abs(DET-1)>tol)) DCM_not_1=TRUE
            # Permuting  DCM
             DCM = array( DCM, c(1, 3, 3)) # 1×3×3
        } else {
        delta <- array(0, dim(DCM))
        d=lapply(1:N, function(cntDCM) delta[,,cntDCM] <- t(matrix(DCM[,,cntDCM],3,3) %*% t(matrix(DCM[,,cntDCM],3,3)) - diag(3)))
        delta <- unlist(d)
        delta <- c(delta)
DET = DCM[1,1,]*DCM[2,2,]*DCM[3,3,] -DCM[1,1,]*DCM[2,3,]*DCM[3,2,]+
DCM[1,2,]*DCM[2,3,]*DCM[3,1,] -DCM[1,2,]*DCM[2,1,]*DCM[3,3,]+
DCM[1,3,]*DCM[2,1,]*DCM[3,2,] -DCM[1,3,]*DCM[2,2,]*DCM[3,1,]
DET = array(DET,c(1,1,N)) # 1×1×N
        if (any(DET<0)) improper=TRUE
        if (ichk && (any(abs(DET-1)>tol))) DCM_not_1=TRUE 
DCMtmp<-lapply(1:N, function(n) t(DCM[,n,]) ) # it works! but must be turned into an array?
DCMtmp<-array(unlist(DCMtmp), dim = c(dim(DCMtmp[[1]]), length(DCMtmp))) 
DCM <- DCMtmp
        }
        # Issuing error messages or warnings
        if (!ignoreAllChk) if (ichk && any(abs(delta)>tol)) warning('Warning: Input DCM is not orthogonal.')
        if (!ignoreAllChk) if (improper) stop('Improper input DCM')
        if (!ignoreAllChk) if (DCM_not_1) warning('Warning: Input DCM determinant off from 1 by more than tolerance.')
        # Denominators for 4 distinct types of equivalent Q equations
denom = cbind(1 +  DCM[,1,1] -  DCM[,2,2] -  DCM[,3,3], 1 -  DCM[,1,1] +  DCM[,2,2] -  DCM[,3,3], 
1 -  DCM[,1,1] -  DCM[,2,2] +  DCM[,3,3], 1 +  DCM[,1,1] +  DCM[,2,2] +  DCM[,3,3])
#denom[which(is.na(denom))]<-0
denom[which(denom<0)]<-0
denom = 2 * sqrt(denom) # N×4
# Choosing for each DCM the equation which uses largest denominator
maxdenom <- apply(denom,1,max)
index <- apply(denom,1,function(x) which(x==max(x)))
if(is.null(dim(maxdenom))) maxdenom <- matrix(maxdenom,ncol=1)    
Q = matrix(NA,N,4) # N×4
if (N==1){
  ii=1
if (index==1) Q <- (cbind( (DCM[ii,2,3]- DCM[ii,3,2]) / maxdenom, 0.25 * maxdenom, ( DCM[ii,1,2]+ DCM[ii,2,1]) / maxdenom,( DCM[ii,1,3]+ DCM[ii,3,1]) / maxdenom) )
if (index==2) Q <- (cbind( (DCM[ii,3,1]- DCM[ii,1,3]) / maxdenom,( DCM[ii,1,2]+ DCM[ii,2,1]) / maxdenom,0.25 * maxdenom,( DCM[ii,2,3]+ DCM[ii,3,2]) / maxdenom) )
if (index==3) Q <- (cbind( (DCM[ii,1,2]- DCM[ii,2,1]) / maxdenom,( DCM[ii,1,3]+ DCM[ii,3,1]) / maxdenom,( DCM[ii,2,3]+ DCM[ii,3,2]) / maxdenom,0.25 * maxdenom) )
if (index==4) Q <- (cbind(0.25 * maxdenom,( DCM[ii,2,3]- DCM[ii,3,2]) / maxdenom,( DCM[ii,3,1]- DCM[ii,1,3]) / maxdenom,( DCM[ii,1,2]- DCM[ii,2,1]) / maxdenom) )
} else {
ii = which(index==1) 
if (length(ii) !=0) Q[ii,] = (cbind( (DCM[ii,2,3]- DCM[ii,3,2]) / maxdenom[ii], 0.25 * maxdenom[ii], ( DCM[ii,1,2]+ DCM[ii,2,1]) / maxdenom[ii],( DCM[ii,1,3]+ DCM[ii,3,1]) / maxdenom[ii]) )
ii = which(index==2) 
if (length(ii) !=0) Q[ii,] = (cbind( (DCM[ii,3,1]- DCM[ii,1,3]) / maxdenom[ii],( DCM[ii,1,2]+ DCM[ii,2,1]) / maxdenom[ii],0.25 * maxdenom[ii],( DCM[ii,2,3]+ DCM[ii,3,2]) / maxdenom[ii]) )
ii = which(index==3) 
if (length(ii) !=0) Q[ii,] = (cbind( (DCM[ii,1,2]- DCM[ii,2,1]) / maxdenom[ii],( DCM[ii,1,3]+ DCM[ii,3,1]) / maxdenom[ii],( DCM[ii,2,3]+ DCM[ii,3,2]) / maxdenom[ii],0.25 * maxdenom[ii]) )
ii = which(index==4) 
if (length(ii) !=0) Q[ii,] = (cbind(0.25 * maxdenom[ii],( DCM[ii,2,3]- DCM[ii,3,2]) / maxdenom[ii],( DCM[ii,3,1]- DCM[ii,1,3]) / maxdenom[ii],( DCM[ii,1,2]- DCM[ii,2,1]) / maxdenom[ii]) )
}
Q
}

Q2EA <- function(Q, EulerOrder='zyx',tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{
# Implementation of quaternion to Euler angles based on D. M. Henderson 1977
# Shuttle Program. Euler Angles, Quaternions, and Transformation Matrices Working Relationships.
# National Aeronautics and Space Administration (NASA), N77-31234/6
# Q - [q1,q2,q3,q4] to EA - [psi,theta,phi]
# Madgwick (zyx) originaly used Q = [phi, theta, psi]
# Jose Gama 2014
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if(is.null(dim(Q))) N <- 1 else N <- dim(Q)[1]
Q<-matrix(Q,N,4)
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
if (EulerOrder=='zyx') { 
EA <- cbind(atan2((2*(Q[,2]*Q[,3] + Q[,1]*Q[,4])),(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)), atan2(-(2*(Q[,2]*Q[,4] - Q[,1]*Q[,3])),sqrt(1-(2*(Q[,2]*Q[,4] - Q[,1]*Q[,3]))^2)),atan2((2*(Q[,3]*Q[,4] + Q[,1]*Q[,2])),(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)))
}
if (EulerOrder=='yxz') { 
EA <- cbind(atan2(2*(Q[,2]*Q[,4] + Q[,1]*Q[,3]), 1-2*(Q[,2]^2 + Q[,3]^2)), atan2(-(2*(Q[,3]*Q[,4] - Q[,1]*Q[,2])),sqrt(1-(2*(Q[,3]*Q[,4] - Q[,1]*Q[,2]))^2)),atan2((2*(Q[,2]*Q[,3] + Q[,1]*Q[,4])),(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)))
}
if (EulerOrder=='xzy') { 
EA <- cbind( - atan2(-(2*(Q[,3]*Q[,4] + Q[,1]*Q[,2])),(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)), atan2(-(2*(Q[,2]*Q[,3] - Q[,1]*Q[,4])),sqrt(1-(2*(Q[,2]*Q[,3] - Q[,1]*Q[,4]))^2)),atan2((2*(Q[,2]*Q[,4] + Q[,1]*Q[,3])),(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)))
}
if (EulerOrder=='zxy') { 
EA <- cbind(atan2(-(2*(Q[,2]*Q[,3] - Q[,1]*Q[,4])),(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)), atan2((2*(Q[,3]*Q[,4] + Q[,1]*Q[,2])),sqrt(1-(2*(Q[,3]*Q[,4] + Q[,1]*Q[,2]))^2)),atan2(-(2*(Q[,2]*Q[,4] - Q[,1]*Q[,3])),(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)))
}
if (EulerOrder=='yzx') { 
EA <- cbind(atan2(-(2*(Q[,2]*Q[,4] - Q[,1]*Q[,3])),(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)), atan2((2*(Q[,2]*Q[,3] + Q[,1]*Q[,4])),sqrt(1-(2*(Q[,2]*Q[,3] + Q[,1]*Q[,4]))^2)),atan2(-(2*(Q[,3]*Q[,4] - Q[,1]*Q[,2])),(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2))) 
}
if (EulerOrder=='xyz') { 
EA <- cbind(atan2(-(2*(Q[,3]*Q[,4] - Q[,1]*Q[,2])),(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)), atan2((2*(Q[,2]*Q[,4] + Q[,1]*Q[,3])),sqrt(1-(2*(Q[,2]*Q[,4] + Q[,1]*Q[,3]))^2)),atan2(-(2*(Q[,2]*Q[,3] - Q[,1]*Q[,4])),(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2))) 
}
if (EulerOrder=='zyz') { 
EA <- cbind(atan2((2*(Q[,3]*Q[,4] - Q[,1]*Q[,2])),(2*(Q[,2]*Q[,4] + Q[,1]*Q[,3]))), atan2(sqrt(1-(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)^2),(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)),atan2((2*(Q[,3]*Q[,4] + Q[,1]*Q[,2])),-(2*(Q[,2]*Q[,4] - Q[,1]*Q[,3]))))
}
if (EulerOrder=='zxz') { 
EA <- cbind(atan2((2*(Q[,2]*Q[,4] + Q[,1]*Q[,3])),-(2*(Q[,3]*Q[,4] - Q[,1]*Q[,2]))), atan2(sqrt(1-(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)^2),(Q[,1]^2 - Q[,2]^2 - Q[,3]^2 + Q[,4]^2)),atan2((2*(Q[,2]*Q[,4] - Q[,1]*Q[,3])),(2*(Q[,3]*Q[,4] + Q[,1]*Q[,2])))) 
}
if (EulerOrder=='yxy') { 
EA <- cbind(atan2((2*(Q[,2]*Q[,3] - Q[,1]*Q[,4])),(2*(Q[,3]*Q[,4] + Q[,1]*Q[,2]))), atan2(sqrt(1-(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)^2),(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)),atan2((2*(Q[,2]*Q[,3] + Q[,1]*Q[,4])),-(2*(Q[,3]*Q[,4] - Q[,1]*Q[,2])))) 
}
if (EulerOrder=='yzy') { 
EA <- cbind(atan2((2*(Q[,3]*Q[,4] + Q[,1]*Q[,2])),-(2*(Q[,2]*Q[,3] - Q[,1]*Q[,4]))), atan2(sqrt(1-(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)^2),(Q[,1]^2 - Q[,2]^2 + Q[,3]^2- Q[,4]^2)),atan2((2*(Q[,3]*Q[,4] - Q[,1]*Q[,2])),(2*(Q[,2]*Q[,3] + Q[,1]*Q[,4])))) 
}
if (EulerOrder=='xzx') { 
EA <- cbind(atan2((2*(Q[,2]*Q[,4] - Q[,1]*Q[,3])),(2*(Q[,2]*Q[,3] + Q[,1]*Q[,4]))), atan2(sqrt(1-(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)^2),(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)),atan2((2*(Q[,2]*Q[,4] + Q[,1]*Q[,3])),-(2*(Q[,2]*Q[,3] - Q[,1]*Q[,4]))))
}
if (EulerOrder=='xyx') { 
EA <- cbind(atan2((2*(Q[,2]*Q[,3] + Q[,1]*Q[,4])),-(2*(Q[,2]*Q[,4] - Q[,1]*Q[,3]))), atan2(sqrt(1-(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)^2),(Q[,1]^2 + Q[,2]^2 - Q[,3]^2 - Q[,4]^2)),atan2((2*(Q[,2]*Q[,3] - Q[,1]*Q[,4])),(2*(Q[,2]*Q[,4] + Q[,1]*Q[,3]))))
}
        #EA = EA * (180/pi) # (N×3) Euler angles in degrees
        theta  = EA[,2]       # (N×1) Angle THETA in degrees
        # Check EA
        if (!ignoreAllChk) if (any(is.complex( EA ))) stop('Unreal\nUnreal Euler EA. Input resides too close to singularity.\nPlease choose different EA type.')
        # Type 1 rotation (rotations about three distinct axes)
        # THETA is computed using ASIN and ranges from -90 to 90 degrees
        if (!ignoreAllChk) {
        if (substr(EulerOrder,1,1) != substr(EulerOrder,3,3)){
	        singularities = abs(theta) > 89.9*(pi/180) # (N×1) Logical index
	        singularities[is.na(singularities)]<-FALSE
	        if (length(singularities)>0) if (any(singularities)) {
                firstsing = which(singularities)[1] # (1×1)
		        stop(paste('Input rotation ', firstsing, ' resides too close to Type 1 Euler singularity.\n',
                       'Type 1 Euler singularity occurs when second angle is -90 or 90 degrees.\n',
                       'Please choose different EA type.',sep=''))
			}
	        } else {
        # Type 2 rotation (1st and 3rd rotation about same axis)
        # THETA is computed using ACOS and ranges from 0 to 180 degrees
	        singularities = (theta<0.1*(pi/180)) | (theta>179.9*(pi/180)) # (N×1) Logical index
	        singularities[is.na(singularities)]<-FALSE
	        if (length(singularities)>0) if (any(singularities)){
                firstsing = which(singularities)[1] # (1×1)
		        stop(paste('Input rotation ', firstsing, ' resides too close to Type 2 Euler singularity.\n',
                       'Type 2 Euler singularity occurs when second angle is 0 or 180 degrees.\n',
                       'Please choose different EA type.',sep=''))
	        }
        }
        }
EA
}

Q2EA.Xiao <- function(Q, EulerOrder='zyx',tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{
# Implementation of quaternion to Euler angles based on:
# J. Xiao, 2013. Princeton Vision Toolkit. Available from: <http://vision.princeton.edu/code.html>
# http://vision.princeton.edu/pvt/GCBreader/quaternion.m
# Q - [q1,q2,q3,q4] to EA - [psi,theta,phi]
# EA in radians
# Madgwick (zyx) originaly used Q = [phi, theta, psi]
# Converted by Jose Gama 2014
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if(is.null(dim(Q))) Q<-matrix(Q,1,4)
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
Q <- Qnormalize(Q)
if (EulerOrder=='zyx') { 
EA <- cbind(atan2(2*(Q[4]*Q[1]- Q[2]*Q[3]),(Q[1]^2+ Q[2]^2-Q[3]^2-Q[4]^2)),
 asin(2*(Q[2]*Q[4]+ Q[3]*Q[1])),
 atan2(2*(Q[2]*Q[1]- Q[3]*Q[4]),(Q[1]^2- Q[2]^2-Q[3]^2+Q[4]^2)))
}
if (EulerOrder=='yxz') { 
EA <- cbind( atan2(2*(Q[3]*Q[1]- Q[4]*Q[2]),(Q[1]^2- Q[2]^2-Q[3]^2+Q[4]^2)),
 asin(2*(Q[2]*Q[1]+ Q[3]*Q[4])),
 atan2(2*(Q[4]*Q[1]- Q[2]*Q[3]),(Q[1]^2- Q[2]^2+Q[3]^2-Q[4]^2)))
}
if (EulerOrder=='xzy') { 
EA <- cbind( atan2(2*(Q[2]*Q[1]- Q[4]*Q[3]),(Q[1]^2- Q[2]^2+Q[3]^2-Q[4]^2)),
 asin(2*(Q[2]*Q[3]+ Q[4]*Q[1])),
 atan2(2*(Q[3]*Q[1]- Q[2]*Q[4]),(Q[1]^2+ Q[2]^2-Q[3]^2-Q[4]^2)))
}
if (EulerOrder=='zxy') { 
EA <- cbind( atan2(2*(Q[2]*Q[3]+ Q[4]*Q[1]),(Q[1]^2- Q[2]^2+Q[3]^2-Q[4]^2)),
 asin(2*(Q[2]*Q[1]- Q[3]*Q[4])),
 atan2(2*(Q[2]*Q[4]+ Q[3]*Q[1]),(Q[1]^2- Q[2]^2-Q[3]^2+Q[4]^2)))
}
if (EulerOrder=='yzx') { 
EA <- cbind( atan2(2*(Q[2]*Q[4]+ Q[3]*Q[1]),(Q[1]^2+ Q[2]^2-Q[3]^2-Q[4]^2)),
 asin(2*(Q[4]*Q[1]- Q[2]*Q[3])),
 atan2(2*(Q[2]*Q[1]+ Q[3]*Q[4]),(Q[1]^2- Q[2]^2+Q[3]^2-Q[4]^2))) 
}
if (EulerOrder=='xyz') { 
EA <- cbind( atan2(2*(Q[2]*Q[1]+ Q[4]*Q[3]),(Q[1]^2- Q[2]^2-Q[3]^2+Q[4]^2)),
 asin(2*(Q[3]*Q[1]- Q[2]*Q[4])),
 atan2(2*(Q[2]*Q[3]+ Q[4]*Q[1]),(Q[1]^2+ Q[2]^2-Q[3]^2-Q[4]^2))) 
}
if (EulerOrder=='zyz') { 
EA <- cbind( atan2((Q[2]*Q[1]+ Q[3]*Q[4]),(Q[3]*Q[1]- Q[2]*Q[4])),
 acos(Q[1]^2-Q[2]^2- Q[3]^2+Q[4]^2),
 atan2((Q[3]*Q[4]- Q[2]*Q[1]),(Q[2]*Q[4]+ Q[3]*Q[1])))
}
if (EulerOrder=='zxz') { 
EA <- cbind( atan2((Q[2]*Q[4]- Q[3]*Q[1]),(Q[2]*Q[1]+ Q[3]*Q[4])),
 acos(Q[1]^2-Q[2]^2- Q[3]^2+Q[4]^2),
 atan2((Q[2]*Q[4]+ Q[3]*Q[1]),(Q[2]*Q[1]- Q[3]*Q[4]))) 
}
if (EulerOrder=='yxy') { 
EA <- cbind( atan2((Q[2]*Q[3]+ Q[4]*Q[1]),(Q[2]*Q[1]- Q[3]*Q[4])),
 acos(Q[1]^2-Q[2]^2+ Q[3]^2-Q[4]^2),
 atan2((Q[2]*Q[3]- Q[4]*Q[1]),(Q[2]*Q[1]+ Q[3]*Q[4]))) 
}
if (EulerOrder=='yzy') { 
EA <- cbind( atan2((Q[3]*Q[4]- Q[2]*Q[1]),(Q[2]*Q[3]+ Q[4]*Q[1])),
 acos(Q[1]^2-Q[2]^2+ Q[3]^2-Q[4]^2),
 atan2((Q[2]*Q[1]+ Q[3]*Q[4]),(Q[4]*Q[1]- Q[2]*Q[3]))) 
}
if (EulerOrder=='xzx') { 
EA <- cbind( atan2((Q[2]*Q[4]+ Q[3]*Q[1]),(Q[4]*Q[1]- Q[2]*Q[3])),
 acos(Q[1]^2+Q[2]^2- Q[3]^2-Q[4]^2),
 atan2((Q[2]*Q[4]- Q[3]*Q[1]),(Q[2]*Q[3]+ Q[4]*Q[1])))
}
if (EulerOrder=='xyx') { 
EA <- cbind( atan2((Q[2]*Q[3]- Q[4]*Q[1]),(Q[2]*Q[4]+ Q[3]*Q[1])),
 acos(Q[1]^2+Q[2]^2- Q[3]^2-Q[4]^2),
 atan2((Q[2]*Q[3]+ Q[4]*Q[1]),(Q[3]*Q[1]- Q[2]*Q[4])))
}
        #EA = - EA  # *(180/pi) (N×3) Euler angles in radians
        theta  = EA[,2]       # (N×1) Angle THETA in radians
        # Check EA
        if (!ignoreAllChk) if (any(is.complex( EA ))) stop('Unreal\nUnreal Euler EA. Input resides too close to singularity.\nPlease choose different EA type.')
        # Type 1 rotation (rotations about three distinct axes)
        # THETA is computed using ASIN and ranges from -90 to 90 degrees
        if (!ignoreAllChk) {
        if (substr(EulerOrder,1,1) != substr(EulerOrder,3,3)){
	        singularities = abs(theta) > 89.9 *(pi/180) # (N×1) Logical index
	        singularities[is.na(singularities)]<-FALSE
	        if (length(singularities)>0) if (any(singularities)) {
                firstsing = which(singularities)[1] # (1×1)
		        stop(paste('Input rotation ', firstsing, ' resides too close to Type 1 Euler singularity.\n',
                       'Type 1 Euler singularity occurs when second angle is -90 or 90 degrees.\n',
                       'Please choose different EA type.',sep=''))
			}
	        } else {
        # Type 2 rotation (1st and 3rd rotation about same axis)
        # THETA is computed using ACOS and ranges from 0 to 180 degrees
	        singularities = (theta < 0.1*(pi/180)) | (theta > 179.9*(pi/180)) # (N×1) Logical index
	        singularities[is.na(singularities)]<-FALSE
	        if (length(singularities)>0) if (any(singularities)){
                firstsing = which(singularities)[1] # (1×1)
		        stop(paste('Input rotation ', firstsing, ' resides too close to Type 2 Euler singularity.\n',
                       'Type 2 Euler singularity occurs when second angle is 0 or 180 degrees.\n',
                       'Please choose different EA type.',sep=''))
	        }
        }
        }
EA
}

Q2EV <- function(Q,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# Q - [q1,q2,q3,q4] to EV - [m1,m2,m3,MU]
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
# Normalize quaternion(s) in case of deviation from unity. 
# User has already been warned of deviation.
v_length=4; isnot_DCM=TRUE;N=dim(Q)[1]
#Qnorms = sqrt(apply(Q^2, 2,sum))

        # Angle MU in radians and sine of MU/2
        Q2 <-matrix(Q[,1:3],ncol=3,byrow=FALSE)
        halfMUrad = matrix(atan2( sqrt(apply(Q2^2,1,sum)), Q[,4] ),ncol=1) # N×1
        #print(halfMUrad)
        SIN = sin(halfMUrad) # N×1
        index = which(SIN==0) # [N×1] Logical index
        
        if (length(index)>0){
            # Initializing
            EV = matrix(0,N,4)
            # Singular cases [MU is zero degrees]
            EV[index, 1] = 1
            # Non-singular cases
            SIN = SIN[-index, 1]
            EV[-index, ] = cbind(Q[-index,1] / SIN, 
                                 Q[-index,2] / SIN, 
                                 Q[-index,3] / SIN, 
                                 halfMUrad * 2 )#* (360/pi)
        } else {
            # Non-singular cases            
            EV = cbind(Q[,1]/SIN, Q[,2]/SIN, Q[,3]/SIN, halfMUrad * 2)#*(360/pi)
        }
        # MU greater than 180 degrees
        index = which(EV[,4] > pi/2) # [N×1] Logical index
        EV[index, ] = cbind(-matrix(EV[index,1:3],ncol=3,byrow=FALSE), 2*pi-EV[index,4])
EV
}

Q2DCM <- function(Q,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# Q - [q1,q2,q3,q4] to DCM - 3x3xN
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
# Normalize quaternion(s) in case of deviation from unity. 
Qn <- Qnormalize(Q)
# User has already been warned of deviation.
N<-dim(Qn)[1]
        Qn  = array(t(Qn), c(1, 4, N))
        DCM= array(0, c(3, 3, N))
       DCM[1,1,] = 1-2*(Qn[1,3,]*Qn[1,3,]+Qn[1,4,]*Qn[1,4,])
       DCM[2,1,] = 2*(Qn[1,2,]*Qn[1,3,]-Qn[1,1,]*Qn[1,4,])
       DCM[3,1,] = 2*(Qn[1,2,]*Qn[1,4,]+Qn[1,1,]*Qn[1,3,])
       DCM[1,2,] = 2*(Qn[1,2,]*Qn[1,3,]+Qn[1,1,]*Qn[1,4,])
       DCM[2,2,] = 1-2*(Qn[1,2,]*Qn[1,2,]+Qn[1,4,]*Qn[1,4,])
       DCM[3,2,] = 2*(Qn[1,3,]*Qn[1,4,]-Qn[1,1,]*Qn[1,2,])
       DCM[1,3,] = 2*(Qn[1,2,]*Qn[1,4,]-Qn[1,1,]*Qn[1,3,])
       DCM[2,3,] = 2*(Qn[1,3,]*Qn[1,4,]+Qn[1,1,]*Qn[1,2,]) 
       DCM[3,3,] = 1-2*(Qn[1,2,]*Qn[1,2,]+Qn[1,3,]*Qn[1,3,])
       if (length(dim(DCM))==3) if(all(dim(DCM)==c(3,3,1))) DCM <- matrix(DCM,3)
DCM
}

# Q2DCM <- function(Q,tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
# {# Q - [q1,q2,q3,q4] to DCM - 3x3xN
# if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
# if (!ignoreAllChk) if (ichk && any(abs(sqrt(apply(Q^2, 2,sum)) - 1) > tol)) warning('Warning: (At least one of the) Input quaternion(s) is not a unit vector')
# # Normalize quaternion(s) in case of deviation from unity. 
# # User has already been warned of deviation.
# N<-dim(Q)[1]
#         Q  = array(t(Q), c(1, 4, N))
#         SQ = Q^2
#         DCM= array(0, c(3, 3, N))
#        DCM[1,1,] = SQ[1,1,]-SQ[1,2,]-SQ[1,3,]+SQ[1,4,]
#        DCM[1,2,] = 2*(Q[1,1,]*Q[1,2,] +Q[1,3,]*Q[1,4,])
#        DCM[1,3,] = 2*(Q[1,1,]*Q[1,3,] -Q[1,2,]*Q[1,4,])
#        DCM[2,1,] = 2*(Q[1,1,]*Q[1,2,] -Q[1,3,]*Q[1,4,])
#        DCM[2,2,] = -SQ[1,1,]+SQ[1,2,]-SQ[1,3,]+SQ[1,4,]
#        DCM[2,3,] = 2*(Q[1,2,]*Q[1,3,] +Q[1,1,]*Q[1,4,])
#        DCM[3,1,] = 2*(Q[1,1,]*Q[1,3,] +Q[1,2,]*Q[1,4,])
#        DCM[3,2,] =  2*(Q[1,2,]*Q[1,3,] -Q[1,1,]*Q[1,4,])
#        DCM[3,3,] = -SQ[1,1,]-SQ[1,2,]+SQ[1,3,]+SQ[1,4,]
#        if (all(dim(DCM)==c(3,3,1))) DCM <- matrix(DCM,3)
# DCM
# }

Q2GL<-function(Q)
{# Q - [q1,q2,q3,q4] to OpenGL translation matrix 4x4xn
#based on
#http://www.tinkerforge.com/doc/Software/Bricks/IMU_Brick_Python.html
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
GL <- lapply(1:N, function(n) 
{
x<-Q[n,1]
y<-Q[n,2]
w<-Q[n,3]
z<-Q[n,4]
tmp <- c(1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y), 0,
 2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x), 0,
 2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y), 0,
 0, 0, 0, 1)
matrix(tmp,nrow=4, ncol=4)
})
if (N==1) GL=array(unlist(GL), dim = c(dim(GL[[1]]) )) else GL=array(unlist(GL), dim = c(dim(GL[[1]]), length(GL))) 
GL
}

EV2EA<-function(EV, EulerOrder='zyx',tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EV - [m1,m2,m3,MU] to EA - [psi,theta,phi]
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if (!is.matrix(EV)) EV <- matrix(unlist(EV),ncol=4,byrow=FALSE)
Q<-EV2Q(EV, tol, ichk, ignoreAllChk)
EA<-Q2EA(Q, EulerOrder, tol, ichk, ignoreAllChk)
EA
}

EA2EV<-function(EA, EulerOrder='zyx',tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EA - [psi,theta,phi] to EV - [m1,m2,m3,MU]
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if (!is.matrix(EA)) EA <- matrix(unlist(EA),ncol=3,byrow=FALSE)
Q<-EA2Q(EA, EulerOrder, ichk, ignoreAllChk)
EV<-Q2EV(Q, tol, ichk, ignoreAllChk)
EV
}

EA2EA<-function(EA, EulerOrder1='zyx',EulerOrder2='zyx',tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{#EA - [psi,theta,phi] to EA - [psi,theta,phi]
if (!is.character(EulerOrder1)) stop('<<EulerOrder1>> must be a string.')
EulerOrder1 <- tolower(EulerOrder1)
if (!(EulerOrder1 %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if (!is.character(EulerOrder2)) stop('<<EulerOrder2>> must be a string.')
EulerOrder2 <- tolower(EulerOrder2)
if (!(EulerOrder2 %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')

if (all(EulerOrder1==EulerOrder2)) return (EA)
Q<- EA2Q(EA, EulerOrder1, ichk, ignoreAllChk)
EA<-Q2EA(Q, EulerOrder2, tol, ichk, ignoreAllChk)
EA
}

EA2DCM<-function(EA, EulerOrder='zyx',tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EA - [psi,theta,phi] to DCM - 3x3xN
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if (!is.matrix(EA)) EA <- matrix(unlist(EA),ncol=3,byrow=FALSE)
Q<-EA2Q(EA, EulerOrder, ichk, ignoreAllChk)
DCM<-Q2DCM(Q, tol, ichk, ignoreAllChk)
if (length(DCM)==3) if (all(dim(DCM)==c(3,3,1))) DCM <- matrix(DCM,3)
DCM
}

EV2DCM<-function(EV, tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# EV - [m1,m2,m3,MU] to DCM - 3x3xN
if (!is.matrix(EV)) EV <- matrix(unlist(EV),ncol=4,byrow=FALSE)
Q<-EV2Q(EV, tol, ichk, ignoreAllChk)
DCM<-Q2DCM(Q, tol,  ichk, ignoreAllChk)
if (length(DCM)==3) if (all(dim(DCM)==c(3,3,1))) DCM <- matrix(DCM,3)
DCM
}

DCM2EV<-function(DCM, tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# DCM - 3x3xN to EV - [m1,m2,m3,MU]
Q<-DCM2Q(DCM, tol, ichk, ignoreAllChk)
EV<-Q2EV(Q, tol, ichk, ignoreAllChk)
EV
}

DCM2EA<-function(DCM, EulerOrder='zyx', tol = 10 * .Machine$double.eps, ichk=FALSE, ignoreAllChk=FALSE)
{# DCM - 3x3xN to EA - [psi,theta,phi]
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
Q<-DCM2Q(DCM, tol, ichk, ignoreAllChk)
EA<-Q2EA(Q, EulerOrder, tol, ichk, ignoreAllChk)
EA
}

vectQrot<-function( Q, rr )
# Rotate a vector by a quaternion
{
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
if (!is.matrix(rr)) rr <-matrix(rr,ncol=3,byrow=FALSE)
N <- dim(Q)[1]
if (!is.matrix(rr)) rr <-matrix( rr,ncol=3,nrow=N,byrow=FALSE)
Qr <- (Qconj(Q) %Q*% cbind(0, rr)) %Q*% Q
Qr[,2:4]
}

Qrot <- function(Q,w,dT)
# Updates current attitude quaternion q
# output - current attitude quaternion
# input - wx, wy, wz - angular rate values
# input - dT - inverse of update rate
{
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
if (!is.matrix(w)) w <-matrix(w,ncol=3,byrow=FALSE)#nrow=N,
Qr <- matrix(0,nrow=N, ncol=4)
Qr<-lapply(1:N, function(n) {
Fx <- w[n,1]*dT
Fy <- w[n,2]*dT
Fz <- w[n,3]*dT
Fm <- sqrt(Fx*Fx+Fy*Fy+Fz*Fz)
sinFm2 <- sin(Fm/2)
cosFm2 <- cos(Fm/2)
if (Fm != 0)
	Qr[n,]<-c(cosFm2, Fx/Fm*sinFm2, Fy/Fm*sinFm2, Fz/Fm*sinFm2)#Qr[n,]
else
    Qr[n,]<-c(1, 0, 0, 0)#Qr[n,]
Qr[n,]<- Q[n,] %Q*% Qr[n,]#Qr[n,]<- Q[n,] %Q*% Qr[n,]
Qr
})
if (N==1) Qr=array(unlist(Qr), dim = c(dim(Qr[[1]]) )) else Qr=array(unlist(Qr), dim = c(dim(Qr[[1]]), length(Qr))) 
if (length(dim(Qr))==3) Qr <- array(Qr,c(N,4))
#if (length(dim(Qr))==3) if (all(dim(Qr)==c(1,4,1))) Qr <- array(Qr,c(1,4))
Qr
}

"%Q+%" <- function(Q1, Q2)
{# quaternion sum
if (!is.matrix(Q1)) Q1 <-matrix(Q1,ncol=4,byrow=FALSE)
if (!is.matrix(Q2)) Q2 <-matrix(Q2,ncol=4,byrow=FALSE)
Q1 + Q2
}

"%Q-%" <- function(Q1, Q2)
{# quaternion difference
if (!is.matrix(Q1)) Q1 <-matrix(Q1,ncol=4,byrow=FALSE)
if (!is.matrix(Q2)) Q2 <-matrix(Q2,ncol=4,byrow=FALSE)
Q1 - Q2
}

Qconj<-function(Q) 
{# quaternion conjugate
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
Q[,2:4] <- -Q[,2:4]
Q
}

Qinv<-function(Q)
{# quaternion inverse
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
if (Qnorm(Q)==1) return (Qconj(Q))
else return (Qconj(Q)/Qnorm(Q))
}

"%Q*%" <- function(Q1, Q2)
{# quaternion product
if (!is.matrix(Q1)) Q1 <-matrix(Q1,ncol=4,byrow=FALSE)
if (!is.matrix(Q2)) Q2 <-matrix(Q2,ncol=4,byrow=FALSE)
N <- dim(Q1)[1]
ab<-matrix(0,ncol=4,nrow=N)
ab<-lapply(1:N, function(n) {
ab[n,1]<-Q1[n,1]*Q2[n,1]-Q1[n,2]*Q2[n,2]-Q1[n,3]*Q2[n,3]-Q1[n,4]*Q2[n,4]
ab[n,2]<-Q1[n,1]*Q2[n,2]+Q1[n,2]*Q2[n,1]+Q1[n,3]*Q2[n,4]-Q1[n,4]*Q2[n,3]
ab[n,3]<-Q1[n,1]*Q2[n,3]-Q1[n,2]*Q2[n,4]+Q1[n,3]*Q2[n,1]+Q1[n,4]*Q2[n,2]
ab[n,4]<-Q1[n,1]*Q2[n,4]+Q1[n,2]*Q2[n,3]-Q1[n,3]*Q2[n,2]+Q1[n,4]*Q2[n,1]
ab
} )
if (N==1) ab=array(unlist(ab), dim = c(dim(ab[[1]]) )) else ab=array(unlist(ab), dim = c(dim(ab[[1]]), length(ab))) 
ab<-matrix(ab,ncol=4,nrow=N)
return(ab)
}

"%Q/%" <- function(Q1, Q2)
{# quaternion division
Q1 %Q*% Qinv(Q2)
}

Qnorm<-function(Q)
{ # norm of a quaternion
if (length(Q)==3) Q<-c(0,Q)
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
if (N==1) sqrt(sum(Q^2))
apply(Q,1,function(n) sqrt(sum(n^2)))
}

Qnormalize<-function(Q)
{ # normalize a quaternion
if (length(Q)==3) Q<-c(0,Q)
if (!is.matrix(Q)) Q <-matrix(Q,ncol=4,byrow=FALSE)
N <- dim(Q)[1]
if (N==1) Qn <- Q/sqrt(sum(Q^2)) else Qn <- Q/apply(Q,1,function(n) sqrt(sum(n^2)))
Qn[which(is.na(Qn))] <- 0
Qn
}

#############################################################################################

isUnitQuaternion<-function(Q)
{# returns TRUE if q is a unit quaternion
# |q|==1 => q is a unit quaternion
# source: Quaternion C++ class by Will Perone, 2003
# http://willperone.net/Code/quaternion.php
if (!is.matrix(Q)) Q <- matrix(Q,ncol=4,byrow=FALSE)
qN <- Qnorm(Q)
qN == 1
}

isPureQuaternion<-function(Q)
{# returns TRUE if q is a pure quaternion
# q==(0,v) => q is a pure quaternion
# source: Quaternion C++ class by Will Perone, 2003
# http://willperone.net/Code/quaternion.php
if (!is.matrix(Q)) Q <- matrix(Q,ncol=4,byrow=FALSE)
(Q[,2] == 0) & (Q[,3] == 0) & (Q[,4] == 0)
}

isRealQuaternion<-function(Q)
{# returns TRUE if q is a real quaternion
# q==(0,v) => q is a pure quaternion
# source: Quaternion C++ class by Will Perone, 2003
# http://willperone.net/Code/quaternion.php
if (!is.matrix(Q)) Q <- matrix(Q,ncol=4,byrow=FALSE)
Q[,1] == 0
}

Qlog<-function(Q)
{# returns the logarithm of a quaternion
# source: Quaternion C++ class by Will Perone, 2003
# http://willperone.net/Code/quaternion.php
if (!is.matrix(Q)) Q <- matrix(Q,ncol=4,byrow=FALSE)
a <- acos(Q[,1])
sina <- sin(a)
Qret <- matrix(0,nrow=dim(Q)[1],4)
sinaZero <- which(sina>0)
Qret[sinaZero,2:4] <- a[sinaZero] * Q[sinaZero,2:4] / sina[sinaZero]
Qret
}

Qexp<-function(Q)
{# returns the exponential of a quaternion
# source: Quaternion C++ class by Will Perone, 2003
# http://willperone.net/Code/quaternion.php
if (!is.matrix(Q)) Q <- matrix(Q,ncol=4,byrow=FALSE)
a <- exp(Q[,1])
sina <- sin(a)
cosa <- cos(a)
Qret <- cbind(cosa,0,0,0)
sinaZero <- which(sina>0)
Qret[sinaZero,2:4] <- sina[sinaZero] * Q[sinaZero,2:4] / a[sinaZero]
Qret
}

"%Q.%" <- function(Q1, Q2)
{# quaternion dot product
apply(Q1 * Q2,1,sum)
}

Qzero <- function(n=NA){
# generate zero-value quaternions
if (is.na(n)) n <- 1
matrix(0,n,4)
}

Qone <- function(n=NA){
# generate one-value quaternions
if (is.na(n)) n <- 1
matrix(c(1,0,0,0),n,4)
}

Qlerp<-function(Q1, Q2, fracT)
{#  linear quaternion interpolation
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
N <- dim(Q1)[1]
if (fracT == 0.0) return (Qzero(N))
if (fracT == 1.0) return (Qone(N))
Qnormalize((Q1*(1-fracT) %Q+% Q2*fracT))
}

Qslerp<-function(Q1, Q2, fracT)
{# spherical linear interpolation
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
N <- dim(Q1)[1]
if (fracT == 0.0) return (Qzero(N))
if (fracT == 1.0) return (Qone(N))
Q3 <- Q2
Qd <- Q1 %Q.% Q2
QdLtz <- which(Qd < 0)
Qd[QdLtz] <- -Qd[QdLtz]
Q3[QdLtz] <- -Q2[QdLtz]
QdLt.95 <- which(Qd < 0.95)
angleQ = acos(Qd[QdLt.95])
Q3[QdLt.95] <- (Q1[QdLt.95] * sin(angleQ * (1-fracT)) + Q3[QdLt.95] * sin(angleQ * fracT))/sin(angleQ)
Q3[-QdLt.95] <- Qlerp(Q1[QdLt.95], Q3[QdLt.95],fracT)
Q3
}

QslerpNoInvert<-function(Q1, Q2, fracT)
{# This version of slerp, used by squad, does not check for theta > 90.
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
N <- dim(Q1)[1]
if (fracT == 0.0) return (Qzero(N))
if (fracT == 1.0) return (Qone(N))
Q3 <- Q2
Qd <- Q1 %Q.% Q2
QdLt.bt95 <- which((Qd < 0.95) & (Qd > -0.95))
angleQ = acos(Qd[QdLt.bt95])
Q3[QdLt.bt95] <- (Q1[QdLt.bt95] * sin(angleQ * (1-fracT)) + Q2[QdLt.bt95] * sin(angleQ * fracT))/sin(angleQ)
Q3[-QdLt.bt95] <- Qlerp(Q1[QdLt.bt95], Q2[QdLt.bt95],fracT)
Q3
}

Qsquad<-function(Q1, Q2, Qa, Qb, fracT)
{# Spherical and Quadrangle linear interpolation
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
Qc <- QslerpNoInvert(Q1,Q2,fracT)
Qd <- QslerpNoInvert(Qa,Qb,fracT)
QslerpNoInvert(Qc,Qd,2*fracT*(1-fracT))
}

Qbezier<-function(Q1, Q2, Qa, Qb, fracT)
{# Shoemake-Bezier interpolation using De Castlejau algorithm
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
Q11 <- QslerpNoInvert(Q1,Qa,fracT)
Q12 <- QslerpNoInvert(Qa,Qb,fracT)
Q13 <- QslerpNoInvert(Qb,Q2,fracT)
QslerpNoInvert(QslerpNoInvert(Q11,Q12,fracT), QslerpNoInvert(Q12,Q13,fracT), fracT)
}

Qspline<-function(Qnm1, Qn, Qnp1)
{# Given 3 quaternions, qn-1,qn and qn+1, calculate a control point to be used in spline interpolation
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
Qni <- Qn
Qni[,2:4] <- -Qni[,2:4]
Qn %Q*% Qexp(( Qlog(Qni %Q*% Qnm1) %Q+% Qlog(Qni %Q*% Qnp1) ) / -4)
}

QangularDifference<-function(Q1, Q2)
{#angular difference between 2 quaternions
acos((Q1 %Q.% Q2)/(Qnorm(Q1)*Qnorm(Q2)))
}

Qrandom <- function(n=NA)
{# generate - uniform random unit quaternions
# 
# http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
# LaValle, Steven M. "Planning algorithms." University of Illinois 2004 (1999).
# pp. 164 - "uniform random unit quaternions leads to uniform random samples over SO(3)"
if (is.na(n)) n <- 1
R <- matrix(stats::runif(n*3),n,3)
r1 <- sqrt(1.0 - R[,1])
r2 <- sqrt(R[,1])
pi2 <- pi * 2.0
t1 <- pi2 * R[,2]
t2 <- pi2 * R[,3]
Q <- cbind(cos(t2)*r2, sin(t1)*r1,cos(t1)*r1, sin(t2)*r2)
Qnormalize(Q)
}

isPureRotationMatrix<-function(DCM, tol=0.01)
{# pure rotation matrix = proper orthogonal matrix = det(m)==1
if (is.matrix(DCM)) return(((det(DCM)-1)<tol) & (det(DCM)>0)) else {
apply(DCM,3,function(x) (((det(x)-1)<tol) & (det(x)>0)) )
}
}

DCMrandom <- function(n=NA, tol = 10 * .Machine$double.eps, ignoreAllChk=FALSE)
{# generate - uniform random DCM
if (is.na(n)) n <- 1
Q <- Qrandom(n)
DCM <- Q2DCM(Q, tol = tol, ignoreAllChk=TRUE)
if (ignoreAllChk) return(DCM)
badQ <- which(isPureRotationMatrix(DCM)==FALSE)
while (length(badQ)>0) {# remove wrong matrices
n <- length(badQ)
Q <- Qrandom(n)
DCM[badQ] <-  Q2DCM(Q, tol = tol, ignoreAllChk=TRUE)
badQ <- which(isPureRotationMatrix(DCM)==FALSE)
}
DCM
}

EArandom <- function(n=NA, EulerOrder='zyx',tol = 10 * .Machine$double.eps, ignoreAllChk=FALSE)
{# generate - uniform random Euler Angles
if (!is.character(EulerOrder)) stop('<<EulerOrder>> must be a string.')
EulerOrder <- tolower(EulerOrder)
if (!(EulerOrder %in% c('zyx','zxy','yxz','xzy','xyz','yzx','zyz','zxz','yxy','yzy','xyx','xzx'))) stop('Invalid input Euler angle order')
if (is.na(n)) n <- 1
Q <- Qrandom(n)
EA <- Q2EA(Q, EulerOrder, tol = tol, ignoreAllChk=TRUE)
if (ignoreAllChk) return(EA)
badQ <- which((abs(EA[,2]) > 89.9 *(pi/180)) & ((EA[,2] < 0.1*(pi/180)) | (EA[,2] > 179.9*(pi/180))))
while (length(badQ)>0) {# fix singularities
n <- length(badQ)
Q <- Qrandom(n)
EA[badQ] <- Q2EA(Q, EulerOrder, tol = tol, ignoreAllChk=TRUE)
badQ <- which((abs(EA[,2]) > 89.9 *(pi/180)) & ((EA[,2] < 0.1*(pi/180)) | (EA[,2] > 179.9*(pi/180))))
}
EA
}

EVrandom <- function(n=NA, tol = 10 * .Machine$double.eps, ignoreAllChk=FALSE)
{# generate - uniform random Euler Vectors
if (is.na(n)) n <- 1
Q <- Qrandom(n)
EV <- Q2EV(Q, tol = tol, ignoreAllChk=TRUE)
if (ignoreAllChk) return(EV)
EV
}

