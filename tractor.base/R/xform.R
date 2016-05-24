extractRotationMatrixFromXform <- function (xformMatrix)
{
    if (!is.matrix(xformMatrix) || !equivalent(dim(xformMatrix),c(4,4)))
        report(OL$Error, "The xform must be a 4x4 matrix")
    
    rotationMatrix <- xformMatrix[1:3,1:3]
    
    columnLengths <- apply(rotationMatrix, 2, vectorLength)
    if (any(columnLengths == 0))
        report(OL$Error, "The specified xform matrix contains one or more zero-length columns")
    
    # Normalise the rotation matrix
    rotationMatrix <- sapply(1:3, function(i) rotationMatrix[,i] / columnLengths[i])
    
    # If the matrix is not orthogonal, find a near orthogonal matrix
    if (!equivalent(rotationMatrix %*% t(rotationMatrix), diag(3), tolerance=1e-6))
    {
        report(OL$Debug, "Rotation matrix is not orthogonal")
        svd <- La.svd(rotationMatrix)
        rotationMatrix <- svd$u %*% svd$vt
    }
    
    return (rotationMatrix)
}

xformToOrientation <- function (xformMatrix, string = TRUE)
{
    rotationMatrix <- extractRotationMatrixFromXform(xformMatrix)
    
    directions <- apply(abs(rotationMatrix) > 0.5, 2, which)
    directions <- directions * sign(rotationMatrix[cbind(directions,1:3)])
    
    if (string)
        return (implode(c("I","P","L","","R","A","S")[directions+4], sep=""))
    else
        return (directions * c(-1,1,1))
}

xformToQuaternion <- function (xformMatrix)
{
    if (!is.matrix(xformMatrix) || !equivalent(dim(xformMatrix),c(4,4)))
        report(OL$Error, "The xform must be a 4x4 matrix")
    
    offset <- xformMatrix[1:3,4]
    
    r <- extractRotationMatrixFromXform(xformMatrix)
    handedness <- sign(det(r))
    
    if (handedness < 0)
        r[,3] <- (-r[,3])
    
    # Compute quaternion parameters (code translated from nifti1_io.c)
    a <- sum(diag(r)) + 1
    if (a > 0.5)
    {
        a <- 0.5 * sqrt(a)
        b <- 0.25 * (r[3,2]-r[2,3]) / a
        c <- 0.25 * (r[1,3]-r[3,1]) / a
        d <- 0.25 * (r[2,1]-r[1,2]) / a
        
        q <- c(a,b,c,d)
    }
    else
    {
        xd <- 1 + r[1,1] - (r[2,2]+r[3,3])
        yd <- 1 + r[2,2] - (r[1,1]+r[3,3])
        zd <- 1 + r[3,3] - (r[1,1]+r[2,2])
        if (xd > 1)
        {
            b <- 0.5 * sqrt(xd)
            c <- 0.25 * (r[1,2]+r[2,1]) / b
            d <- 0.25 * (r[1,3]+r[3,1]) / b
            a <- 0.25 * (r[3,2]-r[2,3]) / b
        }
        else if (yd > 1)
        {
            c <- 0.5 * sqrt(yd)
            b <- 0.25 * (r[1,2]+r[2,1]) / c
            d <- 0.25 * (r[2,3]+r[3,2]) / c
            a <- 0.25 * (r[1,3]-r[3,1]) / c
        }
        else
        {
            d <- 0.5 * sqrt(zd)
            b <- 0.25 * (r[1,3]+r[3,1]) / d
            c <- 0.25 * (r[2,3]+r[3,2]) / d
            a <- 0.25 * (r[2,1]-r[1,2]) / d
        }
        
        q <- c(a,b,c,d)
        if (a < 0)
            q <- (-q)
     }
     
     return (list(q=q, offset=offset, handedness=handedness))
}

quaternionToXform <- function (quaternion)
{
    if (length(quaternion) < 3 || length(quaternion) > 4)
        report(OL$Error, "The quaternion should have length 3 or 4")
    if (length(quaternion) == 4)
        quaternion <- quaternion[2:4]
    
    matrix <- diag(4)
    
    quaternionSumOfSquares <- sum(quaternion^2)
    if (equivalent(quaternionSumOfSquares, 1, tolerance=1e-6))
        q <- c(0, quaternion)
    else if (quaternionSumOfSquares > 1)
        report(OL$Error, "Quaternion parameters are invalid")
    else
        q <- c(sqrt(1 - quaternionSumOfSquares), quaternion)

    matrix[1:3,1:3] <- c(  q[1]*q[1] +   q[2]*q[2] - q[3]*q[3] - q[4]*q[4],
                         2*q[2]*q[3] + 2*q[1]*q[4],
                         2*q[2]*q[4] - 2*q[1]*q[3],
                         2*q[2]*q[3] - 2*q[1]*q[4],
                           q[1]*q[1] +   q[3]*q[3] - q[2]*q[2] - q[4]*q[4],
                         2*q[3]*q[4] + 2*q[1]*q[2],
                         2*q[2]*q[4] + 2*q[1]*q[3],
                         2*q[3]*q[4] - 2*q[1]*q[2],
                           q[1]*q[1] +   q[4]*q[4] - q[3]*q[3] - q[2]*q[2])
    
    return (matrix)
}
