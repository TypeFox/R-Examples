
# compute the combined rotation vector for two input rotation vectors
`rot.comb` <-
function (rotvec1 = c(1, 0, 0, pi/2), rotvec2 = c(0, 0, 1, pi)) 
{
    # compute the quaternion representations
    q1 <- .rot2quaternion(rotvec1)
    q2 <- .rot2quaternion(rotvec2)
    
    # multiply the quaternions
    qprod <- .multquaternion(q1, q2)
    
    return(round(.quaternion2rot(qprod), 4))
}

