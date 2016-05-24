transfMatrix = function(angle=0, offset=c(0,0) ) 
{
#---------------------------------------------------------------------------
#
#   Builds a simple transformation matrix as the concatenation of rotation
#   and translation for 2D at the present. Note that this is an affine
#   transformation so we use homogeneous coordinates.
#
#   Note that translation follows rotation the way this is set up currently.
#
#   Also, we set this up so that we can PRE-multiply by a matrix of
#   coordinates to be transformed with columns in the form; [x, y, 1]
#   One can always transpose to post-multiply coordinates.
#
#   Arguments...
#     angle = the angle for rotation: positive = counter clockwise
#     offset = the amount to translate for (x,y)
#
#   Returns...
#     the transformatin matrix--returns identity by default
#
#Author...									Date: 27-Apr-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   translation matrix...
#
    tm = diag(3)
    tm[3,1:2] = offset

#
#   rotation matrix...
#
    rot = diag(3)
    rot[1,1] = cos(angle)
    rot[2,1] = -sin(angle)
    rot[1,2] = -rot[2,1]
    rot[2,2] = rot[1,1]

#
#   concatenate, translation follows rotation...
#
    transf = rot%*%tm
    
    return(transf)
} #transfMatrix
