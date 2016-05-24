
    
    
    
#   MunsellSpecToHVC()
#
#   convert Munsell notation to numeric Hue,Value,Chroma
#       Hue     is ASTM Hue in the interval (0,100]  (except Hue=0 for neutrals)
#       Value   is in the interval [0,10]
#       Chroma  is positive  (except Chroma=0 for neutrals)
#
#   arguments:
#       MunsellSpecString   a character vector of length n
#
#   return value: an nx3 matrix with the computed HVCs placed in the rows
#
#   author:  Glenn Davis
MunsellSpecToHVC <- function( MunsellSpecString )
    {
    if( ! is.character(MunsellSpecString) )
        {
        #   mess = "MunsellSpecToHVC().  MunsellSpecString is not a string.\n"
        #   cat( mess )
        return(NULL)
        }

    n = length(MunsellSpecString)
    
    out = matrix( as.numeric(NA), n, 3 )
    colnames(out)   = c( "Hue", "Value", "Chroma" )
    rownames(out)   = MunsellSpecString
        
    # Remove whitespace from input Munsell string
    MunsellSpecString <- gsub( ' |\t', '', MunsellSpecString )
    
    # Make all letters in Munsell string upper case
    MunsellSpecString <- toupper(MunsellSpecString)  

    #-----   chromatic colors   ----------#
    pattern = "^([0-9.]+)(R|YR|Y|GY|G|BG|B|PB|P|RP)([0-9.]+)/([0-9.]+)$" 

    #   time_start  = as.double( Sys.time() )
        
    sub1234 = sub( pattern, "\\1 \\2 \\3 \\4", MunsellSpecString )
    
    mask    = (nchar(MunsellSpecString) < nchar(sub1234) )     # a matched string gets longer          #mask = grepl( pattern, MunsellSpecString )   #   ; print( mask )
    
    if( any(mask) )
        {
        #   make a little LUT
        hue = seq( 0, 90, by=10 )
        names(hue)  = c("R","YR","Y","GY","G","BG","B","PB","P","RP")
        
        sub1234    = sub1234[mask]        
        
        dat     = unlist( strsplit( sub1234, ' ', fixed=T ) ) #;  print( dat )
        dat     = matrix( dat, length(sub1234), 4, byrow=T ) #;    print( dat )
        
        out[ mask, 1]   = hue[ dat[ ,2] ]  +  as.numeric(dat[ ,1])
        out[ mask, 2]   = as.numeric( dat[ ,3] )
        out[ mask, 3]   = as.numeric( dat[ ,4] )
        }
        
    #   cat( "Elapsed: ", as.double( Sys.time() ) - time_start, '\n' )
        
    #------   achromatic colors     ----#
    pattern = "^N([0-9.]+)(/|/0)?$"

    sub1    = sub( pattern, "\\1", MunsellSpecString )
        
    mask    = (nchar(sub1) < nchar(MunsellSpecString) ) # a matched string gets shorter.    #mask = grepl( pattern, MunsellSpecString )   #   ; print( mask )
 
    if( any(mask) )
        {
        out[ mask, 1]    = 0
        out[ mask, 2]    = as.numeric( sub1[mask] )
        out[ mask, 3]    = 0
        }
        
    return( out )
    }



#   ColorBlockFromMunsell()
#
#   HVC        a numeric 3-vector, or an nx3 matrix
#
#   if HVC[] is a 3-vector
#       HVC[1]     Munsell hue, on the ASTM D1535 100 point circular scale. All values are valid.
#       HVC[2]     Munsell value, must be between 0 and 10
#       HVC[3]     Munsell chroma, must be non-negative
#
#   if HVC[,] is an nx3 matrix
#       each row is defined as above
#
#   return value:
#       for a 3-vector:     list with: HVC, ISCC-NBS Number, ISCC-NBS Name
#       for an nx3 matrix:  data.frame with columns: HVC, ISCC-NBS Number, ISCC-NBS Name
#
#   the function requires these global data.frames:
#       SystemISCCNBS
#       CentralsISCCNBS
#
#   author:  Glenn Davis
ColorBlockFromMunsell  <-  function( HVC )
    {    
    if( ! is.numeric(HVC) )
        {
        #   mess = "ColorBlockFromMunsell().  HVC is not numeric.\n"
        #   cat( mess )
        return(NULL)
        }

    if( is.matrix(HVC) )
        {
        if( ncol(HVC) != 3 )
            {
            #   mess = "ColorBlockFromMunsell().  HVC does not have 3 columns.\n"
            #   cat( mess )
            return(NULL)
            }        
        
        colnames(HVC)  = c( "H", "V", "C" )        
        class(HVC)     = "model.matrix" 
        
        out = data.frame( HVC=HVC, Number=as.integer(NA), Name=as.character(NA), stringsAsFactors=FALSE )
        
        for( i in 1:nrow(out) )
            {
            result  = ColorBlockFromMunsell( HVC[i, ] )
        
            out$Number[i]   = result$Number
            out$Name[i]     = result$Name
            }
        
        return(out)
        }
    
    if( length(HVC) != 3 )
        {
        #   mess = "ColorBlockFromMunsell().  length of HVC is not 3.\n"
        #   cat( mess )
        return(NULL)
        }        
        

    if( is.null( names(HVC) ) )
        names(HVC) = c( "Hue", "Value", "Chroma" )
    
    out = list( HVC=HVC, Number=as.integer(NA), Name=as.character(NA) )
            
    vmax    = 10            
    valid   = all( is.finite(HVC) )  &&  (0 <= HVC[2])  &&  (HVC[2] <= vmax)  &&  (0 <= HVC[3])
    
    if( ! valid )   return( out )
    
    #   translate hue to the interval [1,101)   (101 is not included)
    HVC[1] = ((HVC[1] - 1 ) %% 100) + 1
        
    if( HVC[2] == vmax )   HVC[2] = vmax - 1.e-6     # because upper comparison below is strict
    
    #   do 6-way comparison.  
    #   Note upper comparisons are strict, and lower comparisons are not strict.
    #   So a point on a boundary is in only 1 block.
    mask.H  = colorscience::SystemISCCNBS$Hmin <= HVC[1]  &  HVC[1] < colorscience::SystemISCCNBS$Hmax
    mask.V  = colorscience::SystemISCCNBS$Vmin <= HVC[2]  &  HVC[2] < colorscience::SystemISCCNBS$Vmax
    mask.C  = colorscience::SystemISCCNBS$Cmin <= HVC[3]  &  HVC[3] < colorscience::SystemISCCNBS$Cmax
        
    theRow  = colorscience::SystemISCCNBS[ mask.H & mask.V & mask.C, ]
    
    if( nrow(theRow) != 1 )
        {
        #   mess = sprintf( "Expected exactly 1 match, but found %d\n", nrow(theRow) )
        #   cat( mess )
        return( out )
        }
        
    out$Number  = theRow$Number
    
    out$Name    = colorscience::CentralsISCCNBS$Name[ out$Number ]
    
    return( out )
    }

#
#   CheckColorLookup()
#   performs lookup with colorBlockFromMunsell() and verifies that the color block number is correct
#
#   DataISCCNBS   data.frame with columns
#       MunsellSpec     a string in Munsell notation
#       Number          the correct ISCC-NBS Number
#
#   return value:  TRUE or FALSE
#
#   author:  Glenn Davis
CheckColorLookup <- function( DataISCCNBS=colorscience::CentralsISCCNBS )
    {
    hvc = MunsellSpecToHVC( DataISCCNBS$MunsellSpec )
    
    block   = ColorBlockFromMunsell( hvc )  #;      print( v )
    
    #   compute errors
    idx = which( block$Number != DataISCCNBS$Number )
    
    out = (length(idx) == 0)
    
    if( out )
        cat( "There are no errors !\n" )
    else
        {
        mess    = sprintf( "There are %d errors\n", length(idx) )
        cat(mess)
        df  = cbind( DataISCCNBS[ idx, ], NumberLookup=block$Number[idx] )
        print( df )
        }

    return( out )
    }

# David H. Brainard
# Cone Contrast and Opponent Modulation Color Spaces
# pp. 563
# PART IV: CONE CONTRAST AND
# OPPONENT MODULATION COLOR SPACES
# DKL Example
# 7/6/95 dhb Wrote it.
# September 2014 Jose Gama converted the code to R
LMS2DKL<-function(bg,diffcone.coords,DKL2LMS=FALSE){
# bg= the background vector for the conversion
# diffcone.coords = the vector we wish to convert
# STEP 1, 2: GET bg and diffcone.coords
# STEP 3: Set M.raw as in equation A.4.9.
# This is found by inserting the background
# values into equation A.4.8. Different
# backgrounds produce different matrices.
# The MATLAB notation below just
# fills the desired 3-by-3 matrix.
M.raw <- matrix(c(1,1,0,1,-bg[1]/bg[2],0,-1,-1,(bg[1]+bg[2])/bg[3]),3,3,byrow=TRUE)
# STEP 4: Compute the inverse of M for
# equation A.4.10. The MATLAB inv() function
# computes the matrix inverse of its argument.
M.raw.inv <- solve(M.raw)
# STEP 5: Find the three isolating stimuli as
# the columns of M.inv.raw. The MATLAB
# notation X[,i] extracts the i-th column
# of the matrix X.
isochrom.raw <- M.raw.inv[,1]
rgisolum.raw <- M.raw.inv[,2]
sisolum.raw <- M.raw.inv[,3]
# STEP 6: Find the pooled cone contrast of each
# of these. The MATLAB norm() function returns
# the vector length of its argument. The MATLAB
# / operation represents entry-by-entry division.
isochrom.raw.pooled <- norm(isochrom.raw / bg,'f')
rgisolum.raw.pooled <- norm(rgisolum.raw / bg,'f')
sisolum.raw.pooled <- norm(sisolum.raw / bg,'f')
# STEP 7: Scale each mechanism isolating
# modulation by its pooled contrast to obtain
# mechanism isolating modulations that have
# unit length.
isochrom.unit <-  matrix(isochrom.raw/ isochrom.raw.pooled)
rgisolum.unit <-  matrix(rgisolum.raw/ rgisolum.raw.pooled)
sisolum.unit <- matrix(sisolum.raw / sisolum.raw.pooled)
# STEP 8: Compute the values of the normalizing
# constants by plugging t~e unit isolating stimuli
# into A.4.9 and seeing what we get. Each vector
# should have only one non-zero entry. The size
# of the entry is the response of the unscaled
# mechanism to the stimulus that should give unit
# response.
lum.resp.raw <- M.raw %*% isochrom.unit
l.minus.m.resp.raw <- M.raw %*% rgisolum.unit
s.minus.lum.resp.raw <- M.raw %*% sisolum.unit
# STEP 9: We need to rescale the rows of M.raw
# so that we get unit response. This means
# mUltiplying each row of M.raw by a constant.
# The easiest way to accomplish the multiplication
# is to form a diagonal matrix with the desired
# scalars on the diagonal. These scalars are just
# the multiplicative inverses of the non-zero
# entries of the vectors obtained in the previous
# step. The resulting matrix M provides the
# entries of A.4.11. The three .resp vectors
# computed should be the three unit vectors
# (and they are).
D.rescale <- matrix(c(1/lum.resp.raw[1],0,0,0,1/l.minus.m.resp.raw[2],0,0,0,1/s.minus.lum.resp.raw[3]),3,3,byrow=TRUE)
M <- D.rescale %*% M.raw
lum.resp <- M %*% isochrom.unit
l.minus.m.resp <- M %*% rgisolum.unit
s.minus.lum.resp <- M %*% sisolum.unit
# STEP 10: Compute the inverse of M to obtain
# the matrix in equation A.4.12.
if (!DKL2LMS) {
M.inv <- solve(M)
M<-M.inv
}
# STEP 11: Multiply the vector we wish to
# convert byM to obtain its DKL coordinates.
DKL.coords <- M %*% diffcone.coords
# STEP 12: convert to spherical coordinates.
# According to the conventions in the original DKL
# paper, azimuth of 0 is along our rgisolum axis,
# azimuth of 90 is along our negative sisolum
# axis. The isochromatic axis has an elevation
# of 90 degrees. To do the conversion. we flip the
# sign of the sisolum coordinate and then do a
# standard conversion to polar coordinates .
RADS.TO.DEGS <- 360/(2*pi)
azimuth.rads <- atan(pracma::mrdivide(-DKL.coords[3], DKL.coords[2]))
isolum.len <- sqrt(DKL.coords[2]^2 + DKL.coords[3]^2)
elevation.rads <- atan(pracma::mrdivide(DKL.coords[1], isolum.len))
azimuth <- RADS.TO.DEGS %*% azimuth.rads
elevation <- RADS.TO.DEGS %*% elevation.rads
list(azimuth.rads=azimuth.rads,isolum.len=isolum.len,elevation.rads=elevation.rads,azimuth=azimuth,elevation=elevation)
}

XYZtoRGB<-function(xc,yc,zc, ColorSystem=c(0.67,0.33,0.21,0.71,0.14,0.08,0.310,0.316)){
# simplified XYZ to RGB conversion, mostly for plots
xr <- ColorSystem[1]
yr <- ColorSystem[2]
zr <- 1.0 - xr - yr
xg <- ColorSystem[3]
yg <- ColorSystem[4]
zg <- 1.0 - xg - yg
xb <- ColorSystem[5]
yb <- ColorSystem[6]
zb <- 1.0 - xb - yb
d <- xr*yg*zb - xg*yr*zb - xr*yb*zg + xb*yr*zg + xg*yb*zr - xb*yg*zr
R <- (-xg*yc*zb + xc*yg*zb + xg*yb*zc - xb*yg*zc - xc*yb*zg + xb*yc*zg) / d
G <- ( xr*yc*zb - xc*yr*zb - xr*yb*zc + xb*yr*zc + xc*yb*zr - xb*yc*zr) / d
B <- ( xr*yg*zc - xg*yr*zc - xr*yc*zg + xc*yr*zg + xg*yc*zr - xc*yg*zr) / d
cbind(R,G,B)
}

chromaticity.diagram<-function(chromaticityCoordinates=colorscience::cccie31, conversionFunction=NULL,...){
# plot the chromaticity diagram AKA "horse shoe"
# conversionFunction CIE1931XYZ2CIE1976uv
pLen<-length(chromaticityCoordinates[["wlnm"]])
n=seq(1,pLen,by=1)
x<-chromaticityCoordinates$x[n]
y<-chromaticityCoordinates$y[n]
if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
z<-conversionFunction(cbind(chromaticityCoordinates$x[n],chromaticityCoordinates$y[n],chromaticityCoordinates$z[n]))
x<-z[,1]
y<-z[,2]
}
dots <-  list(...)
nameDots<-names(dots)
lNameDot<-length(dots)>0
if (!lNameDot) dots <- modifyList(dots, list(x=x,y=y,type='l')) else dots <- modifyList(dots, list(x=x,y=y ))
if (!(any(nameDots=='xlab'))) dots <- modifyList(dots, list(xlab='x'))
if (!(any(nameDots=='ylab'))) dots <- modifyList(dots, list(ylab='y'))
if (!(any(nameDots=='type'))) dots <- modifyList(dots, list(type='l'))
do.call(plot, dots ) # horseshoe
segments(x[1],y[1],x[pLen],y[pLen]) # line
}

Maxwell.triangle<-function(primariesRGB=colorscience::whitepointsRGB, conversionFunction=NULL,...){
# plot the Maxwell triangle
# conversionFunction CIE1931XYZ2CIE1976uv
x<-as.numeric(primariesRGB[1,c('xRed','xGreen','xBlue')])
y<-as.numeric(primariesRGB[1,c('yRed','yGreen','yBlue')])
if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
z<-1-x-y
z<-conversionFunction(cbind(x,y,z))
x<-z[,1]
y<-z[,2]
}
dots <- list(...)
nameDots<-names(dots)
lNameDot<-length(dots)>0
if (!lNameDot) dots <- modifyList(dots, list(x=c(x,x[1]),y=c(y,y[1]),type='l')) else dots <- modifyList(dots, list(x=c(x,x[1]),y=c(y,y[1])))
if (!(any(nameDots=='xlab'))) dots <- modifyList(dots, list(xlab='x'))
if (!(any(nameDots=='ylab'))) dots <- modifyList(dots, list(ylab='y'))
if (!(any(nameDots=='type'))) dots <- modifyList(dots, list(type='l'))
do.call(plot, dots )
}

chromaticity.diagram.color<-function(chromaticityCoordinates=colorscience::cccie31, conversionFunction=NULL,granularity=10,...){
# plot the chromaticity diagram AKA "horse shoe"
# conversionFunction CIE1931XYZ2CIE1976uv
pLen<-length(chromaticityCoordinates[["wlnm"]])
n=seq(1,pLen,by=1)
wlnmVector<-seq(chromaticityCoordinates[1,"wlnm"],chromaticityCoordinates[pLen,"wlnm",],by=1/granularity)
pLen2<-length(wlnmVector)
XYZ<-matrix(unlist(wlnm2xyz(wlnmVector)),ncol=3,byrow=FALSE)
x<-XYZ[,1]
y<-XYZ[,2]
z<-XYZ[,3]
X<-x
Y<-y
Z<-z
if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
az<-conversionFunction(cbind(x,y,z))
x<-az[,1]
y<-az[,2]
z<-1-x-y
}
temp<-XYZtoRGB(X,Y,Z)/255
t1<-temp/max(temp)*255
t1<-round(t1)
t1[which(t1<0)]<-0
t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
dots <- list(...)
nameDots<-names(dots)
lNameDot<-length(dots)>0
if (!lNameDot) dots <- modifyList(dots, list(x=x,y=y,type='p')) else dots <- modifyList(dots, list(x=x,y=y))
if (!(any(nameDots=='xlab'))) dots <- modifyList(dots, list(xlab='x'))
if (!(any(nameDots=='ylab'))) dots <- modifyList(dots, list(ylab='y'))
if (!(any(nameDots=='xlim'))) dots <- modifyList(dots, list(xlim=0:1))
if (!(any(nameDots=='ylim'))) dots <- modifyList(dots, list(ylim=0:1))
if (!(any(nameDots=='col'))) dots <- modifyList(dots, list(col=t2))
if (!(any(nameDots=='pch'))) dots <- modifyList(dots, list(pch=20))
if (!(any(nameDots=='type'))) dots <- modifyList(dots, list(type='p'))
do.call(plot, dots );par(new=TRUE)
#line
x1<-seq(x[1],x[pLen2],length.out=pLen2/7)#
y1<-seq(y[1],y[pLen2],length.out=pLen2/7)#
X1<-seq(X[1],X[pLen2],length.out=pLen2/7)#
Y1<-seq(Y[1],Y[pLen2],length.out=pLen2/7)#
temp<-XYZtoRGB(X1,Y1,1-X1-Y1)/255
t1<-temp/max(temp)*255
t1<-round(t1)
t1[which(t1<0)]<-0
t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
dots <- modifyList(dots, list(x=x1))
dots <- modifyList(dots, list(y=y1))
dots <- modifyList(dots, list(col=t2))
do.call(plot, dots )
}

Maxwell.triangle.color<-function(primariesRGB=colorscience::whitepointsRGB, conversionFunction=NULL,granularity=10,...){
# Maxwell triangle
# conversionFunction CIE1931XYZ2CIE1976uv
pLen<-100
n=seq(1,pLen,by=1)
x<-as.numeric(primariesRGB[1,c('xRed','xGreen','xBlue')])
y<-as.numeric(primariesRGB[1,c('yRed','yGreen','yBlue')])
x1<-seq(x[1],x[2],length.out=pLen)
y1<-seq(y[1],y[2],length.out=pLen)
x2<-seq(x[1],x[3],length.out=pLen)
y2<-seq(y[1],y[3],length.out=pLen)
x3<-seq(x[2],x[3],length.out=pLen)
y3<-seq(y[2],y[3],length.out=pLen)
x<-c(x1,x2,x3)
y<-c(y1,y2,y3)
z<-1-x-y
X<-x
Y<-y
Z<-z
if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
z<-conversionFunction(cbind(x,y,z))
x<-z[,1]
y<-z[,2]
}
temp<-XYZtoRGB(X,Y,Z)/255
t1<-temp/max(temp)*255
t1<-round(t1)
t1[which(t1<0)]<-0
t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
dots <- list(...)
nameDots<-names(dots)
lNameDot<-length(dots)>0
if (!lNameDot) dots <- modifyList(dots, list(x=x,y=y,type='p')) else dots <- modifyList(dots, list(x=x,y=y))
if (!(any(nameDots=='xlab'))) dots <- modifyList(dots, list(xlab='x'))
if (!(any(nameDots=='ylab'))) dots <- modifyList(dots, list(ylab='y'))
if (!(any(nameDots=='xlim'))) dots <- modifyList(dots, list(xlim=0:1))
if (!(any(nameDots=='ylim'))) dots <- modifyList(dots, list(ylim=0:1))
if (!(any(nameDots=='col'))) dots <- modifyList(dots, list(col=t2))
if (!(any(nameDots=='pch'))) dots <- modifyList(dots, list(pch=20))
if (!(any(nameDots=='type'))) dots <- modifyList(dots, list(type='p'))
do.call(plot, dots)
}

Maxwell.triangle.color.fill<-function(chromaticityCoordinates=colorscience::cccie31, conversionFunction=NULL,granularity=10,...){
# plot the Maxwell triangle
# conversionFunction CIE1931XYZ2CIE1976uv
pLen<-length(chromaticityCoordinates[["wlnm"]])
n=seq(1,pLen,by=1)
wlnmVector<-seq(chromaticityCoordinates[1,"wlnm"],chromaticityCoordinates[pLen,"wlnm",],by=1/granularity)
pLen2<-length(wlnmVector)
XYZ<-matrix(unlist(wlnm2xyz(wlnmVector)),ncol=3,byrow=FALSE)
x<-XYZ[,1]
y<-XYZ[,2]
z<-XYZ[,3]
X<-x
Y<-y
Z<-z
if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
az<-conversionFunction(cbind(x,y,z))
x<-az[,1]
y<-az[,2]
}
maxx<-max(x)
maxy<-max(y)
minx<-min(x)
miny<-min(y)
pLen2<-10*granularity;k=1
mx<-seq(minx,maxx, length.out=pLen2)
my<-seq(miny,maxy, length.out=pLen2)
p2<-matrix(0,pLen2^2,2)
p3<-vector(mode='character',pLen2^2)
for (n in 1:pLen2)# triangle
for (m in 1:pLen2){
   XC=mx[m]
   YC=my[n]
   ZC=1-XC-YC
   if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
az<-conversionFunction(cbind(XC,YC,ZC))
XC<-az[,1]
YC<-az[,2]
ZC=1-XC-YC
}
   temp<-XYZtoRGB(mx[m],my[n],1-mx[m]-my[n])/255
t1<-temp/max(temp)*255
t1<-round(t1)
if (any(t1<0)) t2<-NA else t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
   p2[k,1]<- XC
   p2[k,2]<- YC
   p3[k]<- t2
   k<-k+1
}
dots <- list(...)
nameDots<-names(dots)
lNameDot<-length(dots)>0
if (!lNameDot) dots <- modifyList(dots, list(x=p2[,1],y=p2[,2],col=p3,type='p')) else dots <- modifyList(dots, list(x=p2[,1],y=p2[,2],col=p3))
if (!(any(nameDots=='xlab'))) dots <- modifyList(dots, list(xlab='x'))
if (!(any(nameDots=='ylab'))) dots <- modifyList(dots, list(ylab='y'))
if (!(any(nameDots=='xlim'))) dots <- modifyList(dots, list(xlim=0:1))
if (!(any(nameDots=='ylim'))) dots <- modifyList(dots, list(ylim=0:1))
if (!(any(nameDots=='col'))) dots <- modifyList(dots, list(col=p3))
if (!(any(nameDots=='pch'))) dots <- modifyList(dots, list(pch=15))
if (!(any(nameDots=='type'))) dots <- modifyList(dots, list(type='p'))
do.call(plot, dots );par(new=TRUE)
}

chromaticity.diagram.color.fill<-function(chromaticityCoordinates=colorscience::cccie31, conversionFunction=NULL,granularity=10, conversionFunctionInv=NULL,...){
# plot the chromaticity diagram AKA "horse shoe"
# conversionFunction CIE1931XYZ2CIE1976uv
pLen<-length(chromaticityCoordinates[["wlnm"]])
#n=seq(1,pLen,by=1)
wlnmVector<-seq(chromaticityCoordinates[1,"wlnm"],chromaticityCoordinates[pLen,"wlnm",],by=1/granularity)
pLen2<-length(wlnmVector)
XYZ<-matrix(unlist(wlnm2xyz(wlnmVector)),ncol=3,byrow=FALSE)
x<-XYZ[,1]
y<-XYZ[,2]
z<-XYZ[,3]
X<-x
Y<-y
Z<-z
if (is.function(conversionFunction)) {
if (!is.function(conversionFunctionInv)){
mystrfun<-as.character(substitute(conversionFunction))
if (mystrfun=='CIE1931xy2CIE1976uv') conversionFunctionInv<-colorscience::CIE1976uv2CIE1931xy
if (mystrfun=='CIE1931xy2CIE1960uv') conversionFunctionInv<-colorscience::CIE1960UCS2xy # CIE1960uv2CIE1931xy
if (mystrfun=='CIE1976uv2CIE1931xy') conversionFunctionInv<-colorscience::CIE1931xy2CIE1976uv
if (mystrfun=='CIE1960uv2CIE1931xy') conversionFunctionInv<-colorscience::CIE1931xy2CIE1960uv
}
az<-conversionFunction(cbind(x,y,z))
x<-az[,1]
y<-az[,2]
}
boundariesXY<-cbind(x,y)
#line
#x1<-seq(x[1],x[length(x)],length.out=pLen2)
#y1<-seq(y[1],y[length(y)],length.out=pLen2)
#limits
maxx<-max(X)
maxy<-max(Y)
minx<-min(X)
miny<-min(Y)
pLen2<-10*granularity;k=1
dots <- list(...)
nameDots<-names(dots)
lNameDot<-length(dots)>0
boolDots<-FALSE
if ((any(nameDots=='xlim'))) { minx<-dots$xlim[1]
maxx<-dots$xlim[2]
boolDots<-TRUE
}
if ((any(nameDots=='ylim'))) { miny<-dots$ylim[1]
maxy<-dots$ylim[2]
boolDots<-TRUE
}
containedInTheHorseShoe<-FALSE
if (boolDots){
if (is.function(conversionFunctionInv)){
minx<-conversionFunctionInv(c(minx, miny))[1]
maxx<-conversionFunctionInv(c(maxx, maxy))[1]
miny<-conversionFunctionInv(c(minx, miny))[2]
maxy<-conversionFunctionInv(c(maxx, maxy))[2]
}
# are xlim, ylim contained in the horse shoe?
if (all(point.in.polygon(c(minx, maxx, maxx, minx),c(miny, maxy, miny, maxy),boundariesXY[,1],boundariesXY[,2]))) containedInTheHorseShoe<-TRUE
else { # 
maxx<-max(X)
maxy<-max(Y)
minx<-min(X)
miny<-min(Y)
}
}
mx<-seq(minx,maxx, length.out=pLen2)
pLen3<-round(pLen2 /( (maxx-minx)/(maxy-miny) ))
my<-seq(miny,maxy, length.out=pLen3)
# p2<-matrix(0,pLen2*pLen3,2)
# p3<-vector(mode='character',pLen2*pLen3)
# for (n in 1:pLen3)# Horse Shoe
# for (m in 1:pLen2){
#    XC=mx[m]
#    YC=my[n]
#    ZC=1-XC-YC
#    if (!is.null(conversionFunction)) if (is.function(conversionFunction)) {
# az<-conversionFunction(cbind(XC,YC,ZC))
# XC<-az[,1]
# YC<-az[,2]
# ZC=1-XC-YC
# }
# temp<-XYZtoRGB(mx[m],my[n],1-mx[m]-my[n])/255
# t1<-temp/max(temp)*255
# t1<-round(t1)
# 
# #boundariesXY<-rbind(boundariesXY,cbind(x1,y1))
# #boundariesXY<-rbind(boundariesXY,boundariesXY[1,])
# if (any(t1<0)) t1[which(t1<0)]<-1
# #if (!in.out(boundariesXY,cbind(XC,YC))) t2<-NA else t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
# #if (!pnt.in.poly(cbind(XC,YC),boundariesXY)) t2<-NA else t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
# 
# if (containedInTheHorseShoe) {
# t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
# } else {
# if (!point.in.polygon(XC,YC,boundariesXY[,1],boundariesXY[,2])) t2<-NA else t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
# }
# 
#    p2[k,1]<- XC
#    p2[k,2]<- YC
#    p3[k]<- t2
#    k<-k+1
# }
#vectorized version
n  = rep(1:pLen3,each=pLen2)
m  = rep(1:pLen2,pLen3)
   XC=mx[m]
   YC=my[n]
   ZC=1-XC-YC
if (is.function(conversionFunction)) {
az<-conversionFunction(cbind(XC,YC,ZC))
XC<-az[,1]
YC<-az[,2]
ZC=1-XC-YC
}
temp<-XYZtoRGB(mx[m],my[n],1-mx[m]-my[n])/255
t1<-temp/apply(temp,1,max)*255
t1<-round(t1)
if (any(t1<0)) t1[which(t1<0)]<-1
if (containedInTheHorseShoe) {
t2<-apply(t1,1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
} else {
boolInPolygon<-point.in.polygon(XC,YC,boundariesXY[,1],boundariesXY[,2])
t2<-rep(NA,pLen2*pLen3)
t2[which(boolInPolygon==1)]<-apply(t1[which(boolInPolygon==1),],1,function(x) sprintf("#%02X%02X%02X",x[1],x[2],x[3]))
}
p2<-cbind(XC,YC)
p3<-t2
if (!lNameDot) dots <- modifyList(dots, list(x=p2[,1],y=p2[,2],col=p3,type='p')) else dots <- modifyList(dots, list(x=p2[,1],y=p2[,2],col=p3))
if (!(any(nameDots=='xlab'))) dots <- modifyList(dots, list(xlab='x'))
if (!(any(nameDots=='ylab'))) dots <- modifyList(dots, list(ylab='y'))
if (!(any(nameDots=='xlim'))) dots <- modifyList(dots, list(xlim=0:1))
if (!(any(nameDots=='ylim'))) dots <- modifyList(dots, list(ylim=0:1))
#if (!(any(nameDots=='col'))) dots <- modifyList(dots, list(col=p3))
if (!(any(nameDots=='pch'))) dots <- modifyList(dots, list(pch=15))
if (!(any(nameDots=='type'))) dots <- modifyList(dots, list(type='p'))
do.call(plot, dots );par(new=TRUE)
}







LUV2LAB<-function(Luvmatrix) XYZ2Lab(Luv2XYZ(Luvmatrix))
LAB2LUV<-function(Labmatrix) XYZ2Luv(Lab2XYZ(Labmatrix))
XYZ2xyY <- function(XYZmatrix) {
rSum<-apply(XYZmatrix,1,sum)
cbind(XYZmatrix[,1]/ rSum ,XYZmatrix[,2]/ rSum,XYZmatrix[,3])
}

makeChromaticAdaptationMatrix<-function(ChromaticAdaptationAlgorithm='VonKries', illuminantSource='C', illuminantDestination='D65', observer=2, ChromaticAdaptationArray=colorscience::ChromaticAdaptation, referenceWhiteArray=colorscience::XYZperfectreflectingdiffuser)
{# Chromatic Adaptation
which(referenceWhiteArray[,"Illuminant"]==illuminantDestination)
if (observer==2) observerPos<-2:4 else observerPos<-5:7
coneResponseDomainSource<-ChromaticAdaptationArray[,,ChromaticAdaptationAlgorithm,'direct'] %*% 
matrix(unlist(referenceWhiteArray[which(referenceWhiteArray[,"Illuminant"]==illuminantSource),observerPos]),3,1,byrow=TRUE)
coneResponseDomainDestination<-ChromaticAdaptationArray[,,ChromaticAdaptationAlgorithm,'direct'] %*% 
matrix(unlist(referenceWhiteArray[which(referenceWhiteArray[,"Illuminant"]==illuminantDestination),observerPos]),3,1,byrow=TRUE)
d <- diag(3)
diag(d) <- coneResponseDomainDestination/coneResponseDomainSource
ChromaticAdaptationArray[,,ChromaticAdaptationAlgorithm,'inverse'] %*%  d   %*%  ChromaticAdaptationArray[,,ChromaticAdaptationAlgorithm,'direct'] #  adaptation matrix
}

footcandle2lux<-function(ftcl) ftcl * 10.7639104167
# converts foot candle to Lumens/lux
# http://www.translatorscafe.com/cafe/EN/units-converter/illumination
footcandle2watt.sqcentimeter<-function(ftcl) ftcl*0.00000157597517082
# converts foot candle to watts / square centimeter [w/cm^2] (at 555 nm) 
# http://www.translatorscafe.com/cafe/EN/units-converter/illumination
footcandle2candela.steradian.sqmeter<-function(ftcl) ftcl * 10763.9104167 
# converts foot candle to candela steradian / square meter [cd*sr/m^2]
# http://www.translatorscafe.com/cafe/EN/units-converter/illumination

dkl2rgb <- function(dklMatrix, conversionMatrix=NA)
{ # Convert from DKL color space (Derrington, Krauskopf & Lennie) to RGB.
# Package psychopy for Python
# http://www.psychopy.org/epydoc/psychopy.misc-pysrc.html#dkl2rgb
if (!is.matrix(dklMatrix)) dklMatrix <- matrix(dklMatrix, ncol=3,byrow=TRUE)
elev = dklMatrix[,1]
azim = dklMatrix[,2]
radius = dklMatrix[,3]
RGBYLUM <- cbind(radius * sin(pi/180*(elev)), radius * cos(pi/180*(elev))*cos(pi/180*(azim)), radius * cos(pi/180*(elev))*sin(pi/180*(azim)) )

if (is.na(conversionMatrix)) conversionMatrix <- matrix(c(
  #LUMIN      %L-M    %L+M-S  (note that dkl has to be in cartesian coords first!) 
  1.0000, 1.0000, -0.1462,#R 
  1.0000, -0.3900, 0.2094,#G 
  1.0000, 0.0180, -1.0000#B 
),3,3,byrow=TRUE)

rgb <- conversionMatrix %*% t(RGBYLUM)
return (t(rgb))#return in the shape we received it
}

dklCart2rgb <- function(dklMatrix, conversionMatrix=NA)
{ # Like dkl2rgb except that it uses cartesian coords (LM,S,LUM) rather than
# spherical coords for DKL (elev, azim, contr)
if (!is.matrix(dklMatrix)) dklMatrix <- matrix(dklMatrix, ncol=3,byrow=TRUE)
if (is.na(conversionMatrix)) conversionMatrix <- matrix(c(
  #LUMIN      %L-M    %L+M-S  (note that dkl has to be in cartesian coords first!) 
  1.0000, 1.0000, -0.1462,#R 
  1.0000, -0.3900, 0.2094,#G 
  1.0000, 0.0180, -1.0000#B 
),3,3,byrow=TRUE)

rgb <- conversionMatrix %*% t(dklMatrix)
return (t(rgb))#return in the shape we received it
}

rgb2dklCart <- function(rgbMatrix, conversionMatrix=NA)
{ # Converts RGB into Cartesian DKL space
if (!is.matrix(rgbMatrix)) rgbMatrix <- matrix(rgbMatrix, ncol=3,byrow=TRUE)
if (is.na(conversionMatrix)) conversionMatrix <- matrix(c(
0.25145542,  0.64933633,  0.09920825,
0.78737943, -0.55586618, -0.23151325,
0.26562825,  0.63933074, -0.90495899
),3,3,byrow=TRUE)

dkl <- conversionMatrix %*% t(rgbMatrix)
return (t(dkl))#return in the shape we received it
}



rgb2dklV <- function(RGB){
# convert RGB to DKL
# Graph-Based Visual Saliency (MATLAB source code)
# Jonathan Harel
# California Institute of Technology
# http://www.klab.caltech.edu/~harel/share/gbvs.php
# constants used for RGB -> DKL conversion:
lutRGB <- matrix(c(
    0.024935, 0.0076954, 0.042291,
    0.024974, 0.0077395, 0.042346,
    0.025013, 0.0077836, 0.042401,
    0.025052, 0.0078277, 0.042456,
    0.025091, 0.0078717, 0.042511,
    0.02513, 0.0079158, 0.042566,
    0.025234, 0.007992, 0.042621,
    0.025338, 0.0080681, 0.042676,
    0.025442, 0.0081443, 0.042731,
    0.025545, 0.0082204, 0.042786,
    0.025649, 0.0082966, 0.042841,
    0.025747, 0.0084168, 0.042952,
    0.025844, 0.0085371, 0.043062,
    0.025942, 0.0086573, 0.043172,
    0.026039, 0.0087776, 0.043282,
    0.026136, 0.0088978, 0.043392,
    0.026234, 0.0090581, 0.043502,
    0.026331, 0.0092184, 0.043612,
    0.026429, 0.0093788, 0.043722,
    0.026526, 0.0095391, 0.043833,
    0.026623, 0.0096994, 0.043943,
    0.026818, 0.0099198, 0.044141,
    0.027013, 0.01014, 0.044339,
    0.027208, 0.010361, 0.044537,
    0.027403, 0.010581, 0.044736,
    0.027597, 0.010802, 0.044934,
    0.027857, 0.010994, 0.04522,
    0.028117, 0.011186, 0.045507,
    0.028377, 0.011379, 0.045793,
    0.028636, 0.011571, 0.046079,
    0.028896, 0.011764, 0.046366,
    0.029104, 0.012068, 0.046652,
    0.029312, 0.012373, 0.046938,
    0.029519, 0.012677, 0.047225,
    0.029727, 0.012982, 0.047511,
    0.029935, 0.013287, 0.047797,
    0.030273, 0.013663, 0.048326,
    0.03061, 0.01404, 0.048855,
    0.030948, 0.014417, 0.049383,
    0.031286, 0.014794, 0.049912,
    0.031623, 0.01517, 0.050441,
    0.032156, 0.015707, 0.051035,
    0.032688, 0.016244, 0.05163,
    0.033221, 0.016782, 0.052225,
    0.033753, 0.017319, 0.052819,
    0.034286, 0.017856, 0.053414,
    0.034961, 0.018693, 0.054493,
    0.035636, 0.019531, 0.055573,
    0.036312, 0.020369, 0.056652,
    0.036987, 0.021206, 0.057731,
    0.037662, 0.022044, 0.058811,
    0.038623, 0.023246, 0.060044,
    0.039584, 0.024449, 0.061278,
    0.040545, 0.025651, 0.062511,
    0.041506, 0.026854, 0.063744,
    0.042468, 0.028056, 0.064978,
    0.043857, 0.029659, 0.066806,
    0.045247, 0.031263, 0.068634,
    0.046636, 0.032866, 0.070463,
    0.048026, 0.034469, 0.072291,
    0.049416, 0.036072, 0.074119,
    0.051221, 0.038156, 0.076476,
    0.053026, 0.04024, 0.078833,
    0.054831, 0.042325, 0.081189,
    0.056636, 0.044409, 0.083546,
    0.058442, 0.046493, 0.085903,
    0.06039, 0.048737, 0.087996,
    0.062338, 0.050982, 0.090088,
    0.064286, 0.053226, 0.092181,
    0.066234, 0.055471, 0.094273,
    0.068182, 0.057715, 0.096366,
    0.070519, 0.06012, 0.098921,
    0.072857, 0.062525, 0.10148,
    0.075195, 0.06493, 0.10403,
    0.077532, 0.067335, 0.10659,
    0.07987, 0.069739, 0.10914,
    0.082208, 0.072345, 0.11176,
    0.084545, 0.07495, 0.11438,
    0.086883, 0.077555, 0.117,
    0.089221, 0.08016, 0.11963,
    0.091558, 0.082766, 0.12225,
    0.094026, 0.085611, 0.12533,
    0.096494, 0.088457, 0.12841,
    0.098961, 0.091303, 0.1315,
    0.10143, 0.094148, 0.13458,
    0.1039, 0.096994, 0.13767,
    0.10688, 0.10028, 0.14119,
    0.10987, 0.10357, 0.14471,
    0.11286, 0.10685, 0.14824,
    0.11584, 0.11014, 0.15176,
    0.11883, 0.11343, 0.15529,
    0.12208, 0.11695, 0.15903,
    0.12532, 0.12048, 0.16278,
    0.12857, 0.12401, 0.16652,
    0.13182, 0.12754, 0.17026,
    0.13506, 0.13106, 0.17401,
    0.1387, 0.13499, 0.17819,
    0.14234, 0.13892, 0.18238,
    0.14597, 0.14285, 0.18656,
    0.14961, 0.14677, 0.19075,
    0.15325, 0.1507, 0.19493,
    0.15727, 0.15519, 0.19956,
    0.1613, 0.15968, 0.20419,
    0.16532, 0.16417, 0.20881,
    0.16935, 0.16866, 0.21344,
    0.17338, 0.17315, 0.21806,
    0.17805, 0.17796, 0.22291,
    0.18273, 0.18277, 0.22775,
    0.1874, 0.18758, 0.2326,
    0.19208, 0.19238, 0.23744,
    0.19675, 0.19719, 0.24229,
    0.20156, 0.20224, 0.24758,
    0.20636, 0.20729, 0.25286,
    0.21117, 0.21234, 0.25815,
    0.21597, 0.21739, 0.26344,
    0.22078, 0.22244, 0.26872,
    0.2261, 0.22806, 0.27423,
    0.23143, 0.23367, 0.27974,
    0.23675, 0.23928, 0.28524,
    0.24208, 0.24489, 0.29075,
    0.2474, 0.2505, 0.29626,
    0.25299, 0.25651, 0.3022,
    0.25857, 0.26253, 0.30815,
    0.26416, 0.26854, 0.3141,
    0.26974, 0.27455, 0.32004,
    0.27532, 0.28056, 0.32599,
    0.28156, 0.28697, 0.33238,
    0.28779, 0.29339, 0.33877,
    0.29403, 0.2998, 0.34515,
    0.30026, 0.30621, 0.35154,
    0.30649, 0.31263, 0.35793,
    0.3126, 0.31904, 0.36388,
    0.3187, 0.32545, 0.36982,
    0.32481, 0.33186, 0.37577,
    0.33091, 0.33828, 0.38172,
    0.33701, 0.34469, 0.38767,
    0.34325, 0.3511, 0.39361,
    0.34948, 0.35752, 0.39956,
    0.35571, 0.36393, 0.40551,
    0.36195, 0.37034, 0.41145,
    0.36818, 0.37675, 0.4174,
    0.37429, 0.38317, 0.42313,
    0.38039, 0.38958, 0.42885,
    0.38649, 0.39599, 0.43458,
    0.3926, 0.4024, 0.44031,
    0.3987, 0.40882, 0.44604,
    0.40494, 0.41523, 0.45198,
    0.41117, 0.42164, 0.45793,
    0.4174, 0.42806, 0.46388,
    0.42364, 0.43447, 0.46982,
    0.42987, 0.44088, 0.47577,
    0.43623, 0.44689, 0.48128,
    0.4426, 0.45291, 0.48678,
    0.44896, 0.45892, 0.49229,
    0.45532, 0.46493, 0.4978,
    0.46169, 0.47094, 0.5033,
    0.46792, 0.47695, 0.50837,
    0.47416, 0.48297, 0.51344,
    0.48039, 0.48898, 0.5185,
    0.48662, 0.49499, 0.52357,
    0.49286, 0.501, 0.52863,
    0.49805, 0.50701, 0.53392,
    0.50325, 0.51303, 0.53921,
    0.50844, 0.51904, 0.54449,
    0.51364, 0.52505, 0.54978,
    0.51883, 0.53106, 0.55507,
    0.52442, 0.53667, 0.55969,
    0.53, 0.54228, 0.56432,
    0.53558, 0.5479, 0.56894,
    0.54117, 0.55351, 0.57357,
    0.54675, 0.55912, 0.57819,
    0.55182, 0.56433, 0.58304,
    0.55688, 0.56954, 0.58789,
    0.56195, 0.57475, 0.59273,
    0.56701, 0.57996, 0.59758,
    0.57208, 0.58517, 0.60242,
    0.57675, 0.58998, 0.60639,
    0.58143, 0.59479, 0.61035,
    0.5861, 0.5996, 0.61432,
    0.59078, 0.60441, 0.61828,
    0.59545, 0.60922, 0.62225,
    0.60065, 0.61403, 0.62709,
    0.60584, 0.61884, 0.63194,
    0.61104, 0.62365, 0.63678,
    0.61623, 0.62846, 0.64163,
    0.62143, 0.63327, 0.64648,
    0.62584, 0.63808, 0.65088,
    0.63026, 0.64289, 0.65529,
    0.63468, 0.6477, 0.65969,
    0.63909, 0.65251, 0.6641,
    0.64351, 0.65731, 0.6685,
    0.64857, 0.66132, 0.67269,
    0.65364, 0.66533, 0.67687,
    0.6587, 0.66934, 0.68106,
    0.66377, 0.67335, 0.68524,
    0.66883, 0.67735, 0.68943,
    0.67273, 0.68136, 0.69361,
    0.67662, 0.68537, 0.6978,
    0.68052, 0.68938, 0.70198,
    0.68442, 0.69339, 0.70617,
    0.68831, 0.69739, 0.71035,
    0.69221, 0.7022, 0.7141,
    0.6961, 0.70701, 0.71784,
    0.7, 0.71182, 0.72159,
    0.7039, 0.71663, 0.72533,
    0.70779, 0.72144, 0.72907,
    0.71169, 0.72505, 0.73348,
    0.71558, 0.72866, 0.73789,
    0.71948, 0.73226, 0.74229,
    0.72338, 0.73587, 0.7467,
    0.72727, 0.73948, 0.7511,
    0.73247, 0.74349, 0.75507,
    0.73766, 0.74749, 0.75903,
    0.74286, 0.7515, 0.763,
    0.74805, 0.75551, 0.76696,
    0.75325, 0.75952, 0.77093,
    0.75714, 0.76393, 0.77599,
    0.76104, 0.76834, 0.78106,
    0.76494, 0.77275, 0.78612,
    0.76883, 0.77715, 0.79119,
    0.77273, 0.78156, 0.79626,
    0.77792, 0.78677, 0.80132,
    0.78312, 0.79198, 0.80639,
    0.78831, 0.79719, 0.81145,
    0.79351, 0.8024, 0.81652,
    0.7987, 0.80762, 0.82159,
    0.80519, 0.81283, 0.82687,
    0.81169, 0.81804, 0.83216,
    0.81818, 0.82325, 0.83744,
    0.82468, 0.82846, 0.84273,
    0.83117, 0.83367, 0.84802,
    0.83636, 0.83888, 0.85286,
    0.84156, 0.84409, 0.85771,
    0.84675, 0.8493, 0.86256,
    0.85195, 0.85451, 0.8674,
    0.85714, 0.85972, 0.87225,
    0.86364, 0.86613, 0.87819,
    0.87013, 0.87255, 0.88414,
    0.87662, 0.87896, 0.89009,
    0.88312, 0.88537, 0.89604,
    0.88961, 0.89178, 0.90198,
    0.8961, 0.8986, 0.90947,
    0.9026, 0.90541, 0.91696,
    0.90909, 0.91222, 0.92445,
    0.91558, 0.91904, 0.93194,
    0.92208, 0.92585, 0.93943,
    0.92857, 0.93307, 0.94493,
    0.93506, 0.94028, 0.95044,
    0.94156, 0.94749, 0.95595,
    0.94805, 0.95471, 0.96145,
    0.95455, 0.96192, 0.96696,
    0.96364, 0.96954, 0.97357,
    0.97273, 0.97715, 0.98018,
    0.98182, 0.98477, 0.98678,
    0.99091, 0.99238, 0.99339,
    1, 1, 1 ), ncol=3,byrow=TRUE)

lms0 <- c( 34.918538957799996, 19.314796676499999, 0.585610818500000 )

m <- c( 18.32535,  44.60077,   7.46216, 4.09544,  28.20135,   6.66066, 0.02114,   0.10325,   1.05258 )

fac <- 1.0 / (lms0[1] + lms0[2])
mm <- c( sqrt(3.0)*fac, sqrt(3.0)*fac, 0.0, sqrt(lms0[1]*lms0[1]+lms0[2]*lms0[2])/lms0[1]*fac, 
-sqrt(lms0[1]*lms0[1]+lms0[2]*lms0[2])/lms0[2]*fac, 0.0, -fac, -fac, (lms0[1] + lms0[2]) / lms0[3] * fac )

# do a lookup RGB -> rgb:
if ( max(RGB)<=1 ){
  RGB <- ceiling( RGB*255 )
  RGB[which(RGB<1)] <- 1
  RGB[which(RGB>256)] <- 256
} else   RGB <- RGB + 1

aa1 <- lutRGB[ RGB[,1] , 1 ]
aa2 <- lutRGB[ RGB[,2] , 2 ]
aa3 <- lutRGB[ RGB[,3] , 3 ]

# now convert to LMS:
lms1 <- m[1] * aa1 + m[2] * aa2 + m[3] * aa3 - lms0[1]
lms2 <- m[4] * aa1 + m[5] * aa2 + m[6] * aa3 - lms0[2]
lms3 <- m[7] * aa1 + m[8] * aa2 + m[9] * aa3 - lms0[3]

# finally to DKL:
dkl1 <- mm[1] * lms1 + mm[2] * lms2 + mm[3] * lms3
dkl2 <- mm[4] * lms1 + mm[5] * lms2 + mm[6] * lms3
dkl3 <- mm[7] * lms1 + mm[8] * lms2 + mm[9] * lms3

# finally to DKLn:
dkl <- c( dkl1 * 0.5774, dkl2 * 2.7525, dkl3 * 0.4526 )
dkl
}

dkl2dklCart<- function(dklMatrix)
{ # Convert from spherical coords to cartesian coords
if (!is.matrix(dklMatrix)) dklMatrix <- matrix(dklMatrix, ncol=3,byrow=TRUE)
elev = dklMatrix[,1]
azim = dklMatrix[,2]
radius = dklMatrix[,3]
cbind(radius * cos(pi/180*(elev))*sin(pi/180*(azim)), radius * cos(pi/180*(elev))*cos(pi/180*(azim)), radius * sin(pi/180*(elev)) )
}

dklCart2dkl <- function(dklMatrix)
{ # Convert from cartesian coords to spherical coords
if (!is.matrix(dklMatrix)) dklMatrix <- matrix(dklMatrix, ncol=3,byrow=TRUE)
n = dklMatrix[,1]
m = dklMatrix[,2]
p = dklMatrix[,3]
radius <- sqrt(n^2+m^2+p^2)
elevation <- asin(p/radius)
azimuth <- atan2(n,m)
elevation <- elevation *180/pi
azimuth <- azimuth *180/pi
cbind(elevation,azimuth,radius)
}

#CIE 1976 Luv to u', v' CIE 1976
Luv2Yuv<-function(Luvmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser) 
{
if (!is.matrix(Luvmatrix)) Luvmatrix <- matrix(Luvmatrix, ncol=3,byrow=TRUE)
XYZ2Yuv(Luv2XYZ(Luvmatrix,illuminant,observer,RefWhite))
}
#CIE u', v' CIE 1976 to CIE 1976 Luv
Yuv2Luv<-function(Yu.v.matrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser) 
{
if (!is.matrix(Yu.v.matrix)) Yu.v.matrix <- matrix(Yu.v.matrix, ncol=3,byrow=TRUE)
XYZ2Luv(Yuv2XYZ(Yu.v.matrix),illuminant,observer,RefWhite)
}
# chromaticity coordinates u', v' CIE 1976
# source: CIE 1976 color space, From Wikipedia, the free encyclopedia http://en.wikipedia.org/wiki/CIE_1976_color_space
XYZ2Yuv<-function(XYZmatrix) {
if (!is.matrix(XYZmatrix)) XYZmatrix <- matrix(XYZmatrix, ncol=3,byrow=TRUE)
cbind(Y=XYZmatrix[,2],u76=4*XYZmatrix[,1]/(XYZmatrix[,1]+15*XYZmatrix[,2]+3*XYZmatrix[,3]), v76=9*XYZmatrix[,2]/(XYZmatrix[,1]+15*XYZmatrix[,2]+3*XYZmatrix[,3]))
}
Yxy2Yuv<-function(Yxymatrix){
if (!is.matrix(Yxymatrix)) Yxymatrix <- matrix(Yxymatrix, ncol=3,byrow=TRUE)
cbind(Y=Yxymatrix[,1], u76=4*Yxymatrix[,2]/(-2*Yxymatrix[,2]+12*Yxymatrix[,3]+3), v76=9*Yxymatrix[,3]/(-2*Yxymatrix[,2]+12*Yxymatrix[,3]+3)) 
}
Yuv2XYZ<-function(Yu.v.matrix) {
if (!is.matrix(Yu.v.matrix)) Yu.v.matrix <- matrix(Yu.v.matrix, ncol=3,byrow=TRUE)
cbind(X=9*Yu.v.matrix[,1]*Yu.v.matrix[,2]/(4*Yu.v.matrix[,3]), Y=Yu.v.matrix[,1], Z=3*Yu.v.matrix[,1]*(4-Yu.v.matrix[,2])/(4*Yu.v.matrix[,3])-5*Yu.v.matrix[,1])
}
Yuv2xy<-function(Yu.v.matrix){
if (!is.matrix(Yu.v.matrix)) Yu.v.matrix <- matrix(Yu.v.matrix, ncol=3,byrow=TRUE)
cbind(Y=Yu.v.matrix[,1], x=9*Yu.v.matrix[,2]/(6*Yu.v.matrix[,2]-16*Yu.v.matrix[,3]+12), y=4*Yu.v.matrix[,3]/(6*Yu.v.matrix[,2]-16*Yu.v.matrix[,3]+12))
}
# CIE-1976  u', v' to CIE-1931 x, y
# http://www.color-theory-phenomena.nl/10.03.htm
CIE1976uv2CIE1931xy<-function(uvmatrix) {
if (!is.matrix(uvmatrix)) uvmatrix <- matrix(uvmatrix, ncol=2,byrow=TRUE)
cbind( x = ( 27 * uvmatrix[,1] / 4 ) / ( ( 9 * uvmatrix[,1] / 2) - 12 * uvmatrix[,2] + 9 ), y = ( 3 * uvmatrix[,2]) / ( ( 9 * uvmatrix[,1] / 2) - 12 * uvmatrix[,2] + 9 )) 
}
#CIE-1976  u', v' to  CIE-1960 u, v
CIE1976uv2CIE1960uv<-function(uvmatrix){
if (!is.matrix(uvmatrix)) uvmatrix <- matrix(uvmatrix, ncol=2,byrow=TRUE)
cbind(u = uvmatrix[,1], v = ( 2 * uvmatrix[,2] / 3 ))
}
#CIE-1931 X, Y, Z to CIE-1931 x, y, z
CIE1931XYZ2CIE1931xyz<-function(XYZmatrix) {
if (!is.matrix(XYZmatrix)) XYZmatrix <- matrix(XYZmatrix, ncol=3,byrow=TRUE)
cbind(x = XYZmatrix[,1] / ( XYZmatrix[,1] + XYZmatrix[,2] + XYZmatrix[,3] ), 
y = XYZmatrix[,2] / ( XYZmatrix[,1] + XYZmatrix[,2] + XYZmatrix[,3] ), z = 1 - ( XYZmatrix[,1] / ( XYZmatrix[,1] + XYZmatrix[,2] + XYZmatrix[,3] ) + XYZmatrix[,2] / ( XYZmatrix[,1] + XYZmatrix[,2] + XYZmatrix[,3] ) ) )
}

#CIE-1931 X, Y, Z to CIE-1960 u, v
CIE1931XYZ2CIE1960uv<-function(XYZmatrix) {
if (!is.matrix(XYZmatrix)) XYZmatrix <- matrix(XYZmatrix, ncol=3,byrow=TRUE)
cbind(u = 4 * XYZmatrix[,1] / ( XYZmatrix[,1] + 15 * XYZmatrix[,2] + 3 * XYZmatrix[,3] ), 
v = 6 * XYZmatrix[,2] / ( XYZmatrix[,1] + 15 * XYZmatrix[,2] + 3 * XYZmatrix[,3] ))
}
#CIE-1931 x, y to CIE-1960 u, v
CIE1931xy2CIE1960uv<-function(xymatrix){
if (!is.matrix(xymatrix)) xymatrix <- matrix(xymatrix, ncol=2,byrow=TRUE)
cbind(u = 4 * xymatrix[,1] / ( -2 * xymatrix[,1] + 12 * xymatrix[,2] + 3 ), v = 6 * xymatrix[,2] / ( -2 * xymatrix[,1] + 12 * xymatrix[,2] + 3 ))
}
#CIE-1931 X, Y, Z to CIE-1976 u', v'
CIE1931XYZ2CIE1976uv<-function(XYZmatrix){
if (!is.matrix(XYZmatrix)) XYZmatrix <- matrix(XYZmatrix, ncol=3,byrow=TRUE)
cbind(u = 4 * XYZmatrix[,1] / ( XYZmatrix[,1] + 15 * XYZmatrix[,2] + 3 * XYZmatrix[,3] ), 
v = 9 * XYZmatrix[,2] / ( XYZmatrix[,1] + 15 * XYZmatrix[,2] + 3 * XYZmatrix[,3] ))
}
#CIE-1931 x, y to CIE-1976 u', v'
CIE1931xy2CIE1976uv<-function(xymatrix){
if (!is.matrix(xymatrix)) xymatrix <- matrix(xymatrix, ncol=2,byrow=TRUE)
cbind(u = 4 * xymatrix[,1] / ( -2 * xymatrix[,1] + 12 * xymatrix[,2] + 3 ), 
v = 9 * xymatrix[,2] / ( -2 * xymatrix[,1] + 12 * xymatrix[,2] + 3 ))
}
#CIE 1960 color space (CIE 1960 UCS)
# source: CIE 1964 color space, From Wikipedia, the free encyclopedia http://en.wikipedia.org/wiki/CIE_1960_color_space
Yxy2CIE1960UCS<-function(Yxymatrix){
if (!is.matrix(Yxymatrix)) Yxymatrix <- matrix(Yxymatrix, ncol=3,byrow=TRUE)
cbind(u=(0.4661*Yxymatrix[,1]+0.1593*Yxymatrix[,2])/(Yxymatrix[,2]-0.15735*Yxymatrix[,1]+0.2424), v=(0.6581**Yxymatrix[,2])/(Yxymatrix[,2]-0.15735*Yxymatrix[,1]+0.2424))
}
CIE1960UCS2xy<-function(uvMatrix){
if (!is.matrix(uvMatrix)) uvMatrix <- matrix(uvMatrix, ncol=2,byrow=TRUE)
cbind(x=3*uvMatrix[,1]/(2*uvMatrix[,1]-8*uvMatrix[,2]+4), y=2*uvMatrix[,2]/(2*uvMatrix[,1]-8*uvMatrix[,2]+4))
}
# CIE 1964 color space
# source: CIE 1964 color space, From Wikipedia, the free encyclopedia http://en.wikipedia.org/wiki/CIE_1964_color_space
CIE1960UCS2CIE1964<-function(uvYmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser) 
{
if (is.null(dim(uvYmatrix))) if (length(uvYmatrix)>2) uvYmatrix<-matrix(uvYmatrix, ncol=3,byrow=TRUE)
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
xymatrix <- CIE1960UCS2xy(uvYmatrix[,1:2])
XYZmatrix <- xyY2XYZ(cbind(xymatrix,uvYmatrix[,3]))
xr <- XYZmatrix[,1] / Rx * 100
yr <- XYZmatrix[,2] / Ry * 100
zr <- XYZmatrix[,3] / Rz * 100
W <- 25*yr^(1/3)-17
U <- 13*W*(uvYmatrix[,1]-Rx)
V <- 13*uvYmatrix[,3]*(uvYmatrix[,2]-Ry)
cbind(U=U,V=V ,W=W )
}

XYZ2BVR <-function(XYZmatrix){# matrix for XYZ to Landolt BVR transformation 
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
XYZmatrix %*% matrix(c(0.200308,-0.017427,0.012508,0.013093,0.789019,-0.328110,-0.224010,-0.465060,1.381969),3,3,byrow=TRUE)
}

BVR2XYZ <-function(BVRmatrix){# matrix for Landolt BVR to XYZ transformation 
if (is.null(dim(BVRmatrix))) if (length(BVRmatrix)>2) BVRmatrix<-matrix(BVRmatrix, ncol=3,byrow=TRUE)
BVRmatrix %*% matrix(c(4.961444,0.293121,0.902866, 0.096640,1.479324,0.513487,-0.021961,0.348571,0.837347),3,3,byrow=TRUE)
}

RGB2LEF<-function(RGBmatrix)
{# based on: Kang, Henry R, 2006, Computational color technology, Spie Press Bellingham
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-t(matrix(RGBmatrix, ncol=3,byrow=TRUE))
matrix(c(2/3,2/3,2/3,1,-1/2,-1/2,0,sqrt(3)/2,-sqrt(3)/2),3,3,byrow=TRUE) %*% RGBmatrix
}

LEF2RGB<-function(LEFmatrix)
{# based on: Kang, Henry R, 2006, Computational color technology, Spie Press Bellingham
if (is.null(dim(LEFmatrix))) if (length(LEFmatrix)>2) LEFmatrix<-t(matrix(LEFmatrix, ncol=3,byrow=TRUE))
matrix(c(1/2,2/3,0,1/2,-1/3,1/sqrt(3),1/2,-1/3,-1/sqrt(3)),3,3,byrow=TRUE) %*% LEFmatrix
}

RGB2YPbPr<-function(RGBmatrix)
{# based on: ColorObject.pm Graphics/ColorObject version 0.5.0, Copyright 2003-2005 by Alex Izvorski, (Portions Copyright 2001-2003 by Alfred Reibenschuh)
# reference: http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.txt
# Y is [0..1], Pb and Pr are [-0.5..0.5]
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
matrix(c(0.299,0.587,0.114,-0.168736,-0.331264,0.5,0.5,-0.418688,-0.081312),3,3,byrow=TRUE) %*% t(RGBmatrix)
}

YPbPr2RGB<-function(YPbPrmatrix)
{# based on: ColorObject.pm Graphics/ColorObject version 0.5.0, Copyright 2003-2005 by Alex Izvorski, (Portions Copyright 2001-2003 by Alfred Reibenschuh)
# reference: http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.txt
# result is NTSC non-linear rgb
# Pb and Pr are [-0.5..0.5], Y is [0..1]
if (is.null(dim(YPbPrmatrix))) if (length(YPbPrmatrix)>2) YPbPrmatrix<-matrix(YPbPrmatrix, ncol=3,byrow=TRUE)
matrix(c(1.0,0.0,1.402,1.0,-0.344136,-0.714136,1.0,1.772,0.0),3,3,byrow=TRUE) %*% t(YPbPrmatrix)
}

RGB2YCbCr<-function(RGBmatrix)
{# based on: ColorObject.pm Graphics/ColorObject version 0.5.0, Copyright 2003-2005 by Alex Izvorski, (Portions Copyright 2001-2003 by Alfred Reibenschuh)
# reference: http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.txt
# input should be NTSC non-linear rgb
# Y is [16..235], Cb and Cr are [16..239]
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
matrix(c(65.481,128.553,24.966,-37.797,-74.203,112.0 ,112.0,-93.786,-18.214),3,3,byrow=TRUE) %*% t(RGBmatrix)
}

YCbCr2RGB<-function(YPbPrmatrix)
{# based on: ColorObject.pm Graphics/ColorObject version 0.5.0, Copyright 2003-2005 by Alex Izvorski, (Portions Copyright 2001-2003 by Alfred Reibenschuh)
# reference: http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.txt
# Y is [16..235], Cb and Cr are [16..239]
if (is.null(dim(YPbPrmatrix))) if (length(YPbPrmatrix)>2) YPbPrmatrix<-matrix(YPbPrmatrix, ncol=3,byrow=TRUE)
matrix(c(0.00456621,0.0,0.00625893, 0.00456621,-0.00153632,-0.00318811,0.00456621, 0.00791071, 0.0),3,3,byrow=TRUE) %*% t(YPbPrmatrix)
}

RGB2YIQ<-function(RGBmatrix)
{
# assumes RGBmatrix in the range [0,1]
# returns YIQmatrix in the range [-1,1]
# based on Color space conversion by Sophie Kirschner, http://www.blitzbasic.com/codearcs/codearcs.php?code=2953
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
YIQmatrix<-matrix(0,dim(RGBmatrix)[1],3)
YIQmatrix[,1]<- 0.299000*RGBmatrix[,1]+0.587000*RGBmatrix[,2]+0.114000*RGBmatrix[,3]
YIQmatrix[,2]<- 0.595716*RGBmatrix[,1]-0.274453*RGBmatrix[,2]-0.321264*RGBmatrix[,3]
YIQmatrix[,3]<- 0.211456*RGBmatrix[,1]-0.522591*RGBmatrix[,2]+0.311350*RGBmatrix[,3]
YIQmatrix
}

YIQ2RGB<-function(YIQmatrix)
{
# assumes YIQmatrix in the range [-1,1]
# returns RGBmatrix in the range [0,1]
# based on Color space conversion by Sophie Kirschner, http://www.blitzbasic.com/codearcs/codearcs.php?code=2953
if (is.null(dim(YIQmatrix))) if (length(YIQmatrix)>2) YIQmatrix<-matrix(YIQmatrix, ncol=3,byrow=TRUE)
RGBmatrix<-matrix(0,dim(YIQmatrix)[1],3)
RGBmatrix[,1]<- YIQmatrix[,1]+0.9563*YIQmatrix[,2]+0.6210*YIQmatrix[,3]
RGBmatrix[,2]<- YIQmatrix[,1]-0.2721*YIQmatrix[,2]-0.6474*YIQmatrix[,3]
RGBmatrix[,3]<- YIQmatrix[,1]-1.1070*YIQmatrix[,2]+1.7046*YIQmatrix[,3]
RGBmatrix
}

RGB2YUV<-function(RGBmatrix)
{
# assumes RGBmatrix in the range [0,1]
# returns YUVmatrix in the range [-1,1]
# based on Color space conversion by Madk, Sophie Kirschner, http://www.blitzbasic.com/codearcs/codearcs.php?code=2953
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
yuv_uconst <- 1.0/0.436
yuv_vconst <- 1.0/0.615
YUVmatrix<-matrix(0,dim(RGBmatrix)[1],3)
YUVmatrix[,1]<- 0.299*RGBmatrix[,1]+0.587*RGBmatrix[,2]+0.114*RGBmatrix[,3]
YUVmatrix[,2]<- (0.492*(RGBmatrix[,3]-YUVmatrix[,1]))*yuv_uconst
YUVmatrix[,3]<- (0.877*(RGBmatrix[,1]-YUVmatrix[,1]))*yuv_vconst
YUVmatrix[,1]<- YUVmatrix[,1]*2-1
YUVmatrix
}

YUV2RGB<-function(YUVmatrix)
{
# assumes YUVmatrix in the range [-1,1]
# returns RGBmatrix in the range [0,1]
# based on Color space conversion by Madk, Sophie Kirschner, http://www.blitzbasic.com/codearcs/codearcs.php?code=2953
if (is.null(dim(YUVmatrix))) if (length(YUVmatrix)>2) YUVmatrix<-matrix(YUVmatrix, ncol=3,byrow=TRUE)
yuv_uconst <- 1.0/0.436
yuv_vconst <- 1.0/0.615
ly <- (YUVmatrix[,1]+1)/2.0 #YUVmatrix[,2]/yuv_uconst,YUVmatrix[,3]/yuv_vconst)
RGBmatrix<-matrix(0,dim(YUVmatrix)[1],3)
RGBmatrix[,1]<- ly +1.140*YUVmatrix[,3]/yuv_vconst
RGBmatrix[,2]<- ly -0.395*YUVmatrix[,2]/yuv_uconst - 0.581*YUVmatrix[,3]/yuv_vconst
RGBmatrix[,3]<- ly +2.032*YUVmatrix[,2]/yuv_uconst
RGBmatrix
}

RGB2HSL<-function(RGBmatrix)
{ # code adapted from: easyrgb Color conversion math and formulas http://www.easyrgb.com/
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
var.R <- ( RGBmatrix[,1] / 255 )                     #RGB values <- From 0 to 255
var.G <- ( RGBmatrix[,2] / 255 )
var.B <- ( RGBmatrix[,3] / 255 )
var.Min <- apply(RGBmatrix / 255,1,min)    #Min. value of RGB
var.Max <- apply(RGBmatrix / 255,1,max)    #Max. value of RGB
del.Max <- var.Max - var.Min             #Delta RGB value
L <- ( var.Max + var.Min ) / 2
H<-vector('numeric',length(L))
S<-H
LS1.max<-which((del.Max != 0) & (L < 0.5))
LS2.max<-which((del.Max != 0) & (L >= 0.5))
if (length(LS1.max)>0) S[LS1.max] <- del.Max[LS1.max] / ( var.Max[LS1.max] + var.Min[LS1.max] )
if (length(LS2.max)>0) S[LS2.max] <- del.Max[LS2.max] / ( 2 - var.Max[LS2.max] - var.Min[LS2.max] )
del.R <- ( ( ( var.Max - var.R ) / 6 ) + ( del.Max / 2 ) ) / del.Max
del.G <- ( ( ( var.Max - var.G ) / 6 ) + ( del.Max / 2 ) ) / del.Max
del.B <- ( ( ( var.Max - var.B ) / 6 ) + ( del.Max / 2 ) ) / del.Max
R.max<-which((del.Max != 0) & (var.R == var.Max))
G.max<-which((del.Max != 0) & (var.G == var.Max))
B.max<-which((del.Max != 0) & (var.B == var.Max))
if (length(R.max)>0) H[R.max] <- del.B[R.max] - del.G[R.max]
if (length(G.max)>0) H[G.max] <- ( 1 / 3 ) + del.R[G.max] - del.B[G.max]
if (length(B.max)>0) H[B.max] <- ( 2 / 3 ) + del.G[B.max] - del.R[B.max]
Hl0<-which((del.Max != 0) & (H < 0))
if (length(Hl0)>0) H[Hl0] <- H[Hl0] + 1
Hg1<-which((del.Max != 0) & (H > 1))
if (length(Hg1)>0) H[Hg1] <- H[Hg1] - 1
cbind(H=H,S=S,L=L)
}

HSL2RGB<-function(HSLmatrix)
{
if (is.null(dim(HSLmatrix))) if (length(HSLmatrix)>2) HSLmatrix<-matrix(HSLmatrix, ncol=3,byrow=TRUE)
R<-vector('numeric',dim(HSLmatrix)[2])
G<-R
B<-R
S0<-which(HSLmatrix[,2] == 0)
if (length(S0)>0) {
R[S0] <- HSLmatrix[S0,1] * 255
G[S0] <- HSLmatrix[S0,2] * 255
B[S0] <- HSLmatrix[S0,3] * 255
}
S1<-which(HSLmatrix[,2] != 0)
if (length(S1)>0) {
var.2<-vector('numeric',length(S1))
Ll05<-which((HSLmatrix[S1,3] < 0.5))
if (length(Ll05)>0) var.2 <- HSLmatrix[Ll05,3] * ( 1 + HSLmatrix[Ll05,2] )
Lg05<-which((HSLmatrix[S1,3] >= 0.5))
if (length(Lg05)>0) var.2 <- ( HSLmatrix[Lg05,2] + HSLmatrix[Lg05,3] ) - ( HSLmatrix[Lg05,2] * HSLmatrix[Lg05,3] )
var.1 <- 2 * HSLmatrix[S1,3] - var.2
R[S1] <- 255 * Hue.2.RGB( var.1, var.2, HSLmatrix[S1,1] + ( 1 / 3 ) ) 
G[S1] <- 255 * Hue.2.RGB( var.1, var.2, HSLmatrix[S1,1] )
B[S1] <- 255 * Hue.2.RGB( var.1, var.2, HSLmatrix[S1,1] - ( 1 / 3 ) )
}
HSLret<-cbind(R=R,G=G,B=B)
HSLret[1:dim(HSLmatrix)[1],]
}

Hue.2.RGB<-function(v1, v2, vH)
{
S0<-which(vH < 0)
vH[S0] <- vH[S0] + 1
S1<-which(vH > 0)
vH[S1] <- vH[S1] - 1
h6<-which(( 6 * vH ) < 1)
v1[h6] <- ( v1[h6] + ( v2[h6] - v1[h6] ) * 6 * vH[h6] )
h2<-which(( 2 * vH ) < 1)
v1[h2] <- ( v1[h2] + ( v2[h2] - v1[h2] ) * 6 * vH[h2] )
h3<-which(( 3 * vH ) < 2)
v1[h3] <- ( v1[h3] + ( v2[h3] - v1[h3] ) * 6 * vH[h3] )
v1
}

RGB2HSV<-function(RGBmatrix)
{ # code adapted from: easyrgb Color conversion math and formulas http://www.easyrgb.com/
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
var.R <- ( RGBmatrix[,1] / 255 )                     #RGB values <- From 0 to 255
var.G <- ( RGBmatrix[,2] / 255 )
var.B <- ( RGBmatrix[,3] / 255 )
var.Min <- apply(RGBmatrix / 255,1,min)    #Min. value of RGB
var.Max <- apply(RGBmatrix / 255,1,max)    #Max. value of RGB
del.Max <- var.Max - var.Min             #Delta RGB value
V <- var.Max
H<-vector('numeric',length(V))
S<- H
Sn0<-which(del.Max != 0)
if (length(Sn0)>0) S[Sn0]<- del.Max[Sn0] / var.Max[Sn0]
del.R <- ( ( ( var.Max - var.R ) / 6 ) + ( del.Max / 2 ) ) / del.Max
del.G <- ( ( ( var.Max - var.G ) / 6 ) + ( del.Max / 2 ) ) / del.Max
del.B <- ( ( ( var.Max - var.B ) / 6 ) + ( del.Max / 2 ) ) / del.Max
R.max<-which((del.Max != 0) & (var.R == var.Max))
G.max<-which((del.Max != 0) & (var.G == var.Max))
B.max<-which((del.Max != 0) & (var.B == var.Max))
if (length(R.max)>0) H[R.max] <- del.B[R.max] - del.G[R.max]
if (length(G.max)>0) H[G.max] <- ( 1 / 3 ) + del.R[G.max] - del.B[G.max]
if (length(B.max)>0) H[B.max] <- ( 2 / 3 ) + del.G[B.max] - del.R[B.max]
Hl0<-which((del.Max != 0) & (H < 0))
if (length(Hl0)>0) H[Hl0] <- H[Hl0] + 1
Hg1<-which((del.Max != 0) & (H > 1))
if (length(Hg1)>0) H[Hg1] <- H[Hg1] - 1
cbind(H=H,S=S,V=V)
}

HSV2RGB<-function(HSVmatrix)
{
if (is.null(dim(HSVmatrix))) if (length(HSVmatrix)>2) HSVmatrix<-matrix(HSVmatrix, ncol=3,byrow=TRUE)
R<-vector('numeric',dim(HSVmatrix)[2])
G<-R
B<-R
S0<-which(HSVmatrix[,2] == 0)
if (length(S0)>0) {
R[S0] <- HSVmatrix[S0,1] * 255
G[S0] <- HSVmatrix[S0,2] * 255
B[S0] <- HSVmatrix[S0,3] * 255
}
Sn0<-which(HSVmatrix[,2] != 0)
if (length(Sn0)>0) {
var.h<-vector('numeric',dim(HSVmatrix)[2])
var.i<-vector('numeric',dim(HSVmatrix)[2])
var.1<-vector('numeric',dim(HSVmatrix)[2])
var.2<-vector('numeric',dim(HSVmatrix)[2])
var.3<-vector('numeric',dim(HSVmatrix)[2])
var.r<-vector('numeric',dim(HSVmatrix)[2])
var.g<-vector('numeric',dim(HSVmatrix)[2])
var.b<-vector('numeric',dim(HSVmatrix)[2])
var.h[Sn0] <- HSVmatrix[Sn0,1] * 6
var.i[Sn0] <- floor( var.h[Sn0] )
var.1[Sn0] <- HSVmatrix[Sn0,3] * ( 1 - HSVmatrix[Sn0,2] )
var.2[Sn0] <- HSVmatrix[Sn0,3] * ( 1 - HSVmatrix[Sn0,2] * ( var.h[Sn0] - var.i[Sn0] ) )
var.3[Sn0] <- HSVmatrix[Sn0,3] * ( 1 - HSVmatrix[Sn0,2] * ( 1 - ( var.h[Sn0] - var.i[Sn0] ) ) )
ieq0<-which((HSVmatrix[,2] != 0) & (var.i == 0))
if (length(ieq0)>0) {
var.r[ieq0] <- HSVmatrix[ieq0]
var.g[ieq0] <- var.3[ieq0]
var.b[ieq0] <- var.1[ieq0]
}
ieq1<-which((HSVmatrix[,2] != 0) & (var.i == 1))
if (length(ieq1)>0) {
var.r[ieq1] <- var.2[ieq1]
var.g[ieq1] <- HSVmatrix[ieq1]
var.b[ieq1] <- var.1[ieq1]
}
ieq2<-which((HSVmatrix[,2] != 0) & (var.i == 2))
if (length(ieq2)>0) {
var.r[ieq2] <- var.1[ieq2]
var.g[ieq2] <- HSVmatrix[ieq2]
var.b[ieq2] <- var.3[ieq2]
}
ieq3<-which((HSVmatrix[,2] != 0) & (var.i == 3))
if (length(ieq3)>0) {
var.r[ieq3] <- var.1[ieq3]
var.g[ieq3] <- var.2[ieq3]
var.b[ieq3] <- HSVmatrix[ieq3]
}
ieq4<-which((HSVmatrix[,2] != 0) & (var.i == 4))
if (length(ieq4)>0) {
var.r[ieq4] <- var.3[ieq4]
var.g[ieq4] <- var.1[ieq4]
var.b[ieq4] <- HSVmatrix[ieq4]
}
ieq5<-which((HSVmatrix[,2] != 0) & !any(var.i %in% 1:4))
if (length(ieq5)>5) {
var.r[ieq5] <- HSVmatrix[ieq5]
var.g[ieq5] <- var.1[ieq5]
var.b[ieq5] <- var.2[ieq5]
}
R[Sn0] <- var.r[Sn0] * 255                  #RGB results <- From 0 to 255
G[Sn0] <- var.g[Sn0] * 255
B[Sn0] <- var.b[Sn0] * 255
}
HSVret<- cbind(R=R,G=G,B=B)
HSVret[1:dim(HSVmatrix)[1],]
}

RGB2CMY<-function(RGBmatrix)
{#RGB values <- From 0 to 255
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
C <- 1 - ( RGBmatrix[,1] / 255 )
M <- 1 - ( RGBmatrix[,2] / 255 )
Y <- 1 - ( RGBmatrix[,3] / 255 )
cbind(C=C,M=M,Y=Y)
}

CMY2RGB<-function(CMYmatrix)
{#CMY values <- From 0 to 1
if (is.null(dim(CMYmatrix))) if (length(CMYmatrix)>2) CMYmatrix<-matrix(CMYmatrix, ncol=3,byrow=TRUE)
R <- ( 1 - CMYmatrix[,1] ) * 255
G <- ( 1 - CMYmatrix[,2] ) * 255
B <- ( 1 - CMYmatrix[,3] ) * 255
cbind(R=R,G=G,B=B)
}

CMY2CMYK<-function(CMYmatrix)
{#CMY values <- From 0 to 1
if (is.null(dim(CMYmatrix))) if (length(CMYmatrix)>2) CMYmatrix<-matrix(CMYmatrix, ncol=3,byrow=TRUE)
var.K <- rep(1,dim(CMYmatrix)[1])
CK<-which(CMYmatrix[,1] < var.K)
if (length(CK)>0) var.K[CK] <- CMYmatrix[CK,1]
CM<-which(CMYmatrix[,2] < var.K)
if (length(CM)>0) var.K[CM] <- CMYmatrix[CM,2]
CY<-which(CMYmatrix[,3] < var.K)
if (length(CY)>0) var.K[CY] <- CMYmatrix[CY,3]
C <- ( CMYmatrix[,1] - var.K ) / ( 1 - var.K )
M <- ( CMYmatrix[,2] - var.K ) / ( 1 - var.K )
Y <- ( CMYmatrix[,3] - var.K ) / ( 1 - var.K )
cbind(C=C,M=M,Y=Y,K=var.K)
}

CMYK2CMY<-function(CMYKmatrix)
{#CMYK values <- From 0 to 1
if (is.null(dim(CMYKmatrix))) if (length(CMYKmatrix)>2) CMYKmatrix<-matrix(CMYKmatrix, ncol=4,byrow=TRUE)
C <- ( CMYKmatrix[,1] * ( 1 - CMYKmatrix[,4] ) + CMYKmatrix[,4] )
M <- ( CMYKmatrix[,2] * ( 1 - CMYKmatrix[,4] ) + CMYKmatrix[,4] )
Y <- ( CMYKmatrix[,3] * ( 1 - CMYKmatrix[,4] ) + CMYKmatrix[,4] )
cbind(C=C,M=M,Y=Y)
}

XYZ2HunterLab<-function(XYZmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{ # adapted from: easyrgb Color conversion math and formulas http://www.easyrgb.com/
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
xr <- XYZmatrix[,1] / Rx * 100
yr <- XYZmatrix[,2] / Ry * 100
zr <- XYZmatrix[,3] / Rz * 100
L<-10 * sqrt( yr )
a<-17.5 * ( ( ( 1.02 * xr ) - yr ) / L )
b<-7 * ( ( yr - ( 0.847 * zr ) ) / L )
cbind(L=L,a=a,b=b)
}

HunterLab2XYZ<-function(HunterLabmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{ # adapted from: easyrgb Color conversion math and formulas http://www.easyrgb.com/
if (is.null(dim(HunterLabmatrix))) if (length(HunterLabmatrix)>2) HunterLabmatrix<-matrix(HunterLabmatrix, ncol=3,byrow=TRUE)
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
Y<-HunterLabmatrix[,1] / 10
X<-HunterLabmatrix[,2] / 17.5 * Y
Z<-HunterLabmatrix[,3] / 7 * Y
Y<-Y^2
X<- ( X + Y ) / 1.02
Z<- -( Z - Y ) / 0.847
X <- X * Rx
Y <- Y * Ry
Z <- Z * Rz
cbind(X=X,Y=Y,Z=Z)
}

NickersonColorDifference<-function(MunsellHVC1, MunsellHVC2) {
# deltaE = 2/5CdH + 6dV + 3dC
# COLOR TECHNOLOGY in the textile industry Second Edition
m1<-MunsellHVC(MunsellHVC1)
m2<-MunsellHVC(MunsellHVC2)
C<- as.numeric(m1[,'C'])
dC<- C - as.numeric(m2[,'C'])
dV<- as.numeric(m1[,'V']) - as.numeric(m2[,'V'])
dH<-huedegree(m1[,'H']) - huedegree(m2[,'H'])
as.numeric(2*C*dH/5 + 6*dV + 3*dC)
}

CCT2XYZ<-function (CCTmatrix)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
C1 <- 2.0 * pi * 6.626176 * 2.99792458 * 2.99792458
C2 <- (6.626176 * 2.99792458) / 1.380662
nmSeq<-seq(360,830,by=5)
dWavelengthM <- nmSeq * 1.0e-3
dWavelengthM5 <- dWavelengthM * dWavelengthM * dWavelengthM * dWavelengthM * dWavelengthM
blackbody <- C1 / (dWavelengthM5 * 1.0e-12 * (exp(C2 / (CCTmatrix * dWavelengthM * 1.0e-3)) - 1.0))
XYZmatrix<-blackbody * colorscience::ciexyz31[nmSeq-359,c('xbar','ybar','zbar')]
XYZmatrix[,1]<-XYZmatrix[,1] / XYZmatrix[,2]
XYZmatrix[,3]<-XYZmatrix[,3] / XYZmatrix[,2]
XYZmatrix[,2]<-1.0
}

Lab2LCHab<-function(LabMatrix)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(LabMatrix))) if (length(LabMatrix)>2) LabMatrix<-matrix(LabMatrix, ncol=3,byrow=TRUE)
LCHabmatrix<-LabMatrix
LCHabmatrix[,2]<-sqrt(LabMatrix[,2]^2 + LabMatrix[,3]^2)
LCHabmatrix[,3]<-180.0 * atan2(LabMatrix[,3], LabMatrix[,2]) / pi
if (LCHabmatrix[,3] < 0.0) LCHabmatrix[,3]<-LCHabmatrix[,3]+360
LCHabmatrix
}

LCHab2Lab<-function(LCHabmatrix)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(LCHabmatrix))) if (length(LCHabmatrix)>2) LCHabmatrix<-matrix(LCHabmatrix, ncol=3,byrow=TRUE)
LabMatrix<-LCHabmatrix
LabMatrix[,2]<-LCHabmatrix[,2] * cos(LCHabmatrix[,3] * pi / 180.0)
LabMatrix[,3]<-LCHabmatrix[,2] * sin(LCHabmatrix[,3] * pi / 180.0)
LabMatrix
}

Luv2LCHuv<-function(LuvMatrix)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(LuvMatrix))) if (length(LuvMatrix)>2) LuvMatrix<-matrix(LuvMatrix, ncol=3,byrow=TRUE)
LCHuvmatrix<-LuvMatrix
LCHuvmatrix[,2]<-sqrt(LuvMatrix[,2]^2 + LuvMatrix[,3]^2)
LCHuvmatrix[,3]<-180.0 * atan2(LuvMatrix[,3], LuvMatrix[,2]) / pi
if (LCHuvmatrix[,3] < 0.0) LCHuvmatrix[,3]<-LCHuvmatrix[,3]+360
LCHuvmatrix
}

LCHuv2Luv<-function(LCHuvmatrix)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(LCHuvmatrix))) if (length(LCHuvmatrix)>2) LCHuvmatrix<-matrix(LCHuvmatrix, ncol=3,byrow=TRUE)
LuvMatrix<-LCHuvmatrix
LuvMatrix[,2]<-LCHuvmatrix[,2] * cos(LCHuvmatrix[,3] * pi / 180.0)
LuvMatrix[,3]<-LCHuvmatrix[,2] * sin(LCHuvmatrix[,3] * pi / 180.0)
LuvMatrix
}

DIN6167.YellownessIndex<-function(XYZmatrix,illuminant='C',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser) 
{ # source: Basic equations for optical properties, 
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
RxRyRz<-XYZ2RxRyRz(XYZmatrix,illuminant,observer,RefWhite)
100*( RxRyRz[,1] - RxRyRz[,3] ) / RxRyRz[,2]
}

ASTM.D1925.YellownessIndex<-function(XYZmatrix) 
{
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
(100*(1.28*XYZmatrix[,1]-1.06*XYZmatrix[,3]))/XYZmatrix[,2]
}
#ASTM D 1925 Yellowness Index for Plastics
#X, Y and Z are the tri-stimulus values for the calculated for illuminant C
# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012, pp. 6

ASTM.E313.YellownessIndex<-function(XYZmatrix)
{
#ASTM D 1925 Yellowness Index
#X, Y and Z are the tri-stimulus values for the calculated for illuminant C
# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012, pp. 6
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
(100*(1-0.847*XYZmatrix[,3]))/XYZmatrix[,2]
}

ASTM.E313.Whiteness<-function(XYZmatrix)
{
#ASTM D 1925 Whiteness Index
#X, Y and Z are the tri-stimulus values for the calculated for illuminant C
# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012, pp. 6
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
3.388*XYZmatrix[,3]-3*XYZmatrix[,2]
}

CIE.Whiteness<-function(xyYmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# strictly for D65 and 2 or 10 deg observer
#Values bigger than 100 indicate a bluish white
#Values smaller than 100 indicate a yellowish white
# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012
if (is.null(dim(xyYmatrix))) if (length(xyYmatrix)>2) xyYmatrix<-matrix(xyYmatrix, ncol=3,byrow=TRUE)
Rrgbwhitergb<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rrgbwhitex<-unlist(Rrgbwhitergb[paste('X',observer,sep='')])
Rrgbwhitey<-unlist(Rrgbwhitergb[paste('Y',observer,sep='')])
Rrgbwhitez<-unlist(Rrgbwhitergb[paste('Z',observer,sep='')])
xr<-Rrgbwhitex / (Rrgbwhitex + Rrgbwhitey + Rrgbwhitez)
yr<-Rrgbwhitey / (Rrgbwhitex + Rrgbwhitey + Rrgbwhitez)
xyYmatrix[,3]+800*(xr-xyYmatrix[,1])+1700*(yr-xyYmatrix[,2])
}

Berger59.Whiteness<-function(xyYmatrix,illuminant='C',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012
if (is.null(dim(xyYmatrix))) if (length(xyYmatrix)>2) xyYmatrix<-matrix(xyYmatrix, ncol=3,byrow=TRUE)
Rrgbwhitergb<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rrgbwhitex<-unlist(Rrgbwhitergb[paste('X',observer,sep='')])
Rrgbwhitey<-unlist(Rrgbwhitergb[paste('Y',observer,sep='')])
Rrgbwhitez<-unlist(Rrgbwhitergb[paste('Z',observer,sep='')])
xr<-Rrgbwhitex / (Rrgbwhitex + Rrgbwhitey + Rrgbwhitez)
zr<-Rrgbwhitez / (Rrgbwhitex + Rrgbwhitey + Rrgbwhitez)
0.333*xyYmatrix[,2]+125*(xyYmatrix[,3]/zr)-125*(xyYmatrix[,1]/xr)
}

Stensby68.Whiteness<-function(LabHunterMatrix) {
if (is.null(dim(LabHunterMatrix))) if (length(LabHunterMatrix)>2) LabHunterMatrix<-matrix(LabHunterMatrix, ncol=3,byrow=TRUE)
LabHunterMatrix[,1]-3*LabHunterMatrix[,3]+3*LabHunterMatrix[,2]
}
#L, a and b are Hunter Color Coordinates
# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012

Taube60.Whiteness<-function(XYZmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
Rrgbwhitergb<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rrgbwhitex<-unlist(Rrgbwhitergb[paste('X',observer,sep='')])
Rrgbwhitey<-unlist(Rrgbwhitergb[paste('Y',observer,sep='')])
Rrgbwhitez<-unlist(Rrgbwhitergb[paste('Z',observer,sep='')])
400*XYZmatrix[,3]/Rrgbwhitez-3*XYZmatrix[,2]
}

Hunter60.WhitenessIndex<-function(LabHunterMatrix) {
#L, a and b are Hunter Color Coordinates
# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012
if (is.null(dim(LabHunterMatrix))) if (length(LabHunterMatrix)>2) LabHunterMatrix<-matrix(LabHunterMatrix, ncol=3,byrow=TRUE)
LabHunterMatrix[,1]-3*LabHunterMatrix[,3]
}

GanzGrieser.Whiteness<-function(xyYmatrix)
{# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012
if (is.null(dim(xyYmatrix))) if (length(xyYmatrix)>2) xyYmatrix<-matrix(xyYmatrix, ncol=3,byrow=TRUE)
P <- -1868.322
Q <- -3695.690
C <- 1809.441
xyYmatrix[,3]+P*xyYmatrix[,1]+Q*xyYmatrix[,2]+C
}

GanzGrieser.Tint<-function(xyYmatrix)
{# Color iQC and Color iMatch Color Calculations Guide Version 8.0 July 2012
if (is.null(dim(xyYmatrix))) if (length(xyYmatrix)>2) xyYmatrix<-matrix(xyYmatrix, ncol=3,byrow=TRUE)
m <- -1001.223
n <- 748.366
k <- 68.261
m*xyYmatrix[,1]+n*xyYmatrix[,2]+k
}

CIETint<-function(xymatrix,illuminant='D65',observer=2)
{#Tint indices: CIE Tint and ASTM E313 Tint
#Tint E313 = Tint CIE = Tx (x_n - x) - 650 (y_n - y)
#CIE Publication 15:2004, "Colorimetry"
#ASTM E313, "Standard Practice for Calculating Yellowness and Whiteness Indices from Instrumentally Measured Color Coordinates"
if (is.null(dim(xymatrix))) if (length(xymatrix)>2) xymatrix<-matrix(xymatrix, ncol=3,byrow=TRUE)
x<-array(matrix(c(1000,1000,1000,900,900,900,0.3101,0.3457,0.3127,0.3104,0.3477,0.3138,0.3161,0.3585,0.3290,0.3191,0.3595,0.3310),3,6),c(3,2,3))
dimnames(x)<-list(c('C','D50','D65'), c('2','10'), c('T','x','y'))
Txy<-x[illuminant,observer,]
Txy[1]*(Txy[2]-xymatrix[,1])-650*(Txy[3]-xymatrix[,2])
}

DominantWavelength<-function(xyYmatrix, illuminant='D65',observer=2,RefWhiteIllum=colorscience::XYZperfectreflectingdiffuser)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(xyYmatrix))) if (length(xyYmatrix)>2) xyYmatrix<-matrix(xyYmatrix, ncol=3,byrow=TRUE)
Rrgbwhitergb<-RefWhiteIllum[which(RefWhiteIllum[["Illuminant"]]==illuminant ),]
Rrgbwhitex<-unlist(Rrgbwhitergb[paste('X',observer,sep='')])
Rrgbwhitey<-unlist(Rrgbwhitergb[paste('Y',observer,sep='')])
Rrgbwhitez<-unlist(Rrgbwhitergb[paste('Z',observer,sep='')])
xr<-Rrgbwhitex / (Rrgbwhitex + Rrgbwhitey + Rrgbwhitez)
yr<-Rrgbwhitey / (Rrgbwhitex + Rrgbwhitey + Rrgbwhitez)
count <- 1
tArray <- vector('numeric',2)# t
wArray <- vector('numeric',2)# wavelength
cArray <- vector('numeric',2)# cycle
a <- xyYmatrix[,1] - xr
b <- xyYmatrix[,2] - yr
if ((a >= -0.000001) & (a <= 0.000001) & (b >= -0.000001) & (b <= 0.000001)) stop('Error, (x, y) must be different from (xr, yr) at least by 1e-6')
for (nm in seq(360, 830, by=5))
{
i1 <- (nm - 360) / 5
i2 <- ifelse(nm == 830, 0, i1 + 1)
nm2 <- 5 * i2 + 360
x1 <- colorscience::ciexyz31[i1*5 + 1,'xbar'] / (colorscience::ciexyz31[i1*5 + 1,'xbar']+colorscience::ciexyz31[i1*5 + 1,'ybar']+colorscience::ciexyz31[i1*5 + 1,'zbar'])
y1 <- colorscience::ciexyz31[i1*5 + 1,'ybar'] / (colorscience::ciexyz31[i1*5 + 1,'xbar']+colorscience::ciexyz31[i1*5 + 1,'ybar']+colorscience::ciexyz31[i1*5 + 1,'zbar'])
x2 <- colorscience::ciexyz31[i2*5 + 1,'xbar'] / (colorscience::ciexyz31[i2*5 + 1,'xbar']+colorscience::ciexyz31[i2*5 + 1,'ybar']+colorscience::ciexyz31[i2*5 + 1,'zbar'])
y2 <- colorscience::ciexyz31[i2*5 + 1,'ybar'] / (colorscience::ciexyz31[i2*5 + 1,'xbar']+colorscience::ciexyz31[i2*5 + 1,'ybar']+colorscience::ciexyz31[i2*5 + 1,'zbar'])
C <- x1 - xr
d <- y1 - yr
E <- x2 - x1
f <- y2 - y1
s <- (a * d - b * C) / (b * E - a * f)
#cat(nm,i1,i2,nm2,x1,y1,x2,y2,C,d,E,f,s,'\n')
if ((s < 0.0) | (s >= 1.0)) next
Tval <- ifelse(abs(a) >= abs(b),((E * s + C) / a), ((f * s + d) / b))
tArray[count] <- Tval
cArray[count] <- nm
wArray[count] <- (nm2 - nm) * s + nm
count <- count + 1
}
if ((cArray[2] == 830) & (tArray[2] > 0.0)) dominantWavelength <- -wArray[1] else dominantWavelength <- 
ifelse(tArray[1] >= 0.0, wArray[1], wArray[2])
dominantWavelength
}

XYZ2RGB<-function(XYZmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser,RGBModel='sRGB',RefWhiteRGB=colorscience::whitepointsRGB,gamma=NA,RefWhiteIllum=colorscience::XYZperfectreflectingdiffuser,CAT='Bradford',CATarray=colorscience::ChromaticAdaptation) 
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
CATmatrix<-CATarray[, , CAT, 'direct']
Rrgb<-RefWhiteRGB[which(RefWhiteRGB[["description"]]==RGBModel ),]
Rillum<-unlist(Rrgb[ 'whitepointilluminant' ])
Rgamma<-unlist(Rrgb['gamma'])
Rrgbwhitergb<-RefWhiteIllum[which(RefWhiteIllum[["Illuminant"]]==Rillum ),]
Rrgbwhitex<-unlist(Rrgbwhitergb[paste('X',observer,sep='')])
Rrgbwhitey<-unlist(Rrgbwhitergb[paste('Y',observer,sep='')])
Rrgbwhitez<-unlist(Rrgbwhitergb[paste('Z',observer,sep='')])
if (is.na(gamma)) gamma<-Rgamma
if (RGBModel=='sRGB') gamma<- -gamma
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
xr<-unlist(Rrgb['xRed'])
yr<-unlist(Rrgb['yRed'])
xg<-unlist(Rrgb['xGreen'])
yg<-unlist(Rrgb['yGreen'])
xb<-unlist(Rrgb['xBlue'])
yb<-unlist(Rrgb['yBlue'])
m<-matrix(c(xr/yr,xg/yg,xb/yb,1.0,1.0,1.0,(1.0-xr-yr)/yr,(1.0-xg-yg)/yg,(1.0-xb-yb)/yb),3,3,byrow=TRUE)
mi<-solve(m)
sr<- mi %*% matrix(c(Rrgbwhitex,Rrgbwhitey,Rrgbwhitez),3,1)
MtxXYZ2RGB<- t(matrix(sr,3,3,byrow=TRUE) * m)
MtxXYZ2RGB<- solve(MtxXYZ2RGB)
Asz<-CATmatrix %*% matrix(c(Rx,Ry,Rz),3,1,byrow=T)
Adz<-CATmatrix %*% matrix(c(Rrgbwhitex,Rrgbwhitey,Rrgbwhitez),3,1,byrow=T)
XYZ<-XYZmatrix %*% t(CATmatrix)
xyz1<- XYZ * c(Adz/Asz)
xyz2<- xyz1 %*% t(solve(CATmatrix))
xyz3<- xyz2 %*% MtxXYZ2RGB
XYZ<-apply(xyz3,1:2,function(x) Compand(x, gamma))
XYZ
}

RGB2XYZ<-function(RGBmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser,RGBModel='sRGB',RefWhiteRGB=colorscience::whitepointsRGB,gamma=NA,
RefWhiteIllum=colorscience::XYZperfectreflectingdiffuser,CAT='Bradford',CATarray=colorscience::ChromaticAdaptation)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
CATmatrix<-CATarray[, , CAT, 'direct']
Rrgb<-RefWhiteRGB[which(RefWhiteRGB[["description"]]==RGBModel ),]
Rillum<-unlist(Rrgb[ 'whitepointilluminant' ])
Rgamma<-unlist(Rrgb['gamma'])
Rrgbwhitergb<-RefWhiteIllum[which(RefWhiteIllum[["Illuminant"]]==Rillum ),]
Rrgbwhitex<-unlist(Rrgbwhitergb[paste('X',observer,sep='')])
Rrgbwhitey<-unlist(Rrgbwhitergb[paste('Y',observer,sep='')])
Rrgbwhitez<-unlist(Rrgbwhitergb[paste('Z',observer,sep='')])
if (is.na(gamma)) gamma<-Rgamma
if (RGBModel=='sRGB') gamma<- -gamma
RGBmatrix<-apply(RGBmatrix,1:2,function(x) InvCompand(x, gamma))
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
xr<-unlist(Rrgb['xRed'])
yr<-unlist(Rrgb['yRed'])
xg<-unlist(Rrgb['xGreen'])
yg<-unlist(Rrgb['yGreen'])
xb<-unlist(Rrgb['xBlue'])
yb<-unlist(Rrgb['yBlue'])
m<-matrix(c(xr/yr,xg/yg,xb/yb,1.0,1.0,1.0,(1.0-xr-yr)/yr,(1.0-xg-yg)/yg,(1.0-xb-yb)/yb),3,3,byrow=TRUE)
mi<-solve(m)
sr<- mi %*% matrix(c(Rrgbwhitex,Rrgbwhitey,Rrgbwhitez),3,1)
MtxRGB2XYZ<- t(matrix(sr,3,3,byrow=TRUE) * m)
XYZ<-RGBmatrix %*% MtxRGB2XYZ
Adz<-CATmatrix %*% matrix(c(Rx,Ry,Rz),3,1,byrow=T)
Asz<-CATmatrix %*% matrix(c(Rrgbwhitex,Rrgbwhitey,Rrgbwhitez),3,1,byrow=T)
xyz2<- t(apply(XYZ,1, function(x) rowSums( matrix(x,3,3,byrow=T) * CATmatrix)))
xyz2<- xyz2 *  matrix(c(Adz/Asz),dim(XYZ)[1],3,byrow=T)
XYZ<-t(apply(xyz2,1, function(x) c(x) %*% t(solve(CATmatrix))))
XYZ
}

PreucilAngle<-function(RGBmatrix) {
# The Reproduction of Colour  By R. W. G. Hunt
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
apply(RGBmatrix,1,function(x) {
R<-x[1]
G<-x[2]
B<-x[3]
if ((R >= G) & (G >= B)) return(60*(G-B)/(R-B))
if ((G > R) & (R >= B)) return(60*(2-(R-B)/(G-B)))
if ((G >= B) & (B > R)) return(60*(2+(B-R)/(G-R)))
if ((B > G) & (G > R)) return(60*(4-(G-R)/(B-R)))
if ((B > R) & (R >= G)) return(60*(4+(R-G)/(B-G)))
if ((R >= B) & (B > G)) return(60*(6-(B-G)/(R-G)))
})
}

PreucilPercentGreyness<-function(RGBmatrix) { # The Reproduction of Colour  By R. W. G. Hunt pp. 513
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
apply(RGBmatrix,1,function(x) {
M<-max(x)
L<-min(x)
H<-x[which(x != c(M,L))]
100*L/H
})
}

PreucilPercentHueError<-function(RGBmatrix) { # The Reproduction of Colour  By R. W. G. Hunt
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
apply(RGBmatrix,1,function(x) {
M<-max(x)
L<-min(x)
H<-x[which(x != c(M,L))]
100*(M-L)/(H-L)
})
}

RGB2hue<-function(RGBmatrix) {
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
atan2(sqrt(3)*(RGBmatrix[,2]-RGBmatrix[,3]), 2*RGBmatrix[,1]-RGBmatrix[,2]-RGBmatrix[,3])
}

CIE1976chroma<-function(CIELMatrix){
# CIELab and CIELuv input
if (is.null(dim(CIELMatrix))) if (length(CIELMatrix)>2) CIELMatrix<-matrix(CIELMatrix, ncol=3,byrow=TRUE)
sqrt(CIELMatrix[,2]^2+CIELMatrix[,3]^2)
}
CIE1976hueangle<-function(CIELMatrix){
#  CIELab and CIELuv input
if (is.null(dim(CIELMatrix))) if (length(CIELMatrix)>2) CIELMatrix<-matrix(CIELMatrix, ncol=3,byrow=TRUE)
#atan(CIELMatrix[,3]/CIELMatrix[,2])*180/pi
atan2(CIELMatrix[,3],CIELMatrix[,2])*180/pi
}

CIE1976uvSaturation<-function(uvMatrix, whitepoint){
# CIELuv input
if (is.null(dim(uvMatrix))) if (length(uvMatrix)>2) uvMatrix<-matrix(uvMatrix, ncol=3,byrow=TRUE)
13*sqrt((uvMatrix[,2]-whitepoint[1])^2 + (uvMatrix[,3]-whitepoint[2])^2)
}

WestlandBlacknessIndex<-function(CIELabMatrix){
# (Westland, et al., 2006) blackness index
if (is.null(dim(CIELabMatrix))) if (length(CIELabMatrix)>2) CIELabMatrix<- matrix(CIELabMatrix, ncol=3,byrow=TRUE)
9.6523-0.3118*CIELabMatrix[,1]-0.0054*CIELabMatrix[,2]^2-0.0003*CIELabMatrix[,3]^2
}

RGB2LSLM<-function(RGBmatrix) {#Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats Francoise Vienot, Hans Brettel,John D. Mollon
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
RGBmatrix<-RGBmatrix^2.2
cbind(0.209*(RGBmatrix[,1] - 0.5) + 0.715*(RGBmatrix[,2] - 0.5) + 0.076*(RGBmatrix[,3] - 0.5), 0.209*(RGBmatrix[,1] - 0.5) + 0.715*(RGBmatrix[,2] - 0.5) - 0.924*(RGBmatrix[,3] - 0.5),
3.14*(RGBmatrix[,1] - 0.5) - 2.799*(RGBmatrix[,2] - 0.5) - 0.349*(RGBmatrix[,3] - 0.5))
}

LSLM2RGB<-function(LSLMmatrix) {#Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats Francoise Vienot, Hans Brettel,John D. Mollon
if (is.null(dim(LSLMmatrix))) if (length(LSLMmatrix)>2) LSLMmatrix<-matrix(LSLMmatrix, ncol=3,byrow=TRUE)
cbind(LSLMmatrix[,1] - 0.013*LSLMmatrix[,2] + 0.252*LSLMmatrix[,3] + 0.5, LSLMmatrix[,1] + 0.110*LSLMmatrix[,2] - 0.074*LSLMmatrix[,3] + 0.5, LSLMmatrix[,1] - LSLMmatrix[,2] + 0.5)
}

XYZ2LMS<-function(XYZmatrix)   XYZmatrix %*% (matrix(c(0.15514,0.54312,-0.03286,-0.15514,0.45684,0.03286,0,0,0.01608),3,3,TRUE))
#Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats Francoise Vienot, Hans Brettel,John D. Mollon

LMS2XYZ<-function(LMSmatrix)   LMSmatrix %*% matrix(c(2.944813,-3.500978,13.17218,1.000040,1.000040,0.00000,0.000000,0.000000,62.18905),3,3,TRUE)
#Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats Francoise Vienot, Hans Brettel,John D. Mollon

RGB2LMS<-function(RGBmatrix)
{#Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats Francoise Vienot, Hans Brettel,John D. Mollon
if (any(RGBmatrix>1)) RGBmatrix<-RGBmatrix/255
RGBmatrix<-RGBmatrix^2.2
RGBmatrix %*% matrix(c(17.8824, 43.5161, 4.11935, 3.45565, 27.1554, 3.86714, 0.0299566, 0.184309, 1.46709),3,3,TRUE)
}

LMS2RGB<-function(LMSmatrix)
{#Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats Francoise Vienot, Hans Brettel,John D. Mollon
l<-LMSmatrix %*% matrix(c(0.0809444479,-0.130504409,0.116721066,-0.0102485335,0.0540193266,-0.113614708,-0.000365296938,-0.00412161469,0.693511405),3,3,TRUE)
matrix(as.integer(l^(1/2.2)*255),ncol=3,byrow=F)
}

spectra2lux<-function(spectraIn=NA, ciexyzIn=NA,wlIn=NA, wlInterval=NA)
{# convert spectral data to illuminance
#spectraIn = spectra from spectrophotometer
#ciexyzIn = photopic response
#wlIn = range of output wavelengths
#wlInterval = arbitrary wavelength interval to be applied to all series through interpolation
#require(Hmisc)
if (any(is.na(spectraIn))) stop('<<spectraIn>> must be a numeric array nx2')
if (is.null(dim(spectraIn))) stop('<<spectraIn>> must be a numeric array nx2')
if (dim(spectraIn)[2] != 2) stop('<<spectraIn>> must be a numeric array nx2')
if (!is.numeric(spectraIn)) stop('<<spectraIn>> must be a numeric array nx2')
if (any(is.na(wlIn))) wlIn<-c(min(spectraIn[,1]), max(spectraIn[,1]))
if (any(is.na(wlInterval))) wlInterval<-spectraIn[2,1]-spectraIn[1,1]
if (any(is.na(ciexyzIn))) ciexyzIn<-colorscience::ciexyz31
wlMin<- min(wlIn)
wlMax<- max(wlIn)
wlSeq<-seq(wlMin,wlMax,wlInterval)
if ((min(spectraIn[,1])<=wlMin) & (max(spectraIn[,1])>=wlMax)) spectraIn <- spectraIn[which(spectraIn[,1] >=wlMin & spectraIn[,1] <=wlMax ),] else {
# extrapolate 
f <- approxfun(spectraIn[,1], spectraIn[,2])
xtrapolated<-approxExtrap(spectraIn[,1], spectraIn[,2], wlSeq)$y
xtrapolated[which(xtrapolated<0)]<-0
#curve(f(x), wlMin,wlMax, col = "green2")
#points(spectraIn)
spectraIn<-cbind(wlSeq,xtrapolated)
#curve(f(x), wlMin,wlMax, col = "green2",ylim=c(0,40))
#points(spectraIn)
}
if ((min(ciexyzIn[,1])<=wlMin) & (max(ciexyzIn[,1])>=wlMax)) ciexyzIn <- ciexyzIn[which(ciexyzIn[,1] %in% wlSeq),] else {
# extrapolate 
f <- approxfun(ciexyzIn[,1], ciexyzIn[,2])
xtrapolated1<-approxExtrap(ciexyzIn[,1], ciexyzIn[,2], wlSeq)$y
xtrapolated1[which(xtrapolated<0)]<-0
f <- approxfun(ciexyzIn[,1], ciexyzIn[,3])
xtrapolated2<-approxExtrap(ciexyzIn[,1], ciexyzIn[,3], wlSeq)$y
xtrapolated2[which(xtrapolated<0)]<-0
f <- approxfun(ciexyzIn[,1], ciexyzIn[,4])
xtrapolated3<-approxExtrap(ciexyzIn[,1], ciexyzIn[,4], wlSeq)$y
xtrapolated3[which(xtrapolated<0)]<-0
ciexyzIn<-cbind(wlSeq,xtrapolated1,xtrapolated2,xtrapolated3)
}
683 * sum( spectraIn[,2] * ciexyzIn[,3])  # multiply spectral distribution curve * photopic response curve
}

spectra2XYZ<-function(spectraIn=NA, illuminantIn=NA, ciexyzIn=NA,wlIn=NA, wlInterval=NA)
{# convert spectral data to tristimulus values
#spectraIn = spectra from spectrophotometer
#illuminantIn = illuminant
#ciexyzIn = photopic response
#wlIn = range of output wavelengths
#wlInterval = arbitrary wavelength interval to be applied to all series through interpolation
#require(Hmisc)
if (any(is.na(spectraIn))) stop('<<spectraIn>> must be a numeric array nx2')
if (is.null(dim(spectraIn))) stop('<<spectraIn>> must be a numeric array nx2')
if (dim(spectraIn)[2] != 2) stop('<<spectraIn>> must be a numeric array nx2')
if (!is.numeric(spectraIn)) stop('<<spectraIn>> must be a numeric array nx2')
if (any(is.na(wlIn))) wlIn<-c(min(spectraIn[,1]), max(spectraIn[,1]))
if (any(is.na(wlInterval))) wlInterval<-spectraIn[2,1]-spectraIn[1,1]
if (any(is.na(illuminantIn))) illuminantIn<-colorscience::illuminantD65[which(colorscience::illuminantD65[,1] %in% seq(min(colorscience::illuminantD65[,1]), max(colorscience::illuminantD65[,1]), wlInterval)),]
if (any(is.na(ciexyzIn))) ciexyzIn<-colorscience::ciexyz31
wlMin<- min(wlIn)
wlMax<- max(wlIn)
wlSeq<-seq(wlMin,wlMax,wlInterval)
if ((min(illuminantIn[,1])<=wlMin) & (max(illuminantIn[,1])>=wlMax)) illuminantIn <- illuminantIn[which(illuminantIn[,1] %in% wlSeq),] else {
# extrapolate 
f <- approxfun(illuminantIn[,1], illuminantIn[,2])
xtrapolated<-approxExtrap(illuminantIn[,1], illuminantIn[,2], wlSeq)$y
xtrapolated[which(xtrapolated<0)]<-0
illuminantIn<-cbind(wlSeq,xtrapolated)
}
if ((min(spectraIn[,1])<=wlMin) & (max(spectraIn[,1])>=wlMax)) spectraIn <- spectraIn[which(spectraIn[,1] %in% wlSeq),] else {
# extrapolate 
f <- approxfun(spectraIn[,1], spectraIn[,2])
xtrapolated<-approxExtrap(spectraIn[,1], spectraIn[,2], wlSeq)$y
xtrapolated[which(xtrapolated<0)]<-0
#curve(f(x), wlMin,wlMax, col = "green2")
#points(spectraIn)
spectraIn<-cbind(wlSeq,xtrapolated)
#curve(f(x), wlMin,wlMax, col = "green2",ylim=c(0,40))
#points(spectraIn)
}
if ((min(ciexyzIn[,1])<=wlMin) & (max(ciexyzIn[,1])>=wlMax)) ciexyzIn <- ciexyzIn[which(ciexyzIn[,1] %in% wlSeq),] else {
# extrapolate 
f <- approxfun(ciexyzIn[,1], ciexyzIn[,2])
xtrapolated1<-approxExtrap(ciexyzIn[,1], ciexyzIn[,2], wlSeq)$y
xtrapolated1[which(xtrapolated<0)]<-0
f <- approxfun(ciexyzIn[,1], ciexyzIn[,3])
xtrapolated2<-approxExtrap(ciexyzIn[,1], ciexyzIn[,3], wlSeq)$y
xtrapolated2[which(xtrapolated<0)]<-0
f <- approxfun(ciexyzIn[,1], ciexyzIn[,4])
xtrapolated3<-approxExtrap(ciexyzIn[,1], ciexyzIn[,4], wlSeq)$y
xtrapolated3[which(xtrapolated<0)]<-0
ciexyzIn<-cbind(wlSeq,xtrapolated1,xtrapolated2,xtrapolated3)
}
r<-illuminantIn[,2] * spectraIn[,2] # multiply illuminant by spectra
xv<-sum(r * ciexyzIn[,2])/100.0 # sum the product of illuminant, spectra and CIE XYZ Color Matching Functions
yv<-sum(r * ciexyzIn[,3])/100.0
zv<-sum(r * ciexyzIn[,4])/100.0
k<-sum(illuminantIn[,2] * ciexyzIn[,3]) # divide by the photopic response
xv<-xv/k
yv<-yv/k
zv<-zv/k
c(xv,yv,zv)
}

createIsoTempLinesTable <- function(SPD=NA,CIETable = colorscience::ciexyz31, TCS = colorscience::TCSdata){
# generate data for isotemperature lines needed for calculating correlated color temperature
# Light source SPD
# reference data values CIETable
# The CIE 1931 2 degree Standard Colorimetric Observer
# wavelength(nm)   xbar    ybar zbar
# The spectral reflectance data of 14 color test samples for CRI TCS
# wavelength (nm) TCS1 TCS2 TCS3 ... TCS14

  wavelength_spd <- SPD[,1]
  spd <- SPD[,2]

        wavelength <- CIETable[,1]
        xbar <- CIETable[,2]
        ybar <- CIETable[,3]
        zbar <- CIETable[,4]

# Data for isotemperature lines needed for calculating correlated color temperature

# The following provides a table of isotemperature lines for use with the Robertson Method
# (Robertson, 1968) to interpolate isotemperature lines from the CIE 1960 UCS.
# The spacing of the isotemp lines is very small (1 1/MK) so very little
# interpolation is actually needed for determining CCT. The latest (2002)
# recommended values for the physical constants determining blackbody
# radiation spectra are used

dwave = wavelength[2]-wavelength[1] # wavelength increment = 1 nm

ubar = (2/3)*xbar
vbar = ybar
wbar = -0.5*xbar + (3/2)*ybar + 0.5*zbar

# 2002 CODATA recommended values
h = 6.6260693e-34
c2002 = 299792458
k = 1.3806505e-23
c1 = 2*pi*h*c2002^2
c2 = h*c2002/k
MrecpK <- c(0.01, 1:600) # mega reciprocal Kelvin values of isotemperature lines
TisotempLines <- 1/(MrecpK*1e-6)


u=v=sl=m=vector(mode='integer',0)
for (i in 1:length(TisotempLines)){
   spdref = c1 * (1e-9*wavelength)^-5 / (exp(c2/(TisotempLines[i]* 1e-9*wavelength)) - 1)
   spdref = spdref/max(spdref)
   wave = wavelength*1e-9
   
   # Equations from Wyszecki and Sitles, Color Science, 2nd ed. 1982, page
   # 226 and 227
   U = sum(spdref*ubar)
   V = sum(spdref*vbar)
   W = sum(spdref*wbar)
   R = U+V+W
   u[i] = U/R
   v[i] = V/R
   
   Uprime = c1*c2*(TisotempLines[i])^-2*sum(wave^-6*ubar*exp(c2/(wave*TisotempLines[i]))*(exp(c2/(wave*(TisotempLines[i])))-1)^-2)*dwave
   Vprime = sum(c1*c2*TisotempLines[i]^-2*wave^-6*vbar*exp(c2/(wave*TisotempLines[i]))*(exp(c2/(wave*(TisotempLines[i])))-1)^-2)*dwave
   Wprime = sum(c1*c2*TisotempLines[i]^-2*wave^-6*wbar*exp(c2/(wave*TisotempLines[i]))*(exp(c2/(wave*(TisotempLines[i])))-1)^-2)*dwave
   Rprime = Uprime+Vprime+Wprime
   
   sl[i] = (Vprime*R-V*Rprime)/(Uprime*R-U*Rprime)
   m[i] = -1/sl[i]
}
cbind(T=TisotempLines, u=u, v=v, m=m) # isoTempLinesTable
}


spectra2CCT <- function(SPD=NA, isoTempLinesTable=NA,CIETable = colorscience::ciexyz31, TCS = colorscience::TCSdata){
# Correlated Color Temperature CCT
if(any(is.na(isoTempLinesTable))) isoTempLinesTable=createIsoTempLinesTable(SPD)
m <- isoTempLinesTable[,"m"]
TisotempLines <- isoTempLinesTable[,"T"]

wavelength_spd <- SPD[,1]
spd <- SPD[,2]
wavelength <- CIETable[,1]
xbar <- CIETable[,2]
ybar <- CIETable[,3]
zbar <- CIETable[,4]

# Correlated Color Temperature (CCT)
# Interpolate CIE functions to spd increments

xbar <- approx(wavelength,xbar,wavelength_spd)$y
xbar[which(is.na(xbar))] <- 0.0
ybar <- approx(wavelength,ybar,wavelength_spd)$y
ybar[which(is.na(ybar))] <- 0.0
zbar <- approx(wavelength,zbar,wavelength_spd)$y
zbar[which(is.na(zbar))] <- 0.0
# Calculate Chromaticity Coordinates
X <- trapz(wavelength_spd,spd*xbar)
Y <- trapz(wavelength_spd,spd*ybar)
Z <- trapz(wavelength_spd,spd*zbar)
x <- X/(X+Y+Z)
y <- Y/(X+Y+Z)
u <- 4*x/(-2*x+12*y+3)
v <- 6*y/(-2*x+12*y+3)
cat('x <- ',x,'\ty <- ',y,'\n')
ut <- isoTempLinesTable[,"u"]
vt <- isoTempLinesTable[,"v"]
tt <- m #isoTempLinesTable$m

# Find adjacent lines to (us, vs) 
n <- length (TisotempLines)
index <- 0
d1 <- ((v-vt[1]) - tt[1]*(u-ut[1]))/sqrt(1+tt[1]*tt[1])
for (i in 2:n){
    d2 <- ((v-vt[i]) - tt[i]*(u-ut[i]))/sqrt(1+tt[i]*tt[i])
    if (d1/d2 < 0){
        index <- i
        break
   } else d1 <- d2
}
if (index == 0) {
    Tc <- -1 # Not able to calculate CCT, u, v coordinates outside range.
    stop('Not able to calculate CCT, u, v coordinates outside range.')
} else {
    # Calculate CCT by interpolation between isotemperature lines
    Tc <- 1/(1/TisotempLines[index-1]+d1/(d1-d2)*(1/TisotempLines[index]-1/TisotempLines[index-1]))
    return(Tc)
}
}

spectra2CRIGAIFSCI <- function(SPD=NA, isoTempLinesTable=NA, CCT=NA, CIETable = colorscience::ciexyz31, TCS = colorscience::TCSdata){
# CRI, GAI and FSCI
# Color Rendering Index CRI
# Gamut Area Index GAI
# full spectrum index FSCI
if(any(is.na(isoTempLinesTable))) isoTempLinesTable=createIsoTempLinesTable(SPD)
if(any(is.na(CCT))) CCT=spectra2CCT(SPD, isoTempLinesTable)

wavelength_spd <- SPD[,1]
spd <- SPD[,2]
wavelength <- CIETable[,1]
xbar <- CIETable[,2]
ybar <- CIETable[,3]
zbar <- CIETable[,4]
# Interpolate CIE functions to spd increments

xbar <- approx(wavelength,xbar,wavelength_spd)$y
xbar[which(is.na(xbar))] <- 0.0
ybar <- approx(wavelength,ybar,wavelength_spd)$y
ybar[which(is.na(ybar))] <- 0.0
zbar <- approx(wavelength,zbar,wavelength_spd)$y
zbar[which(is.na(zbar))] <- 0.0

# calculate the Color Rendering Indices (CRI and its 14 indices)
# Calculate Reference Source Spectrum, spdref. 
if (CCT < 5000){
    c1 = 3.7418e-16;
    c2 = 1.4388e-2;
    spdref = c1 * (1e-9*wavelength_spd)^-5 / (exp(c2/(CCT* 1e-9*wavelength_spd)) - 1)
} else {
    if (CCT <= 25000){
        #load('CIEDaySn','wavelength','S0','S1','S2');
        S0 <- colorscience::daylightcomponents[["S0"]]
        S1 <- colorscience::daylightcomponents[["S1"]]
        S2 <- colorscience::daylightcomponents[["S2"]]
        if (CCT <= 7000){
            xd = -4.6070e9 / CCT^3 + 2.9678e6 / CCT^2 + 0.09911e3 / CCT + 0.244063
        } else {
            xd = -2.0064e9 / CCT^3 + 1.9018e6 / CCT^2 + 0.24748e3 / CCT + 0.237040
        }
        yd = -3.000*xd*xd + 2.870*xd - 0.275
        M1 = (-1.3515 - 1.7703*xd + 5.9114*yd) / (0.0241 + 0.2562*xd - 0.7341*yd)
        M2 = (0.0300 - 31.4424*xd + 30.0717*yd) / (0.0241 + 0.2562*xd - 0.7341*yd)
        spdref = S0 + M1*S1 + M2*S2
        spdref = approx(wavelength,spdref,wavelength_spd)
        spdref[which(is.na(spdref))] = 0.0
    } else {
        return(NA)
        }
    }

# Interpolate TCS values from 5 nm to spd nm increments
TCS_1 = matrix(0, length(wavelength_spd),14)
for (i in 1:14) TCS_1[,i] = approx(TCS[,1],TCS[,i+1],wavelength_spd,'linear',0)$y
TCS_1[which(is.na(TCS_1))] = 0

# Calculate u, v chromaticity coordinates of samples under test illuminant, uk, vk and
# reference illuminant, ur, vr.
Yki = vector(mode='integer',length=14)
Yri = vector(mode='integer',length=14)
uki = vector(mode='integer',length=14)
vki = vector(mode='integer',length=14)
uri = vector(mode='integer',length=14)
vri = vector(mode='integer',length=14)
X = trapz(wavelength_spd,spd * xbar)
Y = trapz(wavelength_spd,spd * ybar)
Z = trapz(wavelength_spd,spd * zbar)
Yknormal = 100 / Y
Yk = Y*Yknormal
uk = 4*X/(X+15*Y+3*Z)
vk = 6*Y/(X+15*Y+3*Z)
X = trapz(wavelength_spd,spdref * xbar)
Y = trapz(wavelength_spd,spdref * ybar)
Z = trapz(wavelength_spd,spdref * zbar)
Yrnormal = 100 / Y
Yr = Y*Yrnormal
ur = 4*X/(X+15*Y+3*Z)
vr = 6*Y/(X+15*Y+3*Z)
for (i in 1:14){
	X = trapz(wavelength_spd,spd * TCS_1[,i] * xbar)
	Y = trapz(wavelength_spd,spd * TCS_1[,i] * ybar)
	Z = trapz(wavelength_spd,spd * TCS_1[,i] * zbar)
	Yki[i] = Y*Yknormal
	uki[i] = 4*X/(X+15*Y+3*Z)
	vki[i] = 6*Y/(X+15*Y+3*Z)
	X = trapz(wavelength_spd,spdref * TCS_1[,i] * xbar)
	Y = trapz(wavelength_spd,spdref * TCS_1[,i] * ybar)
	Z = trapz(wavelength_spd,spdref * TCS_1[,i] * zbar)
	Yri[i] = Y*Yrnormal
	uri[i] = 4*X/(X+15*Y+3*Z)
	vri[i] = 6*Y/(X+15*Y+3*Z)
}

# Check tolerance for reference illuminant
DC = sqrt((uk-ur)^2 + (vk-vr)^2)

# Apply adaptive (perceived) color shift.
ck = (4 - uk - 10*vk) / vk
dk = (1.708*vk + 0.404 - 1.481*uk) / vk
cr = (4 - ur - 10*vr) / vr
dr = (1.708*vr + 0.404 - 1.481*ur) / vr

ukip = vector(mode='integer',length=14)
vkip  = vector(mode='integer',length=14)
for (i in 1:14){
	cki = (4 - uki[i] - 10*vki[i]) / vki[i]
	dki = (1.708*vki[i] + 0.404 - 1.481*uki[i]) / vki[i]
	ukip[i] = (10.872 + 0.404*cr/ck*cki - 4*dr/dk*dki) / (16.518 + 1.481*cr/ck*cki - dr/dk*dki);
	vkip[i] = 5.520 / (16.518 + 1.481*cr/ck*cki - dr/dk*dki);
}

#  Transformation into 1964 Uniform space coordinates.
Wstarr = vector(mode='integer',length=14)
Ustarr = vector(mode='integer',length=14)
Vstarr = vector(mode='integer',length=14)
Wstark = vector(mode='integer',length=14)
Ustark = vector(mode='integer',length=14)
Vstark = vector(mode='integer',length=14)
for (i in 1:14){
	Wstarr[i] = 25*Yri[i]^.333333 - 17
	Ustarr[i] = 13*Wstarr[i]*(uri[i] - ur)
	Vstarr[i] = 13*Wstarr[i]*(vri[i] - vr)
	
	Wstark[i] = 25*Yki[i]^.333333 - 17
	Ustark[i] = 13*Wstark[i]*(ukip[i] - ur) # after applying the adaptive color shift, u'k = ur
	Vstark[i] = 13*Wstark[i]*(vkip[i] - vr) # after applying the adaptive color shift, v'k = vr
}

# Determination of resultant color shift, delta E.
deltaE = vector(mode='integer',length=14)
R = vector(mode='integer',length=14)
for (i in 1:14){
	deltaE[i] = sqrt((Ustarr[i] - Ustark[i])^2 + (Vstarr[i] - Vstark[i])^2 + (Wstarr[i] - Wstark[i])^2)
	R[i] = 100 - 4.6*deltaE[i]
}
Ra = sum(R[1:8])/8

# fourth, calculate the gamut area formed by the 8 CIE standard color samples
ukii=c(uki[1:8],uki[1])
vkii=1.5*c(vki[1:8],vki[1])
Ga=pracma::polyarea(ukii,vkii)
# Normalize gamut area to equal energy source 
Ga=Ga/0.00728468*100
#fprintf(1,'Gamut Area Index = %.1f\n',Ga)

# Fifth, calculate the FSI (full spectrum index)
# Calculates the Full-spectrum Index

# Interpolate to wavelength interval of 1nm from 380nm to 730nm
numWave = 351
t9=380:730
wavelength_spd <- SPD[,1]
spd <- SPD[,2]
spd=interp1(wavelength_spd,spd,t9,method='spline')
spd[which(is.na(spd))] = 0.0
spd = spd/sum(spd) # Normalize the relative spd so that the total power equals 1
#Equal energy cumulative spd
EEcum=seq(1/numWave,1,by=1/numWave)
#Calculate FSI

sumSqrDiff <- vector(mode='integer',length=numWave)
for (j in 1:numWave ){
cum = cumsum(spd) # A MatLab function for cumulative sums 
sqrDiff = (cum-EEcum)^2
sumSqrDiff[j]=sum(sqrDiff)
spd=circshift(spd,1)
}
FSI=mean(sumSqrDiff)
FSCI=100-5.1*FSI

return(cbind(CRIra = Ra, GamutAreaIndex = Ga, FSCI = FSCI))
}

CIELabtoDIN99<-function(Lab){
# Conversion from CIELAB color space to DIN99 coordinates
# input: c(L,a,b)
# http://de.wikipedia.org/w/index.php?title=Diskussion:DIN99-Farbraum
L <- Lab[1]
a <- Lab[2]
b <- Lab[3]
kE <- 1.0 #brightness factor, 1.0 for CIE reference conditions
kCH <- 1.0 #chroma and hue factor, 1.0 for CIE reference conditions
ang <- 2*pi/360*26 #26 degree angle in radians
L99f <- 100/log(139/100.) #L99 scaling factor - method 1
#L99f <- 303.67 #original L99 factor
L99o <- L99f/kE*log(1+.0039*L) #Lightness correction kE
if ((a==0.0) && (b==0.0)) {
a99o <- 0.0
b99o <- 0.0
} else {
eo <- a*cos(ang) + b*sin(ang)       #a rotation
     fo <- 0.83*(b*cos(ang) - a*sin(ang)) #b rotation/stretching
     Go <- sqrt(eo^2 + fo^2)         #chroma
     C99o <- log(1.0 + 0.075 * Go) / (0.0435 * kCH*kE)  #factor for chroma compression and viewing conditions
    heofo <- atan2(fo,eo) #arctan in four quadrants
     h99o <- heofo+ang  #hue rotation
    a99o <- C99o * cos (h99o)
    b99o <- C99o * sin (h99o)
}
cbind(L99o,a99o,b99o)
}

DIN99toCIELab<-function(Lab99o){
# Conversion from DIN99 coordinates to CIELAB color space
# input: c(L99o,a99o,b99o)
# http://de.wikipedia.org/w/index.php?title=Diskussion:DIN99-Farbraum
L99o <- Lab99o[1]
a99o <- Lab99o[2]
b99o <- Lab99o[3]
kE <- 1.0 #brightness factor, 1.0 for CIE reference conditions
kCH <- 1.0 #chroma and hue factor, 1.0 for CIE reference conditions
ang <- 2*pi/360*26 #26 degrees in radians
L99f <- 100/log(139/100.0) #corrected L99 factor
#L99f <- 303.67 #original L99 factor
L <- (exp(L99o*kE/L99f)-1.0)/0.0039
h99ef <- atan2(b99o,a99o) #arctan in four quadrants
heofo <- h99ef-ang #backwards hue rotation
C99 <- sqrt(a99o^2+b99o^2)     #DIN99 chroma
G <- (exp(0.0435 * kE * kCH * C99)-1.0)/0.075 #factor for chroma decompression and viewing conditions
e <- G*cos(heofo)
f <- G*sin(heofo)
a <- (e * cos(ang) - (f/0.83) * sin(ang)) #rotation by 26 degrees
b <- (e * sin(ang) + (f/0.83) * cos(ang)) #rotation by 26 degrees
cbind(L=L,a=a,b=b)
}

Y2MunsellVtable1D1535<-function(Y){
# CIE XYZ "Y" to Munsell "V"
# NLSQ regression for obtaining similar results to table 1 from ASTM Standard D1535-08
Y*100-9.876522e-06 *(Y-1.749992e+02)^2 + 6.139636e-03*(Y-1.749992e+02) + 1.370407e+00
}

MunsellV2Y<-function(V) 1.1914*V-0.22533*V ^ 2 +0.23352*V ^ 3-0.020484*V ^ 4 +0.00081939*V ^ 5
# Munsell "V" to CIE XYZ "Y" 
# ASTM Standard D1535-08

Y2MunsellV<-function(Y){
# CIE XYZ "Y" to Munsell "V" 
# ASTM Standard D1535-08
a = 2.49268
b = 1.5614
C = 0.985
d = 0.1073
e = 3.084
f = 7.54
g = 0.0133
H = 2.3
j = 0.0084
k = 4.1
m = 0.0221
n = 0.39
P = 0.037
Q = 0.44
S = 1.28
T = 0.53
U = 0.87455
W = 0.9967
if (Y <= 0.9) V <- U * Y ^ W else V <- a * (Y ^ (1 / 3)) - b - C / (((d * Y - e) ^ 2) + f) + g / (Y ^ H)  + j *
sin((k * (Y ^ (1 / 3)) + 1)/ 180 * pi) + m / Y * sin((n * (Y - 2))/ 180 * pi) + P / Q / Y * sin((S * (Y - T))/ 180 * pi)
V
}

XYZ2Lab <- function(XYZmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (!is.matrix(XYZmatrix)) XYZmatrix<-matrix(XYZmatrix,ncol=3,byrow=TRUE)
kE <- 216.0 / 24389.0
kK <- 24389.0 / 27.0
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
xr <- XYZmatrix[,1] / Rx #* 100
yr <- XYZmatrix[,2] / Ry #* 100
zr <- XYZmatrix[,3] / Rz #* 100
fx <- ifelse((xr > kE) , (xr^(1.0 / 3.0)) , ((kK * xr + 16.0) / 116.0))
fy <- ifelse((yr > kE) , (yr^(1.0 / 3.0)) , ((kK * yr + 16.0) / 116.0))
fz <- ifelse((zr > kE) , (zr^(1.0 / 3.0)) , ((kK * zr + 16.0) / 116.0))
L <- 116.0 * fy - 16.0
a <- 500.0 * (fx - fy)
b <- 200.0 * (fy - fz)
cbind(L=L,a=a,b=b)
}

XYZ2Luv <- function(XYZmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (!is.matrix(XYZmatrix)) XYZmatrix<-matrix(XYZmatrix,ncol=3,byrow=TRUE)
kE <- 216.0 / 24389.0
kK <- 24389.0 / 27.0
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])

U = ( 4 * XYZmatrix[,1] ) / ( XYZmatrix[,1] + ( 15 * XYZmatrix[,2] ) + ( 3 * XYZmatrix[,3] ) )
V = ( 9 * XYZmatrix[,2] ) / ( XYZmatrix[,1] + ( 15 * XYZmatrix[,2] ) + ( 3 * XYZmatrix[,3] ) )
vY = XYZmatrix[,2] #/ 100
vY = ifelse (( vY > kE ), vY ^ ( 1/3 ), ( kK * vY  + 16) / 116 )
Ru = ( 4 * Rx ) / ( Rx + ( 15 * Ry ) + ( 3 * Rz ) )
Rv = ( 9 * Ry ) / ( Rx + ( 15 * Ry ) + ( 3 * Rz ) )
L = ( 116 * vY ) - 16
u = 13 * L* ( U - Ru )
v = 13 * L * ( V - Rv )

cbind(L=L,u=u,v=v)
}

xyY2XYZ<-function(xyYmatrix)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (!is.matrix(xyYmatrix)) xyYmatrix<-matrix(xyYmatrix,ncol=3,byrow=TRUE)
if (dim(xyYmatrix)[1]>1) {
z<-apply(xyYmatrix,1,function(r) if (r[1] < 0.000001) return(c(X=0,Y=0,Z=0)) else return(c(X=(r[1] * r[3]) / r[2],Y=r[3],Z= ((1.0 - r[1] - r[2]) * r[3]) / r[2])) )
XYZ<-t(z)
} else if (xyYmatrix[3] < 0.000001) XYZ<-c(X=0,Y=0,Z=0) else XYZ<-cbind(X=(xyYmatrix[1] * xyYmatrix[3]) / xyYmatrix[2],Y=xyYmatrix[3],Z= ((1.0 - xyYmatrix[1] - xyYmatrix[2]) * xyYmatrix[3]) / xyYmatrix[2])
colnames(XYZ)<-c('X','Y','Z')
XYZ
}

XYZ2xyY<-function(XYZmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (!is.matrix(XYZmatrix)) XYZmatrix<-matrix(XYZmatrix,ncol=3,byrow=TRUE)
Den <- rowSums(XYZmatrix)
DenG0<-which(Den > 0.0)
xyYmatrix<-XYZmatrix
xyYmatrix[DenG0,1:2]<- XYZmatrix[DenG0,1:2]/Den
xyYmatrix[DenG0,3]<- XYZmatrix[DenG0,2]
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
x <- Rx / (Rx + Ry + Rz)
y <- Ry / (Rx + Ry + Rz)
xyYmatrix[-DenG0,1]<-x
xyYmatrix[-DenG0,2]<-y
xyYmatrix
}

Lab2XYZ<-function(Labmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (!is.matrix(Labmatrix)) Labmatrix<-matrix(Labmatrix,ncol=3,byrow=TRUE)
L<-Labmatrix[,1]
a<-Labmatrix[,2]
b<-Labmatrix[,3]
kE <- 216.0 / 24389.0
kK <- 24389.0 / 27.0
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
fy <- (L + 16.0) / 116.0
fx <- 0.002 * a + fy
fz <- fy - 0.005 * b
fx3 <- fx * fx * fx
fz3 <- fz * fz * fz
xr <- ifelse((fx3 > kE) , fx3 , ((116.0 * fx - 16.0) / kK))
yr <- ifelse((L > 8.0) , (((L + 16.0) / 116.0)^3.0) , (L / kK))
zr <- ifelse((fz3 > kE) , fz3 , ((116.0 * fz - 16.0) / kK))
X <- xr * Rx
Y <- yr * Ry
Z <- zr * Rz
cbind(X=X,Y=Y,Z=Z)
}

Luv2XYZ<-function(Luvmatrix,illuminant='D65',observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
if (!is.matrix(Luvmatrix)) Luvmatrix<-matrix(Luvmatrix,ncol=3,byrow=TRUE)
L<-Luvmatrix[,1]
u<-Luvmatrix[,2]
v<-Luvmatrix[,3]
kE <- 216.0 / 24389.0
kK <- 24389.0 / 27.0
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
Y <- ifelse((L > 8.0) , (((L + 16.0) / 116.0)^3.0) , (L / kK))
u0 <- (4.0 * Rx) / (Rx + 15.0 * Ry + 3.0 * Rz)
v0 <- (9.0 * Ry) / (Rx + 15.0 * Ry + 3.0 * Rz)
a <- (((52.0 * L) / (u + 13.0 * L * u0)) - 1.0) / 3.0
b <- -5.0 * Y
c <- -1.0 / 3.0
d <- Y * (((39.0 * L) / (v + 13.0 * L * v0)) - 5.0)
X <- (d - b) / (a - c)
Z <- X * a + b
cbind(X=X,Y=Y,Z=Z)
}

Compand<-function(linearV, gammaV)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html

if (gammaV > 0.0) companded <- ifelse((linearV >= 0.0) , (linearV^(1.0 / gammaV)) , -(-linearV^(1.0 / gammaV))) else if (gammaV < 0.0)
{
# sRGB
signV <- 1.0
if (linearV < 0.0)
{
signV <- -1.0
linearV <- -linearV
}
companded <- ifelse((linearV <= 0.0031308), (linearV * 12.92) , (1.055 * (linearV^(1.0 / 2.4)) - 0.055)) 
companded <- companded *signV
}
else
{
# L*
signV <- 1.0
if (linearV < 0.0)
{
signV <- -1.0
linearV <- -linearV
}
companded <- ifelse((linearV <= (216.0 / 24389.0)), (linearV * 24389.0 / 2700.0) , (1.16 * (linearV^(1.0 / 3.0)) - 0.16)) 
companded <- companded *signV
}
return(companded)
}

InvCompand<-function(companded, gammaV)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html

if (gammaV > 0.0) linearV <- ifelse((companded >= 0.0), (companded^gammaV) , -(-companded^gammaV)) else if (gammaV < 0.0) 
{
# sRGB
signV <- 1.0
if (companded < 0.0)
{
signV <- -1.0
companded <- -companded
}
linearV <- ifelse((companded <= 0.04045), (companded / 12.92) , (((companded + 0.055) / 1.055)^2.4)) 
linearV <- linearV *signV
}
else
{
# L*
signV <- 1.0
if (companded < 0.0)
{
signV <- -1.0
companded <- -companded
}
linearV <- ifelse((companded <= 0.08), (2700.0 * companded / 24389.0) , ((((1000000.0 * companded + 480000.0) * companded + 76800.0) * companded + 4096.0) / 1560896.0)) 
linearV <- linearV *signV
}
return(linearV)
}

compuphaseDifferenceRGB<-function(RGB1, RGB2)
{ # difference metric for RGB
# http://www.compuphase.com/cmetric.htm
  rmean <- mean( RGB1[1] + RGB2[1] )
  r <- RGB1[1] - RGB2[1]
  g <- RGB1[2] - RGB2[2]
  b <- RGB1[3] - RGB2[3]
  sqrt((((512+rmean)*r^2)/256) + 4*g^2 + (((767-rmean)*b^2)/256))
}

deltaE1976 <- function(Lab1, Lab2)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html

delL <- Lab1[1] - Lab2[1]
dela <- Lab1[2] - Lab2[2]
delb <- Lab1[3] - Lab2[3]
sqrt(delL * delL + dela * dela + delb * delb)
}

deltaE1994 <- function (Lab1, Lab2, textiles=FALSE)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
k1 <- ifelse((textiles == TRUE),0.048, 0.045)
k2 <- ifelse((textiles == TRUE),0.014, 0.015)
kL <- ifelse((textiles == TRUE),2.0, 1.0)
kC <- 1.0
kH <- 1.0
C1 <- sqrt(Lab1[2] * Lab1[2] + Lab1[3] * Lab1[3])
C2 <- sqrt(Lab2[2] * Lab2[2] + Lab2[3] * Lab2[3])
delA <- Lab1[2] - Lab2[2]
delB <- Lab1[3] - Lab2[3]
delC <- C1 - C2
delH2 <- delA * delA + delB * delB - delC * delC
delH <- ifelse((delH2 > 0.0),sqrt(delH2), 0.0) 
delL <- Lab1[1] - Lab2[1]
sL <- 1.0
sC <- 1.0 + k1 * C1
sH <- 1.0 + k2 * C1
vL <- delL / (kL * sL)
vC <- delC / (kC * sC)
vH <- delH / (kH * sH)
if (textiles == TRUE) return(sqrt(vL * vL + vC * vC + vH * vH)) else return(sqrt(vL * vL + vC * vC + vH * vH))
}

deltaE2000 <- function(Lab1, Lab2)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
kL <- 1.0
kC <- 1.0
kH <- 1.0
lBarPrime <- 0.5 * (Lab1[1] + Lab2[1])
c1 <- sqrt(Lab1[2] * Lab1[2] + Lab1[3] * Lab1[3])
c2 <- sqrt(Lab2[2] * Lab2[2] + Lab2[3] * Lab2[3])
cBar <- 0.5 * (c1 + c2)
cBar7 <- cBar^7
g <- 0.5 * (1.0 - sqrt(cBar7 / (cBar7 + 6103515625.0)))# 6103515625 = 25^7
a1Prime <- Lab1[2] * (1.0 + g)
a2Prime <- Lab2[2] * (1.0 + g)
c1Prime <- sqrt(a1Prime * a1Prime + Lab1[3] * Lab1[3])
c2Prime <- sqrt(a2Prime * a2Prime + Lab2[3] * Lab2[3])
cBarPrime <- 0.5 * (c1Prime + c2Prime)
h1Prime <- (atan2(Lab1[3], a1Prime) * 180.0) / pi
if (h1Prime < 0.0)
h1Prime <- h1Prime +360.0
h2Prime <- (atan2(Lab2[3], a2Prime) * 180.0) / pi
if (h2Prime < 0.0)
h2Prime <- h2Prime +360.0
hBarPrime <- ifelse((abs(h1Prime - h2Prime) > 180.0),(0.5 * (h1Prime + h2Prime + 360.0)), (0.5 * (h1Prime + h2Prime))) 
t <- 1.0 -
0.17 * cos(pi * (      hBarPrime - 30.0) / 180.0) +
0.24 * cos(pi * (2.0 * hBarPrime       ) / 180.0) +
0.32 * cos(pi * (3.0 * hBarPrime +  6.0) / 180.0) -
0.20 * cos(pi * (4.0 * hBarPrime - 63.0) / 180.0)
if (abs(h2Prime - h1Prime) <= 180.0) dhPrime <- h2Prime - h1Prime else dhPrime <- ifelse((h2Prime <= h1Prime),(h2Prime - h1Prime + 360.0), (h2Prime - h1Prime - 360.0)) 
dLPrime <- Lab2[1] - Lab1[1]
dCPrime <- c2Prime - c1Prime
dHPrime <- 2.0 * sqrt(c1Prime * c2Prime) * sin(pi * (0.5 * dhPrime) / 180.0)
sL <- 1.0 + ((0.015 * (lBarPrime - 50.0) * (lBarPrime - 50.0)) / sqrt(20.0 + (lBarPrime - 50.0) * (lBarPrime - 50.0)))
sC <- 1.0 + 0.045 * cBarPrime
sH <- 1.0 + 0.015 * cBarPrime * t
dTheta <- 30.0 * exp(-((hBarPrime - 275.0) / 25.0) * ((hBarPrime - 275.0) / 25.0))
cBarPrime7 <- cBarPrime * cBarPrime * cBarPrime * cBarPrime * cBarPrime * cBarPrime * cBarPrime
rC <- sqrt(cBarPrime7 / (cBarPrime7 + 6103515625.0))
rT <- -2.0 * rC * sin(pi * (2.0 * dTheta) / 180.0)
sqrt(
(dLPrime / (kL * sL)) * (dLPrime / (kL * sL)) +
(dCPrime / (kC * sC)) * (dCPrime / (kC * sC)) +
(dHPrime / (kH * sH)) * (dHPrime / (kH * sH)) +
(dCPrime / (kC * sC)) * (dHPrime / (kH * sH)) * rT)
}

deltaECMC <- function(Lab1, Lab2, L=1, C=1)
{# Based on: Bruce Justin Lindbloom 2013 http:#www.brucelindbloom.com/index.html?ColorCalculator.html
c1 <- sqrt(Lab1[2] * Lab1[2] + Lab1[3] * Lab1[3])
c2 <- sqrt(Lab2[2] * Lab2[2] + Lab2[3] * Lab2[3])
sl <- ifelse((Lab1[1] < 16.0),(0.511) ,((0.040975 * Lab1[1]) / (1.0 + 0.01765 * Lab1[1]))) 
sc <- (0.0638 * c1) / (1.0 + 0.0131 * c1) + 0.638
h1 <- ifelse((c1 < 0.000001),0.0 ,((atan2(Lab1[3], Lab1[2]) * 180.0) / pi)) 
while (h1 < 0.0)
h1 <- h1 +360.0
while (h1 >= 360.0)
h1 <- h1 -360.0
t <- ifelse((h1 >= 164.0) && (h1 <= 345.0),(0.56 + abs(0.2 * cos((pi * (h1 + 168.0)) / 180.0))), (0.36 + abs(0.4 * cos((pi * (h1 + 35.0)) / 180.0)))) 
c4 <- c1^4
f <- sqrt(c4 / (c4 + 1900.0))
sh <- sc * (f * t + 1.0 - f)
delL <- Lab1[1] - Lab2[1]
delC <- c1 - c2
delA <- Lab1[2] - Lab2[2]
delB <- Lab1[3] - Lab2[3]
dH2 <- delA * delA + delB * delB - delC * delC
v1 <- delL / (L * sl)
v2 <- delC / (C * sc)
v3 <- sh
sqrt(v1 * v1 + v2 * v2 + (dH2 / (v3 * v3)))
}

MunsellV2relativeLuminanceY<-function(V) 1.2219*V - 0.23111*V^2 + 0.23951*V^3 - 0.021009*V^4 + 0.0008404*V^5
# convert Munsell value V to relative luminance Y
# Color Appearance Models, Mark D. Fairchild pp. 96

CIEluminanceY2NCSblackness<-function(Y) 100 - 156*Y / (56 + Y )
# approximated NCS blackness s by the CIE luminance factor Y
# Introduction to Color Imaging Science, Lee pp. 366

huedegree<-function(MunIn){
# convert Munsell hue to degree
# based on Takahiro Onodera 2010 Color-Model-Munsell-Util http://annocpan.org/dist/Color-Model-Munsell-Util
# 10.0RP will return 0
# sapply(MunIn[,'H']
m<-sapply(MunIn , function(x) {if (any(is.na(x))) return(NA) else {
sidN<-gregexpr('^(\\d{1,2}\\.\\d{1,2}|\\d{1,2})(RP|YR|Y|GY|G|BG|B|PB|P|R)', x, perl=TRUE)
h1<-substr(x,attr(sidN[[1]], "capture.start")[1,][1],attr(sidN[[1]], "capture.start")[1,][1]+attr(sidN[[1]], "capture.length")[1,][1]-1)
hM<-substr(x,attr(sidN[[1]], "capture.start")[1,][2],attr(sidN[[1]], "capture.start")[1,][2]+attr(sidN[[1]], "capture.length")[1,][2]-1)
h2<-which(hM == colorscience::MunsellHues)
tmp <- (h2-1)*10+as.numeric(h1)
if (tmp>=100) tmp<-tmp %% 100
return(as.numeric(tmp))
} })
m
}

Adjust<-function(myColor, Factor, Gamma,IntensityMax=1.0)# HIDDEN
{# part of heuristic.wlnm2RGB
if (length(myColor)==1) {
if (myColor == 0.0) return(0) else return(round(IntensityMax * (myColor * Factor)^ Gamma))
} else {
w<-which(myColor != 0)
myColor[w] <- round(IntensityMax * (myColor[w] * Factor[w])^ Gamma)
myColor}
}

#Adjust<-function(Color, factorL, Gamma, IntensityMax=1.0){
#if (Color == 0.0) return (0) else return (round(IntensityMax * (Color * factorL)^Gamma))
#}

heuristic.wlnm2RGB<-function(wavelength,Gamma=0.80, IntensityMax=1.0){
#heuristic approximation of RGB values from wavelengths
# original work by Dan Bruton's (www.physics.sfasu.edu/astro/color.html)
# Delphi translation by Earl F. Glynn 2006 (http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm)
# R translation by Jose Gama 2013
if (IntensityMax<0) IntensityMax <- 0
if (IntensityMax>1.0) IntensityMax <- 1.0
wavelength<-round(wavelength)
Red <- Green <- Blue <- factorL <- rep(0,length(wavelength))
w<-which(wavelength %in% 380:439)
if (length(w)>0) {Red[w] <- -(wavelength[w] - 440) / (440 - 380);Green[w] <- 0.0;Blue[w] <- 1.0}
w<-which(wavelength %in% 440:489)
if (length(w)>0) {Red[w] <- 0.0;Green[w] <- (wavelength[w] - 440) / (490 - 440);Blue[w] <- 1.0}
w<-which(wavelength %in% 490:509)
if (length(w)>0) {Red[w] <- 0.0;Green[w] <- 1.0;Blue[w] <- -(wavelength[w] - 510) / (510 - 490)}
w<-which(wavelength %in% 510:579)
if (length(w)>0) {Red[w] <- (wavelength[w] - 510) / (580 - 510);Green[w] <- 1.0;Blue[w] <- 0.0}
w<-which(wavelength %in% 580:644)
if (length(w)>0) {Red[w] <- 1.0;Green[w] <- -(wavelength[w] - 645) / (645 - 580);Blue[w] <- 0.0}
w<-which(wavelength %in% 645:780)
if (length(w)>0) {Red[w] <- 1.0;Green[w] <- 0.0;Blue[w] <- 0.0}
#Let the intensity fall off near the vision limits
w<-which(wavelength %in% 380:419)
if (length(w)>0) factorL[w] <- 0.3 + 0.7*(wavelength[w] - 380) / (420 - 380)
w<-which(wavelength %in% 420:700)
if (length(w)>0) factorL[w] <- 1.0
w<-which(wavelength %in% 701:780)
if (length(w)>0) factorL[w] <- 0.3 + 0.7*(780 - wavelength[w]) / (780 - 700)
R <- Adjust(Red, factorL, Gamma, IntensityMax)
G <- Adjust(Green, factorL, Gamma, IntensityMax)
B <- Adjust(Blue, factorL, Gamma, IntensityMax)
list(R=R,G=G,B=B)
}

xFit_1931<-function(  wave )
{#Chris Wyman Peter-Pike Sloan Peter Shirley Journal of Computer Graphics Techniques Vol. 2, No. 2, 2013
#Simple Analytic Approximations to the CIE XYZ Color Matching Functions
t1 = (wave-442.0)*ifelse((wave<442.0),0.0624,0.0374)
t2 = (wave-599.8)*ifelse((wave<599.8),0.0264,0.0323)
t3 = (wave-501.1)*ifelse((wave<501.1),0.0490,0.0382)
(0.362*exp(-0.5*t1*t1) + 1.056*exp(-0.5*t2*t2) - 0.065*exp(-0.5*t3*t3))
}
yFit_1931<-function( wave )
{#Chris Wyman Peter-Pike Sloan Peter Shirley Journal of Computer Graphics Techniques Vol. 2, No. 2, 2013
#Simple Analytic Approximations to the CIE XYZ Color Matching Functions
t1 = (wave-568.8)*ifelse((wave<568.8),0.0213,0.0247)
t2 = (wave-530.9)*ifelse((wave<530.9),0.0613,0.0322)
0.821*exp(-0.5*t1*t1) + 0.286*exp(-0.5*t2*t2)
}
zFit_1931<-function( wave )
{#Chris Wyman Peter-Pike Sloan Peter Shirley Journal of Computer Graphics Techniques Vol. 2, No. 2, 2013
#Simple Analytic Approximations to the CIE XYZ Color Matching Functions
t1 = (wave-437.0)*ifelse((wave<437.0),0.0845,0.0278)
t2 = (wave-459.0)*ifelse((wave<459.0),0.0385,0.0725)
1.217*exp(-0.5*t1*t1) + 0.681*exp(-0.5*t2*t2)
}

StearnsStearnscorrection<-function(P){
# Stearns and Stearns correction
# based on Computational Colour Science using MATLAB Stephen Westland and Caterina Ripamonti John Wiley & Sons Ltd  2004 pp.35
# correction for spectral bandpass for TRUE reflectance data P (numeric vector)
a <- 0.083
l <- length(P)
cp<-rep(0,l)
for (n in 2:(l-1)) cP[n] = -a*P[n-1] + (1 + 2*a)*P[n] - a*P[n+1]
cP[1] = (1 + a)*P[1] - a*P[2]
cP[l] = (1 + a)*P[l] - a*P[l-1]
cP
}

xy2CCT.McCamy<-function(x,y) {
# convert from chromaticity coordinates to correlated color temperature (approximation)
# C. S. McCamy "Correlated color temperature as an explicit function of chromaticity coordinates" Color Research & Application Volume 17, Issue 2, pages 142-144, April 1992
xe = 0.3320
ye = 0.1858
n = (x - xe)/(y - ye)
-449*n^3 + 3525*n^2 - 6823.3*n + 5520.33
}

xy2CCT.HernandezAndres<-function(x,y) {
# convert from chromaticity coordinates to correlated color temperature (approximation)
# Hernandez-Andres, et al. 1999: "Calculating correlated color temperatures across the entire gamut of daylight and skylight chromaticities"
# http://en.wikipedia.org/wiki/Color_temperature
xe<-0.3366;xe2<-0.3356
ye<-0.1735;ye2<-0.1691
A0<- -949.86315;A02<-36284.48953
A1<-6253.80338;A12<-0.00228
t1<-0.92159;t12<-0.07861
A2<-28.70599;A22<-5.4535*10^-36
t2<-0.20039;t22<-0.01543
A3<-0.00004
t3<-0.07125
n = (x - xe)/(y - ye)
n2 = (x - xe2)/(y - ye2)
CCT<-A0 + A1*exp(-n/t1) + A2*exp(-n/t2) + A3*exp(-n/t3)
if (CCT>50000) CCT<-A02 + A12*exp(-n2/t12) + A22*exp(-n2/t22)
CCT
}

XYZ2CCT.Robertson<-function(X,Y,Z) {
#XYZ to Correlated Color Temperature (method developed by A. R. Robertson)
#Copyright (c) 2003 Bruce Justin Lindbloom. All rights reserved.
#http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_T.html
rt <- c(.Machine$double.eps,  10.0e-6,  20.0e-6,  30.0e-6,  40.0e-6,  50.0e-6,60.0e-6,  70.0e-6,  80.0e-6,  90.0e-6, 100.0e-6, 125.0e-6,
150.0e-6, 175.0e-6, 200.0e-6, 225.0e-6, 250.0e-6, 275.0e-6,300.0e-6, 325.0e-6, 350.0e-6, 375.0e-6, 400.0e-6, 425.0e-6,
450.0e-6, 475.0e-6, 500.0e-6, 525.0e-6, 550.0e-6, 575.0e-6,600.0e-6)# reciprocal temperature (K)
uvt <- matrix(c(0.18006, 0.26352, -0.24341,0.18066, 0.26589, -0.25479,0.18133, 0.26846, -0.26876,0.18208, 0.27119, -0.28539,0.18293, 0.27407, 
-0.30470,0.18388, 0.27709, -0.32675,0.18494, 0.28021, -0.35156,0.18611, 0.28342, -0.37915,0.18740, 0.28668, -0.40955,0.18880, 
0.28997, -0.44278,0.19032, 0.29326, -0.47888,0.19462, 0.30141, -0.58204,0.19962, 0.30921, -0.70471,0.20525, 0.31647, -0.84901,
0.21142, 0.32312, -1.0182,0.21807, 0.32909, -1.2168,0.22511, 0.33439, -1.4512,0.23247, 0.33904, -1.7298,0.24010, 0.34308, 
-2.0637,0.24792, 0.34655, -2.4681,0.25591, 0.34951, -2.9641,0.26400, 0.35200, -3.5814,0.27218, 0.35407, -4.3633,0.28039, 
0.35577, -5.3762,0.28863, 0.35714, -6.7262,0.29685, 0.35823, -8.5955,0.30505, 0.35907, -11.324,0.31320, 0.35968, -15.628,
0.32129, 0.36011, -23.325,0.32931, 0.36038, -40.770,0.33724, 0.36051, -116.45),ncol=3,byrow=TRUE)
if (any(c(X,Y,Z)==0)) stop('X, Y or Z must not be zero')
us = (4.0 * X) / (X + 15.0 * Y + 3.0 * Z)
        vs = (6.0 * Y) / (X + 15.0 * Y + 3.0 * Z)
        dm = 0.0;
        for (i in 1:32) {
                di = (vs - uvt[i,2]) - uvt[i,3] * (us - uvt[i,1])
                if ((i > 0) & (((di < 0) & (dm >= 0)) | ((di >= 0) & (dm < 0)))) break # found lines bounding (us, vs) : i-1 and i
                dm = di
        }
        if (i == 31) stop('Bad input, color temp would be less than minimum of 1666.7 degrees, or too far towards blue')
        di = di / sqrt(1.0 + uvt[i,3] * uvt[i,3])
        dm = dm / sqrt(1.0 + uvt[i - 1,3] * uvt[i - 1,3])
        p = dm / (dm - di) # p = interpolation parameter, 0.0 : i-1, 1.0 : i
        p = 1.0 / (((rt[i]) - (rt[i - 1])) * (p) + (rt[i - 1]))
p
}

kelvin2xy<-function(T) {
# Blackbody radiator color temperature to CIE 1931 x,y chromaticity approximation function.
# Uses approximate formula of Kim et al. 2002: "Design of Advanced Color - Temperature Control System for HDTV Applications" 
# http://fcam.garage.maemo.org/apiDocs/namespace_f_cam.html
# chromaticity x coefficients for T <= 4000K
A.x00 <- -0.2661239
A.x01 <- -0.2343580
A.x02 <- 0.8776956
A.x03 <- 0.179910
# chromaticity x coefficients for T > 4000K
A.x10 <- -3.0258469
A.x11 <- 2.1070379
A.x12 <- 0.2226347
A.x13 <- 0.24039
# chromaticity y coefficients for T <= 2222K
A.y00 <- -1.1063814
A.y01 <- -1.34811020
A.y02 <- 2.18555832
A.y03 <- -0.20219683
# chromaticity y coefficients for 2222K < T <= 4000K
A.y10 <- -0.9549476
A.y11 <- -1.37418593
A.y12 <- 2.09137015
A.y13 <- -0.16748867
# chromaticity y coefficients for T > 4000K
A.y20 <- 3.0817580
A.y21 <- -5.87338670
A.y22 <- 3.75112997
A.y23 <- -0.37001483
invKiloK <- 1000.0/T
if (T <= 4000) {
xc <- A.x00*invKiloK*invKiloK*invKiloK + A.x01*invKiloK*invKiloK + A.x02*invKiloK + A.x03
} else {
xc <- A.x10*invKiloK*invKiloK*invKiloK + A.x11*invKiloK*invKiloK + A.x12*invKiloK + A.x13
}
if (T <= 2222) {
yc <- A.y00*xc*xc*xc +A.y01*xc*xc +A.y02*xc +A.y03
} else if (T <= 4000) {
yc <- A.y10*xc*xc*xc +A.y11*xc*xc +A.y12*xc +A.y13
} else {
yc <- A.y20*xc*xc*xc +A.y21*xc*xc +A.y22*xc +A.y23
}
c(x=xc,y=yc)
}

wlnm2XYZ<-function(wavelength) c(approx(colorscience::ciexyz31[,1],colorscience::ciexyz31[,2],wavelength)$y,approx(colorscience::ciexyz31[,1],colorscience::ciexyz31[,3],wavelength)$y,approx(colorscience::ciexyz31[,1],colorscience::ciexyz31[,4],wavelength)$y)

wlnm2xyz<-function(wavelength) c(approx(colorscience::cccie31[,1],colorscience::cccie31[,2],wavelength)$y,approx(colorscience::cccie31[,1],colorscience::cccie31[,3],wavelength)$y,approx(colorscience::cccie31[,1],colorscience::cccie31[,4],wavelength)$y)

emittanceblackbodyPlanck<-function(wlnm, T){
# emittance of a black body of temperature T at a given wavelength (in metres)
# Planck's radiation law
# http://www.fourmich/documents/specrend/specrend.c
wlm = wlnm * 1e-9
(3.74183e-16 * (wlm^-5.0)) / (exp(1.4388e-2 / (wlm * T)) - 1.0)
}

XYZTannenbaum1974<-function(wavelen) 
c(X=2.37e-5 * (     0.101*(wavelen-380.0)^6.0 ) * ( 6 - (0.101*(wavelen-380.0)*0.5) )^2 * exp(    -0.101*(wavelen-380.0)*0.5 ),
Y=7.585e-4 * (      0.085*(660.0-wavelen)^6.0 ) * exp( -2.0*0.085*(660.0-wavelen)/3.0 ) * ifelse(( 660.0-wavelen )>=0,1,0),
Z=4.167e-1 * ( 0.093 * (wavelen-400.0)^4.0 ) * exp( -0.093 * (wavelen-400.0) ) * ifelse(( 400.0-wavelen )>=0,1,0))
#Chris Wyman Peter-Pike Sloan Peter Shirley Journal of Computer Graphics Techniques Vol. 2, No. 2, 2013
#Simple Analytic Approximations to the CIE XYZ Color Matching Functions

XYZMoonSpencer1945<-function(wavelen)
{#Chris Wyman Peter-Pike Sloan Peter Shirley Journal of Computer Graphics Techniques Vol. 2, No. 2, 2013
#Simple Analytic Approximations to the CIE XYZ Color Matching Functions
wav = wavelen/1000.0
xA = ( 1.237894975e+32 / ( wav^ 365.3388 ) )*exp( -164.9506 / wav )
xB = ( 3.247132422e+54 / ( wav^ 500.0    ) )*exp( -237.5    / wav )
xC = ( 3.933688807e+69 / ( wav^ 336.0385 ) )*exp( -199.1054 / wav )
yA = ( 2.592984973e+32 / ( wav^ 182.1905 ) )*exp( -100.9370 / wav )
yB = ( 0.000010681     / ( wav^ 3.205    ) )*exp( 0.1282    / wav )
yC = ( 1.44280378e+69  / ( wav^ 336.0385 ) )*exp( -199.1054 / wav )
yD = ( 1.07176606e+262 / ( wav^ 1013.5   ) )*exp( -679.045  / wav )
c(X=(xA - xB + xC), Y=ifelse(wav < 0.62, yA, ifelse(wav < 0.67, yA-yB,yC)), Z=(( 5.657180035e+32 / ( wav ^ 365.3388 ) )*exp( -164.9506 / wav )))
}

sprague<-function(spectra,f)
{
# Stephen Westland
# http://www.mathworks.com/matlabcentral/fileexchange/40640-computational-colour-science-using-matlab-2e/content/sprague.m
# Interpolates an n by w matrix of spectra 
# where n is the number of spectra 
# f is an interpolation factor
# e.g. if f=2 the sampling rate is doubled
if (f<2 | ((f-floor(f))>0)) stop('invalid f value - premature termination')
# set the parameters
c1<-c(884,-1960,3033,-2648,1080,-180,508,-540,488,-367,144,-24,-24,144,-367,488,-540,508,-180,1080,-2648,3033,-1960,884)
numSpectra<-dim(spectra)[1]
lengthSpectra<-dim(spectra)[2]
for (i in 1:numSpectra)
{
# select a spectrum
r <- spectra[i,]
# add the extra start and end points
k <- c1[1,]
p1 <- (k * t(r[1:6]))/209
k <- c1[2,]
p2 <- (k * t(r[1:6]))/209
k <- c1[3,]
endV<-length(r)
p3 <- (k * t(r[(endV-5):endV]))/209
k <- c1[4,]
p4 <- (k * t(r[(endV-5):endV]))/209
r <- c(p1, p2, r, p3, p4)
N <- lengthSpectra+4
p <- matrix(0,numSpectra,f*(N-5)+1) 
xx <- seq(1/f,1-1/f,len=(f-1))
for (j in 3:(N-3))
{
a0 <- r[j]
a1 <- (2*r[j-2]-16*r[j-1]+16*r[j+1]-2*r[j+2])/24
a2 <- (-r[j-2]+16*r[j-1]-30*r[j]+16*r[j+1]-r[j+2])/24
a3 <- (-9*r[j-2]+39*r[j-1]-70*r[j]+66*r[j+1]-33*r[j+2]+7*r[j+3])/24
a4 <- (13*r[j-2]-64*r[j-1]+126*r[j]-124*r[j+1]+61*r[j+2]-12*r[j+3])/24
a5 <- (-5*r[j-2]+25*r[j-1]-50*r[j]+50*r[j+1]-25*r[j+2]+5*r[j+3])/24
y <- a0+a1*xx+a2*xx^2+a3*xx^3+a4*xx^4+a5*xx^5
index <- j-2
p[i,(index-1)*f+1] <- r(j)
p[i,((index-1)*f+1+1):((index-1)*f+1+f-1)] <- y
}
p[i,f*(N-5)+1] <- r[N-2]
}
p
}

spectra2ISObrightness<-function(spectraIn=NA, wlIn=NA, RSDmatrix=colorscience::ISObrightnessReflectometerRSD){
# Diffuse blue reflectance factor (ISO brightness), R457,  ISO 2470
# ISO 2470-1 : 2009 PAPER, BOARD AND PULPS - MEASUREMENT OF DIFFUSE BLUE REFLECTANCE FACTOR, PART 1 INDOOR DAYLIGHT CONDITIONS (ISO BRIGHTNESS)
#require(Hmisc)
if (any(is.na(spectraIn))) stop('<<spectraIn>> must be a numeric vector')
if (!is.numeric(spectraIn)) stop('<<spectraIn>> must be a numeric vector')
if (any(is.na(wlIn))) stop('<<wlIn>> must be a numeric vector')
if (!is.numeric(wlIn)) stop('<<wlIn>> must be a numeric vector')
wlMin<- wlIn[1]
wlMax<- wlIn[2]
# extrapolate
f <- approxfun(spectraIn, wlIn)
xtrapolated<-approxExtrap(spectraIn, wlIn, RSDmatrix[,1])$y
xtrapolated[which(xtrapolated<0)]<-0
sum(matrix(xtrapolated,ncol=1) * RSDmatrix[,2]) / sum(RSDmatrix[,2])
}

RxRyRz2XYZ<-function(RxRyRzmatrix=NA,illuminant='C', observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser){
# convert from three filter measurements (reflectance factors) to XYZ
if (is.null(dim(RxRyRzmatrix))) if (length(RxRyRzmatrix)>2) RxRyRzmatrix<-matrix(RxRyRzmatrix, ncol=3,byrow=TRUE)
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
if (illuminant=='C' & observer==2) return(cbind(X = 78.321 * RxRyRzmatrix[,1] + 19.753 * RxRyRzmatrix[,3], Y = 100 * RxRyRzmatrix[,2], Z = 118.232 * RxRyRzmatrix[,3]))
if (illuminant=='D65' & observer==10) return(cbind(X = 76.841 * RxRyRzmatrix[,1] + 17.970 * RxRyRzmatrix[,3], Y = 100 * RxRyRzmatrix[,2], Z = 107.304 * RxRyRzmatrix[,3]))
stop('<<illuminant>> and <<observer>> must be C/2 or D65/10')
}

XYZ2RxRyRz<-function(XYZmatrix=NA,illuminant='C', observer=2,RefWhite=colorscience::XYZperfectreflectingdiffuser){
# convert from XYZ to three filter measurements (reflectance factors)
if (is.null(dim(XYZmatrix))) if (length(XYZmatrix)>2) XYZmatrix<-matrix(XYZmatrix, ncol=3,byrow=TRUE)
R<-RefWhite[which(RefWhite[["Illuminant"]]==illuminant ),]
Rx<-unlist(R[paste('X',observer,sep='')])
Ry<-unlist(R[paste('Y',observer,sep='')])
Rz<-unlist(R[paste('Z',observer,sep='')])
if (illuminant=='C' & observer==2) return(cbind(Rx = ( XYZmatrix[,1] - 0.16707 * XYZmatrix[,3] ) * 78.321, Ry = XYZmatrix[,2]  *100, Rz = XYZmatrix[,3] * 118.232))
if (illuminant=='D65' & observer==10) return(cbind(Rx = ( XYZmatrix[,1] * 10 - 0.16747 * XYZmatrix[,3] ) * 76.841, Ry = XYZmatrix[,2] * 100, Rz = XYZmatrix[,3] * 107.304))
stop('<<illuminant>> and <<observer>> must be C/2 or D65/10')
}

RGB2PhotoYCC<-function(RGBmatrix)
{# based on: ColorObject.pm Graphics/ColorObject version 0.5.0, Copyright 2003-2005 by Alex Izvorski, (Portions Copyright 2001-2003 by Alfred Reibenschuh)
# reference: http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.txt
# reference: http://wwwde.kodak.com/global/en/professional/products/storage/pcd/techInfo/pcd-045.jhtml
# input should be CIE Rec 709 non-linear rgb
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
t(matrix(c(0, 156, 137),dim(RGBmatrix)[1],3,byrow=TRUE)) + matrix(c(1/(255/1.402), 0, 0,0, 1/111.40, 0,0, 0, 1/135.64),3,3,byrow=TRUE) %*% 
matrix(c(0.299,0.587,0.114,-0.299,-0.587,0.866,0.701,-0.587,-0.114),3,3,byrow=TRUE) %*% t(RGBmatrix)
}

PhotoYCC2RGB<-function(PhotoYCCmatrix)
{# based on: ColorObject.pm Graphics/ColorObject version 0.5.0, Copyright 2003-2005 by Alex Izvorski, (Portions Copyright 2001-2003 by Alfred Reibenschuh)
# reference: http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.txt
# reference: http://wwwde.kodak.com/global/en/professional/products/storage/pcd/techInfo/pcd-045.jhtml
# result is CIE 709 non-linear rgb
if (is.null(dim(PhotoYCCmatrix))) if (length(PhotoYCCmatrix)>2) PhotoYCCmatrix<-matrix(PhotoYCCmatrix, ncol=3,byrow=TRUE)
matrix(c(1.0,0.0,1.0,0.99603657,-0.19817126,-0.50936968,1.0204082,1.0204082,0.0),3,3,byrow=TRUE) %*% 
matrix(c(1/(255/1.402), 0, 0,0, 1/111.40, 0,0, 0, 1/135.64),3,3,byrow=TRUE) %*% t(matrix(c(0, -156, -137),dim(PhotoYCCmatrix)[1],3,byrow=TRUE)) + t(PhotoYCCmatrix)
}

xyChromaticitiesVos1978<-function(x,y) cbind((1.0271*x - 0.00008*y - 0.00009)/(0.03845*x + 0.01496*y + 1), (0.00376*x + 1.0072*y + 0.00764)/(0.03845*x + 0.01496*y +1))
# x, y coordinates transformed to Judd (1951) x', y' system

saturationCIELUV<-function(u,v,un,vn) 13*sqrt((u-un)^2+(v-vn)^2)
# CIELUV saturation = chroma normalized by lightness
# u, v = CIELUV coordinates
# un, vn = chromaticity of the white point
# Color by Wikipedians pp. 36

saturationCIELAB<-function(L,a,b) sqrt(a^2+b^2)/L
# CIELAB saturation = chroma normalized by lightness
# L, a, b = CIELAB coordinates
# Color by Wikipedians pp. 36

saturationCIELABEvaLubbe<-function(L,a,b) {
# CIELAB saturation = chroma normalized by lightness
# L, a, b = CIELAB coordinates
# Color by Wikipedians pp. 36
Cab<-sqrt(a^2+b^2)
Cab/sqrt(Cab^2+L^2)
}

saturationCIECAM02<-function(M, Q) sqrt(M/Q)
# square root of the colorfulness divided by the brightness
# Color by Wikipedians pp. 36



