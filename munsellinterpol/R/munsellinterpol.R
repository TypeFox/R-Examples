





hypot<-function(a, b){
# sqrt(a^2 + b^2) without under/overflow. **/
# http://www.java2s.com/Tutorial/Java/0120__Development/sqrta2b2withoutunderoverflow.htm
r<-0.0
if (abs(a) > abs(b)) {
         r <- b/a
         r <- abs(a)*sqrt(1+r^2)
      } else if (b != 0) {
         r <- a/b
         r <- abs(b)*sqrt(1+r^2)
      }
r
}

MunsellToXYZ<-function(MunsellSpec, InterpolateByLuminanceFactor=TRUE)
{#Convert a Munsell specification into XYZ coordinates, by interpolating over the extrapolated Munsell renotation data.
tmpxyY<-MunsellToxyY(MunsellSpec, InterpolateByLuminanceFactor)
x<-tmpxyY$x
y<-tmpxyY$y
Y<-tmpxyY$Y
Status.ind<-tmpxyY$Status.ind
if (Status.ind != 1) return(Status.ind=Status.ind)
XYZ<-xyY2XYZ(cbind(x,y,Y))
attributes(XYZ)$dimnames<-NULL
# Set successful status return code and return
list(X=XYZ[,1],Y=XYZ[,2],Z=XYZ[,3],Status.ind = 1)
}

MunsellToLab<-function(MunsellSpec, InterpolateByLuminanceFactor=TRUE)
{#Convert a Munsell specification into CIE Lab coordinates, by interpolating over the extrapolated Munsell renotation data.
tmpxyY<-MunsellToxyY(MunsellSpec, InterpolateByLuminanceFactor)
x<-tmpxyY$x
y<-tmpxyY$y
Y<-tmpxyY$Y
Status.ind<-tmpxyY$Status.ind
if (Status.ind != 1) return(Status.ind=Status.ind)
Lab<-xyz2lab(cbind(x,y,Y))
# Set successful status return code and return
list(L=Lab[,1],a=Lab[,2],b=Lab[,3],Status.ind = 1)
}

MunsellToLuv<-function(MunsellSpec, InterpolateByLuminanceFactor=TRUE)
{#Convert a Munsell specification into CIE Luv coordinates, by interpolating over the extrapolated Munsell renotation data.
tmpxyY<-MunsellToxyY(MunsellSpec, InterpolateByLuminanceFactor)
x<-tmpxyY$x
y<-tmpxyY$y
Y<-tmpxyY$Y
Status.ind<-tmpxyY$Status.ind
if (Status.ind != 1) return(Status.ind=Status.ind)
Luv<-xyz2luv(cbind(x,y,Y))
# Set successful status return code and return
list(L=Luv[,1],u=Luv[,2],v=Luv[,3],Status.ind = 1)
}

XYZtoMunsell <- function(XYZ)
{
xyY<-XYZ2xyY(XYZ)
tmpMunsell<-xyYtoMunsell(xyY)
Status.ind<-tmpMunsell$Status.ind
if (Status.ind != 1){
# If the routine did not exit from the while loop, then the maximum number of
# iterations has been tried.  Set an error message and return.
Status.dist <- tmpMunsell$Status.dist
Status.num <- tmpMunsell$Status.num
return(list(Status.ind  = Status.ind, Status.dist = Status.dist, Status.num  = Status.num))
}
MunsellVec <- tmpMunsell$MunsellVec
MunsellSpec <- tmpMunsell$MunsellSpec
return(list(MunsellVec  = MunsellVec,MunsellSpec = MunsellSpec,Status.ind  = 1))
}

labtoMunsell <- function(Lab)
{
xyz<-lab2xyz(Lab)
tmpMunsell<-XYZtoMunsell(xyz)
Status.ind<-tmpMunsell$Status.ind
if (Status.ind != 1){
# If the routine did not exit from the while loop, then the maximum number of
# iterations has been tried.  Set an error message and return.
Status.dist <- tmpMunsell$Status.dist
Status.num <- tmpMunsell$Status.num
return(list(Status.ind  = Status.ind, Status.dist = Status.dist, Status.num  = Status.num))
}
MunsellVec <- tmpMunsell$MunsellVec
MunsellSpec <- tmpMunsell$MunsellSpec
return(list(MunsellVec  = MunsellVec,MunsellSpec = MunsellSpec,Status.ind  = 1))
}

luvtoMunsell <- function(Luv)
{
xyz<-luv2xyz(Luv)[1,]
tmpMunsell<-XYZtoMunsell(xyz)
Status.ind<-tmpMunsell$Status.ind
if (Status.ind != 1){
# If the routine did not exit from the while loop, then the maximum number of
# iterations has been tried.  Set an error message and return.
Status.dist <- tmpMunsell$Status.dist
Status.num <- tmpMunsell$Status.num
return(list(Status.ind  = Status.ind, Status.dist = Status.dist, Status.num  = Status.num))
}
MunsellVec <- tmpMunsell$MunsellVec
MunsellSpec <- tmpMunsell$MunsellSpec
return(list(MunsellVec  = MunsellVec,MunsellSpec = MunsellSpec,Status.ind  = 1))
}

XYZ2xyY <- function(XYZ){
# Purpose		Convert XYZ coordinates to xyY coordinates.
if (is.null(dim(XYZ))) XYZ <- matrix(XYZ,1)
w <- which(XYZ[,2] != 0.0)
xyY <- matrix(0,dim(XYZ)[1],3)
if (length(w)>0){
xyY[,3] <- XYZ[w,2]
xyY[w,1] <- XYZ[w,1]  /(XYZ[w,1]+XYZ[w,2]+XYZ[w,3])
xyY[w,2] <- XYZ[w,2]  /(XYZ[w,1]+XYZ[w,2]+XYZ[w,3])
}
xyY
}

xyYtoMunsell <- function(xyY){
if ((length(xyY) %% 3) != 0)  stop('xyY must be n x 3')
if (is.null(dim(xyY))) if (length(xyY)>2) xyY<-matrix(xyY, ncol=3,byrow=TRUE)
# Assign default output values
MunsellVec <- -99
MunsellSpec <- 'ERROR'
Status.ind <- -99
# If the colour stimulus is beyond the MacAdam limits, set an error code and return.
# Determine the MacAdam limits with regard to illuminant C, which is the Munsell standard.
if (!IsWithinMacAdamLimits(xyY, 'C')) return(list(Status.ind = 4))
# The Munsell value is calculated directly from the luminance factor,
# using a lookup table of evaluations of [ASTMD1535-08, Eq. 2].
MunsellValues <- LuminanceFactorToMunsellValue(xyY[,3])
MunsellValue  <- MunsellValues$ASTMTableLookup
# To avoid numerical issues with approximation, set MunsellValue to an
# integer if it is very close to an integer.
if (abs(MunsellValue - round(MunsellValue)) < 0.001) MunsellValue <- round(MunsellValue)
# Express xy chromaticity coordinates as polar coordinates around the achromatic
# point with Munsell value corresponding to Y.
tmp <- MunsellToxyY(MunsellValue)
xcenter <- tmp$x
ycenter  <- tmp$y
Ytemp <- tmp$Y
Status.ind <- tmp$Status.ind
if (Status.ind != 1) return(list(Status.dist = -99,Status.num  = -99,Status.ind = 2))
# rInput and thetaInput are the values the inversion algorithm will attempt to match
x <- xyY[,1];y <- xyY[,2]
thetaInput <- atan2(y-ycenter, x-xcenter)
#rInput <- sqrt(abs(x-xcenter)^2 + abs(y-ycenter)^2)
rInput <- hypot(x-xcenter, y-ycenter)
thetaInput <- ((180/pi) * thetaInput) %% 360 # Convert to degrees
thetaR <- c(thetaInput, rInput)
# Use the following parameters for the inversion algorithm.  ConvergenceThreshold
# is the required Euclidean distance in xy coordinates between the input
# xy, and the xy corresponding to a Munsell sample.  GreyThreshold is the required
# distance from the origin for a colour to be considered grey; it is defined 
# separately from ConvergenceThreshold for numerical reasons.   
ConvergenceThreshold <- 0.0001
GreyThreshold <- 0.001
MaxNumOfTries <- 60
NumOfTries <- 0
# Check for grey
if (rInput < GreyThreshold) return(list(MunsellVec = MunsellValue, MunsellSpec = ColorLabFormatToMunsellSpec(MunsellValue)$MunsellSpecString,Status.ind = 2))
# The input xyY are the CIE coordinates of Illuminant C, when reflected off a Munsell sample
XYZ <- xyY2XYZ(xyY)	
ReflectedXYZ <- XYZ
# Since Y is a percentage reflection, we want to find the CIE coordinates
# for illuminant C, whose Y value is 100.  This is the reference
# illuminant which is reflected off the Munsell sample.
# The chromaticity coordinates for Illuminant C are
   xC <- 0.310062 # whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'x2' ]
   yC <- 0.316148 # whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'y2' ]
# The line above replaces the following line (replaced Jan. 1, 2014, by Paul Centore)
# [xC yC] = IlluminantCWhitePoint();
# Convert to XYZ form, and normalize to have luminance 100.  This will
# be the reference illuminant.
XrYrZr <- xyY2XYZ(cbind(xC, yC, xyY[,3]))
ReferenceXYZ <- cbind((100/XrYrZr[,2])*XrYrZr[,1], 100, (100/XrYrZr[,2])*XrYrZr[,3])
# ReferenceXYZ are the CIE coordinates for a reference light, whose power spectral distribution
# is consistent with Illuminant C, and whose luminance is 100.  The relative
# luminance of the input sample, when illuminated by Illuminant C, is given by Y, which
# varies between 0 and 100, and is the percentage of illuminant C that the sample
# reflects.  Therefore, RefXYZ can be viewed as the reference illuminant for the 
# sample, and XYZ as the CIE coordinates for the reflected illuminant.  These are the
# data needed to apply the CIELAB model.
Lab <- xyz2lab(ReflectedXYZ, ReferenceXYZ)
LhC <- lab2perc(Lab)
# Find an initial Munsell specification for the interpolation algorithm, by using a
# rough transformation from CIELAB to Munsell coordinates.  InitialMunsellSpec is
# not an exact inverse for xyY, but should be close.
L   <- LhC[1]		# Lightness in CIELAB system
hab <- LhC[2]		# Hue angle in CIELAB system
Cab <- LhC[3]		# C*ab in CIELAB system
InitialMunsellSpec <- CIELABtoApproxMunsellSpec(L, Cab, hab)
# Replace value of initial Munsell estimate with known Munsell value
TempCLvec    <- MunsellSpecToColorLabFormat(InitialMunsellSpec)
TempCLvec[2] <- MunsellValue									
# Deliberately underestimate chroma to avoid chromas beyond extrapolated
# renotation values
TempCLvec[3] <- (5/5.5) * TempCLvec[3]

# Use a loop to iterate over different Munsell specifications, whose xyY coordinates
# should progressively approach the input xyY.
CurrentCLVec  <- TempCLvec
EuclideanDifference <- -99
while (NumOfTries <= MaxNumOfTries){
   # After adjusting chroma without adjusting hue, the next iteration adjusts hue
   # without adjusting chroma.
   NumOfTries <- NumOfTries + 1
   # Extract Munsell quantities from Munsell specification
   CurrentCLHueNumber <- CurrentCLVec[1]
   CurrentMunsellChroma <- CurrentCLVec[3]
   CurrentCLHueLetterIndex <- CurrentCLVec[4]
   CurrentChromDiagHueAngle <- MunsellHueToChromDiagHueAngle(CurrentCLHueNumber,CurrentCLHueLetterIndex)
   # Check that current chroma is possible for the hue and value in the Munsell specification.  If the
   # chroma is too high, set it to be the maximum chroma.
   tmp <- MaxChromaForExtrapolatedRenotation(CurrentCLHueNumber, CurrentCLHueLetterIndex, MunsellValue)
   MaxChroma <- tmp$MaxChroma
   StatusCode <- tmp$StatusCode
   if (CurrentMunsellChroma > MaxChroma){
      CurrentMunsellChroma <- MaxChroma
	  CurrentCLVec[3] <- CurrentMunsellChroma
   }
   # Find xy coordinates corresponding to the current Munsell specification.
   
   if (length(CurrentCLVec)==2) stop('164')
   
   tmp = MunsellToxyY(CurrentCLVec)
   xCur <- tmp$x
   yCur <- tmp$y
   YCur <- tmp$Y
   Status.ind <- tmp$Status.ind
      if (Status.ind != 1) return(list(Status.ind = 2,Status.dist = EuclideanDifference,Status.num  = NumOfTries))
   # Hue angles correspond to values of theta.  The new hue angle will differ from the
   # current hue angle approximately as much as the desired theta value (thetaInput) differs from the
   # current theta value.  Call this difference thetaDiff and calculate it.  Other
   # hue angles, with their corresponding thetas and theta differences, will be tried.
   # Make a list, thetaDiffsVec, of the corresponding theta differences.
   thetaCur <- atan2(yCur-ycenter, xCur-xcenter)
#rCur <- sqrt(abs(xCur-xcenter)^2 + abs(yCur-ycenter)^2)
rCur <- hypot((xCur-xcenter), (yCur-ycenter))
thetaCur <- ((180/pi) * thetaCur) %% 360 # Convert to degrees
   thetaDiff <- (360 - thetaInput + thetaCur) %% 360
   if (thetaDiff > 180) thetaDiff <- thetaDiff-360 # Adjust for wraparound if necessary
   thetaDiffsVec <- thetaDiff
   # Start a similar list for hue angles that correspond to the theta differences.
   ChromDiagHueAngles <- CurrentChromDiagHueAngle
   # Also make a list of much the new hue angles differ from CurrentChromDiagHueAngle.  These
   # angles will be near zero, and will avoid potential problems with wraparound.
   ChromDiagHueAngleDiffs <- 0
   # Ideally, thetaDiff will be 0
   # (thetaInput will agree with the theta corresponding to the new hue angle).  
   # Continue constructing the list of theta differences until it contains both
   # negative and positive differences; then find the new hue angle by linear interpolation
   # at the value thetaDiff = 0.
   ctr = 0
   AttemptExtrapolation = FALSE
   while (sign(min(thetaDiffsVec)) == sign(max(thetaDiffsVec))  && AttemptExtrapolation == FALSE){
      ctr = ctr + 1

	  if (ctr > 10) return(list(Status.ind = 3))# Too many attempts.  Return with error message
  
	  # Find another hue angle, by increasing the coefficient of (thetaInput-thetaCur).
	  # Construct a trial Munsell specification by using the new hue angle.
	  TempChromDiagHueAngle <- (CurrentChromDiagHueAngle + ctr * (thetaInput - thetaCur)) %% 360
	  # Record the difference of the trial angle from the current hue angle
	  TempChromDiagHueAngleDiff <- (ctr * (thetaInput - thetaCur)) %% 360
	  if (TempChromDiagHueAngleDiff > 180) TempChromDiagHueAngleDiff <- TempChromDiagHueAngleDiff - 360

      tmp <- ChromDiagHueAngleToMunsellHue(TempChromDiagHueAngle)
      TempHueNumber <- tmp$HueNumber
      TempHueLetterCode <- tmp$HueLetterCode
      TempCLVec <- c(TempHueNumber, MunsellValue, CurrentMunsellChroma, TempHueLetterCode)
	  
	  # Evaluate the trial Munsell specification, convert to polar coordinates, and 
	  # calculate the difference of the resulting trial theta and thetaInput
  
      tmp <- MunsellToxyY(TempCLVec)
      xCurTemp <- tmp$x
      yCurTemp <- tmp$x
      YCurTemp <- tmp$Y
      Status.ind <- tmp$Status.ind
      
	  if (Status.ind != 1){
	# Interpolation is impossible because there are not both positive and negative differences,
	# but extrapolation is possible if there are at least two data points already
	     if (length(thetaDiffsVec) >= 2) AttemptExtrapolation <- TRUE else return(list(Status.dist = EuclideanDifference,Status.num  = NumOfTries,Status.ind=2))
	  }
	  if (AttemptExtrapolation == FALSE){
     
	        thetaCurTemp <- atan2(yCurTemp-ycenter, xCurTemp-xcenter)
#rCurTemp <- sqrt(abs(xCurTemp-xcenter)^2 + abs(yCurTemp-ycenter)^2)
rCurTemp <- hypot((xCurTemp-xcenter), (yCurTemp-ycenter))
         thetaCurTemp           <- ((180/pi)*thetaCurTemp) %% 360 # Express in degrees
         thetaDiff              <- (360 - thetaInput + thetaCurTemp) %% 360
         if (thetaDiff > 180)	thetaDiff = thetaDiff-360	# Adjust for wraparound if necessary

     # Add trial hue angle and theta difference to lists
         thetaDiffsVec          = cbind(thetaDiffsVec, thetaDiff)
 	     ChromDiagHueAngleDiffs <- cbind(ChromDiagHueAngleDiffs, TempChromDiagHueAngleDiff)
	     ChromDiagHueAngles     <- cbind(ChromDiagHueAngles, TempChromDiagHueAngle)
	  }
   }
   # Since the while loop exited successfully, both negative and positive theta
   # differences have been found, or an extrapolation should
   # be attempted.  Interpolate linearly to estimate the hue
   # angle that corresponds to a theta difference of 0
   thetaDiffsVecSort <- sort(thetaDiffsVec)
   Iord <- order(thetaDiffsVec)
   ChromDiagHueAnglesSort <- ChromDiagHueAngles[Iord]
   ChromDiagHueAngleDiffs <- ChromDiagHueAngleDiffs[Iord]
   # The extrapolation option will be used if there is not sufficient data
   # for interpolation
   NewChromDiagHueAngleDiff  <- (approx(thetaDiffsVecSort, ChromDiagHueAngleDiffs, 0, 'linear'))$y %% 360
   
   NewChromDiagHueAngle      <- (CurrentChromDiagHueAngle + NewChromDiagHueAngleDiff) %% 360

   # Adjust the current Munsell specification by replacing the current hue with the
   # new hue.
   tmp <- ChromDiagHueAngleToMunsellHue(NewChromDiagHueAngle)
   NewHueNumber <- tmp$HueNumber
   NewHueLetterCode <- tmp$HueLetterCode
   CurrentCLVec <- cbind(NewHueNumber, MunsellValue, CurrentMunsellChroma, NewHueLetterCode)

   # Calculate the Euclidean distance between the xy coordinates of the newly
   # constructed Munsell specification, and the input xy
   
   if (length(CurrentCLVec)==2) stop('267')
   
   tmp <- MunsellToxyY(CurrentCLVec)
   xCur <- tmp$x
   yCur <- tmp$y
   YCur <- tmp$Y
   Status.ind <- tmp$Status.ind

   if (Status.ind != 1) return(list(Status.ind  = 2,Status.dist = EuclideanDifference,Status.num  = NumOfTries))

   EuclideanDifference  <- sqrt(((x-xCur)*(x-xCur)) + ((y-yCur)*(y-yCur)))

   # If the two xy coordinate pairs are close enough, then exit with a success message
   if (EuclideanDifference < ConvergenceThreshold) return(list(Status.ind  = 1,MunsellVec = CurrentCLVec,MunsellSpec = ColorLabFormatToMunsellSpec(CurrentCLVec)$MunsellSpecString))  # Current Munsell spec is close enough
   
   NumOfTries = NumOfTries + 1
 
   # Extract Munsell quantities from Munsell specification
   CurrentCLHueNumber       = CurrentCLVec[1]
   CurrentMunsellChroma     = CurrentCLVec[3]
   CurrentCLHueLetterIndex  = CurrentCLVec[4]
   CurrentChromDiagHueAngle = MunsellHueToChromDiagHueAngle(CurrentCLHueNumber,CurrentCLHueLetterIndex)
   
   # Check that current chroma is possible for the hue and value in the Munsell specification.  If the
   # chroma is too high, set it to be the maximum chroma.
   tmp <- MaxChromaForExtrapolatedRenotation(CurrentCLHueNumber, CurrentCLHueLetterIndex, MunsellValue)
   MaxChroma <- tmp$MaxChroma
   Status.ind <- tmp$Status.ind
   if (CurrentMunsellChroma > MaxChroma){
      CurrentMunsellChroma <- MaxChroma
	  CurrentCLVec[3] <- CurrentMunsellChroma
   }
   
   # Find xy coordinates corresponding to the current Munsell specification.
   
   if (length(CurrentCLVec)==2) stop('303')
   
   tmp <- MunsellToxyY(CurrentCLVec)
   xCur <- tmp$x
   yCur <- tmp$y
   YCur <- tmp$Y
   Status.ind <- tmp$Status.ind
      
   if (Status.ind != 1) return(list(Status.ind = 2,Status.dist = EuclideanDifference,Status.num  = NumOfTries))

   thetaCur <- atan2(yCur-ycenter, xCur-xcenter)
#rCur <- sqrt(abs(xCur-xcenter)^2 + abs(yCur-ycenter)^2)
rCur <- hypot((xCur-xcenter), (yCur-ycenter))
   thetaCur <- ((180/pi)*thetaCur) %% 360					# Express theta of current point in degrees

   # For this iteration, keep hue (corresponding to thetaCur) constant,
   # and let chroma (corresponding to rCur) vary.
   # Construct a set of chromas, called TempMunsellChromas, and find the (r, theta)
   # values for Munsell specifications with those chromas.
   # Make a list of r values, called rTempValues, that correspond to different chromas
   rTempValues <- rCur
   TempMunsellChromas <- CurrentMunsellChroma
if (NumOfTries >= 2000) OneChroma <- cbind(NumOfTries, rTempValues, TempMunsellChromas) # Use for debugging non-converging cases

   # In order to interpolate, rInput must be within the span of the r values
   # in rTempValues.  Check this condition, and add more chromas and r values
   # until it is satisfied.  Usually, no more than three chromas are necessary.
   ctr = 0
   while (rInput < min(rTempValues) || rInput > max(rTempValues)){
      ctr = ctr + 1
	  if (ctr > 10) return(list(Status.ind = 3))		# Too many attempts to bound rInput.  Return with error message
	  
	  # Try a new chroma, by increasing the exponent on rInput/rCur
	  TempMunsellChroma = ((rInput/rCur)^ctr) * CurrentMunsellChroma
	  # Check that current chroma is possible for the hue and value in the Munsell specification.  If the
      # chroma is too high, set it to be the maximum chroma.
      if (TempMunsellChroma > MaxChroma){
          TempMunsellChroma = MaxChroma			
	  CurrentCLVec[3]   = TempMunsellChroma	
      }

      # Find (r, theta) values for a Munsell specification which is identical to
      # the current Munsell specification, except that the chroma is TempMunsellChroma
	  TempCLVec         = cbind(CurrentCLHueNumber, MunsellValue, TempMunsellChroma, CurrentCLHueLetterIndex)
	  
	  if (length(TempCLVec)==2) stop('348')
	  
      tmp <- MunsellToxyY(TempCLVec)
         xCurTemp <- tmp$x
   yCurTemp <- tmp$y
   YCurTemp <- tmp$Y
   Status.ind <- tmp$Status.ind

	  if (Status.ind != 1) return(list(Status.ind = 2,Status.dist = EuclideanDifference,Status.num  = NumOfTries))

	  thetaCurTemp <- atan2(yCurTemp-ycenter, xCur-xcenter)
#rCurTemp <- sqrt(abs(xCur-xcenter)^2 + abs(yCurTemp-ycenter)^2)
rCurTemp <- hypot((xCur-xcenter), (yCurTemp-ycenter))
      thetaCurTemp  <- ((180/pi)*thetaCurTemp) %% 360 	# Express in degrees

	  # Add r and chroma to lists
      rTempValues        = cbind(rTempValues, rCurTemp)
      TempMunsellChromas = cbind(TempMunsellChromas, TempMunsellChroma)

   }
   
   # Since the while loop has been exited, we have found r values, resulting from
   # different chroma values, that bound rInput.  Linearly interpolate to find a further
   # chroma, which should be an even better approximation.
   # Linear interpolation requires the list of r values to be sorted. 
   rTempSort <- sort(rTempValues)
   Iord <- order(rTempValues)
   TempMunsellChromasSort = TempMunsellChromas[Iord]
   NewMunsellChroma       = approx(rTempSort, TempMunsellChromasSort, rInput)$y
 
   # Adjust the current Munsell specification, by using the new chroma
   CurrentCLVec                = cbind(CurrentCLHueNumber, MunsellValue, NewMunsellChroma, CurrentCLHueLetterIndex)
   
   # Calculate the Euclidean distance between the xy coordinates of the newly
   # constructed Munsell specification, and the input xy
   
   if (length(CurrentCLVec)==2) stop('384')
   
   tmp <- MunsellToxyY(CurrentCLVec)
   xCur <- tmp$x
   yCur <- tmp$y
   YCur <- tmp$Y
   Status.ind <- tmp$Status.ind

   if (Status.ind != 1) return(list(Status.ind = 2,Status.dist = EuclideanDifference,Status.num  = NumOfTries))

   EuclideanDifference = sqrt(((x-xCur)*(x-xCur)) + ((y-yCur)*(y-yCur)))

   # If the two xy coordinate pairs are close enough, then exit with a success message
   if (EuclideanDifference < ConvergenceThreshold) return(list(MunsellVec  = CurrentCLVec,MunsellSpec = ColorLabFormatToMunsellSpec(CurrentCLVec)$MunsellSpecString,Status.ind  = 1))

}

# If the routine did not exit from the while loop, then the maximum number of
# iterations has been tried.  Set an error message and return.
return(list(Status.ind  = 3, Status.dist = EuclideanDifference, Status.num  = NumOfTries))
}



sRGBtoMunsell <- function(sRGB){
# Assign default values
MunsellSpec         <- '-99'
MunsellVec          <- -99
InMacAdamLimitsFlag <- -99
Status.ind          <- -99
# Convert sRGB input into CIE XYZ and xyY coordinates
XYZ <- srgb2xyz(sRGB)
xyY <- XYZ2xyY(XYZ)
# For convenience, return that a blank screen is ideal black
if (max(xyY) == 0) return(list(MunsellSpec = 'N0',MunsellVec=0,InMacAdamLimitsFlag = TRUE,Status.ind = 1))
# Convert CIE coordinates to Munsell coordinates
tmp <- xyYtoMunsell(cbind(xyY[,1], xyY[,2], 100*xyY[,3]))
MunsellSpec  <- tmp$MunsellSpec
MunsellVec  <- tmp$MunsellVec
Status.ind <- tmp$Status.ind
# Check whether CIE coordinates are outside MacAdam limits, so no conversion to Munsell is possible and whether other reasons caused the conversion to fail
if (Status.ind == 4) return(list(Status.ind = 2, InMacAdamLimitsFlag = FALSE)) else if (Status.ind == 2 || Status.ind == 3) return(list(Status.ind = 3))
# Since none of the error conditions hold, the conversion was successful
return(list(MunsellSpec=MunsellSpec, MunsellVec=MunsellVec, Status.ind = 1, InMacAdamLimitsFlag = TRUE))
}

srgb2xyz <- function(RGBmatrix){
if ((length(RGBmatrix) %% 3) != 0)  stop('RGBmatrix must be n x 3')
if (is.null(dim(RGBmatrix))) if (length(RGBmatrix)>2) RGBmatrix<-matrix(RGBmatrix, ncol=3,byrow=TRUE)
XYZmatrix <- matrix(0,dim(RGBmatrix)[1],3)
M <- matrix(c(0.4124, 0.3576, 0.1805, 0.2126, 0.7152, 0.0722, 0.0193, 0.1192, 0.9505),3,3,byrow=TRUE)
DACS <- RGBmatrix/255
RGBmatrix <- matrix(0,dim(DACS)[1],3)
index <- DACS <= 0.04045
index[index==TRUE] <- 1
index[index==FALSE] <- 0
RGBmatrix = RGBmatrix + (index)*(DACS/12.92)
RGBmatrix = RGBmatrix + (1-index)*((DACS+0.055)/1.055)^2.4
XYZmatrix <- t(M %*% t(RGBmatrix))
XYZmatrix
}

MunsellTosRGB <- function(MunsellSpec){
# Assign default values
R              = -99
G              = -99
B              = -99
OutOfGamutFlag = -99
Status.ind     = -99
# Convert from Munsell notation to CIE coordiantes
tmp <- MunsellToxyY(MunsellSpec)
x <- tmp$x
y <- tmp$y
Y <- tmp$Y
Status.ind <- tmp$Status.ind
if (Status.ind != 1) return(list(Status.ind  = Status.ind)) #  Failed to convert Munsell to CIE
# Convert CIE xyY coordinates to CIE XYZ coordinates
XYZ <- xyY2XYZ(cbind(x, y, Y)) * (1/100)
# Convert CIE XYZ coordinates to sRGB coordinates
tmp <- xyz2srgb(XYZ)
sRGB <- tmp$sRGB
OutOfGamutFlag <- tmp$OutOfGamutFlag
# Set successful status
return(list(Status.ind  = 1, sRGB=sRGB, OutOfGamutFlag=OutOfGamutFlag))
}

xyY2XYZ <- function(xyY){
if (is.null(dim(xyY))) xyY <- matrix(xyY,1)
XYZ <- cbind(0,xyY[,3],0)
w <- which(XYZ[,2] !=0)
XYZ[w,1] <- xyY[w,1]*xyY[w,3]/xyY[w,2]
XYZ[w,3] <- (1-xyY[w,1]-xyY[w,2])*xyY[w,3]/xyY[w,2]
XYZ
}

xyz2srgb <- function(XYZ){
if ((length(XYZ) %% 3) != 0)  stop('XYZ matrix must be n x 3')
if (is.null(dim(XYZ))) if (length(XYZ)>2) XYZ<-matrix(XYZ, ncol=3,byrow=TRUE)
M <- matrix(c(3.2406, -1.5372, -0.4986, -0.9689, 1.8758, 0.0415, 0.0557, -0.2040, 1.0570),3,3,byrow=TRUE)
RGB <- t(M %*% t(XYZ))
# START: lines added March 2013 to set out-of-gamut flag.  
# The out-of-gamut flag is a column vector of Boolean true/false values.  Each
# entry corresponds to one row of the input matrix XYZ.
NumberOfInputs <- dim(RGB)[1]
OutOfGamutFlag <- -99 * matrix(1,NumberOfInputs)
for (index in 1:NumberOfInputs){
if (RGB[index,1] < 0 || RGB[index,1] > 1 || RGB[index,2] < 0 || RGB[index,2] > 1 || RGB[index,3] < 0 || RGB[index,3] > 1) OutOfGamutFlag[index] = TRUE else OutOfGamutFlag[index] <- FALSE
}
# END: lines added March 2013 to set out-of-gamut flag
RGB[which(RGB<0)] <- 0
RGB[which(RGB>1)] <- 1
DACS = matrix(0,dim(XYZ)[1],3)
index <- RGB<=0.0031308
index[index==TRUE] <- 1
index[index==FALSE] <- 0
DACS <- DACS+index*(12.92*RGB)
DACS <- DACS+(1-index)*(1.055*RGB^(1/2.4)-0.055)
RGB <- ceiling(DACS*255)
RGB[which(RGB<0)] <- 0
RGB[which(RGB>255)] <- 255
return(list(Status.ind  = 1, sRGB=RGB, OutOfGamutFlag=OutOfGamutFlag))
}

IsWithinMacAdamLimits <- function(xyY, Illuminant){
# Select list of optimal colours from stored lists
data("OptimalColoursForIlluminantC", envir = environment())
data("OptimalColoursForIlluminantD65", envir = environment())
if ((length(xyY) %% 3) != 0)  stop('XYZ matrix must be n x 3')
if (is.null(dim(xyY))) if (length(xyY)>2) xyY<-matrix(xyY, ncol=3,byrow=TRUE)
if (Illuminant=='C') OptColoursInxyY = get("OptimalColoursForIlluminantC", envir  = environment()) else {
if (Illuminant=='D65') OptColoursInxyY = get("OptimalColoursForIlluminantD65", envir  = environment()) else stop('ERROR in data.')
}
# Convert optimal colours from xyY format to XYZ format.  In XYZ format, the
# colour solid is convex, so existing R routines can be used.
OptColoursInXYZ <- matrix(0,dim(OptColoursInxyY)[1],3)
row = dim(OptColoursInxyY)[1]
col = dim(OptColoursInxyY)[2]
for (i in 1:row){
   xOpt      = OptColoursInxyY[i,1]
   yOpt      = OptColoursInxyY[i,2]
   YOpt      = OptColoursInxyY[i,3]
   # Convert colour to XYZ coordinates
   XYZ <- xyY2XYZ(cbind(xOpt, yOpt, YOpt))
   # Store XYZ coordinates in matrix of optimal colours
   OptColoursInXYZ[i,] <- XYZ
}
# Convert input colour to XYZ coordinates
XYZ = xyY2XYZ(xyY)

   # Construct a Delaunay tessellation of the convex XYZ colour solid
   tessellation = delaunayn(OptColoursInXYZ)
   if (Illuminant=='C') CTessellation = tessellation else {
   if (Illuminant=='D65') D65Tessellation = tessellation else stop('ERROR in tessellation data.')
   }
# The routine tsearchn returns the index of the tetrahedron in the Delaunay
# tessellation that contains the input colour.  If no tetrahedron contains
# that colour, then tsearchn returns NaN.
idx = tsearchn(OptColoursInXYZ, tessellation, XYZ)
if (any(is.na(idx))) InsideLimits <- FALSE else InsideLimits <- TRUE
InsideLimits
}

xyz2lab <- function(C,bbb=cbind(95.047,100.000,108.883)){
# The following code made up the original ColorLab routine, before the Dec. 2013 revision
if (is.null(dim(C))) C <- matrix(C,1)
if (is.null(dim(bbb))) bbb <- matrix(bbb,1)
s=dim(C)
ss=dim(bbb)
if (ss[1] != s[1]) bbb=matrix(1,s[1],1) * bbb[1,]
lab <-matrix(0,s[1],3)
for (i in 1:s[1]){
    k=C[i,]
    b=bbb[i,]
    if ((k[2]/b[2])>0.008856){
       F2=(k[2]/b[2])^(1/3)
       L=116*F2-16
    }else{
       L=903.3*(k[2]/b[2])
       F2=(1/116)*(903.3*(k[2]/b[2])+16)
    }
    if ((k[1]/b[1])>0.008856) F1=(k[1]/b[1])^(1/3) else F1=(1/116)*(903.3*(k[1]/b[1])+16)
    if ((k[3]/b[3])>0.008856) F3=(k[3]/b[3])^(1/3) else F3=(1/116)*(903.3*(k[3]/b[3])+16)
    a=500*(F1-F2)
    bb=200*(F2-F3)
    lab[i,]=cbind(L, a, bb)
}
lab
}

lab2xyz <- function(C,bbb=c(95.047,100.000,108.883)){
# Logicol S.r.l., 2014
# EasyRGB color search engine
# http://www.easyrgb.com/
#default Observer= 2deg, Illuminant= D65
if (is.null(dim(C))) C <- matrix(C,1)
var.Y = ( C[,1] + 16 ) / 116
var.X = C[,2] / 500 + var.Y
var.Z = var.Y - C[,3] / 200
if ( var.Y^3 > 0.008856 ) var.Y = var.Y^3 else var.Y = ( var.Y - 16 / 116 ) / 7.787
if ( var.X^3 > 0.008856 ) var.X = var.X^3 else var.X = ( var.X - 16 / 116 ) / 7.787
if ( var.Z^3 > 0.008856 ) var.Z = var.Z^3 else var.Z = ( var.Z - 16 / 116 ) / 7.787
X = bbb[1] * var.X     #ref.X =  95.047     Observer= 2deg, Illuminant= D65
Y = bbb[2] * var.Y     #ref.Y = 100.000
Z = bbb[3] * var.Z     #ref.Z = 108.883
cbind(X,Y,Z)
}

luv2xyz <- function(C,bbb=c(95.047,100.000,108.883)){
# Logicol S.r.l., 2014
# EasyRGB color search engine
# http://www.easyrgb.com/
#default Observer= 2deg, Illuminant= D65
if (is.null(dim(C))) C <- matrix(C,1)
var.Y = ( C[,1] + 16 ) / 116
if ( var.Y^3 > 0.008856 ) var.Y = var.Y^3 else var.Y = ( var.Y - 16 / 116 ) / 7.787
ref.U = ( 4 * bbb[1] ) / ( bbb[1] + ( 15 * bbb[2] ) + ( 3 * bbb[3] ) )
ref.V = ( 9 * bbb[2] ) / ( bbb[1] + ( 15 * bbb[2] ) + ( 3 * bbb[3] ) )
var.U = C[,2] / ( 13 * C[1,] ) + ref.U
var.V = C[,3] / ( 13 * C[1,] ) + ref.V
Y = var.Y * 100
X =  - ( 9 * Y * var.U ) / ( ( var.U - 4 ) * var.V  - var.U * var.V )
Z = ( 9 * Y - ( 15 * var.V * Y ) - ( var.V * X ) ) / ( 3 * var.V )
cbind(X,Y,Z)
}

xyz2luv <- function(XYZ,bbb=c(95.047,100.000,108.883)){
# Logicol S.r.l., 2014
# EasyRGB color search engine
# http://www.easyrgb.com/
#default Observer= 2deg, Illuminant= D65
if (is.null(dim(XYZ))) XYZ <- matrix(XYZ,1)
var.U = ( 4 * XYZ[,1] ) / ( XYZ[,1] + ( 15 * XYZ[,2] ) + ( 3 * XYZ[,3] ) )
var.V = ( 9 * XYZ[,2] ) / ( XYZ[,1] + ( 15 * XYZ[,2] ) + ( 3 * XYZ[,3] ) )
var.Y = XYZ[,2] / 100
if ( var.Y > 0.008856 ) var.Y = var.Y ^ ( 1/3 ) else var.Y = ( 7.787 * var.Y ) + ( 16 / 116 )
ref.U = ( 4 * bbb[1] ) / ( bbb[1] + ( 15 * bbb[2] ) + ( 3 * bbb[3] ) )
ref.V = ( 9 * bbb[2] ) / ( bbb[1] + ( 15 * bbb[2] ) + ( 3 * bbb[3] ) )
L = ( 116 * var.Y ) - 16
u = 13 * L * ( var.U - ref.U )
v = 13 * L * ( var.V - ref.V )
cbind(L,u,v)
}

lab2perc <- function(lab){
# LAB2PERC computes the lightness, chroma and hue angle of a set of colours
# characterized in the CIELAB space.
perc <- cbind(lab[,1], atan2(lab[,3],lab[,2])+2*pi*(+(atan2(lab[,3],lab[,2])<0)), sqrt((lab[,3]^2)+(lab[,2]^2)) )
perc
}

CIELABtoApproxMunsellSpec <- function(L, Cab, hab){
# Convert a CIELAB specification to an approximate Munsell specification.
# Approximate Munsell hue first
# Convert hab from radians to degrees, and make sure it is between 0 and 360
habDegrees <- (180/pi)*hab
habDegrees <- (habDegrees %% 360)
# The Munsell hues are assumed to be evenly spaced on a circle, with 5Y
# at 90 degrees, 5G at 162 degrees, and so on.  Each letter code corresponds
# to a certain sector of the circle.  The following cases extract the letter code.
if (habDegrees == 0) HueLetterCode <- 8 else {
 if (habDegrees <= 36) HueLetterCode <- 7 
else { if (habDegrees <= 72) HueLetterCode <- 6 
else { if (habDegrees <= 108) HueLetterCode <- 5 
else { if (habDegrees <= 144) HueLetterCode <- 4 
else { if (habDegrees <= 180) HueLetterCode <- 3 
else { if (habDegrees <= 216) HueLetterCode <- 2 
else { if (habDegrees <= 252) HueLetterCode <- 1 
else { if (habDegrees <= 288) HueLetterCode <- 10 
else { if (habDegrees <= 324) HueLetterCode <- 9 else HueLetterCode <- 8
}
}
}
}
}
}
}
}
}
# Each letter code is prefixed by a number greater than 0, and less than
# or equal to 10, that further specifies the hue.
HueNumber <- approx(c(0, 36), c(0, 10), habDegrees %% 36)$y
if (HueNumber == 0) HueNumber = 10
# Munsell value can be approximated very accurately with a simple division by 10
Value         = L/10
# This Munsell chroma expression is a very rough approximation, but the best available
Chroma        = Cab/5
# Assemble individual Munsell coordinates into one Munsell specification
ColorLabMunsellVector = c(HueNumber, Value, Chroma, HueLetterCode)
MunsellSpec           = ColorLabFormatToMunsellSpec(ColorLabMunsellVector)$MunsellSpecString
MunsellSpec
}

MaxChromaForExtrapolatedRenotation <- function(HueNumPrefix, CLHueLetterIndex, Value){
#For a given hue and value, find the maximum chroma possible, when interpolating over the extrapolated Munsell renotation data.
# Initialize default return values
MaxChroma = -99
Status.ind = -99
if (Value >= 9.99) return(list(MaxChroma = 0, Status.ind = 1)) # Colour is ideal white, which has no chroma
# Bound Munsell value between two values, ValueMinus and ValuePlus, for which Munsell
# reflectance spectra are available.  Currently, the Munsell values for which data is
# available are 1 or greater.
if (Value < 1 ) return(list(Status.ind = 2)) # Set error and return
if ((Value %% 1) == 0) { # input value is already integer
   ValueMinus = Value
   ValuePlus  = Value
} else { # Input value is not integer
   ValueMinus = floor(Value)	
   ValuePlus  = ValueMinus + 1	
}
# Bound the Munsell hue between two canonical hues.  One of the hues is immediately
# clockwise to the input hue in the CIE diagram, and the other is 
# immediately counterclockwise.  Express the canonical hues
# in terms of the hue list in the description, so that they can be used as indices to
# the extrapolated MacAdam matrix.
tmp = BoundingRenotationHues(HueNumPrefix, CLHueLetterIndex)
ClockwiseHue <- tmp$ClockwiseHue
CtrClockwiseHue <- tmp$CtrClockwiseHue
CLHueNumPrefixCW    = ClockwiseHue[1]
CLHueLetterIndexCW  = ClockwiseHue[2]
CWHueIndex          = (4 * ((7-CLHueLetterIndexCW) %% 10)) + (CLHueNumPrefixCW/2.5)
CLHueNumPrefixCCW   = CtrClockwiseHue[1]
CLHueLetterIndexCCW = CtrClockwiseHue[2]
CCWHueIndex         = (4 * ((7-CLHueLetterIndexCCW) %% 10)) + (CLHueNumPrefixCCW/2.5)
# Combine the two value bounds and two hue bounds to produce four hue-value 
# combinations.  Find the MacAdam limit for each of these combinations.

MaxExtRenChromaMatrix <-matrix(c(14,16,20,24,28,28,24,20,14,12,16,20,24,28,30,24,18,12,12,16,20,24,28,30,26,20,14,12,16,
20,24,26,30,26,26,16,12,16,16,18,20,22,24,26,18,10,10,12,14,18,20,22,26,26,6,8,10,14,16,
20,22,24,26,4,8,10,12,16,18,20,22,24,4,8,10,12,14,18,20,22,24,4,6,8,12,14,16,20,22,24,4,
6,8,12,14,16,18,20,22,4,6,8,12,14,16,18,20,22,6,8,10,12,14,16,18,22,22,6,10,12,14,16,18,
20,22,22,8,14,14,18,20,22,24,28,28,10,18,24,28,32,32,32,32,28,10,18,24,30,32,34,34,32,28,
10,18,24,30,32,34,34,32,28,10,18,24,30,32,34,34,30,24,10,18,24,30,32,34,32,28,24,10,18,
24,28,32,32,32,26,22,10,16,22,26,30,30,30,26,20,8,16,20,24,26,26,26,24,18,8,14,18,22,24,
24,24,22,18,8,12,16,20,20,22,22,20,16,10,12,16,18,20,20,20,18,16,10,14,16,18,20,22,20,16,
12,10,14,18,20,22,24,20,14,10,14,18,22,24,24,24,18,14,10,38,38,38,38,32,26,18,14,10,38,38,
38,38,34,26,20,14,10,38,38,38,38,36,28,22,16,10,38,38,38,38,38,32,24,18,10,32,36,38,38,38,
36,28,22,14,26,30,38,38,38,36,30,26,16,24,28,34,38,38,36,30,26,16,22,26,30,34,36,36,30,26,
16,20,22,26,30,32,32,28,24,16,16,20,22,26,30,30,26,22,14,16,18,20,24,28,28,24,20,14),40,9,byrow=TRUE)

MALimitVMCW  = MaxExtRenChromaMatrix[CWHueIndex,  ValueMinus]
MALimitVMCCW = MaxExtRenChromaMatrix[CCWHueIndex, ValueMinus]
if (ValuePlus <= 9) { # Handle values between 9 and 10 separately
   MALimitVPCW  = MaxExtRenChromaMatrix[CWHueIndex,  ValuePlus]
   MALimitVPCCW = MaxExtRenChromaMatrix[CCWHueIndex, ValuePlus]
   # The maximum chroma at which an interpolating function can be evaluated is the
   # minimum of the maximum chromas for which all neighboring canonical points
   # can be evaluated.
   MaxChroma = min(MALimitVMCW, MALimitVMCCW, MALimitVPCW, MALimitVPCCW)
} else { # Input value is between 9 and 10; use geometric calculation
   FactorList        = MunsellValueToLuminanceFactor(9)	
   LuminanceFactor9  = FactorList$ASTMD153508
   FactorList        = MunsellValueToLuminanceFactor(10)
   LuminanceFactor10 = FactorList$ASTMD153508
   FactorList        = MunsellValueToLuminanceFactor(Value)
   LuminanceFactorV  = FactorList$ASTMD153508
   MaxCWchroma  = approx(c(LuminanceFactor9, LuminanceFactor10), c(MALimitVMCW, 0),  LuminanceFactorV)$y
   MaxCCWchroma = approx(c(LuminanceFactor9, LuminanceFactor10), c(MALimitVMCCW, 0), LuminanceFactorV)$y
   MaxChroma    = min(MaxCWchroma, MaxCCWchroma)
}
# Set success code and return
list(Status.ind = 1,MaxChroma=MaxChroma )
}


LuminanceFactorToMunsellValue <- function(LumFac){
# Convert the luminance factor of an object colour, into a Munsell value
AllMunsVals <- seq(0,10,0.1)
AllLumFaccs = c(0.0000, 0.1201, 0.2370, 0.3521, 0.4666, 0.5818, 0.6990, 0.8193, 0.9439, 1.0738, 1.2101, 
1.3538, 1.5059, 1.6672, 1.8387, 2.0212, 2.2156, 2.4225, 2.6428, 2.8771, 3.1262, 3.3906, 3.6711, 3.9681,
 4.2822, 4.6141, 4.9641, 5.3327, 5.7204, 6.1277, 6.5550, 7.0025, 7.4708, 7.9601, 8.4709,9.0033,9.5577,
10.1344,10.7336,11.3557,12.0007,12.6691,13.3610,14.0766,14.8161,15.5797,16.3676,17.1801,18.0172,18.8791,
19.7661,20.6783,21.6159,22.5791,23.5680,24.5829,25.6239,26.6912,27.7850,28.9055,30.0529,31.2275,32.4294,
33.6590,34.9165,36.2021,37.5161,38.8589,40.2307,41.6319,43.0628,44.5238,46.0153,47.5377,49.0914,50.6769,
52.2946,53.9451,55.6289,57.3466,59.0986,60.8858,62.7086,64.5678,66.4641,68.3982,70.3710,72.3831,74.4355,
76.5291,78.6647,80.8434,83.0661,85.3339,87.6479,90.0092,92.4190,94.8784,97.3889,99.9516,102.5680)
NewhallMunsellValue = approx(AllLumFaccs, AllMunsVals, LumFac)$y
# Calculate lightness (L*) from CIELAB model ([Fairchild2005, Sect. 10.3]),
# and divide by 10 to get the Munsell value.
LumFacAsFrac = LumFac/100
if (LumFacAsFrac > 0.008856) fofomega = LumFacAsFrac^(1/3) else fofomega = 7.787*LumFacAsFrac + (16/116)
Lstar = 116*fofomega - 16
CIELABMunsellValue = Lstar/10
# Calculate Munsell Value using the closed-form approximation from [McCamy1992].
# This approximation has been incorporated in [ASTMD1535-08]
Y = LumFac
if (Y <= 0.9) McCamyMunsellValue = 0.87445 * (Y^0.9967) else McCamyMunsellValue = 2.49268*(Y^(1/3)) - 1.5614  - (0.985/(((0.1073*Y-3.084)^2) + 7.54)) +
(0.0133/(Y^2.3)) + 0.0084*sin(4.1*(Y^(1/3))+1) + (0.0221/Y)*sin(0.39*(Y-2)) - (0.0037/(0.44*Y)) * sin(1.28*(Y-0.53))
AllMunsVals = seq(0,10,0.02)
AllASTMLuminanceFactors <- c(0.000000,0.023740,0.047310,0.070723,0.093989,0.117118,0.140123,0.163012,0.185799,
0.208492,0.231102,0.253641,0.276118,0.298543,0.320928,0.343281,0.365614,0.387936,
0.410257,0.432587,0.454936,0.477314,0.499730,0.522194,0.544715,0.567303,0.589967,
0.612717,0.635561,0.658509,0.681571,0.704754,0.728068,0.751522,0.775125,0.798885,
0.822812,0.846913,0.871197,0.895673,0.920349,0.945234,0.970336,0.995662,1.021222,
1.047023,1.073073,1.099381,1.125954,1.152799,1.179925,1.207340,1.235051,1.263065,
1.291391,1.320035,1.349005,1.378308,1.407952,1.437944,1.468291,1.498999,1.530077,
1.561531,1.593367,1.625594,1.658217,1.691243,1.724679,1.758532,1.792808,1.827513,
1.862655,1.898239,1.934272,1.970759,2.007709,2.045125,2.083016,2.121386,2.160241,
2.199588,2.239433,2.279781,2.320638,2.362010,2.403902,2.446321,2.489271,2.532759,
2.576790,2.621368,2.666500,2.712192,2.758447,2.805272,2.852671,2.900650,2.949214,
2.998368,3.048116,3.098465,3.149418,3.200980,3.253157,3.305953,3.359373,3.413421,
3.468102,3.523421,3.579382,3.635990,3.693248,3.751163,3.809737,3.868975,3.928881,
3.989460,4.050716,4.112653,4.175275,4.238586,4.302590,4.367291,4.432693,4.498800,
4.565616,4.633144,4.701389,4.770353,4.840042,4.910458,4.981605,5.053487,5.126107,
5.199468,5.273575,5.348431,5.424038,5.500401,5.577523,5.655406,5.734055,5.813473,
5.893662,5.974626,6.056368,6.138891,6.222198,6.306293,6.391178,6.476856,6.563330,
6.650603,6.738677,6.827557,6.917244,7.007741,7.099052,7.191178,7.284123,7.377889,
7.472478,7.567894,7.664139,7.761215,7.859126,7.957873,8.057459,8.157887,8.259158,
8.361276,8.464242,8.568059,8.672730,8.778256,8.884640,8.991885,9.099991,9.208963,
9.318801,9.429508,9.541086,9.653537,9.766863,9.881067,9.996150,10.112115,10.228963,
10.346696,10.465317,10.584827,10.705228,10.826523,10.948712,11.071799,11.195785,
11.320671,11.446459,11.573152,11.700751,11.829258,11.958675,12.089003,12.220244,
12.352400,12.485472,12.619463,12.754374,12.890206,13.026961,13.164642,13.303248,
13.442783,13.583248,13.724643,13.866972,14.010235,14.154434,14.299571,14.445646,
14.592662,14.740620,14.889522,15.039369,15.190162,15.341904,15.494595,15.648237,
15.802831,15.958379,16.114883,16.272343,16.430761,16.590140,16.750479,16.911780,
17.074046,17.237277,17.401474,17.566640,17.732774,17.899880,18.067958,18.237010,
18.407037,18.578040,18.750021,18.922981,19.096921,19.271844,19.447749,19.624640,
19.802516,19.981380,20.161233,20.342076,20.523910,20.706737,20.890559,21.075377,
21.261191,21.448004,21.635817,21.824631,22.014448,22.205269,22.397096,22.589929,
22.783771,22.978623,23.174486,23.371361,23.569251,23.768157,23.968079,24.169020,
24.370981,24.573964,24.777969,24.982999,25.189055,25.396139,25.604251,25.813394,
26.023570,26.234779,26.447023,26.660305,26.874624,27.089984,27.306386,27.523831,
27.742321,27.961858,28.182443,28.404078,28.626765,28.850505,29.075300,29.301153,
29.528064,29.756035,29.985069,30.215167,30.446331,30.678562,30.911863,31.146236,
31.381682,31.618204,31.855802,32.094480,32.334239,32.575081,32.817008,33.060023,
33.304126,33.549321,33.795610,34.042994,34.291476,34.541057,34.791741,35.043528,
35.296423,35.550425,35.805539,36.061766,36.319109,36.577570,36.837151,37.097854,
37.359683,37.622640,37.886726,38.151945,38.418299,38.685791,38.954422,39.224197,
39.495117,39.767186,40.040405,40.314777,40.590306,40.866994,41.144844,41.423859,
41.704041,41.985394,42.267920,42.551623,42.836505,43.122569,43.409819,43.698257,
43.987888,44.278713,44.570736,44.863961,45.158390,45.454027,45.750875,46.048938,
46.348219,46.648721,46.950448,47.253403,47.557591,47.863014,48.169676,48.477581,
48.786732,49.097134,49.408790,49.721704,50.035879,50.351320,50.668031,50.986015,
51.305276,51.625820,51.947649,52.270767,52.595180,52.920891,53.247904,53.576225,
53.905856,54.236803,54.569070,54.902662,55.237582,55.573836,55.911428,56.250364,
56.590646,56.932281,57.275273,57.619628,57.965348,58.312441,58.660911,59.010762,
59.362000,59.714631,60.068658,60.424088,60.780926,61.139177,61.498846,61.859939,
62.222462,62.586419,62.951817,63.318661,63.686957,64.056710,64.427926,64.800612,
65.174773,65.550415,65.927544,66.306166,66.686287,67.067913,67.451051,67.835707,
68.221887,68.609598,68.998846,69.389637,69.781979,70.175877,70.571338,70.968370,
71.366978,71.767171,72.168954,72.572334,72.977319,73.383917,73.792132,74.201974,
74.613450,75.026566,75.441330,75.857750,76.275833,76.695586,77.117018,77.540136,
77.964947,78.391460,78.819682,79.249622,79.681288,80.114687,80.549827,80.986718,
81.425366,81.865781,82.307971,82.751945,83.197711,83.645277,84.094652,84.545845,
84.998866,85.453722,85.910422,86.368977,86.829394,87.291683,87.755853,88.221914,
88.689875,89.159745,89.631535,90.105252,90.580908,91.058511,91.538072,92.019601,
92.503107,92.988601,93.476093,93.965592,94.457109,94.950655,95.446239,95.943873,
96.443567,96.945331,97.449176,97.955113,98.463153,98.973307,99.485586,100.000000)
ASTMTableLookupMunsellValue = approx(AllASTMLuminanceFactors, AllMunsVals, LumFac)$y
list(OriginalLuminanceFactor = LumFac, Newhall1943 = NewhallMunsellValue, CIELAB = CIELABMunsellValue, 
McCamy1992 = McCamyMunsellValue, ASTMTableLookup = ASTMTableLookupMunsellValue)
}

ChromDiagHueAngleToMunsellHue <- function(ChromDiagHueAngle){
# Munsell hue as a single number
SingleHueNumber = approx(c(0,45,70,135,160,225,255,315,360), c(0,2,3,4,5,6,8,9,10), ChromDiagHueAngle)$y
# Now the single hue number is converted back to a Munsell hue specification.
if (SingleHueNumber <= 0.5) HueLetterCode = 7
else { if (SingleHueNumber <= 1.5) HueLetterCode = 6
else { if (SingleHueNumber <= 2.5) HueLetterCode = 5
else { if (SingleHueNumber <= 3.5) HueLetterCode = 4
else { if (SingleHueNumber <= 4.5) HueLetterCode = 3
else { if (SingleHueNumber <= 5.5) HueLetterCode = 2
else { if (SingleHueNumber <= 6.5) HueLetterCode = 1
else { if (SingleHueNumber <= 7.5) HueLetterCode = 10
else { if (SingleHueNumber <= 8.5) HueLetterCode = 9
else { if (SingleHueNumber <= 9.5) HueLetterCode = 8 else HueLetterCode = 7
}
}
}
}
}
}
}
}
}
HueNumber = (10*(SingleHueNumber %% 1) + 5) %% 10
if (HueNumber == 0) HueNumber = 10
list(HueNumber=HueNumber,HueLetterCode=HueLetterCode)
}


















getMunsellHueLetterDesignator<-function(munsellSTR)
{# Convert a Munsell specification into ColorLab format specification.
# Based on: MunsellAndKubelkaMunkToolbox by Paul Centore 
if (any(grepl('BG',munsellSTR))) return(2)
if (any(grepl('GY',munsellSTR))) return(4)
if (any(grepl('YR',munsellSTR))) return(6)
if (any(grepl('RP',munsellSTR))) return(8)
if (any(grepl('PB',munsellSTR))) return(10)
if (any(grepl('B',munsellSTR))) return(1)
if (any(grepl('G',munsellSTR))) return(3)
if (any(grepl('Y',munsellSTR))) return(5)
if (any(grepl('R',munsellSTR))) return(7)
if (any(grepl('P',munsellSTR))) return(9)
if (any(grepl('N',munsellSTR))) return(0)
-1
}

ColorLabFormatToMunsellSpec <- function(ColorLabMunsellVector=NA,HueDecimals=NA,ValueDecimals=NA,ChromaDecimals=NA)
{# Convert a ColorLab format specification into a Munsell specification.
# Based on: MunsellAndKubelkaMunkToolbox by Paul Centore 
if (any(is.na(ColorLabMunsellVector))) stop('ColorLabMunsellVector must be a numeric vector')
if (any(!is.numeric(ColorLabMunsellVector))) stop('ColorLabMunsellVector must be a numeric vector')
if ((length(ColorLabMunsellVector)<1) | (length(ColorLabMunsellVector)>4)) stop('ColorLabMunsellVector must be a numeric vector with 4 elements')
if (is.na(HueDecimals)) HueDecimals<-2
if (is.na(ValueDecimals)) ValueDecimals<-2
if (is.na(ChromaDecimals)) ChromaDecimals<-2
if (length(ColorLabMunsellVector)==1){ #Achromatic colour of form [V]
if (ValueDecimals == 0) MunsellSpecString <- sprintf('N%0.0f', ColorLabMunsellVector) else MunsellSpecString <- sprintf(paste('N%',as.character(ValueDecimals+2),
'.',as.character(ValueDecimals),'f',sep=''), ColorLabMunsellVector)
HueString <- 'N'
} else {
if (ColorLabMunsellVector[3] == 0) { #Achromatic colour, of form [H1 V 0 H2]
if (ValueDecimals == 0) MunsellSpecString <- sprintf('N%0.0f', ColorLabMunsellVector[2]) else MunsellSpecString <- sprintf(paste('N%',as.character(ValueDecimals+2),
'.',as.character(ValueDecimals),'f',sep=''), ColorLabMunsellVector[2])
HueString <- 'N'
} else { # Construct the complete output string as a concatenation of a hue string, a value string, and a chroma string.
ColourLetters <- c('B', 'BG', 'G', 'GY', 'Y', 'YR', 'R', 'RP', 'P', 'PB')
   if (HueDecimals == 0){ # Avoid unwanted spaces in final output
       HueString <- sprintf('%0.0f', ColorLabMunsellVector[1])      
   } else {
       HueString <- sprintf(paste('%',as.character(HueDecimals+2),'.',as.character(HueDecimals),'f',sep=''), ColorLabMunsellVector[1])
   }
   HueString <- paste(HueString,ColourLetters[ColorLabMunsellVector[4]],sep='')
   if (ValueDecimals == 0){# Avoid unwanted spaces in final output
       ValueString <- sprintf('%0.0f', ColorLabMunsellVector[2])
  } else {
       ValueString  <- sprintf(paste('%',as.character(ValueDecimals+2),'.',as.character(ValueDecimals),'f',sep=''), ColorLabMunsellVector[2])
   }
   if (ChromaDecimals == 0) # Avoid unwanted spaces in final output
       ChromaString <- sprintf('%0.0f', ColorLabMunsellVector[3])
   else
       ChromaString <- sprintf(paste('%',as.character(ChromaDecimals+2),'.',as.character(ChromaDecimals),'f',sep=''), ColorLabMunsellVector[3])
   }
}
if (any(HueString != 'N')) MunsellSpecString <- paste(HueString,' ',ValueString,'/',ChromaString,sep='')
list(MunsellSpecString=MunsellSpecString, HueString=HueString)
}

MunsellSpecToColorLabFormat <- function(MunsellSpecString)
{ # Convert a Munsell specification into ColorLab format specification.
# Based on: MunsellAndKubelkaMunkToolbox by Paul Centore 
# Remove all spaces from input Munsell string
MunsellString <- gsub(' ','',MunsellSpecString)
# Make all letters in Munsell string upper case
MunsellString <- toupper(MunsellString)
# Read through the Munsell string in order, extracting hue, value, and chroma
ctr <- 1
entry <- substr(MunsellString,ctr,ctr)
# Check for achromatic colour
if (entry=='N'){
   MunsellValue <- as.character(substr(MunsellString,2,nchar(MunsellString)))
   ColorLabMunsellVector <- MunsellValue
   return(ColorLabMunsellVector)#list(MunsellSpecString=MunsellString, HueString=HueString)
}
# Extract hue number from start of Munsell specification.
HueLetterReached <- FALSE
while (grepl('[0-9\\.]',entry)){
   ctr <- ctr + 1
   entry <- substr(MunsellString,ctr,ctr)
}
MunsellHueNumber <- as.numeric(substr(MunsellString,1,(ctr-1)))
# Next, extract hue letter designator 
# entry is 1st letter of Munsell Hue Designator, which might consist of two letters
PossibleNextLetter <- substr(MunsellString,ctr+1,ctr+1)
ColorLabHueLetterDesignator <- -99# Default error value
if (!grepl('[0-9\\.]+',PossibleNextLetter)){
   MunsellHueLetterDesignator <- substr(MunsellString,ctr,(ctr+1))
   if (MunsellHueLetterDesignator=='BG') ColorLabHueLetterDesignator <- 2 else {
   if (MunsellHueLetterDesignator=='GY') ColorLabHueLetterDesignator <- 4 else {
   if (MunsellHueLetterDesignator=='YR') ColorLabHueLetterDesignator <- 6 else {
   if (MunsellHueLetterDesignator=='RP') ColorLabHueLetterDesignator <- 8 else {
   if (MunsellHueLetterDesignator=='PB') ColorLabHueLetterDesignator <- 10 else stop(paste('ERROR1: ',MunsellHueLetterDesignator,' is not a valid hue letter designator'))
   }
   }
   }
   }
   ctr <- ctr + 1
} else {
   MunsellHueLetterDesignator <- substr(MunsellString,ctr,ctr)
   if (MunsellHueLetterDesignator=='B') ColorLabHueLetterDesignator <- 1 else {
   if (MunsellHueLetterDesignator=='G') ColorLabHueLetterDesignator <- 3 else {
   if (MunsellHueLetterDesignator=='Y') ColorLabHueLetterDesignator <- 5 else {
   if (MunsellHueLetterDesignator=='R') ColorLabHueLetterDesignator <- 7 else {
   if (MunsellHueLetterDesignator=='P') ColorLabHueLetterDesignator <- 9 else stop(paste('ERROR2: ',MunsellHueLetterDesignator,' is not a valid hue letter designator'))
   }
   }
   }
   }   
}
if (MunsellHueNumber == 0){
   MunsellHueNumber <- 10
   MunsellHueLetterDesignator <- (MunsellHueLetterDesignator + 1) %% 10
}
# Remainder of Munsell specification string is value and chroma
vec              <- unlist(strsplit(substr(MunsellString,(ctr+1),nchar(MunsellString)), "/", fixed = TRUE))  #sscanf(MunsellString((ctr+1):end), '%f/%f')
MunsellValue     <- vec[1]
MunsellChroma    <- vec[2]
if (MunsellChroma == 0) ColorLabMunsellVector <- MunsellValue else ColorLabMunsellVector <- c(MunsellHueNumber, MunsellValue, MunsellChroma, ColorLabHueLetterDesignator)
as.numeric(ColorLabMunsellVector)
}

FindHueOnRenotationOvoid <- function(MunsellSpec)
{
# The input could be either a Munsell string, such as 4.2R8.1/5.3, or a Munsell vector in ColorLab format
if (is.character(MunsellSpec)) ColorLabMunsellVector <- MunsellSpecToColorLabFormat(MunsellSpec) else ColorLabMunsellVector <- MunsellSpec
# Extract hue, chroma, and value from ColorLab Munsell vector that corresponds to input.
HueNumber <- Value <- Chroma <- HueLetterCode <- 0
if (length(ColorLabMunsellVector) == 1){ # Colour is Munsell grey
# A one-element vector is an achromatic grey, so no interpolation is necessary.  Evaluate directly and return.
   x <- 0.310062 #whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'x2' ]
   y <- 0.316148 #whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'y2' ]
   return(list(x=x, y=y, Status.ind = 1))
} else {
   HueNumber     <- ColorLabMunsellVector[1]
   Value         <- ColorLabMunsellVector[2]
   Chroma        <- ColorLabMunsellVector[3]
   HueLetterCode <- ColorLabMunsellVector[4]
}
# Check that the Munsell value of the input is an integer between 1 and 9
# (If the value is 10, then the colour is ideal white, which was previously assigned the xy coordinates for Illuminant C.)
if (Value < 1 | Value > 9) return(list(Status.ind = 2))# Set error and return
# For numerical convenience, allow Munsell values very close to integers, and round them to integers.
if (abs(Value-round(Value)) > 0.001) return(list(Status.ind = 2))# Set error and return
# Round value to integer, if it is already very close to an integer.
Value <- round(Value)
# Check that the chroma of the input is a positive, even integer.
if (Chroma < 2) return(list(Status.ind = 3))# Set error and return
# For numerical convenience, allow Munsell chromas very close to even integers, and round them to even integers.
if (abs(2*((Chroma/2)-round(Chroma/2))) > 0.001) return(list(Status.ind = 3))# Set error and return
#  Round chroma to positive even integer, if it is already very near a positive even integer
Chroma <- 2*round(Chroma/2)
# Check to see if the input colour is a standard Munsell colour, for which renotation data is available without interpolation.
# f so, make the renotation conversion and return
threshold <- 0.001
if (abs(HueNumber) < threshold | abs(HueNumber-2.5) < threshold | abs(HueNumber-5) < threshold | abs(HueNumber-7.5) < threshold | abs(HueNumber-10) < threshold){
  HueNumber <- 2.5 * round(HueNumber/2.5) # Round to very close standard hue
  tmp <- MunsellToxyYfromExtrapolatedRenotation(ColorLabMunsellVector)
  Status.ind <- tmp$Status.ind
  x <- tmp$x
  y <- tmp$y
  Y <- tmp$Y  
  if (Status.ind != 1) return(list(Status.ind = 5)) else return(list(x=x,y=y,Y=Y,Status.ind = 1))
}
# Find two hues which bound the hue of the input colour, and for which renotation data is available.  Renotation data is available only for hues whose prefix number is 2.5, 5.0, 7.5, or 10.0.
tmp <- BoundingRenotationHues(HueNumber, HueLetterCode)
ClockwiseHue <- tmp$ClockwiseHue
CtrClockwiseHue <- tmp$CtrClockwiseHue
MunsellHueNumberMinus <- ClockwiseHue[1]
CLHueLetterIndexMinus <- ClockwiseHue[2]
MunsellHueNumberPlus <- CtrClockwiseHue[1]
CLHueLetterIndexPlus <- CtrClockwiseHue[2]
# Express the two bounding Munsell colours in ColorLab Munsell vector format, and in ASTM format
CLMVMinus <- c(ClockwiseHue[1], Value, Chroma, ClockwiseHue[2])
CLMVPlus <- c(CtrClockwiseHue[1], Value, Chroma, CtrClockwiseHue[2])
# Find the achromatic point, to be used as the center of polar coordinates
tmp <- MunsellToxyYfromExtrapolatedRenotation(Value)
xGrey <- tmp$x
yGrey <- tmp$y
YGrey <- tmp$Y
Status.ind <- tmp$Status.ind
if (Status.ind != 1) return(list(Status.ind = 5))# Set error and return
# Express the two bounding Munsell colours in polar coordinates, after finding their chromaticity coordinates in the Munsell renotation.  
tmp <- MunsellToxyYfromExtrapolatedRenotation(CLMVPlus)
  Status.ind <- tmp$Status.ind
  xPlus <- tmp$x
  yPlus <- tmp$y
  YPlus <- tmp$Y
if (Status.ind != 1) return(list(Status.ind = 2))# No renotation data available for bounding colour. Set error code and return
THPlus <- atan2(yPlus-yGrey, xPlus-xGrey)
#RPlus <- sqrt(abs(xPlus-xGrey)^2 + abs(yPlus-yGrey)^2)
RPlus <- hypot((xPlus-xGrey), (yPlus-yGrey))
THPlus          <- ((180/pi) * THPlus) %% 360 # Convert to degrees
tmp <- MunsellToxyYfromExtrapolatedRenotation(CLMVMinus)
  Status.ind <- tmp$Status.ind
  xMinus <- tmp$x
  yMinus <- tmp$y
  YPMinus <- tmp$Y
if (Status.ind != 1) return(list(Status.ind = 2)) # No renotation data available for bounding colour.  Set error code and return
THMinus <- atan2(yMinus-yGrey, xMinus-xGrey)
RMinus <- sqrt(abs(xMinus-xGrey)^2 + abs(yMinus-yGrey)^2)
RMinus <- hypot((xMinus-xGrey), (yMinus-yGrey))
THMinus <- ((180/pi) * THMinus) %% 360 # Convert to degrees
LowerTempHueAngle <- MunsellHueToChromDiagHueAngle(MunsellHueNumberMinus, CLHueLetterIndexMinus) ;
TempHueAngle      <- MunsellHueToChromDiagHueAngle(HueNumber,             HueLetterCode)	
UpperTempHueAngle <- MunsellHueToChromDiagHueAngle(MunsellHueNumberPlus,  CLHueLetterIndexPlus)  ;
# Adjust for possible wraparound.  There should be a short arc running counter clockwise from
# the lower hue value to the upper hue value
if (THMinus - THPlus > 180 ) THPlus  <- THPlus + 360
if (LowerTempHueAngle == 0) LowerTempHueAngle <- 360

if (LowerTempHueAngle > UpperTempHueAngle) { # E.g. Lower is 355, Upper is 10
   if (LowerTempHueAngle > TempHueAngle) LowerTempHueAngle <- LowerTempHueAngle - 360 else {
      LowerTempHueAngle <- LowerTempHueAngle - 360
      TempHueAngle      <- TempHueAngle      - 360
   }
}
# The interpolation of the input colour using the two bounding colours will be done 
# by either linear or radial interpolation, as determined by a function call.
tmp <- LinearVsRadialInterpOnRenotationOvoid(MunsellSpec)
InterpStyle.Linear <- tmp$InterpStyle.Linear
InterpStyle.Radial <- tmp$InterpStyle.Radial
InterpStyle.Input <- tmp$InterpStyle.Input
InterpStyle.OnGrid <- tmp$InterpStyle.OnGrid
Status.ind <- tmp$Status.ind
if (Status.ind != 1) return(list(Status.ind = 4)) #Unsuccessful call; return with error message
# Perform interpolation as indicated
if (InterpStyle.Linear == TRUE) { # Use linear interpolation
   x = approx(c(LowerTempHueAngle,UpperTempHueAngle), c(xMinus, xPlus), TempHueAngle)$y
   y = approx(c(LowerTempHueAngle,UpperTempHueAngle), c(yMinus, yPlus), TempHueAngle)$y
   } else { if (InterpStyle.Radial == TRUE) { # Use radial interpolation
   # Interpolate radially along the chroma ovoid. For example, if the input colour is 4B6/7, then 
   # the two bounding points on this ovoid are 2.5B6/6 and 5B6/6.  The new point on the
   # ovoid will have hue angle 60% of the way from the hue angle of 2.5B, to the
   # hue angle at 5B.  The R value of the new point will be between the R values for 2.5B6/6
   # and 5B6/6, in a 60/40 ratio.  
   InterpolatedTheta <- approx(c(LowerTempHueAngle,UpperTempHueAngle), c(THMinus, THPlus), TempHueAngle)$y
   InterpolatedR     <- approx(c(LowerTempHueAngle,UpperTempHueAngle), c(RMinus,  RPlus),  TempHueAngle)$y
   # Find xy chromaticity coordinates for the new point on the chroma ovoid
   x  <- InterpolatedR * cos(InterpolatedTheta/180*pi) + xGrey
   y  <- InterpolatedR * sin(InterpolatedTheta/180*pi) + yGrey
   } else return(list(Status.ind = 4)) # Interpolation style not determined; return with error message
}
   return(list(x=x, y=y, Status.ind = 1))
   }

#ColorLabColour=c(2.5, 5.0, 4.0, 6.0)
#ColorLabColour=c(10,5,4,7)
#ColorLabColour=5
MunsellToxyYfromExtrapolatedRenotation<-function(ColorLabColour){
HueListLetters <- c('R', 'YR', 'Y', 'GY', 'G', 'BG', 'B', 'PB', 'P', 'RP')
ColorLabHueLetters <- c('B', 'BG', 'G', 'GY', 'Y', 'YR', 'R', 'RP', 'P', 'PB')
data("MunsellRenotation", envir = environment())
MunsellRenotation<-get("MunsellRenotation", envir  = environment())
# A one-element colour vector is achromatic grey
if (length(ColorLabColour) == 1){
V  <- ColorLabColour[1]
C <- 0
} else	{# ColorLab vector has four elements
NumericalHuePrefix <- ColorLabColour[1]
HueLetterDesignator <- ColorLabColour[4]
V <- ColorLabColour[2]
C <- ColorLabColour[3]
ColorLabHueLetter <- ColorLabHueLetters[HueLetterDesignator]
}
# Check that value is an integer between 1 and 9.  If not, return with an error message.
if ((V < 1) || (V >9) || ((V %% 1) != 0)) return(list(Status.ind = 3)) else ValueIndex <- V
# Check for a valid chroma input
if (C != 0)   if (((C %% 2) != 0) || (C < 2) || (C > 38)){ # Input chroma is not a multiple of 2
      return(list(Status.ind = 5))
   }
if (length(ColorLabColour) == 1){


  w <- which(MunsellRenotation[["H"]]=='N' & MunsellRenotation[["V"]]==V & MunsellRenotation[["C"]]==C)
  if (length(w)==0) return(list(Status.ind = 2))
  m<-MunsellRenotation[w,] 
  } else {
    w <- which(MunsellRenotation[["H"]]==paste(ColorLabColour[1],ColorLabHueLetter,sep='') & MunsellRenotation[["V"]]==V & MunsellRenotation[["C"]]==C)
    if (length(w)==0) return(list(Status.ind = 2))
    m<-MunsellRenotation[w,]
  }
if (length(m$x)==0) return(list(Status.ind = 2))
return(list(x=m$x,y=m$y,Y=m$Y,Status.ind = 1))
}

#MunsellRenotation[which(MunsellRenotation[["x"]]==0.3879  & MunsellRenotation[["y"]]==0.3398 & MunsellRenotation[["Y"]]==19.7700),]

#ColorLabMunsellVector=c(7.6, 8.9, 2.2, 9.0)
MunsellToxyForIntegerMunsellValue<-function(ColorLabMunsellVector){
MunsellSpecString <- ColorLabFormatToMunsellSpec(ColorLabMunsellVector)$MunsellSpecString
if (length(ColorLabMunsellVector) == 1){
   if (ColorLabMunsellVector[1] == 10){ #Ideal white; set value directly and return
   x <- 0.310062 #whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'x2' ]
   y <- 0.316148 #whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'y2' ]
   return(list(x=x,y=y,Status.ind = TRUE))
   }
   tmp <- MunsellToxyYfromExtrapolatedRenotation(ColorLabMunsellVector)
   Status.ind <- tmp$Status.ind
   x <- tmp$x
   y <- tmp$y
   Y <- tmp$Y
   if (Status.ind == 1) return(list(x=x,y=y,Y=Y,Status.ind = 1)) else return(list(Status.ind = 2))   
}
# Otherwise, ColorLabMunsellVector has four elements with Munsell information
# Extract data from ColorLab version of Munsell specification
MunsellHueNumber    <- ColorLabMunsellVector[1]
IntegerMunsellValue <- ColorLabMunsellVector[2]
MunsellChroma       <- ColorLabMunsellVector[3]
CLHueLetterIndex    <- ColorLabMunsellVector[4]
# Find two chromas which bound the chroma of the input colour, and for which
# renotation data are available.  Renotation data is available only 
# for even chromas.
if (MunsellChroma == 0){ # Munsell grey, for which no interpolation is needed
   tmp <- MunsellToxyYfromExtrapolatedRenotation(IntegerMunsellValue)
   Status.ind <- tmp$Status.ind
   x <- tmp$x
   y <- tmp$y
   Y <- tmp$Y
   if (Status.ind == 1) Status.ind <- 1 else Status.ind <- 2
   return(list(Status.ind = Status.ind))
} else {
   if ((MunsellChroma %% 2) == 0){ # No interpolation needed for Munsell chroma
      MunsellChromaMinus <- MunsellChroma
      MunsellChromaPlus  <- MunsellChroma
   } else {
      MunsellChromaMinus <- 2 * floor(MunsellChroma/2)
      MunsellChromaPlus  <- MunsellChromaMinus + 2
   }
}
if (MunsellChromaMinus == 0) { # Colour within smallest even chroma ovoid
# Smallest chroma ovoid collapses to white point
   xMinus <- 0.310062 # whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'x2' ]
   yMinus <- 0.316148 #whitepointsilluminants[which(whitepointsilluminants[,'illuminant']=='C'), 'y2' ]
} else {
   tmp <- FindHueOnRenotationOvoid(c(MunsellHueNumber,IntegerMunsellValue,MunsellChromaMinus,CLHueLetterIndex))
   xMinus <- tmp$x
   yMinus <- tmp$y
   Status.ind <- tmp$Status.ind
   if (Status.ind != 1)  return(list(Status.ind = 3)) # Unsuccessful call - Return with error code
}
tmp <- FindHueOnRenotationOvoid(c(MunsellHueNumber,IntegerMunsellValue, MunsellChromaPlus, CLHueLetterIndex))
xPlus <- tmp$x
yPlus <- tmp$y
Status.ind <- tmp$Status.ind
if (Status.ind != 1) return(list(Status.ind = 3))# Unsuccessful call - Return with error code
if (MunsellChromaMinus == MunsellChromaPlus){ # Two bounding points are the same
   x <- xMinus # xMinus and xPlus are the same, as are yMinus and yPlus
   y <- yMinus
} else { # Two bounding points are different
   x <- approx(c(MunsellChromaMinus, MunsellChromaPlus), c(xMinus, xPlus), MunsellChroma)$y
   y <- approx(c(MunsellChromaMinus, MunsellChromaPlus), c(yMinus, yPlus), MunsellChroma)$y
}
# Set successful status return code and return
list(x=x, y=y, Status.ind=1)
}

CLHueNumPrefix=.86; CLHueLetterIndex=8
BoundingRenotationHues <- function(CLHueNumPrefix, CLHueLetterIndex)
{
if ((CLHueNumPrefix %% 2.5) == 0) { # No interpolation needed for Munsell hue
   if (CLHueNumPrefix == 0) { # Convert to (0,10] from [0,10)
      CLHueNumPrefixCW <- 10
      CLHueLetterIndexCW <- (CLHueLetterIndex+1) %% 10
   } else {
      CLHueNumPrefixCW    <- CLHueNumPrefix
      CLHueLetterIndexCW  <- CLHueLetterIndex
   }
   CLHueNumPrefixCCW   <- CLHueNumPrefixCW
   CLHueLetterIndexCCW <- CLHueLetterIndexCW
} else { # Interpolate between nearest hues
   CLHueNumPrefixCW    <- 2.5 * (floor( CLHueNumPrefix/2.5))
   CLHueNumPrefixCCW   <- (CLHueNumPrefixCW + 2.5) %% 10
   if (CLHueNumPrefixCCW == 0) CLHueNumPrefixCCW <- 10 # Put in (0 10], not [0 10)
   CLHueLetterIndexCCW <- CLHueLetterIndex
   # Switch hue letter of clockwise point, if necessary.  E.g., go from 0R to 10RP
   if (CLHueNumPrefixCW == 0){
      CLHueNumPrefixCW <- 10
      CLHueLetterIndexCW <- (CLHueLetterIndex+1) %% 10
	  # Wraparound from B to PB, if needed.
	  if (CLHueLetterIndexCW == 0) CLHueLetterIndexCW <- 10 #  Put in (0 10], not [0 10)
  } else CLHueLetterIndexCW <- CLHueLetterIndex # if CLHueNumPrefixCW == 0
   CLHueLetterIndexCCW <- CLHueLetterIndex
}
ClockwiseHue    <- c(CLHueNumPrefixCW, CLHueLetterIndexCW)
CtrClockwiseHue <- c(CLHueNumPrefixCCW, CLHueLetterIndexCCW)
list(ClockwiseHue=ClockwiseHue, CtrClockwiseHue=CtrClockwiseHue)
}

#HueNumber=10;HueLetterCode=7
MunsellHueToChromDiagHueAngle <- function(HueNumber,HueLetterCode)
{
SingleHueNumber <- ((((17-HueLetterCode) %% 10) + (HueNumber/10) - 0.5) %% 10)
ChromDiagHueAngle <- approx(c(0,2,3,4,5,6,8,9,10), c(0,45,70,135,160,225,255,315,360), SingleHueNumber)$y
ChromDiagHueAngle
}

#MunsellSpec=c(7.6, 8.9, 2.0, 9.0)
LinearVsRadialInterpOnRenotationOvoid <- function(MunsellSpec){
# Initialize output variables
InterpStyle.Input  <- NULL
InterpStyle.Linear <- FALSE
InterpStyle.Radial <- FALSE
InterpStyle.OnGrid <- FALSE
# The input could be either a Munsell string, such as 4.2R8.1/5.3,
# or a Munsell vector in ColorLab format.  Determine which, and convert to ColorLab format, if needed.
if (is.character(MunsellSpec)) ColorLabMunsellVector <- MunsellSpecToColorLabFormat(MunsellSpec) else ColorLabMunsellVector <- MunsellSpec
# Extract hue, chroma, and value from ColorLab Munsell vector that corresponds to input.
if (length(ColorLabMunsellVector) == 1) { # Colour is Munsell grey
   Value  <- ColorLabMunsellVector[1]
   Chroma <- 0
} else {
   HueNumber     <- ColorLabMunsellVector[1]
   Value         <- ColorLabMunsellVector[2]
   Chroma        <- ColorLabMunsellVector[3]
   HueLetterCode <- ColorLabMunsellVector[4]
}
# Save input in output structure
InterpStyle.Input  <- ColorLabMunsellVector
# No interpolation needed for greys
if (Chroma == 0) {
   InterpStyle.OnGrid <- TRUE
   Status.ind         <- 1 # Return successfully
}
# Check that the Munsell value of the input is an integer between 1 and 10
if (Value < 1 | Value > 10) return(list(Status.ind = 2)) # Set error and return
# For numerical convenience, allow Munsell values very close to integers, and round them to integers.
if (abs(Value-round(Value)) > 0.001) return(list(Status.ind = 2)) # Set error and return
# Round value to integer, if it is already very close to an integer.
Value <- round(Value)
# If Munsell value is 10, then assume that the input colour is an ideal white, regardless of the other entries in the input vector.  
if (Value == 10){
   InterpStyle.OnGrid <- TRUE
   Status.ind         <- 1 # Return successfully
}
# Check that the chroma of the input is a positive, even integer.
if (Chroma < 2)  return(list(Status.ind = 2)) # Set error and return
# For numerical convenience, allow Munsell chromas very close to even integers, and round them to even integers.
if (abs(2*((Chroma/2)-round(Chroma/2))) > 0.001) return(list(Status.ind = 3)) # Set error and return
Chroma <- 2*round(Chroma/2)
# If the hue is already a standard hue, then no interpolation is needed, so return successfully
if ((HueNumber %% 2.5) == 0){
   InterpStyle.OnGrid <- TRUE
   Status.ind         <- 1 # Return successfully
}
# The ASTM hue is a number between 0 and 100, that uses the same increments as
# HueNumber, and is sometimes easier to work with than Munsell hue
ASTMHue          <- MunsellHueToASTMHue(HueNumber,HueLetterCode)
# Assign ovoid segment interpolation to be linear or radial.  These assignments are based on
# visual inspection of the grid points in [Newhall1943, Figs. 1 through 9], rather than on any mathematical technique. 
if (Value == 1) {
   if (Chroma == 2){
      if ((ASTMHue > 15 && ASTMHue < 30) || (ASTMHue > 60 && ASTMHue < 85)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
  } else { if (Chroma == 4 ){	  
      if ((ASTMHue > 12.5 & ASTMHue < 27.5) || (ASTMHue > 57.5 & ASTMHue < 80)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 6){ 
      if ((ASTMHue > 55 & ASTMHue < 80)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 8 ){
      if ((ASTMHue > 67.5 & ASTMHue < 77.5)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 10){
      if ((ASTMHue > 72.5 & ASTMHue < 77.5)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }}}}
} else {
if (Value == 2){
   if (Chroma == 2){
      if ((ASTMHue > 15 && ASTMHue < 27.5) || (ASTMHue > 77.5 && ASTMHue < 80)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
} else { if (Chroma == 4){
      if ((ASTMHue > 12.5 & ASTMHue < 30) || (ASTMHue > 62.5 & ASTMHue < 80) ) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
  } else { if (Chroma == 6 ){
      if ((ASTMHue > 7.5 & ASTMHue < 22.5) || (ASTMHue > 62.5 & ASTMHue < 80)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
  } else { if (Chroma == 8 ){
      if ((ASTMHue > 7.5 & ASTMHue < 15)  || (ASTMHue > 60 & ASTMHue < 80) ) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
  } else { if (Chroma >= 10 ){
      if (ASTMHue > 65 & ASTMHue < 77.5)  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
  } else InterpStyle.Linear <- TRUE
   }}}}
} else {
if (Value == 3){
   if (Chroma == 2){
      if ((ASTMHue > 10 && ASTMHue < 37.5) || (ASTMHue > 65 && ASTMHue < 85)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 4 ){  
      if ((ASTMHue > 5 & ASTMHue < 37.5) || (ASTMHue > 55 & ASTMHue < 72.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 6 || Chroma == 8 || Chroma == 10){
      if ((ASTMHue > 7.5 & ASTMHue < 37.5) || (ASTMHue > 57.5 & ASTMHue < 82.5)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 12 ){
      if ((ASTMHue > 7.5 & ASTMHue < 42.5) || (ASTMHue > 57.5 & ASTMHue < 80))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
      } else InterpStyle.Linear <- TRUE
   }}}
} else {
if (Value == 4){
   if (Chroma == 2 || Chroma == 4){
      if ((ASTMHue > 7.5 && ASTMHue < 42.5) || (ASTMHue > 57.5 && ASTMHue < 85)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 6 || Chroma == 8 ){
      if ((ASTMHue > 7.5 & ASTMHue < 40) || (ASTMHue > 57.5 & ASTMHue < 82.5) ) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 10 ){
      if ((ASTMHue > 7.5 & ASTMHue < 40) || (ASTMHue > 57.5 & ASTMHue < 80) ) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }}
} else {
if (Value == 5){
   if (Chroma == 2 ){
      if ((ASTMHue > 5 && ASTMHue < 37.5) || (ASTMHue > 55 && ASTMHue < 85)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 4 || Chroma == 6 || Chroma == 8){
      if ((ASTMHue > 2.5 & ASTMHue < 42.5) || (ASTMHue > 55 & ASTMHue < 85))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 10 ){
      if ((ASTMHue > 2.5 & ASTMHue < 42.5) || (ASTMHue > 55 & ASTMHue < 82.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }}
} else {
if (Value == 6){
   if (Chroma == 2 || Chroma == 4){
      if ((ASTMHue > 5 && ASTMHue < 37.5) || (ASTMHue > 55 && ASTMHue < 87.5)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 6 )	 {
      if ((ASTMHue > 5 & ASTMHue < 42.5) || (ASTMHue > 57.5 & ASTMHue < 87.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 8 || Chroma == 10 	){  
      if ((ASTMHue > 5 & ASTMHue < 42.5) || (ASTMHue > 60 & ASTMHue < 85))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 12 || Chroma == 14 )	  {
      if ((ASTMHue > 5 & ASTMHue < 42.5) || (ASTMHue > 60 & ASTMHue < 82.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 16 ){
      if ((ASTMHue > 5 & ASTMHue < 42.5) || (ASTMHue > 60 & ASTMHue < 80))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }}}}
} else {
if (Value == 7	 ) {
   if (Chroma == 2 || Chroma == 4 || Chroma == 6){
      if ((ASTMHue > 5 && ASTMHue < 42.5) || (ASTMHue > 60 && ASTMHue < 85)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 8 	  ){
      if ((ASTMHue > 5 & ASTMHue < 42.5) || (ASTMHue > 60 & ASTMHue < 82.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 10 	)  {
      if ((ASTMHue > 30 & ASTMHue < 42.5) || (ASTMHue > 5 & ASTMHue < 25) || (ASTMHue > 60 & ASTMHue < 82.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 12  )	  {
      if ((ASTMHue > 30 & ASTMHue < 42.5) || (ASTMHue > 7.5 & ASTMHue < 27.5) || (ASTMHue > 80 & ASTMHue < 82.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 14 ){
      if ((ASTMHue > 32.5 & ASTMHue < 40) || (ASTMHue > 7.5 & ASTMHue < 15) || (ASTMHue > 80 & ASTMHue < 82.5))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }}}}
} else {
if (Value == 8){
   if ((Chroma == 2 || Chroma == 4 || Chroma == 6 || Chroma == 8 || Chroma == 10 || Chroma == 12)) {
      if ((ASTMHue > 5 && ASTMHue < 40) || (ASTMHue > 60 && ASTMHue < 85)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 14 ){
      if ((ASTMHue > 32.5 & ASTMHue < 40) || (ASTMHue > 5 & ASTMHue < 15) || (ASTMHue > 60 & ASTMHue < 85))  InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }	  
} else {
if (Value == 9){
   if (Chroma == 2 || Chroma == 4){
      if ((ASTMHue > 5 && ASTMHue < 40) || (ASTMHue > 55 && ASTMHue < 80)) InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma == 6 || Chroma == 8 || Chroma == 10 || Chroma == 12 || Chroma == 14	  ){
      if ((ASTMHue > 5 && ASTMHue < 42.5))   InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else { if (Chroma >= 16 ){
      if ((ASTMHue > 35 & ASTMHue < 42.5))   InterpStyle.Radial <- TRUE else InterpStyle.Linear <- TRUE
   } else InterpStyle.Linear <- TRUE
   }}}
}
}
}
}
}
}
}
}
# Set successful status return code and return
list(InterpStyle.Input=InterpStyle.Input,InterpStyle.Linear=InterpStyle.Linear,InterpStyle.Radial=InterpStyle.Radial,InterpStyle.OnGrid=InterpStyle.OnGrid, Status.ind = 1)
}


MunsellHueToASTMHue <- function(HueNumber,HueLetterCode)
{
ASTMHue <- 10*((7-HueLetterCode) %% 10) + HueNumber
if (ASTMHue == 0) ASTMHue = 100
ASTMHue
}


MunsellValueToLuminanceFactor <- function(MunsellValue)
{# Calculate luminance factor from quintic polynomial in ([Newhall1943, p. 417]).
Term1 <-  1.2219    * MunsellValue
Term2 <- -0.23111   * (MunsellValue^2)
Term3 <-  0.23951   * (MunsellValue^3)
Term4 <- -0.021009  * (MunsellValue^4)
Term5 <-  0.0008404 * (MunsellValue^5)
NewhallLuminanceFactor <- Term1 + Term2 + Term3 + Term4 + Term5
# Calculate luminance factor in terms of CIELAB model as given in ([Fairchild2005, Sect. 10.3]).
Lstar <- 10 * MunsellValue
fy    <- (Lstar + 16)/116
delta <- 6/29
if (fy > delta) CIELABLuminanceFactor <- fy^3 else CIELABLuminanceFactor <- 3*(fy - (16/116))*(delta^2)
CIELABLuminanceFactor <- 100 * CIELABLuminanceFactor
# Calculate luminance factor from expression in ([ASTMD1535-08, p. 4]).  This
# luminance factor is just 0.975 times the luminance factor in [Newhall1943].
Term1 <-  1.1914    * MunsellValue	
Term2 <- -0.22533   * (MunsellValue^2)
Term3 <-  0.23352   * (MunsellValue^3)
Term4 <- -0.020484  * (MunsellValue^4)
Term5 <-  0.00081939 * (MunsellValue^5)
ASTMLuminanceFactor <- Term1 + Term2 + Term3 + Term4 + Term5
list(OriginalMunsellValue = MunsellValue, Newhall1943 = NewhallLuminanceFactor, CIELAB = CIELABLuminanceFactor, ASTMD153508 = ASTMLuminanceFactor)
}

MunsellToxyY<-function(MunsellSpec, InterpolateByLuminanceFactor=TRUE)
{#Convert a Munsell specification into xyY coordinates, by interpolating over the extrapolated Munsell renotation data.
if (is.character(MunsellSpec)) ColorLabMunsellVector <- MunsellSpecToColorLabFormat(MunsellSpec) else ColorLabMunsellVector <- MunsellSpec
if (length(ColorLabMunsellVector)>1){
  H <- ColorLabMunsellVector[1]
  V <- ColorLabMunsellVector[2]
  C <- ColorLabMunsellVector[3]
  HueLetterCode <- ColorLabMunsellVector[4]
} else {
H <- 'N'
V <- as.numeric(ColorLabMunsellVector)
C <- 0
HueLetterCode<-getMunsellHueLetterDesignator(MunsellSpec)
}
if ((V < 1) | (V > 10)) stop('Munsell value must be within the limits of the renotation data (between 1 and 10)')
LuminanceFactors <- MunsellValueToLuminanceFactor(V)
Y <- LuminanceFactors$ASTMD153508  #MunsellV2Y(V)
if ((abs(V - round(V)) < 0.001)){
ValueMinus <- round(V)
ValuePlus  <- round(V)
} else {
ValueMinus <- floor(V)
ValuePlus  <- ValueMinus + 1
}
if (H=='N') ColorLabVecMinus <-ValueMinus else ColorLabVecMinus <- c(H,ValueMinus,C,HueLetterCode)
mtemp <- MunsellToxyForIntegerMunsellValue(ColorLabVecMinus)
xminus<-mtemp$x
yminus<-mtemp$y
Status.ind<-mtemp$Status.ind
if (Status.ind != 1)  return(list(Status.ind = 2))
if (length(ColorLabMunsellVector) == 1) ColorLabVecPlus <- ValuePlus  else {   # Colour is Munsell grey
if (ValuePlus == 10) ColorLabVecPlus <- ValuePlus  else ColorLabVecPlus <- c(H,ValuePlus,C,HueLetterCode)
}
mtemp <-  MunsellToxyForIntegerMunsellValue(ColorLabVecPlus)
xplus<-mtemp$x
yplus<-mtemp$y
Status.ind <- mtemp$Status.ind
if (Status.ind != 1)  return(list(Status.ind = 2)) # No reflectance data available for bounding colour.  Set error code and return
# Interpolate between the xy coordinates for the lighter and darker Munsell samples
if (ValueMinus == ValuePlus){
   x <- xminus
   y <- yminus
} else {
   LuminanceFactors <- MunsellValueToLuminanceFactor(ValueMinus)
   YMinus           <- LuminanceFactors$ASTMD153508
   LuminanceFactors <- MunsellValueToLuminanceFactor(ValuePlus)	
   YPlus            <- LuminanceFactors$ASTMD153508
   if (InterpolateByLuminanceFactor) { # Interpolate over luminance factor
   x <- approx(c(YMinus, YPlus), c(xminus, xplus), Y)$y
   y <- approx(c(YMinus, YPlus), c(yminus, yplus), Y)$y
   } else {
   # The following option can be used to test interpolating over value vs interpolating over luminance factor
   # Interpolate over value instead of over luminance factor
      x <- approx(c(ValueMinus, ValuePlus), c(xminus, xplus), V)$y
      y <- approx(c(ValueMinus, ValuePlus), c(yminus, yplus), V)$y
   }
}
# Set successful status return code and return
list(x=x,y=y,Y=Y,Status.ind = 1)
}

MunsellHVC<-function(MunIn)
{# Convert a Munsell string to HVC
m<-sapply(MunIn, function(x) {if (is.na(x)) return(cbind(H=NA,V=NA,C=NA)) else {
if (substr(x,1,1)=='N'){
if (nchar(x)==1) return(cbind(H=NA,V=NA,C=NA)) else {
H1<-''
C<-'0'
H2<-'N'
V<-substr(x,2,nchar(x))
if (!(grepl('^(\\d{1,2}|\\d{1,2}\\.\\d{1,2})$',V))) return(cbind(H=NA,V=NA,C=NA))
}
} else {
sidN <- gregexpr('^(\\d{1,2}\\.\\d{1,2}|\\d{1,2})(RP|YR|Y|GY|G|BG|B|PB|P|R) +(\\d{1,2}\\.\\d{1,2}|\\d{1,2})/(\\d{1,2}\\.\\d{1,2}|\\d{1,2})', x, perl=TRUE)
H1<-substr(x,attr(sidN[[1]], "capture.start")[1,][1],attr(sidN[[1]], "capture.start")[1,][1]+attr(sidN[[1]], "capture.length")[1,][1]-1)
H2<-substr(x,attr(sidN[[1]], "capture.start")[1,][2],attr(sidN[[1]], "capture.start")[1,][2]+attr(sidN[[1]], "capture.length")[1,][2]-1)
V<-substr(x,attr(sidN[[1]], "capture.start")[1,][3],attr(sidN[[1]], "capture.start")[1,][3]+attr(sidN[[1]], "capture.length")[1,][3]-1)
C<-substr(x,attr(sidN[[1]], "capture.start")[1,][4],attr(sidN[[1]], "capture.start")[1,][4]+attr(sidN[[1]], "capture.length")[1,][4]-1)
if (C=='') return(cbind(H=NA,V=NA,C=NA))
}
cbind(H=paste(H1,H2,sep=''),V=V,C=C)
} })
m<-t(m)
colnames(m)<-c('H','V','C')
rownames(m)<-1:length(MunIn)
m
}

