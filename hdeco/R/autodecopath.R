"autodecopath" <-
function () {

  #############################################################
  # 
  # TITLE:  autodecopath()
  # AUTHOR: SANDOR KABOS
  # DATE:   20 AUG, 2003
  # CALLS:  
  # NEEDS:  
  # NOTES:  BUILD VFONAL DECOMPOSITION PATH MATRIX BASED ON
  # 		DEFAULT REQUIREMENTS FOR KABOS TYPE DATA
  #
  #############################################################

  BE <- .QND
  BE <- names(BE)
  H <- length(BE)
  
  # CHECK HOW MANY IMAGES
  if(dim(.QKEP)[1]==1) {
    EGYKE <- T
  }
  else {
    EGYKE <- F
  }
 
  # CREATE VECTORS OF TRUE/FALSE FOR X, Y, AND Z VARIABLES
  YEK <- "Y"==substring(BE,1,1)
  ZEK <- "Z"==substring(BE,1,1)
  XEK <- "X"==substring(BE,1,1)
  # STORE COUNTS OF HOW MANY OF EACH TYPE OF VARIABLE THERE ARE
  NX <- sum(XEK)
  NY <- sum(YEK)
  NZ <- sum(ZEK)
  
  # DIFFERENT APPROACHES DEPENDING ON NUMBER OF IMAGES
  if(EGYKE){  # ONE IMAGE
    # BUILD A MATRIX OF -1s
    VFONAL <- matrix(-1,nrow=NY,ncol=H,dimnames=list(NULL,BE))
    # CONVERT ALL Z VARIABLES TO 0s
    VFONAL[,ZEK] <- 0
    # CONVERT Y VARIABLES TO 1s ON LOWER TRIANGLE
    for (i in 1:NY) {
      VFONAL[i,(NX+1):(NX+i)] <- 1
    }
    # CONVERT DIAGONALS FOR STEPS >=2 TO 2s
    if(NY > 1) {
      for (i in 2:NY) {
        VFONAL[i,(NX+i)] <- 2
      }
    } 
  }
  else {  # TWO IMAGES
    VFONAL <- matrix(-1,nrow=(NY+1),ncol=H,dimnames=list(NULL,BE))
    VFONAL[,XEK] <- 0
    VFONAL[,ZEK] <- 1
    for (i in 1:NY) VFONAL[i+1,(NX+1):(NX+i)] <- 1
    if(NY > 1) {
      for (i in 2:NY) {
        VFONAL[i+1,(NX+i)] <- 2
      }
    } 
  }

  # APPLY DEFAULT TITLE
  attr(VFONAL,"cim") <- "Default Decomposition Path"

  # RETURN THE DECOMPOSITION PATH
  return(VFONAL)

}

