#########################################################################
#
# Copyright Weierstrass Institute for Applied Analysis and 
#           Stochastics (WIAS) & Humboldt Universitaet zu Berlin, 
#           Institut fuer Mathematik, Germany 2006
# *********************************************************
#
# Name:          phantom.r
#                ---------------
# Author:        Joern Schulz
# Stand:         16.08.2006
#
#########################################################################

phantom <- function(n=257, design="A", addIm="none", DebugLevel="Normal")
{
# This function create two dimensional phantom data. It's generate a big 
# ellipse, that represents the head and several smaller ellipses, that
# represents  pathological areas to be located in the space of the bigger 
# ellipse.
#
# Syntax of function:
# -------------------
#  phantom(nameP1=valueP1, nameP2=valueP2,...)
#  with:
#                nameP1  = name of parameter P1 as a string enclosed in ' '
#                valueP1 = value of parameter P1
#
# for P2,...etc. correspondingly.
# 
# parameter default     meaning
# name      value
# ------------------------------------------------------------------------
# n         257         n is the number of columns and rows in the generated
#                       phantom. It is assumed that the number of cloumns equal
#                       to the number of rows.
# design    'A'         The design characterized the phantom data. It is
#                       possible to define different ellipse, with different
#                       intensity. There are three default-designs, these are
#                       'A','B' and 'C'. To define a own design of the phantom
#                       you have to note following conditions:
# 
#                       - P is a (n,5) or (n,6) matrix, whereas n > 0
#                       - all elements x of P have to be between -1 and 1
#                       
#                            (:,1)  (:,2)  (:,3)  (:,4)  (:,5)  (:,6)
#
#                       P = [ ...    ...    ...    ...    ...    ...     (1,:)
#                                                .
#                                                .
#                                                .
#                             ...    ...    ...    ...    ...    ...  ]  (n,:)
#
#                       at what each row define a ellipse and where
#                       - (:,1) is the the additive intensity of the corresponding
#                               ellipse
#                       - (:,2) the x-coordinate of the center of the ellipse
#                       - (:,3) the y-coordinate of the center of the ellipse
#                       - (:,4) the half length of the horizontal axis
#                       - (:,5) the half length of the vertical axis
#                       - (:,6) the angle in degree between the x-axis of the 
#                               ellipse and the x-axis of the grid.
#                               This parameter is optional, but the number of 
#                               columns have to be the same in all rows.
# addIm     "blurred1"  There are six default-designs, these are "blurred1", "blurred2", 
#                       "keen1", "keen2", "simple1" and "simple2". The default add-Image
#                       are generated with the function 'partEllipse'. If addIm="none", no
#                       image is add to the phantom. A further possibility is a matrix A
#                       with the the same size as the phantom, addIm=A.
#                       Also possible is to set addIm="None".
#
# Example:
# --------
# A <- phantom()
# B <- phantom(design="B")
# C <- phantom(design="C")
# viewData(list(A,"design='A'",
#               B, "design='A'",
#	          C, "design='A'"))
#

# variante mit switch
# design <- switch(design, "A"=modelA, "B"=modelB, "C"=modelB)

  DL1 <- logDebug(DebugLevel)
  DebugLevel <- DL1[[3]]
  DL2 <- DL1[[2]]
  DL1 <- DL1[[1]]

  if(is.character(design)){
  
   # ==========================================================
   # given default models
   # ==========================================================
    model <- function(design){
    if (design == "A"){
      #
      #   Shepp-Logan head phantom with values from
      #   Peter Toft: "The Radon Transform - Theory and Implementation", 
      #   Ph.D. thesis. Department of Mathematical Modelling, 
      #   Technical University of Denmark, June 1996,
      #   Web: http://pto.linux.dk/PhD/.
      #   
      #      
      #                  A    x0      y0     a       b      alpha
      #                 ------------------------------------------
        par <- matrix(c(  1,   0 ,    0   , .69  , .92,        0,   
                        -.7,   0 ,  -.0184, .6624, .8740,      0,
                        -.2,  .22,    0   , .1100, .3100,    -18,
                        -.2, -.22,    0   , .1600, .4100,     18, 
                         .3,   0 ,   .35  , .2100, .2500,      0,
                         .2,   0 ,   .1   , .0460, .0460,      0,
                         .1,   0 ,  -.1   , .0460, .0460,      0,
                         .2, -.08,  -.605 , .0460, .0230,      0,
                         .1,   0 ,  -.606 , .0230, .0230,      0,
                        -.3,  .06,  -.605 , .0230, .0460,      0 ),
                      nrow=10 , ncol=6 , byrow=TRUE)
    } 
    else if(design == "B"){
      #
      #  Self constructed data.
      #
      #                   A     x0     y0     a      b    
      #                 ---------------------------------
        par <- matrix(c( .5 ,   0 ,    0 ,   .9 ,   .9 ,
                         .1 , -.6 ,    0 ,   .2 ,   .2 ,
                         .3 ,   0 ,    0 ,   .2 ,   .2 ,
                        -.1 ,  .6 ,    0 ,   .2 ,   .2 ,
                        -.2 ,   0 ,   .55,   .2 ,   .2 ,
                         .2 ,  .0 ,  -.55,   .2 ,   .2 ,
                         .4 , -.55,   .55,   .05,   .05,
                        -.3 ,  .55,   .55,   .05,   .05,
                         .5, -.55,  -.55,   .05,   .05,
                        -.4,  .55,  -.55,   .05,   .05 ),
                      nrow=10 , ncol=5 , byrow=TRUE)
    }
    else if(design == "C"){
      #
      #  Self constructed data.
      #
      #                    A     x0     y0     a      b    alpha  
      #                  ---------------------------------------
        par <- matrix(c( .3 ,   0 ,    0 ,   .9 ,   .6 ,    0,
                         .7 ,   0 ,    0 ,   .2 ,   .06,  115,
                        -.1 ,   0 ,    0 ,   .2 ,   .4 ,    0,
                        -.1 ,   0 ,    0 ,   .2 ,   .1 ,   45,
                         .15, -.6 ,   .05,   .1,    .03,    0,
                         .25,  .65,  -.05,   .05,   .1 ,   15 ),
                      nrow=6 , ncol=6 , byrow=TRUE)
    } 
    else if(design == "D"){
      #
      #  Self constructed data.
      #
      #                  A    x0      y0     a       b      alpha
      #                 ------------------------------------------
        par <- matrix(c(  1,   0 ,    0   , .69  , .92, ,      0,   
                        -.8,   0 ,  -.0184, .6624, .8740,      0,),
                      nrow=10 , ncol=6 , byrow=TRUE)
    }
    else stop("The default design ", design, " is not defined")
    return(par)
    } # End of function model
  
    design <- model(design)
    ddesign <- dim(design)
    nd1 <- ddesign[1]
    nd2 <- ddesign[2]
  
  } else if(is.vector(design)){
    if(length(design) != 5 && length(design) != 6  ){
       stop("The size ", length(design), " doesn't be allowed." )
    }
    nd1 <- 1
    nd2 <- length(design)
    ddesign <- c(nd1,nd2)
    design <- matrix(design, nrow=nd1, ncol=nd2, byrow=TRUE)
  
  } else if(is.matrix(design)){
    ddesign <- dim(design)
    nd1 <- ddesign[1]
    nd2 <- ddesign[2]
    if (is.null(ddesign) || (nd2 != 5 & nd2 !=6)){
       stop("The size of the matrix doesn't be allowed.")
    }
  
  } else stop("The design doesn't correspond to conditions.")

  if (is.character(addIm)){
    if (!(addIm %in% c("blurred1", "blurred2", "keen1", "keen2", "simple1", "simple2", "none")))
        stop("The value from addim ", addIm, " isn't possile.")
  } else if(is.matrix(addIm)){
    if (dim(addIm)[1] != n || dim(addIm)[2] != n){
        stop("The dimensions of 'addIm' and 'design' have to be the same.")
    }
  } else stop("The data-typ of 'addIm' has to be of typ 'character' or 'matrix'.")

 # ========================================================================
 # Creation and start the example

  if(DL1) cat("Creation of the source-data --> ")
  data <- matrix(0,n,n)
  Ax <- matrix(seq(-1,1,length.out=n), nrow=n, ncol=n, byrow=TRUE)
 # Matrix Ax 90° gedreht
  rot90Ax <- matrix(seq(1,-1,length.out=n), nrow=n, ncol=n, byrow=FALSE)

  for (i in 1:nd1){
	intensity <- design[i,1]
	x0 <- design[i,2]
	y0 <- design[i,3]
	a2 <- (design[i,4])^2
	b2 <- (design[i,5])^2
	if (nd2 == 5){
		Axnew <- Ax - x0
		Aynew <- rot90Ax - y0
		index <- which( ((Axnew)^2) / a2 + ((Aynew)^2) / b2  <= 1 )
		data[index] <- data[index] + intensity
		}
	else if (nd2 == 6){
		alpha <- design[i,6]
		Axnew <- Ax-x0
		Aynew <- rot90Ax - y0                        # Artificial PET data
		Axnew <- (Axnew)*cos((alpha*pi)/180) + (Aynew)*sin((alpha*pi)/180)
		Aynew <- (Aynew)*cos((alpha*pi)/180) + (Axnew)*sin((alpha*pi)/180)
		index <- which( ((Axnew)^2) / a2 + ((Aynew)^2) / b2  <= 1 )
		data[index] <- data[index] + intensity
		}
	else  stop("Incorrect dimension.")
  }

 # add a additional figure to source-data
  if (is.character(addIm)){
    if (addIm == "blurred1"){
        data <- data + partEllipse(mod="hcirc3", x=0, y=0.25, intensity=0.3, n=n,
                       re1=0.01, re2=0.07,ring.wide=0.2,
                       no.xax1=-1, no.xax2=1, no.yax1=0, no.yax2=1,in.r=0.0001,
                       DebugLevel)
    } else if(addIm == "blurred2"){
        data <- data + partEllipse(mod="hcirc3", x=.0, y=-0.4, intensity=0.4, n=n,
                       re1=0.2, re2=0.05, ring.wide=0.18,
                       no.xax1=-1, no.xax2=1, no.yax1=-0.18, no.yax2=0.1,in.r=0.0005,
			     DebugLevel)
    } else if(addIm == "keen1"){
        data <- data + partEllipse(mod="hcirc2", x=0, y=0.25, intensity=0.25, n=n,
                       re1=0.01, re2=0.07,ring.wide=0.05,
                       no.xax1=-1, no.xax2=1, no.yax1=0, no.yax2=1,in.r=0.001)
    } else if(addIm == "keen2"){
        data <- data + partEllipse(mod="hcirc2", x=0.35, y=-0.4, intensity=0.1, n=n, 
                       re1=0.001,re2=0.005,ring.wide=0.1,
                       no.xax1=-1, no.xax2=1, no.yax1=-1, no.yax2=1,in.r=0.0008)
        data <- data + partEllipse(mod="hcirc2", x=-0.35, y=-0.38, intensity=0.05,
                       n=n, re1=0.005,re2=0.009,ring.wide=0.2,
                       no.xax1=-1, no.xax2=0.01, no.yax1=-1, no.yax2=0.04,in.r=0.0001)
    } else if(addIm == "simple1"){
        data <- data + partEllipse(mod="hcirc1", x=0, y=0, intensity=0.1, n=n, 
                       re1=0.25,re2=0.25,ring.wide=0.07,
                       no.xax1=-1, no.xax2=1, no.yax1=-1, no.yax2=1)
        data <- data + partEllipse(mod="hcirc1", x=0, y=0, intensity=-0.1, n=n, 
                       re1=0.25,re2=0.5,ring.wide=0.07,
                       no.xax1=-1, no.xax2=1, no.yax1=-1, no.yax2=1)
    } else if(addIm == "simple2"){
        data <- data + partEllipse(mod="hcirc1", x=.0, y=-0.4, intensity=0.4, n=n,
                       re1=0.2, re2=0.05, ring.wide=0.18,
                       no.xax1=-1, no.xax2=1, no.yax1=-0.18, no.yax2=0.1)
    } else if(addIm == "simple3"){
        data <- data + partEllipse(mod="hcirc1", x=0, y=0.0, intensity=0.6, n=n,
                       re1=0.1, re2=0.17,ring.wide=0.5,
                       no.xax1=-1, no.xax2=1, no.yax1=0, no.yax2=1,in.r=0.0001)
        data <- data + partEllipse(mod="hcirc1", x=0.265, y=0.0, intensity=0.6, n=n,
                       re1=0.0024, re2=0.0024,ring.wide=1, 
                       no.xax1=-1, no.xax2=1, no.yax1=-1, no.yax2=0)
        data <- data + partEllipse(mod="hcirc1", x=-0.265, y=0.0, intensity=0.6, n=n,
                       re1=0.0024, re2=0.0024,ring.wide=1, 
                       no.xax1=-1, no.xax2=1, no.yax1=-1, no.yax2=0)
        data <- data + partEllipse(mod="hcirc1", x=0, y=-0.4, intensity=0.4, n=n,
                       re1=0.005, re2=0.01,ring.wide=1 )

    } else if(addIm!="none") stop("A value of addim ", addIm, " isn't possile.")
  } else if(is.numeric(addIm)){
    data <- data + addIm
  }

  if(DL1) cat("complete. \n")

  data <- scaleImage(data)
  return(data)
}





