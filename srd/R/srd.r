#roxygen embedded tags:

#' Draws scaled rectangle diagrams to represent categorical data.
#'
#' @param  data  A binary data matrix of ones and zeroes (or TRUE and FALSE) with k columns, with 1 or TRUE 
#' representing presence of an attribute that represents a rectangle.
#' @param weight  A frequency count vector attaches a frequency weight to each row of \code{d}.
#' @param labelcell Attaches labels to cells. 
#' One of "none", "frequency", "area", "error", or "residual". Use "frequency" for cell frequencies (default), "area" for areas, "error" for percentage error, "residual" for Pearson residuals
#' @param  new TRUE or FALSE. Fits a new configuration, new=FALSE  assumes an existing configuration already exists
#' @param criterion Criterion to determine fitting optimisation.  One of "lsqs", "labs", "logl", or "chi2". 
#' Use "lsqs" for least squares (default), "labs" for least absolute difference, "logl" for minus log-likelihood, 
#' "chi2" for chi-square
#' @param colour Determines colour of configuration. One of "earthy", "bright", or  "mono". Use "earthy" for muted
#'  colours, "bright" for rainbow colours (default), "mono" is monochrome (shades of grey). 
#'  If \code{cells=TRUE},  "earthy" and "bright" are interpreted as shades of yellow and blue respectively.
#' @param labelrectangle To attach labels to each of the k rectangles. 
#' One of "none", "inside", or "arrow". Use "none" for no labels,  "inside" (default) for inside rectangle, 
#' "arrow" is an arrow pointer from outside of the unit square.
#' @param y A vector of data (with same number of rows as \code{d}) for the third axis in a 3D projection . 
#' Default NA is equivalent to a vector of zeroes. 
#' Projected values are  mean values of \code{y}, or \eqn{\sum y / \sum w} if \code{w} is specified. Projected values are of
#' each rectangle, or if \code{cells=TRUE} of each cell (see \code{cells} parameter)
#' @param  w  A numeric vector associated with  \code{y}. 
#' Statistics of \code{y} for cells/rectangles (see cells argument) represented by \eqn{\sum y / \sum w}. 
#' Default NA gives 3D projection as mean values of \code{y} (equivalent to \code{w=1}). 
#' \code{w} can be used, for example, to obtain 3D projection axis as a "incidence rate" if \code{y} is 
#' a binary indicator of an outcome and \code{w} is person-time.
#' @param cells TRUE or FALSE. True for a diagram shaded by mean of y in each cell of 
#' intersecting rectangles and 3D projection showing each cell. Default (FALSE) is  transparent represention of each of k rectangles, 
#' and 3D projection the mean of \code{y}
#' of each rectangle. 
#' @param motion Determines whether and how configuration is put into 3D mode. 
#' One of values "none", "frames", "movie".  "none"  for a flat single image, "frames" for broken series 
#' of rotated images, "movie"  for an animated rotating movie. 
#' @param thin Parameter to attempt to avoid thin and elongated rectangles. The 
#' \code{criterion} is penalised by multiplying by \eqn{x^thin} where is \eqn{x}
#'  is maximal thinness (length/breadth) of the \eqn{q} rectangles. 
#'  Recommended value \eqn{<=1}. Default NA is no penalising, equivalent to \eqn{thin=0}.
#'  
#' @return Error message, or indication of successful completion.
#' @export
#' @author Roger Marshall, <rj.marshall@@auckland.ac.nz>, The University of Auckland, New Zealand
#' 
#' @details A function to visualise  k attributes as a Venn-like  diagram 
#' with areas scaled, as best as possible, to be proportional to  frequency. 
#' It is a representation is of a \eqn{2^q} contingency table. 
#' Up to six attributes (rectangles) can be displayed.
#' Fitting is done by optimising the congruence between cell area and cell frequency by the 
#' chosen criterion. A flat diagram can be thrown into an animated 3D projection, showing
#' mean/rate values of the associated \code{y} with \code{motion} equal to "frames" or "movie".
#' @import plyr colorspace animation
#' @references 
#' Marshall, RJ.  Displaying  categorical data relationships by scaled rectangle diagrams. Statistics in Medicine 20, 1077-1088, 2001.
#' Marshall, RJ. Scaled rectangle diagrams can be used to visualise clinical and epidemiological data. 
#' Journal of Clinical Epidemiology 58, 974-981, 2005.
#' Marshall, RJ. Determining and visualising at-risk groups in case-control data. 
#' Journal of Epidemiology and Biostatistics  6,  343-348, 2001

#' @examples  
#' library(datasets)
#' data(Titanic)
#' x <- data.frame(Titanic)
#' x$women<-  ifelse(x$Sex=="Female",1,0)
#' x$child <- ifelse(x$Age=="Child",1,0)
#' x$first <- ifelse(x$Class=="1st",1,0)
#' x$survived <- ifelse(x$Survived=="yes",1,0)
#' d <- x[c("women","child","first")]
#' srd(data=d,y=x$survived,weight=x$Freq,motion="frames")
#' 
#'
#' library(survival)
#' data(pbc) 
#' ## create srd and show projection axis as incidence rate (by specifying w as time)
#' pbc$male <- ifelse(pbc$sex=="m",1,0)
#' d  <- with(pbc,cbind(ascites,hepato,spiders,male))
#' pbc$event <- ifelse(pbc$status==2,1,0)
#' srd(d, y=pbc$event, w=pbc$time, colour="earthy", 
#'    labelrectangle="arrow",cells=FALSE,motion="frames", criterion="chi2", thin=0.5)
#' 
srd <- function(data=stop("'data' must be specified"),weight=NA,labelcell="frequency", new=TRUE,
  criterion="logl",colour="earthy",labelrectangle="inside",y=NA,w=NA,cells=FALSE,motion="none",thin=NA)
{
 


 #library(stringr)
 #library(colorspace)
 #library(plyr)
 #library(animation) 
 
 #  change character attributes to numeric for subsequent code and check values
 lst <- list("earthy","bright","mono")
  xxxx <- colour
  colour <- which(lst==colour)
   if(is.na(colour[1])) {
     return(message(paste("Invalid colour",xxxx)))
   }
  
  lst <- list("lsqs","labs","logl","chi2")
  xxxx <- criterion
  crit <- which(lst==criterion)
  if(is.na(crit[1])) {
    return(message(paste("Invalid criterion",xxxx)))
  }
  
  lst <- list("none","inside","arrow")
  xxxx <- labelrectangle
  labelrectangle <- which(lst==labelrectangle) - 1
  if(is.na(labelrectangle[1])) {
   return(message(paste("Invalid labelrectangle",xxxx)))
  }
  
  lst <- list("none","frequency","area","error","residual")
  xxxx <- labelcell
  labelcell <- which(lst==labelcell) - 1
  if(is.na(labelcell[1])) {
    return(message(paste("Invalid labelcell",xxxx)))
  }
  lst <- list("none","frames","movie")
  xxxx <- motion
  labelmotion <- which(lst==motion)
  if(is.na(labelmotion[1])) {
    return(message(paste("Invalid motion",xxxx)))
  }
  
  ## can only accommodate up to 6 rectangles 
  if(ncol(data)> 6){
    message("srd only possible for up to six rectangles")
    return(message("Terminated"))
  }
 
## thinness parameter if NA is equaivalent to zero
 if(is.na(thin)==TRUE) thin <- 0
  
  #-----------------------------------------------------------------------
  ## create flat representation, single image
  #d <- as.matrix(d)
  
## make, weight=1, w=1, and y=0 values if NAs
  if(length(weight)==1){ 
  if(is.na(weight)){weight <- matrix(data=1, ncol=1,nrow=nrow(data))}
  }
if(length(w)==1){ 
 if(is.na(w)){w <- matrix(data=1, ncol=1,nrow=nrow(data))}
}
if(length(y)==1){ 
 if(is.na(y)){y <- matrix(data=0, ncol=1,nrow=nrow(data))}
 
}


print(paste("Creating srd for ", ncol(data)," rectangles"))
print(paste("Optimising criterion ", criterion))
if(thin !=0) print(paste("Penalised for thinness, parameter=",thin))


  if(motion=="none"){

  drawn <- draw.srd(d=data,weight, labelcell=labelcell,new=new,colour=colour, crit, 
          labelrectangle=labelrectangle,Y=y,project=c(0,0,0),cells=cells, w, thin)
  if(drawn$fail){return(message("srd failed"))}
  else{return("SRD created")}
  }
  ##---------------------------------------------------------------------------------------------------------
  ## create  animated representation, requiring package(animation)
  
  if(motion=="movie"){
    
    drawn <- draw.srd(d=data,weight, labelcell=labelcell,new=new,colour=colour, crit, 
             labelrectangle=labelrectangle,Y=y,project=c(0,0,0), cells=cells,w,thin)                                                                                                                                                                                                                                                                                                                                                                                                  
    if(drawn$fail){return(message("srd failed"))}
#    library(animation)
    
    note <- c(paste(" Creating animation frames. Please wait. "))  
    print(note)
    
    
    saveGIF({
      
      ani.options(interval=0.2)
      
      renew <- FALSE 
      for (i in 2:90){
        
        ya <-  2*(i-1)
        xa <-  2*(i-1)
        za <-  2*(i-1)
        
        ##  note in rotation beyond 0,0,0  labelrectangle=0, turns off labelcells
        lr<- labelrectangle
        lc <- labelcell
        if(cells==TRUE & labelrectangle>0){lr <- 0}
        if(cells==FALSE & labelcell>0){ lc <- 0}
        
        draw.srd(d=data,weight, labelcell=lc,new=FALSE,colour=colour, crit, 
                 labelrectangle=lr,Y=y,project=c(xa,ya,za),cells=cells,w,thin)
    
      }
      
      ## restore the options
      ##ani.options(oopt)
    })
    return("Movie created")
  }
##------------------------------------------------------------------
## create series of frames in 3D
 if(motion=="frames") {
   
   drawn <- draw.srd(d=data,weight, labelcell=labelcell,new=new,colour=colour, crit, 
            labelrectangle=labelrectangle,Y=y,project=c(0,0,0),cells=cells, w,thin)
     if(drawn$fail){return(message("srd failed"))}
    #  manually step rotations, pausing with browser()
    renew <- FALSE 
    for (i in 2:90){
      
      ya <-  2*(i-1)
      xa <-  2*(i-1)
      za <-  2*(i-1)
      
      
      devAskNewPage(ask=TRUE)
      lr<- labelrectangle
      lc <- labelcell
      if(cells==TRUE & labelrectangle>0){lr <- 0}
      if(cells==FALSE & labelcell>0){ lc <- 0}
      
      draw.srd(d=data,weight, labelcell=lc,new=FALSE,colour=colour, crit, 
               labelrectangle=lr,Y=y,project=c(xa,ya,za),cells=cells,w,thin)

    }
 
 
  draw.srd(d=data,weight, labelcell=labelcell,new=renew,colour=colour, crit, 
           labelrectangle=labelrectangle,Y=y,project=c(0,0,0),cells=cells,w,thin)
  return("Frames created")    
}

}
##-----------------------------------------------------------
## function to  create a Scaled Rectangle Diagram
draw.srd<- function(d,f=matrix(data=1, ncol=1,nrow=nrow(d)),labelcell=0,new=TRUE, 
          crit=3,colour=2,labelrectangle=1, Y=matrix(data=0, ncol=1,nrow=nrow(d)),
          project=c(0,0,0),cells=FALSE,w=matrix(data=1, ncol=1,nrow=nrow(d)),thin=0)
  {
  ##------------------------------------------------------------------------------------------------------
  ## d is data matrix of ones and zeroes
  ## number of data rows n and columns for q rectangles
  ## f is frequency count vector (default vector of 1's)
  ## labelcell=0, no labelcells (default), 1= cell frequencies, 2= area frequencies, 3= % error
  ## new  re-fits a new SRD  (default TRUE). Otherwise the existing configuration is used 
  ## crit=1 least squares (default), 2=least absolute difference, 3= minus log-likelihood
  ## colour=1  pale blue/brown colours, 2=rainbow colours (default), 3=monochrome  (if cells=FALSE)
  ## colour=1,2,3,4 for pale yellow, pale blue, black and white intensity shading if cells=TRUE
  ##       = 4 for full non-see-through rectangles
  ## labelrectangle=0 no rectangle labelcell, =1 (default) inside rectangle, 2=arrow pointer
  ## Y vector of  variable for third axis in a 3-D projection, default  zeroes
  ## w is weight associated with  3-D projection axis.
  ## project - vector of viewing projection angles (degrees) for 3D projection. Default no projection (0,0,0)
  ## cells=FALSE. If TRUE gives representation of cells of the configuation and projected to mean cell value
  ##         
  ##---------------------------------------------------------------------------------------------------------
  if(is.matrix(d)==FALSE) {d <- as.matrix(d)}
  if(is.matrix(Y)==FALSE) {Y <- as.matrix(Y)}
  if(is.matrix(f)==FALSE) {f <- as.matrix(f)}
  if(is.matrix(w)==FALSE) {w <- as.matrix(w)} 
  
  fail <- FALSE

  n_in <- nrow(d)
  n_out <- n_in
  q <- ncol(d)

## Delete rows of d  and corresponding elements of Y if incomplete

rnames <- colnames(d)
yvar <- colnames(Y)

ok <- complete.cases(d)
if(sum(!ok) ){


if(q==1){
##seems to need fix up when q=1 to ensure d is a matrix  
  
  d <-  d[ok]
  
d <- matrix(d,ncol=1,dimnames=list(c(NULL),colnames(d)[1]))
}

else{
  
  d <-  d[ok,]
  
}
n_out <- nrow(d)


##yvar <- colnames(Y)
 Y <- Y[ok] 
 f <- f[ok]

w <- w[ok]

}
##  if new=TRUE  create a new srd configuration by invoking FORTRAN
##  subroutine srd.
  if(new){
   if(q>6){
      message("Cannot create an SRD for q>6 rectangles")
      
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
   }
    ##test for strictly binary  variables
    for (i in 1:q){
      
      if(q>1){
    test <-unique(d[,i])}
    else
      {test <- unique(d[,])}

    if(length(test)!=2){
      v <-colnames(d)[i]
      message(paste("Variable ",v, " not binary"))
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
    }
    
    else
    {if((test[1]!=0 | test[2] !=1)&(test[1]!=1 | test[2] !=0)){
      v <-colnames(d)[i]
      message(paste("Variable ",v," is binary but not coded 0 1"))
     fail <- TRUE 
      output <- list(fail=fail)
      return(output)
      
    }}   
    }
    ncells <- 2**q
    
   ## load Intel FORTRAN created DLL to  fit rectangle diagram 
   ##  dyn.load commented out in  CRAN package release
    
  ##dyn.load("C:/temp/q6000036408/Dll_srd/x64/Release/Dll_srd.dll")

    
    
    ## invoke associated srd subroutine  to fit rectangle coordinates and % error
    x<- matrix(nrow=4,ncol=7)
    y<- matrix(nrow=4,ncol=7)
    cstat<- matrix(nrow=5,ncol=64)

#   required random number seed for FORTRAN 77 srd subroutine to be
#   passed as argument.    
    seed <- round(100000*runif(1))
#    seed <- 12345


print(paste("Rectangle labels:"))
for(i in 1:q){
message(paste("     ", colnames(d)[i]))
}
doptm <- 0
thinx <- 0
    run_srd <- .Fortran("srd",d=as.integer(d),f=as.integer(f),n_out=as.integer(n_out),
                q=as.integer(q),x=as.double(x),y=as.double(y),e=as.single(0),cstat=as.single(cstat),
                crit=as.integer(crit),seed=as.integer(seed),thin=as.single(thin), 
                doptm=as.double(doptm),thinx=as.single(thinx),NAOK=TRUE )

# save data and configuration as global objects, for subsequent new=FALSE
dsave <- NULL
run_srd_save <- NULL
dsave <<- d
run_srd_save <<- run_srd 
if(n_out<n_in) {
    note <- c(paste("Missing data (d), n reduced from ",n_in, "to ", n_out))  
   message(note)
       
    }

if(thin==0) {print(paste("Optimised criterion=",round(run_srd$doptm,digits=6)))}
else
{print(paste("Optimised (penalised) criterion=",round(run_srd$doptm,digits=6)))
}
print(paste("Maximal thinness=",round(run_srd$thinx,digits=2)))
  
  }

##  if new=FALSE  use the previously saved objects run_srd_save and dsave
  else #  for if(new)
    
  {

     if(is.na(run_srd_save[1])==TRUE){
      note <- c(paste("No previous saved new=TRUE configuration"))  
      message(note)
      
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
     }
    if(identical(d,dsave)==FALSE){
      note <- c(paste("previous data with new=TRUE data does not match current"))  
      message(note)
      
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
    }
    
    run_srd <- run_srd_save   
  }    
##  Configuration of  rectangles on flat plane in run_srd$x and run_srd$y

## now consider  the 3D Y axis 

rot_x <- run_srd$x
rot_y <- run_srd$y
rot_x_ws <- run_srd$x
rot_y_ws <- run_srd$y
front <- FALSE

z<-  matrix( ncol=1,nrow=28)  
   
   ## account for missing Y
   okY <- complete.cases(Y)
   Y <- Y[okY]
   nY <- length(Y)
   f <- f[okY]

w <-  w[okY]

   if(q>1){ 
   d <- d[okY,]}
   else{
     d <- d[okY]
   }

##   obtain z values of the q rectangles, as mean of Y
##   same at each of 4 rectangle corners
## use dsum as indicator as whether inside one of the rectangles.
  dsum <- vector(length=length(Y))

dsum <- 0
   for(i in 1:q){
  if(q>1){   
  YFi <- Y*f*d[,i]  
  Fi  <- f*d[,i]*w
  dsum  <-  dsum +d[,i]}
  else{ YFi <- Y*f*d  
        Fi  <- f*d*w
  dsum  <-  dsum + d}
  
  sumn=sum(YFi)
  sumd=sum(Fi)
  mean<-sumn/sumd
  z[1+(i-1)*4]=mean
  z[2+(i-1)*4]=mean
  z[3+(i-1)*4]=mean
  z[4+(i-1)*4]=mean
  


   }  #for(i 1:q)
  
  ##browser()
  ## q+1 th is the 4 corners of unit square
  #overall mean 
  Yf <- Y*f
  fw <- f*w
   omean <- sum(Yf)/sum(fw)
#  mean of whitespace,  1-dsum is 0,1 indicator of being in white space

dac <- ifelse((dsum >0) == TRUE,0,1)
YFc <- Y*f*dac  
Fc  <- f*dac*w

if(sum(Fc)==0) {
  ws_mean <- omean
  nowhitespace <- TRUE}
else
{nowhitespace <- FALSE
ws_mean <- sum(YFc)/sum(Fc)}
  z[1+(q)*4]=ws_mean
  z[2+(q)*4]=ws_mean
  z[3+(q)*4]=ws_mean
  z[4+(q)*4]=ws_mean


  zmax<-max(z,na.rm=TRUE)
  zmin<-min(z,na.rm=TRUE)
#rescale z
   if(zmax == zmin){ 
                     zmin_scaled <- zmax
                     zmax_scaled <- zmin
                     ws_scaled <-   zmax    }
else{

## rescale Z values (overwrite original z) 
 z<- (z/(zmax-zmin))*0.5
 

zmin_scaled <-min(z,na.rm=TRUE)
zmax_scaled <-max(z,na.rm=TRUE)
ws_scaled <- (ws_mean/(zmax-zmin))*0.5
}
# 28 ?   = 7 * 4 corners, 7= 6 rectangles + 1 border
##browser()
  for (i in 1:28){ 
 #   angles <- c(130,330,45)
    xyz <- c(run_srd$x[i],run_srd$y[i],z[i])
  test <- rotate(xyz,project)
#  run_srd$x[i]<- test[1]
rot_y[i]<- test[2]
rot_x[i]<- test[1]
##  z[i]<- test[3] 
## now for back projection onto whitespace to make "hole" 
xyz <- c(run_srd$x[i],run_srd$y[i],ws_scaled)
test <- rotate(xyz,project)
#  run_srd$x[i]<- test[1]
rot_y_ws[i]<- test[2]
rot_x_ws[i]<- test[1]

  }


 
 ncell <-2**q


## variable to attempt to deal with direction of view and looking through "hole"

##browser()
front <- FALSE
if(nowhitespace==FALSE){
front <-  ifelse(omean < ws_mean,TRUE,FALSE)
}
##}



#------------------------------------------------------------------
## collapse to get cell combination and frequencies
 ##browser()
## Try fix to deal with possible "!" in colnames(d) which upsets ddply. 
## Rename columns as 1, 2,3.... 
## This still outputs correct labels for cells=FALSE (I think, keep an eye)
##  but bombs for cells=TRUE
##for(j in 1:q){
##  colnames(d)[j] <- paste(j)
##}


  df <- as.data.frame( cbind(d,f) )
  ## if cells  need to establish mean values in each cell  
 ##browser() 
#  library(plyr)
#  library(stringr)
  colnames(df)[ncol(df)] <- c("f")
#ddply fussy about ! character in name. Replace with X
colnames(df) <- gsub(pattern="!",replacement="X",x=colnames(df))


#ddply fussy about " " blanks character in name. Replace with _
colnames(df) <- gsub(pattern=" ",replacement="_",x=colnames(df))

 cellf <- ddply(df, colnames(df)[1:(ncol(df)-1)], function(xxx){sum(xxx$f,na.rm=TRUE)}  )

## and group frequencies
##browser() 

groupf <- vector(length=q) 
for(i in 1:q){ groupf[i] <- sum(cellf[,i]*cellf$V1)}

#browser()
##--------------------------------------------------------------

#  if cell represented, need to establish  Y values for each cell

if(cells==TRUE){
  Yf <- Y*f
  fw <- f*w
  dYf <- as.data.frame( cbind(d,fw,Yf) )
## if cells  need to establish mean values in each cell  

 # library(plyr)
# library(stringr)
  colnames(dYf)[ncol(dYf)] <- c("Yf")
  colnames(dYf)[ncol(dYf)-1] <- c("fw")
#ddply fussy about ! character in name. Replace with X

colnames(dYf) <- gsub(pattern="!",replacement="X",x=colnames(dYf))
colnames(dYf) <- gsub(pattern=" ",replacement="_",x=colnames(dYf))

##  browser()
  testYf <- ddply(dYf, colnames(dYf)[1:(ncol(dYf)-2)], function(xx){sum(xx$Yf,na.rm=TRUE)/sum(xx$fw,na.rm=TRUE)}  )
 ##browser()
 ## creates object testYf$V1 of means for each binary combination
 ## rescale Z values (overwrite original z) 
 if(zmax != zmin){testYf$V1<- (testYf$V1/(zmax-zmin))*0.5

 ## order cells for correct 3D view layering
 testYf <- testYf[order(testYf$V1,decreasing=FALSE),] 
 }
 
}
  
#-------------------------------------------------------------------------------------
## call the drawing function

  draw.srd.enc(rot_x,rot_y, rot_x_ws, rot_y_ws, run_srd$e, run_srd$q,rnames,run_srd$cstat,labelcell,colour,labelrectangle,
               run_srd$x,run_srd$y,z, ws_scaled,project,front,nowhitespace,cells,testYf,yvar)
#--------------------------------------------------------------------------------------
##  add z axis in 3D projection

if( project[3]!=0 | project[1]!=0 | project[2]!=0) {
  if( project[3]!=180 | project[1]!=180 | project[2]!=180) 

   xyz<- c(0,0,zmin_scaled)
  origin <- rotate(xyz,project)
  xyz<- c(0,0,zmax_scaled)
  
  end <- rotate(xyz,project)
  if((origin[1] != end[1]) | (origin[2] != end[2]) ) {
  arrows(origin[1],origin[2], end[1], end[2],length = 0.1, angle = 25)
  }
  text(origin[1]+0.8*(end[1]-origin[1]),origin[2] +0.8*(end[2]-origin[2]),yvar,cex=1.1, c(1,0) )
  text(end[1],end[2],paste(signif(zmax,4)),cex=1.0,c(0,0))
  text(origin[1],origin[2],paste(signif(zmin,4)),cex=1.0,c(0,0))
  
  omean <-  (omean/(zmax-zmin))*0.5
## mark axis with an X  for overall mean 
  xyz <- c(-0.02,-.02,omean)
  e <- rotate(xyz,project)
  xyz <- c(0.02,0.02,omean)
  o <- rotate(xyz,project)
  segments(e[1],e[2],o[1],o[2],col="blue")

xyz <- c(-0.02,0.02,omean)
e <- rotate(xyz,project)
xyz <- c(0.02,-0.02,omean)
o <- rotate(xyz,project)
segments(e[1],e[2],o[1],o[2],col="blue")

##  text(0,0,paste(xa,ya,za))
}
fail <- FALSE
output <- list(fail=fail)
return(output)

}
## end of draw.srd
#########################################################################
## Drawing and labelcellling function
## note input x,y: 
##  x,y  rotated  projected coordinates of full rectangles  to mean value of Yvar
##  x_ws, y_ws rotated coordinates of back-projected rectangles onto whitespace
##  x_raw,y_raw unrotated x and y coordinates for no projection.
##  

draw.srd.enc<- function(x, y, x_ws, y_ws, e, q, rnames,cstat,labelcell,colour,labelrectangle,
               x_raw,y_raw,z, ws_scaled,project,front,nowhitespace,cells, testYf,yvar) { 
  par(bg = "#EEEEEE", lwd=1)
  plot.new()                 ## new plotting area
  xmin <- min(x, na.rm=TRUE)
  xmax <- max(x, na.rm=TRUE)
  ymin <- min(y, na.rm=TRUE)
  ymax <- max(y, na.rm=TRUE)
  
  ##browser()

  plot.window(xlim=c(xmin-.05,xmax+.4),ylim=c(ymin-.2,ymax+.2))
 
  ##  plot.window(xlim=c(-0.15,1.0),ylim=c(0,1))
  heading <- " "
  if(labelcell==1){
    heading <- "Cell frequencies n_i"
  }
  if(labelcell==3){
    heading <- "Cell error 100(n_i-a_i)/n"
  }
  if(labelcell==4){
    heading <- "Pearson residuals (n_i-a_i)/sqrt a_i"
  }  
  
  if(labelcell==2){
    heading <- "Cell area a_i "
  }
  ##browser()
  if(cells & is.null(yvar)==FALSE) {heading <-paste(heading, ": Shading intensity of ",yvar) }
  
  m.list <- as.list(1:q)     ##set up list of matrices with x, y co-ordinates
#------------------------------------------------------------
##  if   cells are to be displayed 
if(cells == TRUE)
  {
  ncells <- nrow(testYf)
  
  xpos <- vector(length=ncells)
  ypos <- vector(length=ncells)
for (i in 1:ncells ){
  
  bin <- as.numeric(testYf[i,1:q])

 #browser()
# library(colorspace)
## possibility that no variation in Y
if(testYf$V1[ncells]==testYf$V1[1]) {delta <-  0.8}
else{
 delta <-  (1 - (1-0.8^ncells)*((testYf$V1[i]-testYf$V1[1]) /
                                  (1.5*(testYf$V1[ncells]-testYf$V1[1]))) )
}## browser()
if(colour <= 1) {cellcol <- rgb(red=delta, blue=delta*0.8, green=delta)}
if(colour == 2) {cellcol <- rgb(red=0.8*delta, blue=delta, green=delta)}
if(colour >= 3) {cellcol <- rgb(red=delta, blue=delta, green=delta)}

#firstcall <- (project ==0)
# need to ensure that  this cell combination bin corresponds to
# ordering of cells in cstat[]. ibin is binary ordering in cstat[]

ibin <- sum(2^(which(bin==1)-1)) +1 

if(cstat[4+5*(ibin-1)]>0){
drawcell(x_raw,y_raw,q,project, testYf$V1[i],bin,colour=cellcol) 
}

## elements 1 and 2 of cstat are x,y coordinates of  cell midpoints for labelcells
##  transformed to rotated positions 
  Z<- testYf$V1[i]
  ##browser()
   xyz <- c(cstat[1+5*(ibin-1)],cstat[2+5*(ibin-1)],Z)
   test <- rotate(xyz,project)
xpos[ibin] <- test[1]
ypos[ibin] <- test[2]  




}  #for(i in 1:ncells
  
}

else 
#------------------------------------------------------------
##  if rectangles (subgroups)  are to be displayed 

{  
  
  ## draw whitespace 
  bin <- vector(length=q)
  bin <- rep(0,times=q)
  #bin[1] <- 1
  #bin[2] <- 1

  if(nowhitespace==FALSE & front ==FALSE){
#browser()
    drawcell(x_raw,y_raw,q,project, ws_scaled,bin)}
  
  
  ##browser()
  ## colour==1 SPAN earthy colours 
  if(colour==1) {col<-as.vector(c("#58D6B880","#A7C8CC80","#A9D85E80","#6B785B80","#D8AD7580","#BFFFD680"))}
  
## colours 2 =bright colour
  if(colour==2) {col<-as.vector(c("#0000FF80","#FF000080","#00FF0080","#F0FF0080","#0FF0F080","#F0B0E080"))  ## create list of colours
  }
##colour=3 black and white
if(colour==3) {col<-as.vector(c("#fafafa80","#CCCCCC80","#AAAAAA90","#88888890","#66655590","#00000090"))}
#3 default ordering  o  is as listed 
o <- c(1:q)

zm <- z[ seq(1,4*q,by=4)]
o <- order(zm,decreasing=FALSE)
## colour =4  no see-thoughness
if(colour==4)  {
## ensure  rectangle with largest Z is on top
  
  col<-as.vector(c("#BFFFD6","#58D6B8","#A7C8CC","#A9D85E","#6B785B","#D8AD75"))}
  
  ##  Draw the q rectangles 
  for (ii in 1:q){ 
    i <- o[ii]

    m.list[[i]] <- as.matrix(cbind(x[(4*(i)-3):(4*i)], y[(4*(i)-3):(4*i)]))   ## extract components of x and y matrices into m.list which is list of matrices for each rectangle
  }
  
  for (ii in 1:q){
   i <- o[ii] 
    polygon(m.list[[i]], col=col[i])  #### draw a rectangle with each i_th matrix and choose i_th colour in list col.
  }

if(colour==4) {  
## put  rectangle boundaries  
for (ii in 1:q){
  i <- o[ii] 
  polygon(m.list[[i]], col=NA, lty=c("dashed"))  
}
}


ncells <- 2**q

xpos <- vector(length=ncells)
ypos <- vector(length=ncells)

for (i in 1:ncells ){
  
  
  ## elements 1 and 2 of cstat are x,y coordinates of  cell midpoints
  ##  transformed to rotated positions 
  xm <- cstat[1+5*(i-1)]
  ym <- cstat[2+5*(i-1)]
  if(is.na(ym)==FALSE & is.na(xm)==FALSE){
 ##browser()
    Z <- NA
  for(j in 1:(q)){
    
    j4 <- 4*j
    maxx <- max(c(x_raw[j4-3],x_raw[j4]))
    minx <- min(c(x_raw[j4-3],x_raw[j4]))
    maxy <- max(c(y_raw[j4-3],y_raw[j4-2]))
    miny <- min(c(y_raw[j4-3],y_raw[j4-2]))
    if(xm <   maxx & xm > minx ) {
      if(ym < maxy & ym > miny ) {
        
       Z <- zm[j] 
        
      }}
  
  }
   if(is.na(Z)) { Z <- ws_scaled}
  ##browser()
  xyz <- c(cstat[1+5*(i-1)],cstat[2+5*(i-1)],Z)
  test <- rotate(xyz,project)
  xpos[i]<- test[1]
  ypos[i]<- test[2] 
  
  }
} 

  

## ensure  whitspace in front if ws_mean > omean

if(nowhitespace==FALSE & front==TRUE){
  bin <- vector(length=q)
  bin <- rep(0,times=q)
  drawcell(x_raw,y_raw,q,project, ws_scaled,bin)
}
}  #else for  if(cells)

  
  nabsent <- 0 
  
  
  ##output cell labelcells 

#    cstat is 5 x 2**q matrix of cells statistics
#         1,2  are x,y coordinates of midpoint
#         3  is cell frequency    
#         4  is cell area of fitted configuration
#         5  is cell error   

  cstat<- matrix(cstat,nrow=5,ncol=64)
  ncell<- 2**q
  
  ## compute chi-square
  maxr=0
  resid<- matrix(nrow=ncell,ncol=1)
  absentlabelcell<- matrix(nrow=ncell,ncol=1) 

  for (i in 1:ncell){
    resid[i]=((cstat[3,i]-cstat[4,i]))/sqrt(max(cstat[4,i],0.5))
    if(resid[i] >maxr) {maxr=resid[i]}

  }
  r2 <- resid*resid

  chi<- sum(r2)

  pvalue <-  1- pchisq(chi,ncell-1)

if(cells & is.null(yvar)==FALSE) {txtcol="red"}
else {txtcol="black"}

  
  for (i in 1:ncell)

  {
   # put on labelcells, but only if  cell area >0 
   
    if(cstat[4,i] >0){
    ## add  Error% labelcells 
    if(labelcell == 3){
      text( col=txtcol,xpos[i],ypos[i], round(cstat[5,i],2), adj = c(0.5,0),cex=0.75)
    }
    ##  add cell frequency labelcells 
    if(labelcell == 1){  
      text( col=txtcol,xpos[i],ypos[i], cstat[3,i], adj = c(0.5,0),cex=0.75)
    }
    ##  add cell areas labelcells 
    if(labelcell == 2){  
      text(col=txtcol, xpos[i],ypos[i], round(cstat[4,i],1), adj = c(0.5,0), cex=0.75)
    }
    ## add  residual labelcells 
    if(labelcell == 4){
      text( col=txtcol,xpos[i],ypos[i], round(resid[i],2), adj = c(0.5,0),cex=0.75)
    }
    
#     if(cstat[1,i]>1.0 & cstat[3,i]>0){ 
#       nabsent=nabsent+1
#       if(labelcell == 3){absentlabelcell[nabsent]<-round(cstat[5,i],2)}
#       if(labelcell == 2){absentlabelcell[nabsent]<-round(cstat[4,i],1)}
#       if(labelcell == 1){absentlabelcell[nabsent]<-cstat[3,i]}
#       if(labelcell == 4){absentlabelcell[nabsent]<-round(resid[i],2)}
#       
#       
#     } 
  }

else
  {  if(cstat[3,i]>0){ 
        nabsent=nabsent+1
        if(labelcell == 3){absentlabelcell[nabsent]<-round(cstat[5,i],2)}
        if(labelcell == 2){absentlabelcell[nabsent]<-round(cstat[4,i],1)}
        if(labelcell == 1){absentlabelcell[nabsent]<-cstat[3,i]}
        if(labelcell == 4){absentlabelcell[nabsent]<-round(resid[i],2)}
  }
  }
  }

##browser()

  ## output note that there are absent cells  
  if(nabsent>0){
    text( 0,-0.08, paste("No. of Absent cells=",nabsent,":"), adj = c(0,1), cex=0.75)
    if(labelcell !=0) {
      for (i in 1:nabsent)
      {
        text( 0.25+0.1*i,-0.08, absentlabelcell[i], adj = c(0,1), cex=0.75) 
      }
    }
  }


if(labelcell==4){  
  
  title(heading,paste("E%=",round(e,2), "Chi-sq P=",round(pvalue,4)))}
else
{title(heading,paste("E%=",round(e,2) ) )}


#---------------------------------------------------------------
## output  rectangle labelcells and arrows
##browser()

if(labelrectangle > 0){
  if(labelrectangle== 1){
    # label inside  rectangle
    for (i in 1:q){
 
      if(cells) {txtcol="white"}
      else {txtcol="black"}
 #  put label in bottom left corner     
      k <- 3
      offset <- 0.01
      
      text(col=txtcol, x[4*i-k]+0.01, y[4*i-k]+offset, rnames[i], adj = c(0,0),cex=1.0)
    }}
  if(labelrectangle== 2){
    for (i in 1:q){
      k <- 1
      offset <- 0.04*i
      if(i- round((i/2))*2==0) {
        k <- 0
        offset <- -0.04*i}
      ## arrow on left 
      ##    arrows(xmin-0.1,y[4*i-k]+offset, x[4*i-k]+0.01, y[4*i-k]+offset,length = 0.1, angle = 25)
      ##    text( xmin-.15, y[4*i-k]+offset, rnames[i], adj = c(0,1),cex=1.0)
      arrows(xmax+0.01,y[4*i-k]-offset, x[4*i-k]-0.01, y[4*i-k],length = 0.1, angle = 25)
      text( xmax+0.01, y[4*i-k]-offset, rnames[i], adj = c(0,1),cex=1.0)
    }}
  
}


   
return("Configuration fitted")
}
## end of draw.srd.enc -------------------------------------------------------------------
####################################################################
rotate<- function(xyz,angles) { 
  
  xa=angles[1]
  ya=angles[2]
  za=angles[3]
  
  X=xyz[1]
  Y=xyz[2]
  Z=xyz[3]

COSya=cos(ya*0.017453292)
SINya=sin(ya*0.017453292)
COSxa=cos(xa*0.017453292)
SINxa=sin(xa*0.017453292)
COSza=cos(za*0.017453292)
SINza=sin(za*0.017453292)

XD <- X*(COSya*COSza) + Y*(-COSya*SINza) + Z*(SINya) 

YD  <- X*(SINxa*SINya*COSza+COSxa*SINza) + 
     Y*(-SINxa*SINya*SINza+COSxa*COSza) + 
     Z*(-SINxa*COSya) 

ZD <- X*(-COSxa*SINya*COSza+SINxa*SINza) + 
       Y*(COSxa*SINya*SINza+SINxa*COSza) + 
       Z*(COSxa*COSya)
output <- c(XD,YD,ZD)


}
####################################################################


rotate.srd <- function(d,f=matrix(data=1, ncol=1,nrow=nrow(d)),labelcell=0,new=TRUE, 
              crit=3,colour=1,labelrectangle=1, Y=matrix(data=0, ncol=1,nrow=nrow(d)),movie=FALSE,cells=FALSE,w=matrix(data=1, ncol=1,nrow=nrow(d)))
{
  ##------------------------------------------------------------------------------------------------------
  ## d is data matrix of ones and zeroes
  ## f is frequency count vector (default =1)
  ## labelcell=0, no labelcells (default), 1= cell frequencies, 2= area frequencies, 3= % error
  ## new  re-fits a new SRD  (default TRUE). Otherwise the existing configuration is used 
  ## crit=1 least squares (default), 2=least absolute difference, 3= minus log-likelihood
  ## colour=1  earthy muted colours, 2=bright rainbow colours (default), 3=monochrome
  ## labelrectangle=0 no rectangle labelcell, =1 (default) insided rectangle, 2=arrow pointer
  ## Y vector of  variable for third axis in a 3-D projection, default  zeroes
  ## movie=TRUE for an animation. Default =FALSE for manual rotation with "return" button
  ##---------------------------------------------------------------------------------------------------------
  
  draw.srd(d,f, labelcell=labelcell,new,colour=colour, crit, labelrectangle=labelrectangle,Y,project=c(0,0,0),cells=cells, w)
  
if(movie){

#  library(animation)
  
  note <- c(paste(" Creating animation frames. Please wait. "))  
  print(note)
  

saveGIF({
  
  ani.options(interval=0.2)
  
  renew <- FALSE 
  for (i in 2:90){
  
    ya <-  2*(i-1)
    xa <-  2*(i-1)
    za <-  2*(i-1)

 ##  note in rotation beyond 0,0,0  labelrectangle=0, turns off labelcells   
    draw.srd(d,f,labelcell=labelcell,new=renew,colour=colour, crit, labelrectangle=labelrectangle,Y,project=c(xa,ya,za), cells=cells,w)
    
  }
  
  ## restore the options
  ##ani.options(oopt)
})
}

else {
  
  #  manually step rotations, pausing with brwser()
renew <- FALSE 
for (i in 2:90){
  
  ya <-  2*(i-1)
  xa <-  2*(i-1)
  za <-  2*(i-1)

  
devAskNewPage(ask=TRUE)
draw.srd(d,f, labelcell=labelcell,new=renew,colour=colour, crit, labelrectangle=0,Y,project=c(xa,ya,za),cells=cells,w)
  
}
}



draw.srd(d,f, labelcell=labelcell,new=renew,colour=colour, crit, labelrectangle=labelrectangle,Y,project=c(0,0,0),cells=cells)

}

######################################################################
drawcell <- function(x,y,q,project,ws_scaled,bin,colour="#FFFFFF")
{  

  bintest <- vector(length=q)
  nn <- 2*(q+1)
  xx <- vector(length=nn)
  yy <- vector(length=nn)
  xrot <- vector(length=4)
  yrot <- vector(length=4)
  nn1 <- nn-1 
for(k in 1:(q+1)){
  k4 <- 4*k
  k2 <- k*2
  xx[k2-1] <- x[k4 -3]
  xx[k2] <- x[k4]
  yy[k2-1] <- y[k4-3]
  yy[k2] <- y[k4-2]
}


xx <- sort(xx, method="quick")



yy <- sort(yy, method="quick")



#browser()
for (k in 1:(nn1)){
  
  xxk <- xx[k]
  xxk1 <- xx[k+1]
  xm <- (xxk1 +xxk)*0.5
  
  for (l in 1:(nn1)){
    yyl <- yy[l]
    yyl1 <- yy[l+1]  
    ym <- (yyl1 +yyl)*0.5  

 
  bintest <- NA
  
   for(j in 1:q){
##     browser()
     j4 <- 4*j
     bintest[j] <- 0
     if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
      
        bintest[j] <- 1
        
       }}
     
 
   }
  
  if(identical(bin,bintest) == TRUE) {
 
## filling a  subcell. Extract the unrotated subcell rectangle corners
    ## at z value corresponding 
    toosmalltodraw <- (abs(xxk-xxk1)<0.001) |  (abs(yyl-yyl1)<0.001)
    
##    if(toosmalltodraw == FALSE){

      p1 <- c(xxk1,yyl1,ws_scaled)
      p2 <- c(xxk1,yyl,ws_scaled)
      p3 <- c(xxk,yyl,ws_scaled)
      p4 <- c(xxk,yyl1,ws_scaled)
##  rotate       
      t1 <- rotate(p1,project)
      t2 <- rotate(p2,project)
      t3 <- rotate(p3,project)
      t4 <- rotate(p4,project)
      
      xrot[1] <- t1[1]
      yrot[1] <- t1[2]
      
      xrot[2] <- t2[1]
      yrot[2] <- t2[2]
      
      xrot[3] <- t3[1]
      yrot[3] <- t3[2]
      
      xrot[4] <- t4[1]
      yrot[4] <- t4[2]
#if(l==5){ browser()}      
polygon(xrot, yrot, border=NA, col=colour)    ## fill subcell colour   
#if(l==5){ browser()}
#}

## determine edges, depending on whether adjacent above, below, left, right

#----------------------------------------------------------------
##  on LEFT 
edge_left <- FALSE

if(k == 1){
  # must be subcell on left boundary
  edge_left <- TRUE}
else
{
  #centre of cell on the left
  xm <- (xxk +xx[k-1])*0.5
  ym <- (yyl1 +yyl)*0.5
 bintest <- NA
 for(j in 1:q){
   j4 <- 4*j
   ##     browser()
   bintest[j] <- 0
   if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
       bintest[j] <- 1
     }}
 } 
 
if(identical(bin,bintest)==FALSE){edge_left <- TRUE}
} 

#if(l==5){ browser()}

#----------------------------------------------------------------
##  on RIGHT 

edge_right <- FALSE

if(k == nn1){
  # must be subcell on right boundary
  edge_right <- TRUE}
else
{xm <- (xx[k+2] +xxk1)*0.5 
 ym <- (yyl1 +yyl)*0.5 
 bintest <- NA
 for(j in 1:q){
   j4 <- 4*j
   ##browser()
   bintest[j] <- 0
   if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
       bintest[j] <- 1
     }}
 } 
 if(identical(bin,bintest)==FALSE){edge_right <- TRUE}
} 

#----------------------------------------------------------------
 ##  above 

 edge_above <- FALSE
 
 if(l == nn1){
   # must be subcell on above boundary
   edge_above <- TRUE}
 else
 {xm <- (xxk1 +xxk)*0.5 
  ym <- (yy[l+2] +yyl1)*0.5 
  bintest <- NA
  for(j in 1:q){
    j4 <- 4*j
    ##     browser()
    bintest[j] <- 0
    if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
      if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
        bintest[j] <- 1
      }}
  } 
  if(identical(bin,bintest)==FALSE){edge_above <- TRUE}
 }  
#----------------------------------------------------------------
##  below 
edge_below <- FALSE

if(l == 1){
  # must be subcell on above boundary
  edge_below <- TRUE}
else
{xm <- (xxk1 +xxk)*0.5
 ym <- (yyl +yy[l-1])*0.5
 bintest <- NA
 for(j in 1:q){
   j4 <- 4*j
   ##     browser()
   bintest[j] <- 0
   if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
       bintest[j] <- 1
     }}
 } 
 if(identical(bin,bintest)==FALSE){edge_below <- TRUE}
}  

## Need to now draw the edges 



if(edge_right){
  
  p1 <- c(xxk1,yyl1,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk1,yyl,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col='black')}


if(edge_left){
  
  p1 <- c(xxk,yyl1,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk,yyl,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col='black')}
 


if(edge_above){
  
  p1 <- c(xxk1,yyl1,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk,yyl1,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col='black')}



if(edge_below){
  
  p1 <- c(xxk1,yyl,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk,yyl,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col='black')}

#browser()
##}  #if(toosmalltodraw==F)
} #if(identical(bin,bintest)

}}   #over k,l
if(identical(bin,rep(0,times=q)) ){
##   is occasional bug drawing left edge, attempt to fix
##   by always putting in unit square
p1 <- c(0,0,ws_scaled)
p2 <- c(0,1,ws_scaled)
p3 <- c(1,1,ws_scaled)
p4 <- c(1,0,ws_scaled)
##  rotate       
t1 <- rotate(p1,project)
t2 <- rotate(p2,project)
t3 <- rotate(p3,project)
t4 <- rotate(p4,project)

xrot[1] <- t1[1]
yrot[1] <- t1[2]

xrot[2] <- t2[1]
yrot[2] <- t2[2]

xrot[3] <- t3[1]
yrot[3] <- t3[2]

xrot[4] <- t4[1]
yrot[4] <- t4[2]
#if(l==5){ browser()}      
polygon(xrot, yrot, border="black", col=NA)    ## fill subcell colour   
#if(l==5){ browser()}
}
#browser()
}
##################################################
##end of drawcell
