#Functions

#For Tribometer
# Use direct data from tribometer
# Delete all lines till header data. Change Friction Force to FrictionForce and delete [m] and [N]

#For Profilometer
# All files should be csv
# Header should be Distance and Profile
# Seperation is , and text is ""
# Sample File
# "Distance","Profile"
# "0,00050","-0,08"
# "0,00100","-0,09"

#setwd(paste(path.package("trib"),"/exdata/profile",sep=""))


#' \code{trib.names} returns the names of .csv
#' files in working directory and stores table values to a global variabla All,
#' filenames to global variabla filenames and returns truncated names of files.
#'
#' @param x Numeric value length of file name for truncation. Truncated file
#'   name is given at output for ease of use.
#' @param DataFF logical scalar, checking whether to use CSV style for
#'   Tribometer of Profilometer
#' @param DEncode is a string that specifies file encoding for import
#'
#' @return returns truncated file names, also generates list All and filenames
#'   vector
#'
#' @section Warning:
#' If wrong type is used at DataFF \code{read.csv} will have problems
#'
#' @examples
#' setwd(paste(path.package("trib"),"/exdata/tribometer",sep=""))
#'
#' trib.names(2)
#'
#' setwd(paste(path.package("trib"),"/exdata/profile",sep=""))
#'
#' trib.names(2,DataFF=FALSE)
#'



#Get sample names and import CSV

trib.names <- function(x,DataFF=TRUE,DEncode="UTF-8") {   #give file name length returns names of files generates Global var All, l_All and filenames
  filenames <- list.files(pattern="*.csv")

  if (DataFF){
  All <- lapply(filenames,function(i){
    read.csv(i, sep = "",quote = "\"", dec = ",", fill = TRUE, comment.char = "",encoding=DEncode)
  })
  l_All <- length(All)

}
  else {
    filenames <- list.files(pattern="*.csv")
    All <- lapply(filenames,function(i){
      read.csv(i, sep = ",",quote = "\"", dec = ",", fill = TRUE, comment.char = "",encoding=DEncode)
    })
    l_All<-length(All)
  }
  namesAll<-substr(filenames,start=1,stop=x)      #set column names
tribOut <- list(data=All,names=namesAll)

  return(tribOut)
}

#' \code{trib.CoF} gets file name length and returns coeficient of friction in
#' bar plot generates global var CF
#'
#' @section Warning: Before using this function please set your working
#'   directory to the directory that contains .csv files from tribometer. File
#'   structure should be as follows
#'   \tabular{ccccc}{
#'
#'   Time  \tab  Distance   \tab  Laps  \tab  mu  \tab   FrictionForce \cr
#'   0,100 \tab 0     \tab      0,00000  \tab 0,18478 \tab      0,55435
#'   }
#'
#'   You can check example datasets at package files exdata/tribometer
#'
#' @param x Numeric value length of file name for truncation. Truncated file
#'   name is given at output for ease of use.
#' @param DATAOUT logical scalar checking whether give data out or not
#' @return Plots a bar plot of Coeficient of Friction generated from all files
#'   in working directory and generates a vector called CF that contains Average
#'   Coefficient of Friction values
#'
#'   Set work directory to file containing example data for Tribometer directory
#'   should be
#' @examples
#' setwd(paste(path.package("trib"),"/exdata/tribometer",sep=""))
#' trib.CoF(2)
#' tribDATA<-trib.CoF(2,DATAOUT = TRUE)


#plot Average Coeficient of Friction
trib.CoF <- function(x,DATAOUT=FALSE) { #gets file name length input returns coeficient of friction in bar plot generates global var CF

  tribAll<-trib.names(x)
  All<-tribAll$data
  namesAll<-tribAll$names
  l_All<-length(All)
  CF <- vector(mode = "numeric")     #generate vector for values
  for (i in 1:l_All){

    CF[i]<-mean(All[[i]][[4]])}             #get all coefficient of friction data into CF

  names(CF)<-namesAll
  barplot(CF,main="Average Friction Coefficient",
          xlab="Samples",ylab="CoF",xlim=c(0,length(CF)*1.25+1), ylim = c(0,1.2*max(CF)),col=c(1:l_All))    #generate bar plot
  legend("topright", legend=signif(CF,6),col=c(1:l_All),pch=16,cex=0.8)                                 #add legends

  if (DATAOUT){
    return(CF)
  }

  }



#' \code{trib.FF} gets file name length, filterValue and Friction vs Distance
#' combined graph and/or single graphs for data
#'
#'
#' @section Warning: Before using this function please set your working
#'   directory to the directory that contains .csv files from tribometer. File
#'   structure should be as follows
#'   \tabular{ccccc}{
#'
#'   Time  \tab  Distance   \tab  Laps  \tab  mu  \tab   FrictionForce \cr
#'   0,100 \tab 0     \tab      0,00000  \tab 0,18478 \tab      0,55435
#'   }
#'
#'   You can check example datasets at package files exdata/tribometer
#'
#' @param x Numeric value length of file name for truncation. Truncated file
#'   name is given at output for ease of use.
#' @param filterValue Numeric value for setting for rolling average from package
#'   \code{zoo}
#' @param CombinedGraph logical scalar checking whether to plot combined graph
#'   of all data or not
#' @param SingleGraph logical scalar checking whether to plot single graphs of
#'   all data or not
#'
#' @return Plots Friction vs Distance scatter plots with given filter value if
#'   CombinedGraph and or SingleGraph options are set, else will return nothing
#'
#' Set work directory to file containing example data
#' for Tribometer directory should be
#' @examples
#' setwd(paste(path.package("trib"),"/exdata/tribometer",sep=""))
#'
#' trib.FF(2,100)
#'
#' trib.FF(2, 100, CombinedGraph=FALSE)
#'
#' trib.FF(2, 100, SingleGraph=FALSE)

#plot tribometer plots
trib.FF <- function(x,filterValue,CombinedGraph=TRUE,SingleGraph=TRUE) { #file name length x and gets filter value y , Combined Graph and Single Graph

  tribAll<-trib.names(x)
  All<-tribAll$data
  namesAll<-tribAll$names
  l_All<-length(All)

  fltN <- filterValue                           #filter value
  flt1 <- rep(1/fltN,fltN)
  flt2 <- rep(1/(fltN+1),fltN+1)

  maxD<-vector(mode = "numeric")
  for (i in 1:l_All){maxD[i]<-max(All[[i]]$Distance)}
  xlim<-c(0,max(maxD)*1.25)



  #All
  if (CombinedGraph){
  plot(100,xlab='Distance (m)',ylab = 'Friction Force (N)',xlim=xlim,ylim=c(0,2))

  for (i in 1:l_All){

    points(filter(All[[i]]$Distance,flt1, sides=1),filter(All[[i]]$FrictionForce, flt1, sides=1),pch='.',col=i)
    points(filter(All[[i]]$Distance,flt1, sides=1)[length(filter(All[[i]]$Distance,flt1, sides=1))],filter(All[[i]]$FrictionForce, flt1, sides=1)[length(filter(All[[i]]$FrictionForce, flt1, sides=1))],pch=i,col=i,cex=1)  #Added last points for black and white

  }
  legend("topright", legend= namesAll,col=c(1:l_All),pch=c(1:l_All),cex=.5)
  }
  #Plot each alone
  if (SingleGraph){
  for (i in 1:l_All){

    plot(filter(All[[i]]$Distance,flt1, sides=1),filter(All[[i]]$FrictionForce, flt1, sides=1),pch='.',col=i,xlab='Distance (m)',ylab = 'Friction Force (N)',xlim=xlim,ylim=c(0,2))
    legend("topright", legend= namesAll[i],col=c(1:l_All),pch=16,cex=.5)

  }
  }
}


#'\code{trib.profile} gets file name length, filterValue and combined graph,
#'single graphs, Wear Track Volume bar plot and full profile for data and
#'generates WTArea global variable
#'
#'
#'@section Warning:
#'
#'  Before using this function please set your working directory to the
#'  directory that contains .csv files from profilometer File structure should
#'  be as follows: Seperation is \code{,} and text is \code{""}
#'  \tabular{ccc}{
#'  "Distance" \tab , \tab "Profile" \cr
#'  "0,00050"\tab , \tab "-0,08"
#'  }
#'
#'  You can check example datasets at package files exdata/profilometer
#'
#'  This version doesn't support linear wear tracks. You can still use
#'  \code{trib.profile} but you have to set r to any numeric value and use
#'  WTArea vector to calculate Wear Track Volume
#'
#'@param x Numeric value length of file name for truncation. Truncated file name
#'  is given at output for ease of use.
#'@param filterValue Numeric value for setting for rolling average from package
#'  \code{zoo}
#'@param r Numeric value for setting radius (in mm) of wear track that will be used
#'  for generating Wear Track Volume. This version doesn't support reciprocating
#'  tests so to generate Wear Track Volume please use WTArea vector instead and
#'  give 1 for r value
#'@param CombinedTrack logical scalar checking whether to plot combined graph of
#'  all data or not
#'@param SingleTrack logical scalar checking whether to plot single graphs of
#'  all data or not
#'@param FullProfile logical scalar checking whether to plot full profiles of
#'  each data in a seperate file or not
#'@param WTVBar logical scalar checking whether to plot bar graph of Wear Track
#'  Volume of all data
#'
#'@return Plots profile scatter plots with given filter value if CombinedTrack,
#'  SingleTrack, FullProfile or a bar plot if WTVBar options are set, else will
#'  return nothing. Generates Wear Track Area (WTArea) global variable
#'
#' Set work directory to file containing example data
#' for profilometer directory should be
#' @examples
#' setwd(paste(path.package("trib"),"/exdata/profile",sep=""))
#'
#' trib.profile(2,100,3)
#'
#' trib.profile(2, 100, 3, CombinedTrack=FALSE)
#'
#' trib.profile(2, 100, 3, SingleTrack=FALSE)
#'
#' trib.profile(2, 100, 3, FullProfile=TRUE )
#'
#' trib.profile(2, 100, 3, WTVBar=FALSE)


#Profilometer

trib.profile <- function(x,filterValue,r,CombinedTrack=TRUE,SingleTrack=TRUE,FullProfile=FALSE,WTVBar=TRUE) { #get file name length x and filter value y, plot combined wear track, plot single wear track, plot Full profile, plot bar graph of Wear Track Volume

tribAll<-trib.names(x,DataFF = FALSE)
All<-tribAll$data
namesAll<-tribAll$names
l_All<-length(All)
filenames<-namesAll

PosProfile <- data.frame(P1 = numeric(l_All),P2 = numeric(l_All))         #Wear Track dataframe for plot
NicePosProfile <- data.frame(P1 = numeric(l_All),P2 = numeric(l_All))      #Averaged out
WTArea <- matrix(,nrow = 1, ncol=l_All)         #Wear Track matrix for plot
fltN <- filterValue                               #lag
flt1 <- rep(1/fltN,fltN)          # Filter

for (i in 1: length(filenames)){
  EndV<-length(All[[i]]$Profile)-1        # Set End point
  ProfileFiltered<- filter(All[[i]]$Profile,flt1,sides=1) #filter profile

  MinF <- min(ProfileFiltered,na.rm=TRUE)   #Find filtered min
  PMinF <- match(MinF,ProfileFiltered)      # Get position of filtered min
  RL <-PMinF-500 # Left Range
  # Get left maximum from minimum
  MaxF1 <-max(ProfileFiltered[RL:PMinF],na.rm=TRUE)  # Find max value on filtered
  PMaxF1 <- match(MaxF1,ProfileFiltered[RL:PMinF])+RL  #Find position of maximum value in filtered vector

  Max1<-max(All[[i]]$Profile[RL:PMinF],na.rm=TRUE)      #Find maximum value on actual
  PMax1 <- match(Max1,All[[i]]$Profile[RL:PMinF])+RL


  PosProfile$P1[i]<-PMax1

  # Getting maximum value near filterd maximum value
  # Check for PMaxF1-PMax1 if it's 0.03 near use that value for Maximum

  while (All[[i]]$Distance[PMaxF1]-All[[i]]$Distance[PMax1]>0.03) {
    PMax1<-PMax1+1;
    Max1 <-max(All[[i]]$Profile[PMax1:PMinF],na.rm=TRUE);
    PMax1<- match(Max1,All[[i]]$Profile[PMax1:PMinF])+PMax1}

  # Get right maximum from minimum
  RR<- PMinF+400 #Right Range
  MaxF2 <-max(ProfileFiltered[PMinF:RR])  # ,na.rm=TRUE Find max value on filtered

  PMaxF2 <- match(MaxF2,ProfileFiltered[PMinF:RR])+PMinF  #Find position of maximum value in filtered vector

  Max2<-max(All[[i]]$Profile[PMinF:RR],na.rm=TRUE)      #Find maximum value on actual
  PMax2<- match(Max2,All[[i]]$Profile[PMinF:RR])+PMinF                # Get positon of maximum value on actual
  PosProfile$P2[i]<-PMax2


  # Getting maximum value near filterd maximum value
  # Check for PMaxF1-PMax1 if it's 0.03 near use that value for Maximum

  while (All[[i]]$Distance[PMax2]-All[[i]]$Distance[PMaxF2]>0.03) {
    PMax2<-PMax2+1;
    Max2 <-max(All[[i]]$Profile[PMinF:PMax2],na.rm=TRUE);
    PMax2<- match(Max2,All[[i]]$Profile[PMinF:PMax2])+PMinF

  }

  Min<-min(All[[i]]$Profile[PMax1:PMax2])      #Find minimum value of curve
  AUCc <- sum(diff(All[[i]]$Distance[PMax1:PMax2])*rollmean(All[[i]]$Profile[PMax1:PMax2]-Min,2))     #Calculate area under curve to 0 use min to have all values greater than zero

  xL <- c(All[[i]]$Distance[PMax1],All[[i]]$Distance[PMax2])   # A vector with start and end values for x
  yL <- c(All[[i]]$Profile[PMax1],All[[i]]$Profile[PMax2])     # A vector with start and end values for y
  AUCl <- sum(diff(xL)*rollmean(yL-Min,2))                  #Calculate area under line
  AUC <- (AUCl-AUCc)/1000                                         #Calculate area between in mm


  if (FullProfile){
    plot(All[[i]]$Distance, All[[i]]$Profile,pch='.',xlab='Distance mm',ylab='Profile um') # Plot Full Data
    legend("topright", substr(filenames[i], start = 1, stop = 2),cex=.5)    # Add filename legend
    points(filter(All[[i]]$Distance,flt1,sides=1),filter(All[[i]]$Profile,flt1,sides=1),col='blue',pch='.')
    points(All[[i]]$Distance[PMax1],All[[i]]$Profile[PMax1],col='red')
    points(All[[i]]$Distance[PMax2],All[[i]]$Profile[PMax2],col='red')
  }

  if (SingleTrack){
    plot(All[[i]]$Distance[PMax1:PMax2], All[[i]]$Profile[PMax1:PMax2],pch='.', xlab='Distance (mm)',ylab='Y?kseklik (um)')
    polygon(All[[i]]$Distance[PMax1:PMax2],All[[i]]$Profile[PMax1:PMax2],col="green")
    legend(All[[i]]$Distance[(PMax2-PMax1)/2+PMax1],All[[i]]$Profile[PMax1/2],legend=c("Area (mm)",AUC), cex=.5)    #Write value
    legend("bottomright", legend=substr(filenames[i], start = 1, stop = 2),cex=.5)
  }

  WTArea[i]<- AUC

}

r <- r #mm radius
WTVol <- WTArea*pi*r*2 # Wear Track Formula
  if (WTVBar){
    barplot(t(as.matrix(WTVol)), beside = TRUE, main = "Wear Track Volume", xlab = "Sample", ylab = "Volume (mm3)",col=c(1:l_All), ylim = range(0,max(WTVol)*1.4),names.arg = colnames(WTVol) )

    for (i in 1: length(filenames)){
      text(i+0.5,WTVol[1,i]*1.3,signif(WTVol[1,i],4))
    }
  }

  if (CombinedTrack){

  Nicer <-mean(All[[1]]$Distance[PosProfile$P1]+All[[1]]$Distance[PosProfile$P2])

  plot(1,xlab='Distance (mm)',ylab='Profile (?m)',xlim=c(Nicer/2-0.2,Nicer/2+0.2),ylim=c(-2,0.5),type = "l")  #

  for (i in 1: length(filenames)){
    Avg<-(All[[1]]$Distance[PosProfile[i,1]]+All[[1]]$Distance[PosProfile[i,2]])/2
    Nicer2<-Nicer/2-Avg

    lines(All[[i]]$Distance[PosProfile[i,1]:PosProfile[i,2]]+Nicer2, All[[i]]$Profile[PosProfile[i,1]:PosProfile[i,2]]-max(All[[i]]$Profile[PosProfile[i,1]:PosProfile[i,2]]),col=i, pch='.',xlab='Distance mm',ylab='Profile ?m') # Plot Full Data

    points(All[[i]]$Distance[PosProfile[i,1]:PosProfile[i,2]]+Nicer2, All[[i]]$Profile[PosProfile[i,1]:PosProfile[i,2]]-max(All[[i]]$Profile[PosProfile[i,1]:PosProfile[i,2]]),col=i, pch='.',xlab='Distance mm',ylab='Profile ?m') # Plot Full Data
    points(All[[i]]$Distance[PosProfile[i,2]]+Nicer2, All[[i]]$Profile[PosProfile[i,2]]-max(All[[i]]$Profile[PosProfile[i,1]:PosProfile[i,2]]),col=i, pch=i,xlab='Distance mm',ylab='Profile ?m') # Plot last point for black and white
    }

  legend("topright", legend= substr(filenames, start = 1, stop = 2),col=c(1:l_All),pch=c(1:l_All),cex=.5)    # Add filename legend
  }


}


