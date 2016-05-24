# Individuals-Methods.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' @title Wrapper function SimulateIndividuals
#' @name simulateIndividuals
#' @description This function simulates individuals as an Individuals object.
#' 
#' Will simulate \code{nb} individuals in receptors fields of a \code{landscape}.
#' @details The Individuals object output includes for each individual the coordinates, the date of birth, the life duration and the intern toxic concentration.
#' @rdname simulateIndividuals-constructor-class
#' @param objectL A Landscape object
#' @param n Number of individuals to simulate
#' @param mintime Start simulation time
#' @param maxtime End simulation time
#' @param dob A vector for the Date Of Birth of each individual 
#' @param life_duration A vector for the life duration of each individual
#' @param toxic_threshold A vector for the intern toxic threshold value leading to death for each individual 
#' @return A S4 \code{Individuals} object
#' @seealso \link{loadIndividuals}
#' @include Class-Landscape.R Class-Individuals.R
#' @export
simulateIndividuals <- function(objectL,n=200,mintime=1,maxtime,dob,life_duration,toxic_threshold) {
  
  res <- new('Individuals')
  res@n = n
  res@mintime=mintime
  res@maxtime=maxtime
  
  # Poisson (aléatoire)
  res@coordinate<-spsample(getSPReceptors(objectL),n,"random",iter=4,bb=objectL@thelandscape@bbox)
  res@coordinate@bbox<-objectL@thelandscape@bbox
  res@xmin<-objectL@xmin
  res@xmax<-objectL@xmax
  res@ymin<-objectL@ymin
  res@ymax<-objectL@ymax
  res@dob=dob
  #res@adult=res@dob
  res@life_duration=life_duration
  res@intern_toxic=matrix(data=0,nrow=n,ncol=maxtime)
  res@toxic_threshold=toxic_threshold
  return (res)
}

#' @title Plot method for Individuals
#' @name Individuals plot 
#' @description Will plot individuals spatial positions.
#' 
# @name plot
#' @param x An Individuals object
# @param y missing (not use)
#' @param add if True the new plot will overlap an already plot image (default False) 
#' @param ... further graphical parameters (\code{par})
#' @param plot.legend plot legend (default TRUE)
#' @rdname Individuals-plot-methods
#' @aliases plot,Individuals,ANY-method
#' @export
setMethod(f="plot",
          signature=c("Individuals",NULL),
          definition=function(x,y,add=F,...,plot.legend=TRUE) {
            if(add==F){
              plot(c(x@xmin,x@xmax),c(x@ymin,x@ymax),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
              rect(x@xmin,x@ymax,x@xmax,x@ymin,lwd=1) 
            }
            plot(x@coordinate,pch=4,lwd=3,col="blue",add=T,...)
            if(plot.legend){ legend(x=x@coordinate@bbox[1,1]+(x@coordinate@bbox[1,2]-x@coordinate@bbox[1,1])*0.3,y=x@coordinate@bbox[2,1],c("individuals"),pch=4,pt.lwd=3,col=c("blue")) }
          }
)

#' @title Plot method for Individuals
#' 
#' @description Will plot individuals positions and state at a time of the simulation.
#' @details "Red" cross means that the individual is dead becouse of toxic exposition.
#' "Green" cross means that the individual is dead by natural death.
#' "Green" to "Red" points give the gradient of toxic concentration before the threshold.
#' 
# @param x An Individuals object
#' @param y time of the simulation to display individuals
# @param add if True this plot will plot overlap an already plot image (default False) 
# @param ... further graphical parameters (\code{par})
#' @rdname Individuals-plot-methods
#' @aliases plot,Individuals,num-method
#' @export
setMethod(f="plot",
          signature=c("Individuals","numeric"),
          definition=function(x,y,add=F,...,plot.legend=TRUE) {
            if(add==F) {
              plot(c(x@xmin,x@xmax),c(x@ymin,x@ymax),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
              rect(x@xmin,x@ymax,x@xmax,x@ymin,lwd=1) 
            }
            p<-colorRampPalette(c("green","red"))
            rainbow<-p(10)
            #rect(min(x@coordinate$x)-50,min(x@coordinate$y)-50,max(x@coordinate$x)+50,max(x@coordinate$y)+50,col="white")
#             if( length(ind<-which(x@dob<=y & x@adult>y & x@intern_toxic[,y]<x@toxic_threshold)) != 0) {
#               plot(x@coordinate[ind],pch=1,col="blue",add=T,...)
#             }
#             if( length(ind<-which(x@dob<=y & x@intern_toxic[,y]<x@toxic_threshold & x@dob+x@life_duration>y)) != 0) {
#               #plot(x@coordinate[ind],pch=16,col="green",add=T,...)
#               for(i in ind) {
#                 points(x@coordinate[i],pch=16,col=rainbow[round((x@intern_toxic[i,y]*length(rainbow))/x@toxic_threshold[ind])+1])
#               }
#             }
#             if( length(ind<-which(x@dob<=y & x@intern_toxic[,y]>=x@toxic_threshold & x@dob+x@life_duration>y)) != 0) {
#               plot(x@coordinate[ind],pch=16,col="red",add=T,...)
#             }
#             if( length(ind<-which(x@dob+x@life_duration<=y & max(x@intern_toxic[,x@mintime:x@maxtime])>=x@toxic_threshold)) != 0) {
#               plot(x@coordinate[ind],pch=4,col="red",add=T,...)
#             }
#             # may show only natural dead
#             if( length(ind<-which(x@dob+x@life_duration<=y & max(x@intern_toxic[,x@mintime:x@maxtime])<x@toxic_threshold)) != 0) {
#               plot(x@coordinate[ind],pch=3,col="white",add=T,...)
#             }
            lt<-getIndividualsLife(x)
            
            if( length(ind<-which(lt[,y] == -1)) != 0 ) {
              plot(x@coordinate[ind],pch=3,cex=1.5,lwd=3,col=rgb(0,180,0,maxColorValue = 255),add=T,...)
            }
            
            if( length(ind<-which(lt[,y] == -2)) != 0 ) {
              plot(x@coordinate[ind],pch=4,cex=1.5,lwd=3,col="red",add=T,...)
            }
            
            if( length(ind<-which(lt[,y] > 0)) != 0 ) {
              for(i in ind) {
                points(x@coordinate[i],pch=16,cex=1.5,lwd=3,col=rainbow[round((x@intern_toxic[i,y]*length(rainbow))/x@toxic_threshold[ind])+1])
              }
            }
            
            if(plot.legend) {
              legend(x=x@coordinate@bbox[1,1]+(x@coordinate@bbox[1,2]-x@coordinate@bbox[1,1])*0.3,y=x@coordinate@bbox[2,1],c("natural death","toxic death","alive"),pch=c(3,4,16),pt.cex=1.5,pt.lwd=3,col=c(rgb(0,180,0,maxColorValue = 255),"red","green"),title="Individuals",bg=rgb(255,255,255,alpha=100,maxColorValue=255))
              legend(x=x@coordinate@bbox[1,1]+(x@coordinate@bbox[1,2]-x@coordinate@bbox[1,1])*0.6,y=x@coordinate@bbox[2,1],c("low","medium","threshold"),pch=c(16,16,16),pt.cex=1.5,pt.lwd=3,col=c(rainbow[1],rainbow[5],rainbow[10]),title="Internal concentration",bg=rgb(255,255,255,alpha=100,maxColorValue=255))
            }
          }
)


#' @title Show a summary of Individuals information
#' @description print a summary of Individuals information
# @name show Individuals
#' @param object An Individuals object
#' @rdname Individuals-show-method
#' @aliases show,Individuals-method
setMethod(f="show",
          signature="Individuals",
          definition=function(object) {
            cat("*** Class Individuals, method Show ***\n")
            cat(sprintf("* Numbers = %i\n",object@n))
            
            nbIndMin <- min(10,length(object@dob))
            if(nbIndMin!=0) {
              cat(sprintf("* %s first Individuals\n",nbIndMin))
              cat("Coordinate :\n")
              for(i in 1:nbIndMin) { print(object@coordinate[i])}
              cat("Date of Birth :\n")
              print(object@dob[1:nbIndMin])
              
              cat("Intern concentrations (date of birth + 10) :\n")
              for(i in 1:nbIndMin) {
                indtimemax <- min(object@dob[i]+10,object@maxtime)
                cat(i," : ", object@intern_toxic[i,object@dob[i]:indtimemax],"\n")
              }
              #print(object@intern_toxic[1:nbIndMin,object@dob[1:nbIndMin]:(object@dob[1:nbIndMin]+10)])
            }
            cat("*** End Show(Individuals) ***\n")  
          }          
)

# @title Print part of Individuals info
# @name print
#' print Individuals information
#' @param x An Individuals object
#' @param ... further arguments passed to or from other methods.
#' @rdname Individuals-print-class
#' @aliases print,Individuals-method 
setMethod(f="print",
          signature="Individuals",
          function(x,...) {
            cat("*** Class Individuals, method print ***\n")
            cat(sprintf("* Numbers = %i\n",x@n))
            
            nbIndMax <- max(10,length(x@dob))
            if(nbIndMax!=0) {
              cat(sprintf("* %s Individuals\n",nbIndMax))
              cat("Coordinate :\n")
              for(i in 1:nbIndMax) { print(x@coordinate[i])}
              cat("Date of Birth :\n")
              print(x@dob[1:nbIndMax])
              
              cat("Intern concentrations (date of birth + 10) :\n")
              print(x@intern_toxic[1:nbIndMax,])
            }
            cat("*** End print(Individuals) ***\n")  
          }
)

# getIndividuals[i] retourne info individue i 
# [i,t] retourne info a temps t
# @name [
#' @title Get an individual infomation
#' @description Get an Individual information
#' @param x An Individuals object
#' @param i individual index
#' @param j time of information
#' @param ... further arguments passed to or from others methods.
#' @param drop logical value (default = TRUE)
#' @return a data.frame
#' @rdname Individuals-get-methods
#' @aliases [,Individuals,numeric,ANY,ANY-method
# @usage x[i,j]
#' @export
setMethod(
  f="[",
  signature=c(x="Individuals",i="numeric",j="numeric",drop="logical"),
  definition=function(x,i,j,...,drop=T) {
   if(i<=x@n && missing(j)) {
     df<-data.frame(id=c(i),coodinate=c(x@coordinate[i]),dob=c(x@dob[i]),life_expectencies=c(x@life_duration[i]),toxic_threshold=c(x@toxic_threshold[i]))
     res<-cbind(df,"intern_toxic"=matrix(x@intern_toxic[i,],ncol=length(x@intern_toxic[i,])))
     return(res)
     #return(list("coodinate"=x@coordinate[i],"dob"=x@dob[i],"concentration"=x@intern_toxic[i,]))
   }
   else {
     if(i<=x@n && !missing(j) && j<=x@maxtime) {
       df<-data.frame(id=c(i),coodinate=c(x@coordinate[i]),dob=c(x@dob[i]),life_expectencies=c(x@life_duration[i]),toxic_threshold=c(x@toxic_threshold[i]),intern_toxic=c(x@intern_toxic[i,j]))
       return(df)
       #return(list("coodinate"=x@coordinate[i],"dob"=x@dob[i],"concentration"=x@intern_toxic[i,j]))
     }
     else {return(data.frame())}
   }
 }
)

#' @title Method to get Individuals Life information
#' @name getIndividualsLife
# @param object An Inidividuals object
#' @param ... other parameters
#' @rdname Individuals-getIndividualsLife-method
#' @exportMethod getIndividualsLife
setGeneric(name="getIndividualsLife",
           def=function(object,...)
             standardGeneric("getIndividualsLife")
)

#' getIndividualsLife
#'
#' @description Get individuals toxic concentration over the simulation time.
#' If intern concentration overtakes the toxic threshold value is "-2", that means the individual is dead because of higher toxic concentration.
#' Otherwise value is "-1" means the individual is dead in natural way. The value "0" means that the individual is not alive yet.
#' @param object An Individuals object
#' @return a matrix indexed by individual ID in rows and by time in columns.
#' @aliases getIndividualsLife,Individuals-method
#' @rdname Individuals-getIndividualsLife-method
setMethod(f="getIndividualsLife",
          signature="Individuals",
          definition=function(object) {
            res<-object@intern_toxic
            for(ind in 1:object@n)
            {
              for(t in object@dob[ind]:min(object@maxtime,(object@dob[ind]+object@life_duration[ind])))
              {
                if(res[ind,t] >= object@toxic_threshold[ind]) { res[ind,t:object@maxtime]<- -2;break}
              }
              if( (object@dob[ind]+object@life_duration[ind]) <= object@maxtime) {if(res[ind,object@dob[ind]+object@life_duration[ind]] != -2) {res[ind,(object@dob[ind]+object@life_duration[ind]):object@maxtime]<- -1}}
            }
            return(res)
          }
)

#' @title Wrapper function loadIndividuals
#' @name loadIndividuals
#' @description Wrapper function to create an Individuals object using SpatialPoints and dataframe.
#' 
#' The SpatialPoints object and the data.frame have to contain the same number of coordinates and rows.
#' @rdname load-Individuals-class
#' @param objectL a Landscape object
#' @param sp a SpatialPoint object (individuals coordinates)
#' @param data a data.frame containing individuals attributes. Rows numbers as individuals ID, columns names as dob (date of birth) | life_duration | toxic_threshold
#' @param mintime Start simulation time
#' @param maxtime End simulation time
#' @examples  
#' \dontrun{
#' # simulate individuals coordinates (SpatialPoints object):
#' coordinates <- spsample(getSPReceptors(land),n=2)
#' df <- data.frame("dob"=c(1,8),"life_duration"=c(20,20),
#'            "toxic_threshold"=c(15,15),row.names = c(1,2))
#' ind <- loadIndividuals(objetL=land,sp=coordinates,data=df,mintime=1,maxtime=60)
#' }
#' @return an \code{\link{Individuals-class}} object
#' @export
loadIndividuals<-function(objectL,sp,data,mintime,maxtime) {
  
  res <- new("Individuals")
  
  res@n=length(sp)
  res@coordinate=sp
  res@xmin=objectL@xmin
  res@xmax=objectL@xmax
  res@ymin=objectL@ymin
  res@ymax=objectL@ymax
  res@dob=data$dob
  res@life_duration=data@life_duration
  #res@adult=data$adult
  res@intern_toxic=matrix(data=0,nrow=res@n,ncol=maxtime)
  res@toxic_threshold=data$toxic_threshold
  res@mintime=mintime
  res@maxtime=maxtime
  
  return(res)
}