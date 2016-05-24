ibdhap.barplot <-
function( x,
                            data.type = c("h", "g", "r"),
                            col = c("grey", "red", "white"),
                            position = NA,
                            top = 1,
                            bottom = 0,
                            density = 50,
                            ... )
{

#x is a single column from ibdhap.make.states output

#define parameters based on haplotype/genotype/or reduced data
if(length(data.type)>1){ print("data.type improperly indicated, please see ?ibdhap.summary")
                         return(0)}
else if( is.element("h", data.type)){ no.ibd.ind = 15}
else if( is.element("g", data.type)){no.ibd.ind = 9}
else if( is.element("r", data.type)){no.ibd.ind = 4}
else{ print("data.type improperly indicated, please see ?ibdhap.summary")
      return(0)}


#change x into only holding information about ibd, no ibd, or no call

  x[((x>0)&(x<no.ibd.ind))] <- 1 #some ibd
  x[(is.element(x,no.ibd.ind))] <- 2 #no ibd


#get the segment lengths

lengths<-ibdhap.seg.lengths(x, position)

colors.vec<-make.col.vec(lengths[,1], col)

start.points<-cumul.sum(lengths[,2])
start.points<-c(0, start.points)

end.points<-start.points[2:length(start.points)]     #everything but the first point
start.points<-start.points[-length(start.points)]    #get rid of the last point


plot(NULL, ylim = c(bottom, top), xlim =c(start.points[1], end.points[length(end.points)]), yaxt = 'n', ...)

rect(xleft = start.points,
     xright = end.points,
     ytop = top,
     ybottom = bottom,
     density = density,
     col = colors.vec )

}

