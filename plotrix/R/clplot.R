# CLPLOT
# Jan 2014 - Fixed a semibug in the way min/max cut levels are
# 	generated, also cleaned up some random stuff.
# 20 Oct.  2008 Now plots even pathological data sets without gaps
#  Parameters:
#	levels--	vector of desired cutpoints. Program will sort them.
#	cols--	vector of desired color sequence (strings or numeric references)
#	x, y--		The data, of  course.  
#	'...' 		is intended for arguments to pass to PLOT or LINES calls
#	lty--		plot parameter lty (see par() )
#	showcuts--	Set to TRUE to plot gridlines at color cut levels  
#
clplot<-function(x,y, ylab=deparse(substitute(y)), xlab=deparse(substitute(x)),
levels=seq(min(y)+(max(y)-min(y))/5, max(y)-(max(y)-min(y))/5, length.out=4),
cols=c("black", "blue", "green", "orange", "red"), lty=1, showcuts=FALSE, ...) 
{
if(missing(y)){
	ylab=xlab # no reason to recalc "deparse(substitute(x))"
	y<-as.numeric(x)
	x<-seq(1,length(y))
# should have done this long ago: fix xlabel
	xlab<-'index'
	}
xx<-as.numeric(x)
yy<-as.numeric(y)
levels<-sort(as.numeric(levels))

#  Jan 2014:  Bugfix here: if levels has values outside min or max of data,
# we do NOT want to creat cuts at the min,max values! And this needs to be done
# PRIOR to checking lengths of levels vs colors; and check against cuts,
# not against levels. 
#### new code
cuts<-sort(levels)
if(min(yy)< min(levels) ) cuts<- c(min(yy),cuts)
if(max(yy)> max(levels) ) cuts <- c(cuts, max(yy))
#  set cols length to fit with number of cuts. (number of colors to be used)
if (length(cuts)> length(cols)+1) {
	cat("Warning: not enough colors. Will repeat.\n")
	}
#notice this will also truncate cols if it's longer than number of slices.
cols<-rep(cols,length.out=length(cuts)-1)
#build 'empty' graph 
plot(xx,yy,type='n', xlab=xlab, ylab=ylab, ...)
# initialize the list variable which will hold all modded slices
modslice<-list()
newxx<-xx
newyy<-yy
# newxx/yy will add new points as generated for each slice
for(jj in 1 : (length(cuts)-1))	
		{
		botsl<-(cuts[jj]>newyy)
		topsl<-(cuts[jj+1]<newyy)
		botrl<-rle(botsl)
		botlist<-cumsum(botrl$lengths[-length(botrl$lengths)])
		toprl<-rle(topsl)
		toplist<-cumsum(toprl$lengths[-length(toprl$lengths)])
			# get the matching indices
		thematch<-na.omit(match(toplist,botlist))   #finds indices of botlist I like
			# so the *list vars have the points to use for interpolation
		slxx<-newxx
		slyy<-newyy
			# if thematch is empty, avoid crashing FOR loop by using BREAK
		for(kk in 1:length(thematch))
			{
			if(length(thematch)==0) break
				# use botlist[thematch..] to get interpolation points
				# since that's how thematch was created
				# ny1/2 will be part of new points at cut levels
			ny1<-cuts[jj]
			ny2<-cuts[jj+1]
			# stick in a couple vars for readability
			bot1<-botlist[thematch[kk]]
			bot2<-bot1+1
			newx<-approx(c(newyy[bot1],newyy[bot2]),c(newxx[bot1],newxx[bot2]),ny1)
			nx1<-newx$y
			newx<-approx(c(newyy[bot1],newyy[bot2]),c(newxx[bot1],newxx[bot2]),ny2)
			nx2<-newx$y
			slxx<-c(slxx[1:(bot1+(kk-1)*2)],nx1,nx2,slxx[(bot2+(kk-1)*2):length(slxx)])
			slyy<-c(slyy[1:(bot1+(kk-1)*2)],ny1,ny2,slyy[(bot2+(kk-1)*2):length(slyy)])
			} #end of thematch loop
		# Now replace newxx, newyy with slxx/yy
		newxx<-slxx
		newyy<-slyy			
		} # end of FOR (jj )
for(j in 1: (length(cuts)-1))
	 {
	 slice <-rep(0,length(newxx))
#get cut of data desired -- slicelog used to build NA strings outside cutlevels
	 slicelog<-as.logical(cuts[j]<=newyy & newyy<=cuts[j+1])
	 is.na(slice)<-!slicelog
	 sly<-slice+newyy
	 slx<-newxx
	# identify lengths of data and NAs
	 runs<-rle(!is.na(as.logical(sly)))
	#note: this starts thepos at the END of the first 'length' so I
	#never bang into the x[1],y[1] data point.  
	# Similarly, it stops at beginning of last 'length'.
	 thepos<-0   #restart run length count 
# generally if runs$lengths=1, means all NA, no data in slice
	 if ( length(runs$lengths)>1 )
			{
			for (i in 1:(length(runs$lengths)-1)) 
				{
				thepos<-thepos+runs$lengths[i]
			# "add" interp point at the end of each run
			# whichway chooses which 'side' of interval to interpolate towards
				whichway<-slicelog[thepos]
			#pick correct cut value  - is subslice is going up or down
			#whichcut chooses the cut to use
				whichcut<-as.logical(newyy[thepos]>=cuts[j]&newyy[thepos+1]>=cuts[j])
			#note the interpolation is TO x FROM y
				xint3<-approx(c(newyy[thepos:(thepos+1)]),c(newxx[thepos:(thepos+1)]),cuts[j+whichcut])
				newx <-xint3$y
				newy <-cuts[j+whichcut]
			# adding in "i-1" to the splitpoint to get correct location
				slx<-c(slx[1:(thepos+(i-1))],newx,slx[(thepos+1+(i-1)):length(slx)])
				sly<-c(sly[1:(thepos+(i-1))],newy,sly[(thepos+1+(i-1)):length(sly)])
			}  #end of for (runs) loop; all necessary points have been added
		} # end of if (length(runs)) 
	  modslice<-c(modslice,list(cbind(slx,sly)))
	  names(modslice)<-c(names(modslice)[1:(length(names(modslice))-1)],paste("newsl",j,sep=""))
	} # end of j loop.

# modslice will have one list element for each slice, empty or not
mapply(function(x,y) suppressWarnings(lines(x[,1],x[,2], col=y, type='l', lty=lty,...)), modslice,cols)
#suppression lets me pass '...' to both plot and lines calls, e.g. titles to plot() above
if(showcuts==TRUE) 
	{
	lnx<-c(min(newxx),max(newxx))
	lncut<-cuts[2:(length(cuts)-1)] 
	mapply(function(rrr,sss) lines(lnx,c(rrr,sss),lty='dotted', col='black'),lncut,lncut)
	} #end of showcuts IF
#save interesting stuff 
stuff<-list(xin=x,yin=y,cuts=cuts)
stuff<-c(stuff,modslice)
names(stuff)<-c('xin','yin','cuts',names(modslice))
return(invisible(stuff))
}
