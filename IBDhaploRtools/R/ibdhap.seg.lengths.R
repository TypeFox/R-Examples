ibdhap.seg.lengths <-
function( x, position=NA ){
#x is a single column from ibdhap.make.states output
## NOT INCLUDING THE haplotype/genotype names!,
# thus it is the ibd.states from one pair of individuals (genotypes)
# or set of 4 haplotypes


n.marker<- length(x)  #number of markers

#if positions are given, we use them, otherwise "length" refers to
# number of SNPs
if(is.element(position[1],NA)){position <- 1:(n.marker) }



# obtain a vector of ibd state change points --index where ibd state changes
	change.points<-c(1)

	for(imarker in 2:n.marker)
	{
		prev.val<-x[imarker-1]
		val <- x[imarker]

		if( prev.val!=val)
		{
		   change.points=c(change.points, imarker)
		}

	}

	# tidy up the end of change.points
  if(change.points[length(change.points)]!= n.marker){change.points=c(change.points, n.marker)}

	change.points.pos<-position[change.points]
	seg.lengths<-diff(change.points.pos)
	ibd.state<-x[change.points[1:length(seg.lengths)] ]

  return( as.data.frame(cbind(ibd.state = ibd.state, seg.lengths = seg.lengths)))

	}

