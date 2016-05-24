# 1) find nearby peaks in the peak_data data frame close to the i-th feature in 
# the target data frame
# 2) the peak_data object is a data frame with chromosome name (chr) and
# starting/ending position of feature (start/end) 
# 3) the target object has the same format requirement as peak_data
# 4) n_limit is the maximum number of nearby features in peak_data defined as 
# close to the i-th feature in the target data frame
# 5) extend defines the scope of the genomic regions to look for nearby
# features within a window centered on the i-th feature in target and 
# expands extend nt to each direction
get_nearby<-function(i,peak_data,target,n_limit=20,extend=100000) 
{
  pos=target[i,]
  repeat
  {
    ids=which(peak_data$chr==pos$chr & peak_data$start>pos$start-extend & 
                peak_data$end<pos$end+extend) 
    if (length(ids)>n_limit & extend>100) {extend=extend/2} else {break}
  }
  ids
}
