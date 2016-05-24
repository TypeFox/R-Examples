
DefineSubseries =function(x,wdw =min(length(x),30), Overlap =wdw-1){      
#===============================================================================
# FUNCTION : DefineSubseries
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 14/10/2010 11:39:08
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
## Defining a matrix of subseries (column by column) from a series  with a certain window wdw and overlap.   
## Dimensions of this matrix wdw*(NbSubseries_overlap)
## Subseries are encoded row by row.        
#-------------------------------------------------------------------------------
# USAGE: 
# DefineSubseries(x,wdw =min(length(x),30), Overlap =wdw-1)
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# x : a series
# wdw: In case distances are measured between "long" series, it could be advantageous to make use of a windowed semimetric. wdw encode the length of the window in which the semimetric will be computed between the subseries. By default there is no windowing if the length of the series ( =ncol(x1) ) is smaller than or equal to 30, and a windows length of 30 otherwise.
# Overlap : In case a windowing of the series is applied, Overlap determines how the subseries overlap each other. By default, a one-step-sliding distance is computed.
#-------------------------------------------------------------------------------
# DETAILS:

#-------------------------------------------------------------------------------
# VALUE:
  
#-------------------------------------------------------------------------------
# WARNING/NOTE :

#-------------------------------------------------------------------------------
# REFERENCES: 

#-------------------------------------------------------------------------------
# SEE ALSO:

#-------------------------------------------------------------------------------
# EXAMPLES:

#===============================================================================
#===============================================================================

x= as.numeric(x) # if x comes from a dataframe, this ensure it is dealt with as a numeric vetor and not a dataframe with one observation

Decalage = wdw-Overlap       # Default: =1.

NbSubseries= (length(x) - wdw +1) # Nb of subseries in case of Decalage = 1.
NbSubseries_overlap= floor((length(x) - wdw)/(wdw-Overlap))+1 #ceiling(NbSubseries/Decalage ) # Nb of subseries for which the semidistance will be calculated. If Overlap is set to its default value, NbSubseries_overlap = NbSubseries.

if (NbSubseries_overlap !=  (length(x) - wdw)/(wdw-Overlap) +1) {
warning('WARNING: the length of the series (ncol(x)) is not a multiple of the Overlap parameter. The series are thus truncated to the closer multiple of Overlap. ')
}

                       
### Transforming  the series j in subseries
x.Windowed_tmp = matrix(NA,ncol=wdw, nrow=NbSubseries)
                                                                        
for (k in 1:NbSubseries){x.Windowed_tmp[k,]= x[k:(k+wdw-1)]}  
      
### Select subseries depending on the Ovelap parameter
x.Windowed = x.Windowed_tmp[seq(1,nrow(x.Windowed_tmp),by =Decalage),]

## Subseries are returned as rows of the matrix x.Windowed.
return(x.Windowed) 

}


        