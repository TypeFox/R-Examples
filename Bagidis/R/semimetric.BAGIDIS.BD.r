semimetric.BAGIDIS.BD= function(Details1, 
                                Breakpoints1,
                                Details2=Details1, 
                                Breakpoints2=Breakpoints1,
                                NbSubseries =1,
                                p = 2, 
                                wk=NULL, 
                                Param=0.5,                             
                                Evol =FALSE
                             
                             ){  
#===============================================================================
# FUNCTION : semimetric.BAGIDIS
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial
# Sciences, UCLouvain. 
# Last update: 05/05/2012 
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# This function computes the Bagidis semidistance between curves.
# If several curves are provided, it returns a matrix of semidistances.
# A user-friendly interface for this function is given by semimetric.BAGIDIS with the option method='BD'.
#-------------------------------------------------------------------------------
# ARGUMENTS: 
# Details1, Breakpoints 1:  matrices containing the details and breakpoints of series out of a dataset DATA1 containing a set of series of identical length.  
# Details2, Breakpoints 2:  matrices containing the details and breakpoints of series out of a dataset DATA2 containing a set of series of identical length as in DATA1.      
# NbSubseries : in case an evolving (windowed) semidistance must be computed, Nbsubseries gives the  number of data measurements in a windowed segment.  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# VALUE:
# dissimilarity.matrix :  Matrix of semidistances between the nrow(DATA1) series of DATA1 and the nrow(DATA2) series of DATA2. Dimensions: nrow(DATA1) x nrow(DATA2) . 
# dissimilarity.evol :  Array of local matrices of semidistances between the windowed series of DATA1 and DATA2. Dimensions: nrow(DATA1) x nrow(DATA2) x Nb_subseries. Nb_subseries is determined by the three quantities nrow(DATA1), wdw and Overlap.   
#-------------------------------------------------------------------------------
# REFERENCES: 
# Timmermans C., von Sachs R., 2010, BAGIDIS,a New Method for Statistical Analysis of Differences Between Curves with Sharp Discontinuities, discussion paper 1030 of Institute of Statistics, Biostatistics and Actuarial Sciences. 
#-------------------------------------------------------------------------------
# SEE ALSO:
# BUUHWE , BAGIDIS.dist ,  semimetric.BAGIDIS, semimetric.BAGIDIS.TS ...
#-------------------------------------------------------------------------------


 
#-------------------------------------------------------------------------------
# 1. Dimensions
#-------------------------------------------------------------------------------

NbSeries1 =   ncol(Details1)/NbSubseries
NbSeries2 =   ncol(Details2)/NbSubseries

#-------------------------------------------------------------------------------
# 2. Computes partial dissimilarities 
#-------------------------------------------------------------------------------



Partial_ij = numeric(length = NbSubseries*NbSeries2*NbSeries1)
k=1

for (i in 1:NbSeries1){
for (j in 1:NbSeries2){
for (l in 1:NbSubseries){
Partial_ij[k] = BAGIDIS.dist.BD(Details1[,(i-1)*NbSubseries+l],Breakpoints1[,(i-1)*NbSubseries+l],Details2[,(j-1)*NbSubseries+l],Breakpoints2[,(j-1)*NbSubseries+l],p, wk, Param)

k=k+1
} #EO FOR l 
} #EO FOR j
} #EO FOR i

dissimilarity.evol_tmp = array(Partial_ij, dim=c(NbSubseries,NbSeries2,NbSeries1) )
dissimilarity.matrix=t(apply(dissimilarity.evol_tmp,c(2,3),sum))
if (Evol) { 
dissimilarity.evol=aperm(dissimilarity.evol_tmp,c(3,2,1))  }


#-------------------------------------------------------------------------------
#  4. Return values of the function
#-------------------------------------------------------------------------------

if (!Evol){
return(dissimilarity.matrix)}
else{
return(list('dissimilarity.matrix' =dissimilarity.matrix,
              'dissimilarity.evol' =dissimilarity.evol))}

}#EOF

#===============================================================================




