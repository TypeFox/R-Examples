# simulate.r #############################################################################################################
# FUNCTION:       DESCRIPTION:     
#  rHAC			  Returns simulated vectors from general HAC.
##########################################################################################################################

rHAC = function(n, hac){
    Sample = rnacopula(n, hac2nacopula(hac))
    colnames(Sample) = .get.leaves(hac$tree)
    Sample
}