############################################################################
# generics to  "usual"  methods
############################################################################

#if(!isGeneric("plotCLT")){
  setGeneric("plotCLT", function(Tn, ...) standardGeneric("plotCLT"))
#}

## The original code was not working under R 2.8.0 (r45868)
## Error in setMethod("plotCLT", "DiscreteDistribution", function(Tn, k,  : 
##  no existing definition for function "plotCLT"
## Error: unable to load R code in package 'distrTeach'
## Execution halted
## ERROR: lazy loading failed for package 'distrTeach'
