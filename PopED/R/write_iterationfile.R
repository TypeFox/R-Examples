## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

write_iterationfile <- function(strSearch,iteration,xt,a,x,ni,groupsize,fim,ofv,poped.db){

## #Make sure that this doesn't crash the search
## try 
##     fn = file(poped.db$settings$strIterationFileName, 'w')
##     if((fn == -1)){
##         stop(sprintf('iteration file could not be open'))
##     }

##     strReturn = 'optimalDesign'

## source(poped.db$settings$strIterationFileName) 
##      returnArgs <-  fileparts(poped.db$settings$strIterationFileName) 
## strPath <- returnArgs[[1]]
##  strFunctionName <- returnArgs[[2]]

##     fprintf(fn,'%s <- function()\n\n',strReturn,strFunctionName){
##     fprintf(fn,'%s%s --\n','% -- Auto generated iteration file for the run ',poped.db$settings$modtit)
##     fprintf(fn,'%s #s, #s\n','% -- Created with',getPopedVersionString(),datestr_poped(poped.db$settings$Engine$Type))
##     fprintf(fn,'%s%s at iteration #d\n\n','% -- Created from the search ',strSearch,iteration)


##     fprintf(fn,'%s ----\n\n','% ---- Current optimal design')

##     fprintf(fn,'%s --\n','% -- The current optimal sampling schedule')
##     write_matlab_matrix(fn,'optimalDesign$xt',xt,size(xt,2),size(xt,1))
##     fprintf(fn,'%s --\n','% -- The current optimal covariates')
##     write_matlab_matrix(fn,'optimalDesign$a',a,size(a,2),size(a,1))
##     fprintf(fn,'%s --\n','% -- The current optimal discrete variables')
##     write_matlab_matrix(fn,'optimalDesign$x',x,size(x,2),size(x,1))
##     fprintf(fn,'%s --\n','% -- The current optimal groupsizes')
##     write_matlab_matrix(fn,'optimalDesign$groupsize',groupsize,size(groupsize,2),size(groupsize,1))
##     fprintf(fn,'%s --\n','% -- The current optimal sampling pattern')
##     write_matlab_matrix(fn,'optimalDesign$ni',ni,size(ni,2),size(ni,1))
##     fprintf(fn,'\n%s --\n','% -- The current optimal Fisher Information Matrix')
##     write_matlab_matrix(fn,'optimalDesign$fim',fim,size(fim,2),size(fim,1))
##     fprintf(fn,'%s --\n','% -- The current optimal objective function value')
##     write_matlab_matrix(fn,'optimalDesign$ofv',ofv,size(ofv,2),size(ofv,1))

##     fprintf(fn,'%s #s', strReturn,'}')

##     fclose(fn)
## catch
##     fclose(fn)
## }

## return(list(    fprintf(fn=    fprintf(fn,' %s =' %s )) 
}
