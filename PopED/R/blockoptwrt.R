## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

blockoptwrt <- function(fn,optsw,opt_xt=optsw[2],
                        opt_a=optsw[4],opt_x=optsw[4],
                        opt_samps=optsw[1],opt_inds=optsw[5]){ #Write which design parameters that is optimized with respect to

fprintf(fn,'==============================================================================\n')
fprintf(fn,'Optimization of design parameters\n\n')
if((opt_samps)){
    fprintf(fn,'* Optimize Samples per subject\n')
}
if((opt_xt)){
    fprintf(fn,'* Optimize Sampling Schedule\n')
}
if((opt_x)){
    fprintf(fn,'* Optimize Discrete variables\n')
}
if((opt_a)){
    fprintf(fn,'* Optimize Covariates\n')
}
if((opt_inds)){
    fprintf(fn,'* Optimize Number of individuals per group\n')
}
fprintf(fn,'\n')
return( ) 
}
