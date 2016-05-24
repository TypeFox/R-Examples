print.RAMP=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("Important main effects:", x$mainInd,'\n');
      cat("Coefficient estimates for main effects:", signif(x$beta.m,digits),'\n');
      cat("Important interaction effects:", x$interInd,'\n');
      cat("Coefficient estimates for interaction effects:", signif(x$beta.i,digits),'\n'); 
      cat("Intercept estimate:", signif(x$a0[x$cri.loc],digits),'\n');
     }
