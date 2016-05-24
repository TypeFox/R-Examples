print.nPV <-
function(x, ...){
cat("Asymptotic experimental design for predictive values:\n")

cat(x$NSETS, "setting(s) generated from given input parameters:
 Sensitivity (se), specifity (sp), prevalence (prev),
 thresholds of negative and positive predictive values (NPV0, PPV0), 
 intended power for rejecting H0: NPV<=NPV0; H0: PPV<=PPV0 (power),
 true values of NPV and PPV resulting from given sp, se, prev.\n")
print(x$outDAT)
checknPV(x)

cat("Proportion of true positives (propP) and sample sizes (n.) for the\n",
x$NSETS, "sets for NPV and PPV:\n")

print(as.data.frame(x), ...)
invisible(x)
}

