checkInput <- function(errors = NULL, strata = NULL, sampframe = NULL) {
    if (!is.null(strata)) {
        colnames(strata) <- toupper(colnames(strata))
        if (sum(grepl("N", colnames(strata))) < 1) 
            stop("In strata dataframe the indication of population (N) is missing")
        if (sum(grepl("X", colnames(strata))) < 1) 
            stop("In strata dataframe the indication of at least one auxiliary variable (X) is missing")
        if (sum(grepl("STRAT", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of stratum (STRATUM) is missing")
        if (sum(grepl("DOM1", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of at least Ine domain (DOM1) is missing")
        if (sum(grepl("CENS", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of strata to be sampled or censused (CENS) is missing")
        if (sum(grepl("COST", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of interviewing cost in strata (COST) is missing")
        if (sum(grepl("M1", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of at least one mean (M1) is missing")
        if (sum(grepl("S1", toupper(colnames(strata)), fixed = TRUE)) < 
            1) 
            stop("In strata dataframe the indication of at least one standard deviation (S1) is missing")
        if (sum(grepl("M+[0123456789]", toupper(colnames(strata)), 
            perl = TRUE)) != sum(grepl("S+[0123456789]", toupper(colnames(strata)), 
            perl = TRUE)) + 1) 
            stop("In strata dataframe the number of means (Mx) differs from the number of standard deviations (Sx)")
    }
    if (!is.null(errors)) {
        colnames(errors) <- toupper(colnames(errors))
        if (sum(grepl("DOM", toupper(colnames(errors)), fixed = TRUE)) < 
            1) 
            stop("In errors dataframe the indication of domain (DOM) is missing")
        if (sum(grepl("CV", toupper(colnames(errors)), fixed = TRUE)) < 
            1) 
            stop("In errors dataframe the indication of at least one constraint (CV) is missing")
    }
    if (!is.null(errors) && !is.null(strata)) {
        if (sum(grepl("S+[0123456789]", toupper(colnames(strata)), 
            perl = TRUE)) != sum(grepl("CV", toupper(colnames(errors)), 
            fixed = TRUE))) 
            stop("In strata dataframe the number of means and std deviations differs from the number of coefficient of variations in errors dataframe")
		if (nrow(errors) != length(levels(as.factor(strata$DOM1))))
			stop("In the 'errors' dataframe the number of domain values if different from the number of domain values in 'strata' dataframe")
	}
    if (!is.null(sampframe)) {
        colnames(sampframe) <- toupper(colnames(sampframe))
        if (sum(grepl("X", colnames(sampframe))) < 1) 
            stop("In frame dataframe the indication of at least one auxiliary variable (X) is missing")
        if (sum(grepl("Y", colnames(sampframe))) < 1) 
            stop("In frame dataframe the indication of at least one target variable (Y) is missing")
        if (sum(grepl("DOMAINVALUE", colnames(sampframe))) < 
            1) 
            stop("In frame dataframe the indication of the domain (DOMAINVALUE) is missing")
    }
    if (!is.null(sampframe) && !is.null(strata)) {
        if (sum(grepl("S+[0123456789]", toupper(colnames(strata)), 
            perl = TRUE)) != sum(grepl("Y", toupper(colnames(sampframe)), 
            fixed = TRUE))) 
            stop("In frame dataframe the number of target variables differ from the number of means and std deviations in strata dataframe")
        if (sum(grepl("X", toupper(colnames(strata)), fixed = TRUE)) != 
            sum(grepl("X", toupper(colnames(sampframe)), fixed = TRUE))) 
            stop("In frame dataframe the number of auxiliary variables (X) differ from the number of auxiliary variables in strata dataframe")
    }
    if (!is.null(strata) || !is.null(sampframe) || !is.null(errors)) 
        cat("\nInput data have been checked and are compliant with requirements\n") else cat("\nNo input data indicated\n")
}
