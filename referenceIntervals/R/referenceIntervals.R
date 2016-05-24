data(sysdata, envir=environment())

horn.outliers = function (data)
{
#   This function implements Horn's algorithm for outlier detection using
#   Tukey's interquartile fences.

	boxcox = car::powerTransform(data);
	lambda = boxcox$lambda;
	transData = data^lambda;
    descriptives = summary(transData);
    Q1 = descriptives[[2]];
    Q3 = descriptives[[5]];
    IQR = Q3 - Q1;

	out = transData[transData <= (Q1 - 1.5*IQR) | transData >= (Q3 + 1.5*IQR)];
	sub = transData[transData > (Q1 - 1.5*IQR) & transData < (Q3 + 1.5*IQR)];

    return(list(outliers = out^(1/lambda), subset = sub^(1/lambda)));
}

dixon.outliers = function (data)
{
	# This outlier detection method implements Dixon and Dean's algorithm to find
	# only a single outlier, if it exists.

	d = sort(data);
	dixResult = outliers::dixon.test(data);
	pResult = dixResult[[3]];
	result = strsplit(dixResult[[2]], " ");
	if(pResult <= 0.05) {
		out = result[[1]][3];
		if(result[[1]][1] == "highest") {
			sub = data[data < d[length(d)]];
		}
		else {
			sub = data[data > d[1]];
		}
	}
	else {
		out = as.numeric(c());
		sub = data;
	}
	return(list(outliers = out, subset = sub));
}

cook.outliers = function (data)
{
#	This function identifies outliers based on their Cook Distance from a fitted
#	regression.  In this case, the regression is truly just the mean.

	fit = lm(data ~ 1);
	cooks_dist = cooks.distance(fit);

	out = data[as.numeric(names(cooks_dist[cooks_dist > (4/length(cooks_dist)) |
		cooks_dist > 1]))];
	sub = data[as.numeric(names(cooks_dist[cooks_dist <= (4/length(cooks_dist)) &
		cooks_dist <= 1]))];

	return(list(outliers = out, subset = sub));
}

vanderLoo.outliers = function (data)
{
#	This function identifies outliers using Mark van der Loo's Method I algorithm for
#	outlier detection in the extremevalues package.

	result = extremevalues::getOutliers(data, method = "I");
	indices = c(result$iLeft, result$iRight);

	out = data[indices];
	sub = data[!data %in% out];

	return(list(outliers = out, subset = sub));
}

robust = function (data, indices = c(1:length(data)), refConf = 0.95)
{
#	This implements the robust algorithm as outlined in the CLSI document C28-A3c.

    data = sort(data[indices]);
    n = length(data);
    median = summary(data)[[3]];
    Tbi = median;
    TbiNew = 10000;
    c = 3.7;
    MAD = summary(abs(data - median))[[3]];
    MAD = MAD / 0.6745;
    smallDiff = FALSE;
    repeat {
        ui = (data - Tbi) / (c * MAD);
        ui[ui < -1] = 1;
        ui[ui > 1] = 1;
        wi = (1 - ui^2)^2;
        TbiNew = (sum(data * wi) / sum(wi));
        if((abs(TbiNew - Tbi)) < 0.000001){
            break;
        }
        Tbi = TbiNew;
    };
    ui = NULL;
    ui = (data - median) / (205.6 * MAD);
    sbi205.6 = 205.6 * MAD * sqrt((n * sum(((1-ui[ui>-1 & ui<1]^2)^4)*ui[ui>-1 & ui<1]^2)) / 
                (sum((1-ui[ui>-1 & ui<1]^2)*(1-5*ui[ui>-1 & ui<1]^2)) * 
                max(c(1, -1 + sum((1-ui[ui>-1 & ui<1]^2)*(1-5*ui[ui>-1 & ui<1]^2))))));
    
    ui = NULL;
    ui = (data - median) / (3.7 * MAD);
    sbi3.7 = 3.7 * MAD * sqrt((n * sum(((1-ui[ui>-1 & ui<1]^2)^4)*ui[ui>-1 & ui<1]^2)) / 
                (sum((1-ui[ui>-1 & ui<1]^2)*(1-5*ui[ui>-1 & ui<1]^2)) * 
                max(c(1, -1 + sum((1-ui[ui>-1 & ui<1]^2)*(1-5*ui[ui>-1 & ui<1]^2))))));
    
    ui = NULL;
    ui = (data - Tbi) / (3.7 * sbi3.7);
    St3.7 = 3.7 * sbi3.7 * sqrt((sum(((1-ui[ui>-1 & ui<1]^2)^4)*ui[ui>-1 & ui<1]^2)) / 
                (sum((1-ui[ui>-1 & ui<1]^2)*(1-5*ui[ui>-1 & ui<1]^2)) * 
                max(c(1, -1 + sum((1-ui[ui>-1 & ui<1]^2)*(1-5*ui[ui>-1 & ui<1]^2))))));
    
    tStatistic = qt(1 - ((1 - refConf)/2), (n-1));
    margin = tStatistic * sqrt(sbi205.6^2 + St3.7^2);
	robustLower = Tbi - margin;
    robustUpper = Tbi + margin;
    RefInterval = c(robustLower, robustUpper);
    
    return (RefInterval);
}

nonparRI = function(data, indices = 1:length(data), refConf = 0.95)
{
#	Nonparametric calculation of reference interval, written to be a compatible
#	boot statistic function.

	d = data[indices];
    results = c(quantile(d, (1 - refConf)/2, type = 6), quantile(d, 1-((1 - refConf)/2), type = 6));
    
    return (results);
}

refLimit = function(data, out.method = "horn", out.rm = FALSE, RI = "p", CI = "p", 
					refConf = 0.95, limitConf = 0.90){

	cl = class(data);
	if(cl == "data.frame"){
		frameLabels = colnames(data);
		dname = deparse(substitute(data));
		result = lapply(data, singleRefLimit, dname, out.method, out.rm, RI, CI, refConf, limitConf);
		for(i in 1:length(data)){
			result[[i]]$dname = frameLabels[i];
		}
		class(result) = "interval";
	}
	else{
		frameLabels = NULL;
		dname = deparse(substitute(data));
		result = singleRefLimit(data, dname, out.method, out.rm, RI, CI, refConf, limitConf);
	}
	
	return(result);
}

singleRefLimit = function(data, dname = "default", out.method = "horn", out.rm = FALSE, 
						RI = "p", CI = "p", refConf = 0.95, limitConf = 0.90)
{
#	This function determines a reference interval from a vector of data samples.
#	The default is a parametric calculation, but other options include a non-parametric
#	calculation of reference interval with bootstrapped confidence intervals around the
#	limits, and also the robust algorithm for calculating the reference interval with
#	bootstrapped confidence intervals of the limits.

	if(out.method == "dixon"){
		output = dixon.outliers(data);
	}
	else if(out.method == "cook"){
		output = cook.outliers(data);
	}
	else if(out.method == "vanderLoo"){
		output = vanderLoo.outliers(data);
	}
	else{
		output = horn.outliers(data);
	}
	if(out.rm == TRUE){
		data = output$subset;
	}
	outliers = output$outliers;
    n = length(data);
    mean = mean(data, na.rm = TRUE);
    sd = sd(data, na.rm = TRUE);
    norm = NULL;
    
#	Calculate a nonparametric reference interval.    
    if(RI == "n"){
    	
    	methodRI = "Reference Interval calculated nonparametrically";

        data = sort(data);
		holder = nonparRI(data, indices = 1:length(data), refConf);
        lowerRefLimit = holder[1];
        upperRefLimit = holder[2];
        if(CI == "p"){
        	CI = "n";
        }
    }
    
#	Calculate a reference interval using the robust algorithm method.    
    if(RI == "r"){
    
    	methodRI = "Reference Interval calculated using Robust algorithm";
    	
        holder = robust(data, 1:length(data), refConf);
        lowerRefLimit = holder[1];
        upperRefLimit = holder[2];
        CI = "boot";  
    }
    
#	Calculate a reference interval parametrically, with parametric confidence interval
#	around the limits.    
    if(RI == "p"){
	
#		http://www.statsdirect.com/help/parametric_methods/reference_range.htm
#		https://en.wikipedia.org/wiki/Reference_range#Confidence_interval_of_limit

		methodRI = "Reference Interval calculated parametrically";
		methodCI = "Confidence Intervals calculated parametrically";
		
		refZ = qnorm(1 - ((1 - refConf) / 2));
		limitZ = qnorm(1 - ((1 - limitConf) / 2));
		
		lowerRefLimit = mean - refZ * sd;
		upperRefLimit = mean + refZ * sd;
        se = sqrt(((sd^2)/n) + (((refZ^2)*(sd^2))/(2*n)));
        lowerRefLowLimit = lowerRefLimit - limitZ * se;
        lowerRefUpperLimit = lowerRefLimit + limitZ * se;
        upperRefLowLimit = upperRefLimit - limitZ * se;
        upperRefUpperLimit = upperRefLimit + limitZ * se;
        
        shap_normalcy = shapiro.test(data);
        shap_output = paste(c("Shapiro-Wilk: W = ", format(shap_normalcy$statistic, 
        					digits = 6), ", p-value = ", format(shap_normalcy$p.value,
        					digits = 6)), collapse = "");
        ks_normalcy = suppressWarnings(ks.test(data, "pnorm", m = mean, sd = sd));
        ks_output = paste(c("Kolmorgorov-Smirnov: D = ", format(ks_normalcy$statistic,
        					digits = 6), ", p-value = ", format(ks_normalcy$p.value,
        					digits = 6)), collapse = "");
        if(shap_normalcy$p.value < 0.05 | ks_normalcy$p.value < 0.05){
        	norm = list(shap_output, ks_output);
        				
        }
        else{
        	norm = list(shap_output, ks_output);
        }
    }
    
#	Calculate confidence interval around limits nonparametrically.    
    if(CI == "n"){
    	
    	if(n < 120){
    		cat("\nSample size too small for non-parametric confidence intervals, 
    		bootstrapping instead\n");
    		CI = "boot";
    	}
    	else{
    	
    		methodCI = "Confidence Intervals calculated nonparametrically";
    		
    		ranks = nonparRanks[which(nonparRanks$SampleSize == n),];
  		  	lowerRefLowLimit = data[ranks$Lower];
    		lowerRefUpperLimit = data[ranks$Upper];
    		upperRefLowLimit = data[(n+1) - ranks$Upper];
    		upperRefUpperLimit = data[(n+1) - ranks$Lower];
		}
    }
    
#	Calculate bootstrapped confidence intervals around limits.    
	if(CI == "boot" & (RI == "n" | RI == "r")){
	
		methodCI = "Confidence Intervals calculated by bootstrapping, R = 5000";
		
		if(RI == "n"){
			bootresult = boot::boot(data = data, statistic = nonparRI, refConf = refConf, R = 5000);
		}
		if(RI == "r"){
			bootresult = boot::boot(data = data, statistic = robust, refConf = refConf, R = 5000);
		}
    	bootresultlower = boot::boot.ci(bootresult, conf = limitConf, type="basic", index = 1);
    	bootresultupper = boot::boot.ci(bootresult, conf = limitConf, type="basic", index = 2);
    	lowerRefLowLimit = bootresultlower$basic[4];
    	lowerRefUpperLimit = bootresultlower$basic[5];
    	upperRefLowLimit = bootresultupper$basic[4];
    	upperRefUpperLimit = bootresultupper$basic[5];
    }
    
    RVAL = list(size = n, dname = dname, out.method = out.method, out.rm = out.rm, 
    			outliers = outliers, methodRI = methodRI, methodCI = methodCI, 
    			norm = norm, refConf = refConf, limitConf = limitConf,
    			Ref_Int = c(lowerRefLimit = lowerRefLimit, upperRefLimit = upperRefLimit), 
    			Conf_Int = c(lowerRefLowLimit = lowerRefLowLimit, 
    						lowerRefUpperLimit = lowerRefUpperLimit, 
        					upperRefLowLimit = upperRefLowLimit, 
        					upperRefUpperLimit = upperRefUpperLimit));
    class(RVAL) = "interval";
    return(RVAL);
}

print.interval = function (x, digits = 4L, quote = TRUE, prefix = "", ...)
{
	if(class(x[[1]]) == "interval"){
		lapply(x, print.interval.sub);
	}
	else{
		print.interval.sub(x);
	}
}

print.interval.sub = function (x, digits = 4L, quote = TRUE, prefix = "", ...) 
{
    cat("\n");
    cat(strwrap(x$methodRI, prefix = "\t"), sep = "\n");
    cat(strwrap(x$methodCI, prefix = "\t"), sep = "\n");
    cat("\n");
    cat("data:  ", x$dname, "\n", sep = "");
    cat("N: ", x$size, "\n", sep = "");
    if (!is.null(x$refConf)) {
        cat(format(100 * x$refConf), "% Reference Interval", "\n", sep = "");
    }
    if (!is.null(x$limitConf)) {
        cat(format(100 * x$limitConf), "% Confidence Intervals\n", "\n", sep = "");
    }
    if (!is.null(x$out.method)) {
    	cat("Outlier detection method:", x$out.method, "\n");
    }
    if (!is.null(x$outliers) & length(x$outliers) > 0) {
    	if(x$out.rm){
			cat("Removed outliers: ");
		}
		else{
			cat("Suspected outliers: ");
		}
		cat(strwrap(paste(format(x$outliers, digits = 6), collapse = ", ")), sep = "\n");
    }
    else {
    	cat("No outliers detected\n");
    }
    if (!is.null(x$norm)) {
        cat(strwrap(x$norm, prefix = "\n"), "\n\n");
    }
    if (!is.null(x$Ref_Int)) {
        cat("\nReference Interval: ");
        cat(strwrap(paste(format(x$Ref_Int, digits = 6), collapse = ", ")), sep = "\n");
    }
    if (!is.null(x$Conf_Int)) {
        cat("Lower Confidence Interval: ");
        cat(strwrap(paste(format(x$Conf_Int[1:2], digits = 6), collapse = ", ")), sep = "\n");
        cat("Upper Confidence Interval: ");
        cat(strwrap(paste(format(x$Conf_Int[3:4], digits = 6), collapse = ", ")), sep = "\n");
    }
    cat("\n");
    invisible(x);
}

plot.interval = function (x, main = NULL, ...)
{
	original.parameters = par();
	
	if(class(x[[1]]) != "interval"){
		range = max(x[["Conf_Int"]][[4]]) - min(x[["Conf_Int"]][[1]]);
		y_low = min(x[["Conf_Int"]][[1]]) - 0.05 * range;
		y_high = max(x[["Conf_Int"]][[4]]) + 0.05 * range;
	
		plot.new();
		plot.window(xlim=c(0,2), ylim=c(y_low,y_high));
	
		segments(1, min(x[["Ref_Int"]]), 1, max(x[["Ref_Int"]]), col = "red");
		segments(1-0.05, x[["Conf_Int"]][[1]], 1+0.05, x[["Conf_Int"]][[1]], col = "blue");
		segments(1-0.05, x[["Conf_Int"]][[2]], 1+0.05, x[["Conf_Int"]][[2]], col = "blue");
		segments(1-0.05, x[["Conf_Int"]][[3]], 1+0.05, x[["Conf_Int"]][[3]], col = "blue");
		segments(1-0.05, x[["Conf_Int"]][[4]], 1+0.05, x[["Conf_Int"]][[4]], col = "blue");
		axis(1, at=1:1, labels=x[["dname"]]);
		axis(2);
		if(!is.null(main)){
			title(main = main);
		}
		else {
			title(main="Reference Range");
		}
		title(xlab="Parameter");
		title(ylab="Units");
		legend(x="topright", col = c("red", "blue"), lty = 1, inset = c(0, -0.09),
				legend = c("Reference Interval", "Confidence Intervals"), cex = 0.75,
				xpd = TRUE);
		box();
	}
	if(class(x[[1]]) == "interval"){
		numRanges = length(x);
		intervals = unlist(sapply(x, "[", "Conf_Int"));
		labels = unlist(sapply(x, "[", "dname"));
		range = max(intervals) - min(intervals);
		y_low = min(intervals) - 0.05 * range;
		y_high = max(intervals) + 0.05 * range;

		plot.new();
		plot.window(xlim=c(0,numRanges + 1), ylim=c(y_low,y_high));

		for(i in 1:numRanges){
			segments(i, x[[i]]$Ref_Int[1], i, x[[i]]$Ref_Int[2], col = "red");
			segments(i-0.05, x[[i]]$Conf_Int[[1]], i+0.05, x[[i]]$Conf_Int[[1]], col = "blue");
			segments(i-0.05, x[[i]]$Conf_Int[[2]], i+0.05, x[[i]]$Conf_Int[[2]], col = "blue");
			segments(i-0.05, x[[i]]$Conf_Int[[3]], i+0.05, x[[i]]$Conf_Int[[3]], col = "blue");
			segments(i-0.05, x[[i]]$Conf_Int[[4]], i+0.05, x[[i]]$Conf_Int[[4]], col = "blue");
		}
		
		axis(1, at=1:numRanges, labels = FALSE);
		text(x = seq(1, numRanges, by=1), par("usr")[3] - 0.2, labels = labels, 
			cex = 0.75, srt = 90, pos = 1, offset = 2, xpd = TRUE);
		axis(2);
		if(!is.null(main)){
			title(main = main);
		}
		else {
			title(main="Reference Range");
		}
		title(xlab="Parameter");
		title(ylab="Units");
		legend(x="topright", col = c("red", "blue"), lty = 1, inset = c(0, -0.09),
				legend = c("Reference Interval", "Confidence Intervals"), cex = 0.75,
				xpd = TRUE);
		box();
	}
	par(original.parameters[-c(13, 19, 21:23)]);
}