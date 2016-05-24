#plot spectral density
spectR <- function(phylo,method=c("standard")){

##kurtosis
kurtosis.sub <-
    function (x, na.rm = FALSE, method = c("moment"), ...)
{
    
    method = match.arg(method)

    stopifnot(NCOL(x) == 1)

    # Warnings:
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("argument is not numeric or logical: returning NA")
        return(as.numeric(NA))}

    # Remove NAs:
    if (na.rm) x = x[!is.na(x)]

    # Kurtosis:
    n = length(x)
    if (is.integer(x)) x = as.numeric(x)
    if (method == "moment") {
        kurtosis = sum((x-mean(x))^4/as.numeric(var(x))^2)/length(x)
    }
     if (method == "excess") {
        kurtosis = sum((x-mean(x))^4/var(x)^2)/length(x) - 3
    }

    if (method == "fisher") {
        kurtosis = ((n+1)*(n-1)*((sum(x^4)/n)/(sum(x^2)/n)^2 -
            (3*(n-1))/(n+1)))/((n-2)*(n-3))
    }

    # Return Value:
    kurtosis
}


#skewness
skewness <- function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
	if (na.rm) x <- x[!is.na(x)] 
	n <- length(x)
     (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
	}
    else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}

#peak height
peak_height <- function(x)
{
	kernelG <- function(x, mean=0, sd=1) dnorm(x, mean = mean, sd = sd)
	dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096,
                from = min(x) - 3*sd, to = max(x) + 3*sd, adjust = 1,
                ...) {
 	if(has.na <- any(is.na(x))) {
    	x <- na.omit(x)
    if(length(x) == 0)
        stop("no finite or non-missing data!")
  }
  	sd <- (if(is.numeric(bw)) bw[1] else bw(x)) * adjust
		X <- seq(from, to, len = n)
			M <- outer(X, x, kernel, sd = sd, ...)
		structure(list(x = X, y = rowMeans(M), bw = sd,
                 call = match.call(), n = length(x),
                 data.name = deparse(substitute(x)),
                 has.na = has.na), class =  "density")
}
integr <- function(x, f)
{
       
       # var is numeric
       if (!is.numeric(x))
       {
              stop('The variable of integration "x" is not numeric.')
       }

       # integrand is numeric
       if (!is.numeric(f))
       {
              stop('The integrand "f" is not numeric.')
       }

       # length(var)=length(int)
       if (length(x) != length(f))
       {
              stop('The lengths of the variable of integration and the integrand do not match.')
       }

      # get lengths of var and integrand
       n = length(x)

       # trapezoidal integration
       integral = 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))

       # print definite integral
       return(integral)
}
	d <- dens(log(m))
	dsc <- max(d$y/integr(d$x,d$y))
	
	#print peak height
	return(dsc)
}
		
	if(method=="standard"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=F)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=1)

	#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
		
		#get principal eigenvalue
		max(m) -> principal_eigenvalue
					
		#get kurtosis
		kurtosis.sub(m) -> kurtosis
		
		#get skewness
		skewness(m) -> skewness
		
		#get peak height
		peak_height(m) -> peak_height
		
		#output
		res<-list(eigenvalues=e$values, principal_eigenvalue=principal_eigenvalue,asymmetry=skewness,peakedness1=kurtosis,peakedness2=peak_height, eigengap=eigenGap[,1])
	}
	
	if(method=="normal"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=T)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=0)
	
		#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
				
		#get principal eigenvalue
		max(m) -> principal_eigenvalue
					
		#get kurtosis
		kurtosis.sub(m) -> kurtosis
		
		#get skewness
		skewness(m) -> skewness
		
		#get peak height
		peak_height(m) -> peak_height
		
		#output
		res<-list(eigenvalues=e$values, principal_eigenvalue=principal_eigenvalue,asymmetry=skewness,peakedness1=kurtosis,peakedness2=peak_height,eigengap=eigenGap[,1])			
			}

	class(res)	<- "spectR"
	return(res)					
}


