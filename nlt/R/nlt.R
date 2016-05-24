`nlt` <-
function (x, f, J, Pred = AdaptPred, neighbours = 1, closest = FALSE, 
    intercept = TRUE, nkeep = 2, trule = "median",verbose=TRUE,do.orig=FALSE,returnall=FALSE) 
{

n <- length(x)
vec <- matrix(0, J, n - nkeep)
deni <- df <- NULL
aveghat <- matrix(0, 1, n)
ghatnat <- NULL
for (i in 1:J) {
	if(verbose){
		cat(i,"...\n")
	}
        v <- sample(1:n, (n - nkeep), FALSE)
        vec[i, ] <- as.row(v)
        deni <- denoiseperm(x, f, pred = Pred, neigh = neighbours, 
            int = intercept, clo = closest, keep = nkeep, rule = trule, 
            per = v,returnall=FALSE)
        aveghat <- aveghat+ as.row(deni)
}

aveghat <- aveghat/J

if(do.orig){
	df <- denoise(x, f, pred = Pred, neigh = neighbours, int = intercept, 
        clo = closest, keep = nkeep, rule = trule,returnall=FALSE)
	ghatnat <- as.row(df)
}

if(returnall){
	return(list(vec = vec, ghatnat = ghatnat, aveghat = aveghat))
}
else{
	return(aveghat)
}

}

