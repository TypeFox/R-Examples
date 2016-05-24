H.sampler <-
function (x = "community matrix (spp=col,obs=row)", n = "sample size vector", 
    nit = "number of iterations to use", base = exp(1), corr = FALSE, p = NULL, method = "Shannon") 
{

    if (n == "sample size vector"){n=1}else{}
    if (nit == "number of iterations to use"){nit=1000}else{}
    if (is.vector(x)==TRUE){
    				if(length(p)==0)  p <- x/sum(x)
    				nrowx<-1
    				sv <- 1:length(x)
    				} else {
    						if(length(p)==0) {p = x * 0
    						for (i in 1:nrow(x)) {
        					p[i, ] = x[i, ]/apply(x, 1, sum)[i]}}
        					nrowx<-nrow(x)
        					sv = 1:ncol(x)}
    
    out = array(NA, c(nrowx, length(n)))
    
    

    for (h in 1:nrowx) {				# loop over communities (rows)
        iout = numeric(length(n))
        for (i in 1:length(n)) {
            jout = numeric(nit)
            if (n[i] == 0) {
                iout[i] = 0
            }
            else {
                for (j in 1:nit) {
                  if(is.vector(x)==TRUE) {obs = sample(sv, n[i], replace = TRUE, prob = p)} else {
                  	obs = sample(sv, n[i], replace = TRUE, prob = p[h,])}
                  	
                  obs = count(sv, obs)
                  
                  if (method=="Shannon") jout[j] = Hs(obs)
                  if (method=="Gene diversity") jout[j] = Gd(obs)
                }
                iout[i] = mean(jout)
            }
        }
        out[h, ] <- iout
    }
    names = numeric(length(n))
    for (i in 1:length(n)) {
        names[i] = paste("N", n[i], sep = "")
    }
    colnames(out) <- names
    return(out)
}

