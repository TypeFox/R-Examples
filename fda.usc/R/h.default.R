h.default=function (fdataobj, prob=c(0.025,0.25),len=51, metric = metric.lp,
                                         Ker = "AKer.norm", type.S ="S.NW",...)
{
    #type.S<-deparse(substitute(type.S))
    if (!is.fdata(fdataobj))   fdataobj = fdata(fdataobj)
    if (is.matrix(metric)) {mdist=metric}
    else {mdist=metric(fdataobj,fdataobj,...)}
    n=nrow(fdataobj)
    if (is.function(Ker))    if (body(Ker)==body(AKer.norm))  Ker= "AKer.norm"

    if (type.S=="S.NW") {
            diag(mdist) = Inf
            h0 <- apply(mdist, 1, min, na.rm = TRUE)
            h.max = max(h0)
            h.med = median(h0)
            class(mdist)<-"matrix"
            q.min <- quantile(mdist, probs = prob[1], na.rm = TRUE,type=4)
            q.max = quantile(mdist, probs = prob[2], na.rm = TRUE,type=4)
            h.min = max(c(drop(q.min),h.max))#h.med
            h.max = max(c(drop(q.max), h.max))
            if (Ker== "AKer.norm" ) {
                    h.max = min(q.max, h.max) #antes min
                    h.min = min(q.min, h.med) #antes min
                  }
            h = unique(seq(h.min, h.max, len = len))
            h<-h+h*0.0001
        }
        else if (type.S=="S.KNN") {
            if (len>n) {len=n-2;print("len=nrow-2")}
            h.min = floor(quantile(1:n, probs = prob[1], na.rm = TRUE,
                type = 4))
            h.max = floor(quantile(1:n, probs = prob[2], na.rm = TRUE,
                type = 4))
            h.min = max(2,h.min)
            h.max = max(h.min,h.max)
#            h =seq.int(h.min,h.max,len=len)
            h =unique(floor(seq(h.min,h.max,len=len)))
        }
        else {
        #           stop("Error in type.S argument
                }
        return(h)
 }



