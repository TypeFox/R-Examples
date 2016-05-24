skewness.ade <-
function (x, na.rm = FALSE, w=NULL)
{
#################################################################
wmean<- function (x, weights = NULL)
{
    if (!length(weights))
        return(mean(x, na.rm = na.rm))
    if (na.rm) {
        s <- !is.na(x + weights)
        x <- x[s]
        weights <- weights[s]
    }
    sum(weights * x)/sum(weights)
}
#################################################################
    if (is.matrix(x))
        apply(x, 2, skewness.ade, na.rm = na.rm)
    else if (is.vector(x)) {
            if(!is.null(w) & na.rm){
        w <- w[!is.na(x)]
        }
        if(na.rm)  x <- x[!is.na(x)]

        n <- length(x)
        
        if(is.null(w)) {
        out<- (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
        return(out)
         }
        if(!is.null(w)){

         out<- ( sum(((x - wmean(x, w))*w)^3)/sum(w, na.rm=T) )      /        ( sum(((x - wmean(x,w))*w)^2)/sum(w, na.rm=T) )^(3/2)
         #w<- w/(sum(w)/length(x))
         #out<-   (1/sum(w))*(sum(w^(3/2) * ((x-mean(x))/sd(x))^3  ))

         
         return(out)

        }
    }
    else if (is.data.frame(x))
        sapply(x, skewness.ade, na.rm = na.rm)
    else skewness.ade(as.vector(x), na.rm = na.rm)
}
