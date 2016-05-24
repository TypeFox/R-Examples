lis.test <-
function (x, y, alternative = c("two.sided", "less", "greater"),method=c("JLMn","Ln","JLn")) 
{
    n <- length(x)
    m <- length(y)
    if (n != m) stop("x and y should have the same size")
    if (n < 4) stop("not enough data")
    if (n >= 201) { stop("sample size must be less than 201")}
    if (min(length(unique(x)),length(unique(y))) < n)  stop("cannot compute p-values with ties")
    test <- match.arg(method)
    if(test=="Ln"){
        ln<-Ln(x,y)
        NN<-as.numeric(TLN[[n]])
        aa<-as.numeric(names(unlist(TLN[[n]])))
    }
    if(test=="JLn"){
        ln <- JLn(x,y)
        NN<-as.numeric(TJLN[[n]])
        aa<-as.numeric(names(unlist(TJLN[[n]])))
    }
    if(test=="JLMn"){
        ln<-JLMn(x,y)
        NN<-as.numeric(TJLMN[[n]])
        aa<-as.numeric(names(unlist(TJLMN[[n]])))
    }
    FF<-sum(NN[aa<ln])/sum(NN)     
    alternative <- match.arg(alternative)   
    if(alternative=="greater"){
        pvalue <- 1-FF
    }
    if(alternative=="less"){
        pvalue <- FF
    }
    if(alternative=="two.sided"){
        if(FF < 0.5 ){ pvalue <- 2*FF } else{ pvalue= 2*(1-FF)}
    }      
    res<-list(p.value=pvalue,sample.estimate=ln,method=test,altenative=alternative)
    print(res)
}
