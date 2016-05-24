`arandi` <-
function(cl1,cl2, adjust=TRUE){
    if(length(cl1)!=length(cl2)) stop("cl1 and cl2 must have same length")
    tab.1 <- table(cl1)
    tab.2 <- table(cl2)
    tab.12 <- table(cl1,cl2)
    if(adjust){
                correc <- sum(choose(tab.1,2))*sum(choose(tab.2,2))/choose(length(cl2),2)
                return((sum(choose(tab.12,2))-correc)/(0.5*sum(choose(tab.1,2))+0.5*sum(choose(tab.2,2))-correc) )}
    else{ 
        1+(sum(tab.12^2)-0.5*sum(tab.1^2)-0.5*sum(tab.2^2))/choose(length(cl2),2)
    }
}

