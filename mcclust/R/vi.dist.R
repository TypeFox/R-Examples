`vi.dist` <-
function(cl1,cl2,parts=FALSE, base=2){ # wenn parts=TRUE, werden die Komponenten der VI ebenfalls berechnet
    if(length(cl1) != length(cl2)) stop("cl1 and cl2 must have same length")
    
   # entropy 
   ent <- function(cl){
        n <- length(cl)
        p <- table(cl)/n
         -sum(p*log(p, base=base))
   }
   # mutual information
   mi <- function(cl1,cl2){
        p12 <- table(cl1,cl2)/length(cl1)
        p1p2 <- outer(table(cl1)/length(cl1),table(cl2)/length(cl2))
        sum(p12[p12>0]*log(p12[p12>0]/p1p2[p12>0], base=base))
   }
   
   if(!parts) return(ent(cl1) + ent(cl2) -2*mi(cl1,cl2))
   ent1 <- ent(cl1)
   ent2 <- ent(cl2)
   mi12 <- mi(cl1,cl2)
   c("vi"=ent1+ent2-2*mi12, "H(1|2)" =ent1-mi12, "H(2|1)"=ent2 -mi12)
}
