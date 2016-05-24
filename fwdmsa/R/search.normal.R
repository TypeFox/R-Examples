# Aangepast op 15 februari 2008
"search.normal" <-
function(X, lowerbound =.3, alpha = .05, verbose = T){

   # Internal functions
   
   any.neg <- function(x){if(any(x < 0))T else F}

   adjusted.alpha <- function(alpha, K) alpha/(K[1]*(K[1]-1)*.5 + sum(K[-1]))

   fitstring <- function(string.arg,length.arg) substr(paste(string.arg,"                        "),1,length.arg)
   
   newH <- function(j,in.this.set, x, lowerbound, Z.c){

     newX <- cbind(x[,in.this.set==1],x[,j])
     H.list <- coefH(newX)
     if (H.list$Hi[length(H.list$Hi)] < lowerbound) return(-98) # less than lower bound
     Zi <- coefZ(newX)$Zi
     if (Zi[length(Zi)] < Z.c) return(-97)                      # not significant
     return(H.list$H)
   }

   # initial calculations
   
   X <- check.data(X)
   item.label <- dimnames(X)[[2]]
   N <- nrow(X)
   S <- var(X)

   if(any(is.na(diag(S/S)))) stop("At least one item has no variance")

   Smax <- var(apply(X,2,sort))
   Hij <- S/Smax
   Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
   Zij <- (S * sqrt(N-1))/sqrt(Sij)

   J <- nrow(Hij)
   result <- rep(-99,J);
   j <- 0
   InSet <- rep(0,J)
   scale <- 0

   # start scaling
   repeat{
     scale <- scale + 1
     step <- 1
     K <- rep(0,J)

     if(verbose){ 
       cat("",fill=T)
       cat("SCALE",scale,fill=T)
       cat("",fill=T)
     }  

     # Are there two items left?
     if(length(InSet[InSet==0]) < 2){
       if(verbose) cat("Less than two items left. PROCEDURE STOPS",fill=T)
       break
     }

     # Compute the critical value for Zij 
     
     K[step] <- length(InSet[InSet == 0])
     Z.c <- abs(qnorm(adjusted.alpha(alpha,K)))

     # Select the first two items

     Hselect <- Hij
     Hselect[abs(Zij) < Z.c] <- -99
     Hselect[InSet > 0 & InSet < scale,] <- -99
     Hselect[,InSet > 0 & InSet < scale] <- -99
     Hselect[col(Hselect) >= row(Hselect)] <- -99

     # Check if there are any feasible values left
     if(max(Hselect) == -99){
       if(verbose) cat("Scale ", scale," could not be formed due to H < ",lowerbound,". PROCEDURE STOPS",fill=T)
       break
     }

     first.item <- row(Hselect)[Hselect==max(Hselect)]
     second.item <- col(Hselect)[Hselect==max(Hselect)]
     maxHij <- Hij[first.item,second.item]



     # Check if H of two item-scale is greater than c
     if(maxHij < lowerbound){
       if(verbose) cat("Scale ", scale," could not be formed due to H < ",lowerbound,". PROCEDURE STOPS",fill=T)
       break
     }


     # Add the first two items to the scale
     if(verbose){
       cat("Item: ",fitstring(item.label[first.item],20)," Scale", scale," H = ",round(maxHij,2),fill=T)
       cat("Item: ",fitstring(item.label[second.item],20)," Scale", scale," H = ",round(maxHij,2),fill=T)
     }  
     InSet[first.item] <- scale
     InSet[second.item] <- scale

# Wat te doen als er meerdere maximale waarden van Hij zijn tijdens het selectieproces.

     # Adding new items
     repeat{
       step <- step + 1

       # exclude items from previous scales
       in.this.set <- InSet
       in.this.set <- ifelse(InSet == scale, 1,0)
       in.this.set <- ifelse(InSet <  scale & InSet > 0,-1,in.this.set)

       # exclude items having a negative covariance with the already selected items
       neg1 <- apply(Hij[in.this.set==1,],2,any.neg)
       neg2 <- apply(Hij[,in.this.set==1],1,any.neg)
       in.this.set[neg1|neg2 & in.this.set==0] <- -1

       # Are there items left after the exclusion?
       available.items <- which(in.this.set==0)
       if(length(available.items)==0){
         if(verbose) cat("Scale ", scale," is completed. No items left with Hij => 0",fill=T)
         break
       }

       # Compute H and Hi of potentially new items
       result[in.this.set!=0] <- -99  # items already selected in other scales
       K[step] <- length(available.items)
       Z.c <- abs(qnorm(adjusted.alpha(alpha,K)))
       for (j in available.items) result[j] <- newH(j,in.this.set, X, lowerbound, Z.c)


       # Is maximum value Hi greater than c?
       if(max(result) < lowerbound){
         if(verbose) cat("Scale ", scale," is completed. No items left such that Hi > ",lowerbound,".",fill=T)
         break
       }

       # Add the newly selected item to the scale
       new.item <- row(as.matrix(result))[result==max(result)]
       InSet[new.item] <- scale
       if(verbose) cat("Item: ",fitstring(item.label[new.item],20)," Scale", scale," H = ",round(max(result),2),fill=T)
     }
  # start with next scale
  }
  InSet <- as.matrix(InSet)
  dimnames(InSet) <- list(item.label,"Scale")
  return(InSet)
}
