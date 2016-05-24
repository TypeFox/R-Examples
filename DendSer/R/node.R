






r0 <- function(i,o,se,h=NULL) 
{
 #reflects the current node
  r <- se[1,i]:se[4,i]
  o[r] <- o[rev(r)]
  o
}

t0 <- function(i,o,se,h=NULL) 
{ 
  #translates the current node
   r <- se[1,i]:se[4,i]
   r1 <- c(se[3,i]:se[4,i],se[1,i]:se[2,i])
    o[r] <- o[r1]

  o
}

c0 <- function(i,o,se,h=NULL) 
{
 #reflects the current node
  r <- se[1,i]:se[4,i]
  o1 <- o
  o1[r] <- o[rev(r)]
  if (length(r)==2)
   return(o1)
   else{
   o2 <- o
   r1 <- c(se[3,i]:se[4,i],se[1,i]:se[2,i])
   o2[r] <- o2[r1]
  cbind(o1,o2,deparse.level=0)
}}


r1 <-function(i,o,se,h=NULL){

 #looks at all combinations of the reflection of sub-nodes of current node 
 # all returned columns are unique
 a <- se[1,i]
 b <- se[2,i]
 c <- se[3,i]
 d <- se[4,i]
 if (a!=b && c!=d){
 	o2 <- o
 	o2[a:b] <- o2[b:a]
 	o3 <- o
 	x <- c:d
 	o3[x] <- o3[d:c]
 	o4 <- o2
 	o4[x] <- o3[x]
 	return(cbind(o2,o3,o4,deparse.level=0))
 }
 else if (a!=b){
 	o2 <- o
 	o2[a:b] <- o2[b:a]
 	return(o2)
 }
 else if (c!=d){
 	o3 <- o
 	o3[c:d] <- o3[d:c]
 	return(o3)
 	
 }
 else return(NULL)
 
 }



r01 <-function(i,o,se,h=NULL){

 #looks at all combinations of the reflection of sub-nodes of current node # all returned columns are unique
 a <- se[1,i]
 b <- se[2,i]
 c <- se[3,i]
 d <- se[4,i]
 if (a!=b && c!=d){
 	o2 <- o
 	o2[a:b] <- o2[b:a]
 	o3 <- o
 	cd <- c:d
 	o3[cd] <- o3[d:c]
 	o4 <- o2
 	o4[cd] <- o3[cd]
 	or <- cbind(o,o2,o3,o4,deparse.level=0)
 	or[a:d,] <- or[d:a,]
 	return(cbind(or,o2,o3,o4,deparse.level=0))
 }
 else if (a!=b){
 	o2 <- o
 	o2[a:b] <- o2[b:a]
 	or<- cbind(o,o2,deparse.level=0)
 	or[a:d,] <- or[d:a,]
 	return(cbind(or,o2,deparse.level=0))
 }
 else if (c!=d){
 	o3 <- o
 	o3[c:d] <- o3[d:c]
 	or<- cbind(o,o3,deparse.level=0)
 	or[a:d,] <- or[d:a,]
 	return(cbind(or,o3,deparse.level=0))
 	
 }
 else {
 	or <- o
 	or[a:d] <- or[d:a]
 	return(or)
 	}
 }


t1 <-function(i,o,se,h){

 #looks at all combinations of the translation of sub-nodes of current node 
 # all returned columns are unique
 m <- h$merge
 a <- se[1,i]
 b <- se[2,i]
 c <- se[3,i]
 d <- se[4,i]
 if (a==b && c==d) return(NULL)
 else if (a==b || c==d){
  j <- max(m[i,])
  return(t0(j,o,se,NULL))
 }
 else {
   j <- m[i,1]
   k <- m[i,2]
   if (se[1,j] != se[1,i]) {
   	  x <- j
   	  j <- k
   	  k <- x
    }
    tl <- t0(j,o,se,NULL)
    tr <- t0(k,o,se,NULL)
    tboth <- tl
    tboth[se[3,i]:se[4,i]] <- tr[se[3,i]:se[4,i]]
   return(cbind(tl,tr,tboth,deparse.level=0))
 }
 }
 
 
t01 <-function(i,o,se,h){

 #looks at all combinations of the translation of sub-nodes of current node 
 # and translating self
 a <- t0(i,o,se,h)
# se <- preph(m,o)
 b <- t1(i,o,se,h)	
 if (is.null(b))
	  return(a)
 else {
 	s1 <- se[1,i]
	s2 <- se[3,i]
	n1 <- se[2,i]-s1+1
	n2 <- se[4,i]-s2+1
 	b <- cbind(b,deparse.level=0)
 	c <- b
	c[s1:(s1+n2-1),] <- b[s2:(s2+n2-1),]	
	c[(s1+n2):(s1+n2+n1-1),] <- b[s1:(s1+n1-1),]
	return(cbind(a,b,c,deparse.level=0))	 	
 }
 }


r0 <- cmpfun(r0)
t0 <- cmpfun(t0)
c0 <- cmpfun(c0)
r1 <- cmpfun(r1)
t1 <- cmpfun(t1)
t01 <- cmpfun(t01)
r01 <- cmpfun(r01)