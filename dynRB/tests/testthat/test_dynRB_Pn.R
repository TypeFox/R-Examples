
t1a <- runif(100,min=0,max=1)
t2a <- runif(100,min=0.7,max=1)
S1 <- data.frame(Species="A", t1a, t2a)
t2a<-c(t2a[1:100]-0.31)
S2 <- data.frame(Species="B", t1a, t2a)
A <- rbind(S1, S2)
  
r<-dynRB_Pn(A, correlogram = FALSE)
  
expect_equal(c(1,1,1,1),r$result[,3])
expect_equal(c(1,0,0,1),r$result[,4])


###############################
 
t1a <- runif(100,min=0,max=1)
t2a <- runif(100,min=0,max=1)
S1 <- data.frame(Species="A", t1a, t2a)
S2 <- data.frame(Species="B", t1a, t2a)
A <- rbind(S1, S2)

r<-dynRB_Pn(A, correlogram = FALSE)

expect_equal(c(1,1,1,1),r$result[,3])
expect_equal(c(1,1,1,1),r$result[,4])
  
  
###############################  
  

t1a <- c(runif(100,min=0,max=1),runif(0,min=0,max=0.1))
t2a <- c(runif(100,min=0.7,max=1),runif(0,min=0,max=0.1))
S1 <- data.frame(Species="A", t1a, t2a)
t2a<-c(t2a[1:100]-0.3)
S2 <- data.frame(Species="B", t1a, t2a)
A <- rbind(S1, S2)

#---

t3a <- c(runif(100,min=0,max=1),runif(0,min=0,max=0.1))
t4a <- c(runif(100,min=0,max=1),runif(0,min=0,max=0.1))
S1 <- data.frame(Species="A", t3a, t4a)
S2 <- data.frame(Species="B", t3a, t4a)
B <- rbind(S1, S2)
C<-cbind(A,B[,2:3])
A<-C

r<-dynRB_Pn(A, correlogram = FALSE)
x <- as.numeric(r$result[2,3:6])
expect_equal(c(1,0,1,1),x)
x <- as.numeric(r$result[3,3:6])
expect_equal(c(1,0,1,1),x)

  













  
  
  