
t1a <- rep(seq(0,1,length=100),100)
t2a <- sort(rep(seq(0,1,length=100),100))

S1 <- data.frame(Species="A", t1a, t2a)

t1a <- rep(seq(0,1,length=10),10)
t2a <- sort(rep(seq(0,1,length=10),10))
S2 <- data.frame(Species="B", t1a, t2a)

A <- rbind(S1, S2)

r<-dynRB_Vn(A, correlogram = FALSE)
x <- as.numeric(r$result[1,2:3] > 0.98)
expect_equal(c(1,1),x)



#############################


t1a <- rep(seq(0,1,length=100),100)
t2a <- sort(rep(seq(0.4,1,length=100),100))

t3a <- t1a[t2a<0.7]
t4a <- t2a[t2a<0.7]
t1a <- t1a[t2a>=0.7]
t2a <- t2a[t2a>=0.7]

S1 <- data.frame(Species="A", t1a=c(t1a,t3a), t2a=c(t2a,t4a))
S2 <- data.frame(Species="B", t1a=t3a, t2a=t4a)

A <- rbind(S1, S2)

r<-dynRB_Vn(A, correlogram = FALSE)
x <- as.numeric(r$result[2,2:3] > c(0.97,0.48) & r$result[2,2:3] <= c(1,0.52))

expect_equal(c(1,1),x)


##############################


t1a <- rep(seq(0,1,length=100),100)
t2a <- sort(rep(seq(0,1,length=100),100))

t3a <- t1a[t2a<0.7]
t4a <- t2a[t2a<0.7]-0.02
t1a <- t1a[t2a>=0.7]
t2a <- t2a[t2a>=0.7]
S1 <- data.frame(Species="A", t1a=c(t1a,t3a), t2a=c(t2a,t4a))
S2 <- data.frame(Species="B", t1a=t3a, t2a=t4a)
A <- rbind(S1, S2)

r<-dynRB_Vn(A, correlogram = FALSE)
x <- as.numeric(c(r$result[1,2:3],r$result[2,2:3]) > c(0.97,0.97,0.97,0.67) & c(r$result[1,2:3],r$result[2,2:3]) <= c(1,1,1,0.73))

expect_equal(c(1,1,1,1),x)







