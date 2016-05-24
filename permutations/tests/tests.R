require("permutations")

id <- cycle(list(list()))
stopifnot(is.id(id))


a1 <- c(3L, 14L, 9L, 18L, 10L, 19L, 4L, 17L, 12L, 7L, 5L, 15L, 1L, 
        16L, 8L, 13L, 11L, 2L, 6L, 20L)

problem <- as.word(a1)
problem*problem


j1 <- c(12L, 11L, 1L, 9L, 7L, 2L, 3L, 8L, 17L, 15L, 18L, 5L, 6L, 13L, 16L, 4L, 14L, 10L, 19L, 20L)
j2 <- c(11L, 12L, 14L, 3L, 10L, 8L, 9L, 13L, 20L, 16L, 6L, 2L, 17L, 15L, 7L, 1L, 18L, 5L, 4L, 19L)
v1 <- as.word(j1)
v2 <- as.word(j2)
word_prod_single(v1,v2)

# to create an identity cyclist, use list(c()).

b1 <- cycle(list(list(1:2,3:5,6:10,11:15)))
b2 <- b1 + as.cycle(100:119)

w <- list(c(1),c(2,5,3),c(4))
x <- list(c(1),c(2,5,3),c(4,7))
y <- list(c(4,6),c(2,5,1),c(8,3))
z <- list(c(4,6,8,2,3,1,7,9,5))


# following returns cyclists:
nicify_cyclist(x)
nicify_cyclist(y)
nicify_cyclist(z)

thing1 <- cycle(list(nicify_cyclist(x),nicify_cyclist(y),nicify_cyclist(z),list(c())))
thing2 <- cycle(list(x,y,z,list(c())))


ww <- word(matrix(replicate(12,sample(9)),ncol=9,byrow=TRUE))
xx <- word(matrix(replicate(12,sample(9)),ncol=9,byrow=TRUE))
jj <- ww
size(jj) <- 10

stopifnot(all(sgn(ww)*sgn(xx) == sgn(ww*xx)))
identical(as.word(as.cycle(xx)),xx)



stopifnot(inverse(thing1)*thing1 == id)
stopifnot(thing1*inverse(thing1) == id)


a <- as.word(thing1)
a[1] <- a[2]


as.cycle(1:9)


# combine cyclists using cycle(list(x,y,z))

z1 <- list(c(3, 17),39,c(10, 14, 5, 8, 11),c(18, 12), c(16, 7, 19, 13, 15, 4, 20, 9, 6, 1),32,22)
z2 <- list(1:4,10:11,20:33)

fish <- cycle(list(x,y,z)) 


x2 <- c(7,9,3,5,8,6,1,4,2)

x3 <- rbind(
    a = c(5L, 8L, 4L, 6L, 3L, 9L, 1L, 7L, 2L),
    b = c(2L, 5L, 3L, 9L, 1L, 7L, 6L, 4L, 8L),
    c = c(9L, 5L, 7L, 1L, 4L, 3L, 8L, 6L, 2L),
    d = 1:9,
    e = c(1L, 2L, 3L, 5L, 4L, 6L, 7L, 8L, 9L)
)

x4 <- cbind(x3,matrix(11:20,5,2))

id <- cycle(list(list()))

xjj <- cycle(list(
    list(c(1,2,4),c(3,6)),
    list(c(1,2),c(3,4,5,6,7)),
    list(c(1,2,4),c(3,6)),
    list(c(1,2,4),c(3,16)),
    list(c(1,2,4),c(3,6)),
    list(c(1,2,4),c(3,6)),
    list(),   # empty!
    list(c(1,2,4),c(3,56)),
    list(c(1,2,4),c(3,56),20:21)
    ))


cl3 <- list(c(3,25,1),5:10,c(17,19))


x3 <- word(x3)
#x3^6


shape_cyclist(list(list()))


m <- megaminx
p <- megaminx
m[1] <- id
m[2] <- as.cycle(1:10)
m[3] <- id
m[4] <- as.cycle(1:40)

m[6]*m[7]


stopifnot(all(p*inverse(p) == id))
stopifnot(all(inverse(p)*p == id))

stopifnot(all(p   ==   p^1))
stopifnot(all(p^2 == p*p^1))
stopifnot(all(p^3 == p*p^2))
stopifnot(all(p^4 == p*p^3))
stopifnot(all(p^5 == p*p^4))


rep(m,0)
sgn(thing1)



o1 <- rperm(31,20)
o2 <- rperm(31,20)
o3 <- rperm(31,20)

stopifnot(o1^(o2*o3) == (o1^o2)^o3)

trim(as.word(id))
trim(word(col(matrix(0,4,5))))

nothing <- word(col(matrix(0,4,5)))
trim(nothing)


o1[1] <- id   #former bug



## check Hall-Witt:

x <- rperm(100,7)
y <- rperm(100,8)
z <- rperm(100,9)

uu <-    # should be the identity
commutator(commutator(x,y),z^x) *
commutator(commutator(z,x),y^z) *
commutator(commutator(y,z),x^y) 
