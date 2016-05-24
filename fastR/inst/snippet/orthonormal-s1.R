x <- c(1,1,1)
y <- c(1,1,-2)
w <- y - project(y,x)
dot(x,w)                                    # confirm normality
# these two column vectors are orthogonal and have correct span
cbind( x / vlength(x), w / vlength(w) )
