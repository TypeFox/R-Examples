x = expand.grid(1:3, 1:2, 3:1)
x = x[sample(seq1(1, nrow(x)), nrow(x)),]
x

ord = roworder(x)
ord

x[ord,]
