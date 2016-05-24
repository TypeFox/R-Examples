a = expand.grid(2:3, 3:6)
a = a[sample(seq1(1, nrow(a)), nrow(a)),]
a

b = expand.grid(3:4, 2:5)
b = b[sample(seq1(1, nrow(b)), nrow(b)),]
b

i = rowmatch(a, b)
i
b[na.omit(i),] # matching rows
a[is.na(i),] # non matching rows
