# 5) permutations: without replacement: distinct items

I = iterpc(5, 2, ordered=TRUE)
getall(I)


# 6) permutations: with replacement: distinct items

I = iterpc(5, 2, replace=TRUE, ordered=TRUE)
getall(I)
    
# 7) permutations: without replacement: non distinct items

x = c("a", "a", "b", "c")
I = iterpc(table(x), 2, ordered=TRUE)
# or I = iterpc(c(2,1,1), 2, label=c("a", "b", "c"), ordered=TRUE)
getall(I)


# 8) permutations: with replacement: non distinct items

x = c("a", "a", "b", "c")
I = iterpc(table(x), 2, replace=TRUE, ordered=TRUE)
# or I = iterpc(c(2,1,1), 2, label=c("a", "b", "c"), replace=TRUE, ordered=TRUE)
getall(I)
