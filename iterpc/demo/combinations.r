# 1) combinations: without replacement: distinct items

I = iterpc(5, 2)
getall(I)

# 2) combinations: with replacement: distinct items

I = iterpc(5, 2, replace=TRUE)
getall(I)
 
    
# 3) combinations: without replacement: non distinct items

x = c("a", "a", "b", "c")
I = iterpc(table(x), 2)
# or I = iterpc(c(2,1,1), 2, label=c("a", "b", "c"))
getall(I)

    
# 4) combinations: with replacement: non distinct items

x = c("a", "a", "b", "c")
I = iterpc(table(x), 2, replace=TRUE)
# or I = iterpc(c(2,1,1), 2, label=c("a", "b", "c"), replace=TRUE)
getall(I)
        
    
