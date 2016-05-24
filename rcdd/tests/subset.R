
 library(rcdd)

 .Call("test_my_subset", integer(0), as.integer(1:3), as.integer(3))

 .Call("test_my_subset", as.integer(1:3), integer(0), as.integer(4))

 .Call("test_my_subset", as.integer(1:3), as.integer(1), as.integer(5))

 .Call("test_my_subset", as.integer(1:3), as.integer(1:3), as.integer(6))

 .Call("test_my_subset", as.integer(1:2), as.integer(1:3), as.integer(7))

 #### need to test for collisions

 pow2 <- 5
 hashsize <- 2^pow2

 i <- seq(1:1000)
 hashi <- (2654435761 * i) %% hashsize
 hash1 <- i[hashi == hashi[1]]

 set2 <- hash1[1:10]
 print(set2)

 set1 <- sample(set2, 4)

 .Call("test_my_subset", as.integer(set1), as.integer(set2), as.integer(pow2))
 
 set1 <- sample(set2, 4)

 .Call("test_my_subset", as.integer(set1), as.integer(set2), as.integer(pow2))

 set1 <- union(set1, 1001)

 .Call("test_my_subset", as.integer(set1), as.integer(set2), as.integer(pow2))

 options(error=dump.frames)

 .Call("test_my_subset", as.integer(set1), as.integer(hash1), as.integer(pow2))

 .Call("test_my_subset", as.integer(set1), as.integer(hash1),
      as.integer(pow2 + 3))

 options(error=NULL)

 sets <- list(1, 2, c(1, 2), 3, 4, 5, c(4, 5), c(4, 5, 6), c(1, 3, 5))
 sets <- lapply(sets, as.integer)

 .Call("nonred", sets, as.integer(pow2))

 maximal(sets)

 ##### intersections #####

 .Call("all_intersect", sets, as.integer(pow2))

 sets <- list(c(1, 2, 3), c(2, 3, 4), c(3, 4, 5), c(4, 5, 1), c(5, 1, 2))
 sets <- lapply(sets, as.integer)

 .Call("all_intersect", sets, as.integer(pow2))
 
 allIntersect(sets)

 ##### unions #####

 .Call("all_union", sets, as.integer(pow2))

 sets <- list(1, 2, c(3, 4), c(1, 3), 5)
 sets <- lapply(sets, as.integer)

 .Call("all_union", sets, as.integer(pow2))

 allUnion(sets)

