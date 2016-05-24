x <- seq(0,30,by=0.25)
density <- 0.3 * dnorm(x,8,2) + 0.7 * dnorm(x,16,3)
myplot <- xyplot(density~x,type="l", 
    main="pdf of a mixture of normals"
    )
