ball <- function(n = 4000, r = 1, dim = 2){ 
    M <- matrix(0,n,dim)
    count <- 0;
    rr <- r ^ 2
    while (count < n) 
    {
        p <- runif (dim, -r, r)        
        if (sum(p ^ 2) <= rr) M[count <- count + 1,] <- p      
    }
    M
}

nballvolume <- function (r, n) 
{
    k <- as.integer(n / 2)
    if (n %% 2 == 0)
    {
        pi ^ k / factorial(k) * r ^ n
    }
    else
    {
        2 * factorial(k) * (4 * pi) ^ k / factorial(n) * r ^ n
    }
}

mouseset <- function (n=4000, r.head=2, r.left.ear=1.1, r.right.ear=1.1, left.ear.dist=2.5, right.ear.dist=2.5, dim=2)
{
    vh <- nballvolume (r.head, dim)
    vl <- nballvolume (r.left.ear, dim)
    vr <- nballvolume (r.right.ear, dim)
    
    if (dim < 2) stop ("Illegal argument: dim < 2")
    pos.h <- rep(0, dim)
    
    pos.l <- pos.h
    pos.r <- pos.h
    
    l.offset <- left.ear.dist / sqrt(2)
    r.offset <- right.ear.dist / sqrt(2)
    
    pos.l[1] <- pos.l[1] - l.offset
    pos.l[2] <- pos.l[2] + l.offset
    
    pos.r[1] <- pos.r[1] + r.offset
    pos.r[2] <- pos.r[2] + r.offset
    
    hh = r.head ^ 2
    ll = r.left.ear ^ 2
    rr = r.right.ear ^ 2
    
    centers <- rbind(pos.h, pos.l, pos.r)
    rs      <- c(r.head, r.left.ear, r.right.ear)
    rrs     <- c(hh, ll, rr)
    
    M = matrix(0, n, dim)  
    
    count <- 0
    
    while(count < n)
    {
        gen <- min (1000, n - count)
        hits <- 0
        s <- sample(x = c(1, 2, 3), size=gen, prob=c(vh, vl, vr), replace=T)
        for (i in s)
        {
            r <- rs[i]
            random.p <- runif(dim, -r, +r)
            p <- centers[i,] + random.p
            if (sum(random.p ^ 2) < rrs[i])       
                if (i == 1) M[count <- count + 1,] <- p                  
            else if (i == 2) 
            {
                if (sum((p - pos.h) ^ 2) > hh && sum((p - pos.r) ^ 2) > rr) M[count <- count + 1,] <- p
            }
            else if (i == 3) 
            {
                if (sum((p - pos.h) ^ 2) > hh && sum((p - pos.l) ^ 2) > ll) M[count <- count + 1,] <- p                      
            }
        }       
    }
    M  
}
