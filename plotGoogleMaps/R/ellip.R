ellip <-
function(a,b, 
         alpha = 0, # degrees
         loc = c(0,0), n = 100,  Id="1"){
    B <- min(a,b)
    A <- max(a,b)
    ## B <= A
    d2 <- (A-B)*(A+B)                   #= A^2 - B^2
    phi <- 2*pi*seq(0,1, len = n)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt(B^2 + d2 * sp^2)
    xy <- r * cbind(cp, sp)
    ## xy are the ellipse points for alpha = 0 and loc = (0,0)
    al <- alpha * pi/180
    ca <- cos(al)
    sa <- sin(al)
   pts<- xy %*% rbind(c(ca, sa), c(-sa, ca)) + cbind(rep(loc[1],n),
                                                rep(loc[2],n))
   pts[n,]<-pts[1,]

   pol<-Polygon(pts,hole=FALSE)
   pol<-Polygons(list(pol),ID=Id)
   # Returnd Polygons obj. it is a list of Polygon objs.
   return(pol)
}
