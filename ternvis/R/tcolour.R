tcolour <-
function( p=cbind(1,1,1)/3, 
                     q=cbind(1,1,1)/3,
                     m=0.7,
                     flip=FALSE,
                     dich="none",
                     theta0=0             # theta0 = 60 (degrees) gives CMY colours
                     ) {   
   n <- nrow(p)
   tcolour <- rbind(NA * p[,1])

#
# first define some utility functions
#



dalton <- function(th){
dalton <- NA * th
for (i in 1:length(th))
{
if (!is.na(th[i])){
if (th[i]/(2.*pi) >  1/3 && th[i]/(2.*pi) <= 1/2) {dalton[i] <- pi/3 + 2*(th[i]-2*pi/3)}
if (th[i]/(2.*pi) <= 1/3)                         {dalton[i] <- th[i]/2}
if (th[i]/(2.*pi) >  1/2)                         {dalton[i] <- th[i]}
}
}
dalton
}                                # angle transformation for daltonism



theta_of_pq <- function(p=cbind(1,1,1)/3,q=cbind(1,1,1)/3)
{   
   x  <- xf(p)
   y  <- yf(p)     
      thout <- atan2(y-as.numeric(yf(q)),x-as.numeric(xf(q)))         # find wrt climatology 
      thout <- (7*pi/6 - thout) %% (2*pi)                             # make bottom left have angle 0 and go clockwise
    thout 
}   

KL <- function(p,q) {              # function used to calculate relative entropy
      if (p > 0)   {KL <- p * log(p/q) }
      if (p == 0)  {KL <- 0}
      if (q == 0)  {KL <- 0}
      KL
   }

H3 <- function(p=c(1,1,1)/3,q=c(1,1,1)/3) { # the subjective certainty of p (relative to q)
   Hout <- NA
   if (! is.na(p[1])) 
    if (! is.na(p[2])) 
     if (! is.na(p[3])) 
      Hout <-KL(p[1],q[1])+
             KL(p[2],q[2])+
	     KL(p[3],q[3]) 
   Hout <- Hout / log(1/min(q)) # normalise so that maximum is 1 in corner
   Hout
   }


   ptemp <- p
   qtemp <- q
   if (flip)                 # optionally transpose categories B and A
      {
       p[,1] <- ptemp[,3]
       p[,3] <- ptemp[,1]
       q[,1] <- qtemp[,3]
       q[,3] <- qtemp[,1]
      }
   
   for (i in 1:nrow(p))
   {
      p[i,] <- tscale(p[i,])
   }
   
	 
         th <- theta_of_pq(p=p,q=q)   # angle
         H  <- H3(p=p,q=q)            # entropy (relative to climatology)


th <- (th - (pi * theta0 / 180)) %% (2*pi)  # colour rotation, usually theta0=0
                                          # theta0=60 gives CMY
                                          # theta0 = 90 good for protanopia
                                      

th <- dalton(th) # apply nonlinear transformation to angles to allow for deutanopia (Daltonism) 
 
for (j in 1:n)
{
if(!is.na(th[j]) && !is.na(H[j])){

      tcolour[j] = hsv( h = th[j]/(2.*pi),
                        s = H[j]^m,
		        v = 1 )                 
       }
}
if (dich != "none")      tcolour <- dichromat(tcolour,type=dich) # simulate effect of colour blindness
for (i in 1:n)
{
if (is.na(th[i])) tcolour[i] <- NA
}
     tcolour
}
