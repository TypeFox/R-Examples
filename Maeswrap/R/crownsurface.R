
#----------------------------------------------------------------------------------------------------#

# Find SA with CW and CL

# Paraboloid
S.paraboloid <- function(CW,CL){
  
  a <- CW/2
  h <- CL
  
  pi*a/(6*h^2) * ((a^2 + 4*h^2)^1.5 - a^3) + pi*a^2
} 

# Ellipsoid surface area (actually prolate/oblate spheroid)
# From Maximilian Posch.
S.ellipsoid <- function(CW,CL){
  
  r1 <- CW/2  # rotational axis radius
  r2 <- CL/2
  
  omega <- function(x){        
    if(identical(x,1))x<-1.0001
    if(x < 1)res <- 1+asin(sqrt(1-x^2))/(x*sqrt(1-x^2))
    if(x > 1)res <- 1+log(x + sqrt(x^2-1))/(x*sqrt(x^2-1))
    res
  }
  omega_V <- Vectorize(omega)
  
  S <- 2*pi*r1^2*omega_V(r1/r2) 
  
  return(S)
}


# Cone, includes bottom circle.
S.cone <- function(CW,CL)pi*(CW/2)*sqrt((CW/2)^2 + CL^2) + pi*(CW/2)^2    

# Box
S.box <- function(CW,CL)2*CW^2 + 4*CW*CL

# Cylinder
S.cyl <- function(CW,CL)2*pi*(CW/2)*CL + 2*pi*(CW/2)^2



SAfun <- function(CW,CL, crownform){
  
  firstlet <- strsplit(crownform,"")[[1]][1]
  
  res <- switch(firstlet,
                c = S.cone(CW,CL),
                e = S.ellipsoid(CW,CL),
                p = S.paraboloid(CW,CL))
  
  return(res)
}




#-----------------------------------------------------------------------------------------------#

# find CW if CL and SA given


# Radius of paraboloid given surface area and ratio of height to radius.
# Includes lower circle.
# b is h/r1 (height divided by radius)
CW.paraboloid <- function(SA, CWCL){
  
  b <- 2/CWCL
  
  2*sqrt(SA/ ( pi/(6*b^2) * ((1+4*b^2)^1.5 - 1) + pi )  )
  
}


CW.ellipsoid <- function(SA, CWCL){
  
  # From Posch:
  omega <- function(x){
    
    res <- rep(NA, length(x))
    
    res[x<1] <- 1+asin(sqrt(1-x[x<1]^2))/(x[x<1]*sqrt(1-x[x<1]^2))
    res[x>1] <- 1+log(x[x>1]+sqrt(x[x>1]^2-1))/(x[x>1]*sqrt(x[x>1]^2-1))
    res
  }
  
  CW <- 2*sqrt(SA/(2*pi*omega(CWCL)))
  
  return(CW)
}



# Cone
CW.cone <- function(SA, CWCL){
  
  b <- 2/CWCL
  
  2*sqrt(SA / (pi*(1+sqrt(1+b^2))))
  
}

# Cylinder
CW.cyl <- function(SA, CWCL){
  
  gam <- 2/CWCL
  
  r <- sqrt(SA/(2*pi*(1+gam)))
  2*r
}

# Box
CW.box <- function(SA, CWCL){
  
  b <- 2/CWCL
  
  2*sqrt(SA/(8+8*b))
  
}


CWfun <- function(shape=c("BOX","CONE","PARA","ELIP","CYL"), 
                  SA, CWCL){
  
  shape <- match.arg(shape)
  cw <- switch(shape,
               BOX = CW.box(SA, CWCL),
               CONE = CW.cone(SA, CWCL),
               PARA = CW.paraboloid(SA, CWCL),
               ELIP = CW.ellipsoid(SA, CWCL),
               CYL = CW.cyl(SA, CWCL)
  )
  
  return(cw)  
}


# 
# ## Test: 
# CW.cone(S.cone(CW=2,CL=10), CWCL=2/10)
# CW.ellipsoid(S.ellipsoid(CW=2,CL=10), CWCL=2/10)
# CW.box(S.box(CW=2,CL=10), CWCL=2/10)
# CW.paraboloid(S.paraboloid(CW=2,CL=10), CWCL=2/10)
# CW.cyl(S.cyl(CW=2,CL=10), CWCL=2/10)
