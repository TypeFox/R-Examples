require(onion)

test <- function(x, TOL= 1e-10){
  stopifnot(Mod(x)<TOL)
  return(TRUE)
}

f <- function(...){
  ## First the quaternions:
  stopifnot(Hi*Hj ==  Hk)
  stopifnot(Hj*Hi == -Hk)
  stopifnot(Hj*Hk ==  Hi)
  stopifnot(Hk*Hj == -Hi)
  stopifnot(Hk*Hi ==  Hj)
  stopifnot(Hi*Hk == -Hj)

  stopifnot(Hi*Hi == -H1)
  stopifnot(Hj*Hj == -H1)
  stopifnot(Hk*Hk == -H1)

  stopifnot(H1*H1 == H1)
  stopifnot(H1*Hi == Hi)
  stopifnot(H1*Hj == Hj)
  stopifnot(H1*Hk == Hk)

  stopifnot(H1*H1 == H1)
  stopifnot(Hi*H1 == Hi)
  stopifnot(Hj*H1 == Hj)
  stopifnot(Hk*H1 == Hk)
  
  stopifnot(Hi*Hj*Hk == -H1)

  ## Quaternion zero times table:
  stopifnot(H0*H1 == H0)
  stopifnot(H0*Hi == H0)
  stopifnot(H0*Hj == H0)
  stopifnot(H0*Hk == H0)

  stopifnot(H1*H0 == H0)
  stopifnot(Hi*H0 == H0)
  stopifnot(Hj*H0 == H0)
  stopifnot(Hk*H0 == H0)

  ## And some quaternion additions:
  stopifnot(H1 + Him == Hall)
  stopifnot(Hi + Hj + Hk == Him)
  stopifnot(H1 + Hi + Hj + Hk == Hall)

  ## And some quaternion subtractions:
  stopifnot(Hi - Hi == H0)
  stopifnot(Hall - Hi - Hj - Hk == H1)
  stopifnot(Hall - Him == H1)

  ## Now all 64 of the octonions:
  stopifnot(O1*O1  == O1 )
  stopifnot(O1*Oi  == Oi )
  stopifnot(O1*Oj  == Oj )
  stopifnot(O1*Ok  == Ok )
  stopifnot(O1*Ol  == Ol )
  stopifnot(O1*Oil == Oil)
  stopifnot(O1*Ojl == Ojl)
  stopifnot(O1*Okl == Okl)
  
  stopifnot(Oi*O1  ==  Oi )
  stopifnot(Oi*Oi  == -O1 )
  stopifnot(Oi*Oj  ==  Ok )
  stopifnot(Oi*Ok  == -Oj )
  stopifnot(Oi*Ol  ==  Oil)
  stopifnot(Oi*Oil == -Ol )
  stopifnot(Oi*Ojl == -Okl)
  stopifnot(Oi*Okl ==  Ojl)

  stopifnot(Oj*O1  ==  Oj )
  stopifnot(Oj*Oi  == -Ok )
  stopifnot(Oj*Oj  == -O1 )
  stopifnot(Oj*Ok  ==  Oi )
  stopifnot(Oj*Ol  ==  Ojl)
  stopifnot(Oj*Oil ==  Okl)
  stopifnot(Oj*Ojl == -Ol )
  stopifnot(Oj*Okl == -Oil)

  stopifnot(Ok*O1  ==  Ok )
  stopifnot(Ok*Oi  ==  Oj )
  stopifnot(Ok*Oj  == -Oi )
  stopifnot(Ok*Ok  == -O1 )
  stopifnot(Ok*Ol  ==  Okl)
  stopifnot(Ok*Oil == -Ojl)
  stopifnot(Ok*Ojl ==  Oil)
  stopifnot(Ok*Okl == -Ol )

  stopifnot(Ol*O1  ==  Ol )
  stopifnot(Ol*Oi  == -Oil)
  stopifnot(Ol*Oj  == -Ojl)
  stopifnot(Ol*Ok  == -Okl)
  stopifnot(Ol*Ol  == -O1 )
  stopifnot(Ol*Oil ==  Oi )
  stopifnot(Ol*Ojl ==  Oj )
  stopifnot(Ol*Okl ==  Ok )

  stopifnot(Oil*O1  ==  Oil)
  stopifnot(Oil*Oi  ==  Ol )
  stopifnot(Oil*Oj  == -Okl)
  stopifnot(Oil*Ok  ==  Ojl)
  stopifnot(Oil*Ol  == -Oi )
  stopifnot(Oil*Oil == -O1 )
  stopifnot(Oil*Ojl == -Ok )
  stopifnot(Oil*Okl ==  Oj )

  stopifnot(Ojl*O1  ==  Ojl)
  stopifnot(Ojl*Oi  ==  Okl)
  stopifnot(Ojl*Oj  ==  Ol )
  stopifnot(Ojl*Ok  == -Oil)
  stopifnot(Ojl*Ol  == -Oj )
  stopifnot(Ojl*Oil ==  Ok )
  stopifnot(Ojl*Ojl == -O1 )
  stopifnot(Ojl*Okl == -Oi )

  stopifnot(Okl*O1  ==  Okl)
  stopifnot(Okl*Oi  == -Ojl)
  stopifnot(Okl*Oj  ==  Oil)
  stopifnot(Okl*Ok  ==  Ol )
  stopifnot(Okl*Ol  == -Ok )
  stopifnot(Okl*Oil == -Oj )
  stopifnot(Okl*Ojl ==  Oi )
  stopifnot(Okl*Okl == -O1 )



  ## And the zero octonion times table:
  stopifnot(O0*O0  == O0)

  stopifnot(O0*O1  == O0)
  stopifnot(O0*Oi  == O0)
  stopifnot(O0*Oj  == O0)
  stopifnot(O0*Ok  == O0)
  stopifnot(O0*Ol  == O0)
  stopifnot(O0*Oil == O0)
  stopifnot(O0*Ojl == O0)
  stopifnot(O0*Okl == O0)

  stopifnot(O1*O0  == O0)
  stopifnot(Oi*O0  == O0)
  stopifnot(Oj*O0  == O0)
  stopifnot(Ok*O0  == O0)
  stopifnot(Ol*O0  == O0)
  stopifnot(Oil*O0 == O0)
  stopifnot(Ojl*O0 == O0)
  stopifnot(Okl*O0 == O0)

  ## And some octonion additions:
  stopifnot(O1 + Oim == Oall)
  stopifnot(Oi + Oj + Ok + Ol + Oil + Ojl + Okl == Oim)
  stopifnot(H1 + Oi + Oj + Ok + Ol + Oil + Ojl + Okl == Oall)

  ## And some subtractions:
  stopifnot(Oil - Oil == O0)
  stopifnot(Oall - Oim == O1)
  
  ## Dummy return value:
  return(TRUE)
}

g <- function(...){
  ## Just pick some random quaternions
  x <- as.quaternion(c(pi,sqrt(2),-3,10.1),single=TRUE)
  y <- as.quaternion(c(exp(1),-2.22222,1/4,-1),single=TRUE)
  z <- as.quaternion(c(exp(-0.1), 0.1122, -2, -0.001),single=TRUE)

  ## Verify associativity:
  test(associator(x,y,z))

  ## And distributivity:
  test(x*(y+z) - (x*y+x*z)) 

  ## And *power* associativity of the octonions:
  jj1 <- x + Oil*y + Oj*z
  test( jj1*(jj1*jj1) - (jj1*jj1)*jj1) 
 
  ## And distributivity of octonions:
  jj2 <- as.octonion(pi+1:8,single=TRUE)
  jj3 <- as.octonion(1.123^(1:8) ,single=TRUE)
  test(jj1*(jj2+jj3) - (jj1*jj2+jj1*jj3))

  ## And alternativity of octonions:
  test(jj1*(jj1*jj2) - (jj1*jj1)*jj2 )
  test(jj1*(jj2*jj1) - (jj1*jj2)*jj1 )

  ## Dummy return value
  return(TRUE)
}

options(use.R=TRUE)
f()
g()
options(use.R=FALSE)
f()
g()


x <- as.octonion(c(1,4,sqrt(2),pi,pi/3, 1e-2, -4,1-pi),single=TRUE)
y <- as.octonion(1:8,single=TRUE)
z <- as.octonion(sqrt(17:10),single=TRUE)

options(use.R = TRUE)
jj.T <- associator(x,y,z)

options(use.R = FALSE)
jj.F <- associator(x,y,z)

test(jj.T-jj.F)



# Now some randomish checks that verify vectorized addition:


h <- function(a){
  test(a-a)
  test(a + (-1)*a)
  test((-1)*a + a)

  test( (a+a  )-2*a)
  test( (a+a  )-a*2)
  test( (a+a+a)-3*a)
  test( (a+a+a)-a*3)

  test(a+1-a-1)

  test(a+a[1]-a-a[1])
  
  test(a/a - 1)
  test(a^2/a - a)
  test(a^3/a^2 - a)
  test(a^4/a^2 - a^2)

  test( (a+a)/a - 2)
  test( (a+a)/(a*a) - 2/a)

  test(a*a       - a^2)
  test(a*a*a     - a^3)    #recall that octonions are *power* associative
  test(a*a*a*a   - a^4)

  test(1/a - a^(-1))
  test(1/a^2 - a^(-2))

  test(1/(a^2) - (1/a)^2)
  test(1/(a^3) - (1/a)^3)

  test( (a/a[1])*a[1] - a)
  
  if(is.quaternion(a)){
    test(associator(a,a+1,a+Hi))
  } else if (is.octonion(a)){
    test(associator(a*(3*Ok+4*Oil),a*(Ok+Oil),a*(Ok+2*Oil))) 
  } else {
    stop("a must be quaternion or octonion")
  } 
}

h(as.octonion(matrix(1:5,nrow=8,ncol=10)))
h(as.quaternion(matrix(1:5,nrow=4,ncol=20)))
