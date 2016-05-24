require(elliptic)
require(MASS)

safety.factor <- 100
"test" <- function(x,abs.error=1e-6){stopifnot(abs(x)<abs.error*safety.factor)}



# equations 16.20.1 to 16.20.3, Jacobi's imaginary transform:
u <-  seq(from=1,to=4+1i,len=10)
m <- 0.1+0.1123312i
test(sn(1i*u,m=m) - 1i*sc(u,m=1-m))
test(cn(1i*u,m=m) -    nc(u,m=1-m))
test(dn(1i*u,m=m) -    dc(u,m=1-m))


# eqs 16.28.1-16.28.5, p576:
test(abs(e16.28.1(z=1:600,m=0.234+0.1i)),abs.error=2e-15)
test(abs(e16.28.2(z=1:600,m=0.234+0.1i)),abs.error=2e-15)
test(abs(e16.28.3(z=1:600,m=0.234+0.1i)),abs.error=2e-15)
test(abs(e16.28.4(z=1:600,m=0.234+0.1i)),abs.error=2e-15)

test(abs(e16.28.5(m=seq(from=-0.234+0.1i, to=0.44-1i,len=100))),abs.error=2e-15)

#Now, try and get 16.28.6, p576: theta1dash=theta2*theta3*theta4:
m <- 0.5
derivative <- function(small){(theta1(small,m=m)-theta1(0,m=m))/small}
right.hand.side <-  theta2(0,m=m)*theta3(0,m=m)*theta4(0,m=m)
test(derivative(1e-7)-right.hand.side)
#looks fine.




#now compare e16.37.1 (product form for Neville's theta.s) with
#eq 16.36.6:
test( e16.37.1(0.5,0.5) - theta.s(0.5,0.5))



#now compare   e16.37.2 with eq 16.36.6, part 2:
test(e16.37.2(0.75,0.5) - theta.c(0.75,0.5))

# An identity of the form pq= pn/qn:
test(theta.c(0.5,0.5)/theta.d(0.5,0.5)- cd(0.5,0.5),abs.error=1e-9)



#Now check Laurent series for equianharmonic case, table on p656:
page.656 <- ck(g=c(0,1),n=12)-as.vector(t(cbind(0,0,c(1/28,1/10192,1/5422144,3/(5*13^2*19*28^4)))))
test(abs(page.656),abs.error=1e-19)



# Example 2, p579:
test(dn(0.20,0.19)-0.996253, abs.error=1e-6)

# Example 3, p579:
test(dn(0.2,0.81)-0.98406, abs.error=1e-5)

# Example 4, p580:
test(cn(0.2,0.81)-0.980278, abs.error=1e-6)

# Example 5, p580:
test(dc(0.672,0.36)-1.174,abs.error=1e-4)

# Example 6, p580:
test(Theta(0.6,m=0.36)-0.97357,abs.error=1e-5)

# Example 7, p581:
test(cs(0.5360162,0.09)-1.6918083,abs.error=1e-7)

# Example 8, p581:
test(sn(0.61802,0.5)-0.56458,abs.error=1e-5)

#Example 9, p581:
test(sc(0.61802,m=0.5)-0.68402,abs.error=1e-5)

#Example 11, p581:
test(cs(0.99391,m=0.5)-0.75,abs.error=1e-5)

# Example 8, p666, LHS
test(P(z=0.07 + 0.1i, g=c(10,2)) - (-22.97450010 - 63.0532328i),abs.error=1e-7)

# Now check sigma() against some Maple arguments:
test(sigma(1+0.4i,g=c(2+0.3i,1-0.99i)) - (1.006555817+0.3865197102i),abs.error=1e-9)
test(sigma(10-8i,g=c(1-0.4i,2.1-0.7i))-(-1.033893831e18 + 6.898810975e17i),1e11)
test(sigma(4,g=c(2,3)) - (-80.74922381),abs.error=1e-7)

#Now verify that g2.fun() and g3.fun() are in fact unimodular:
o <- c(1,1i)
 test(abs(unimodularity(7,o,FUN=g2.fun, maxiter=100)-g2.fun(o)),abs.error=1e-9)
 test(abs(unimodularity(7,o,FUN=g3.fun, maxiter=100)-g3.fun(o)),abs.error=2e-9)

M <- congruence(c(4,9))
test(abs(g.fun(o) - g.fun(M %*% o,maxiter=840)),2e-13)

# Verify Jacobi's formula numerically:
test(theta1dash(z=0,q=0.1+0.2i) -  theta1.dash.zero.q(0.1+0.2i),abs.error=3e-16)

#Now verify theta1.dashdash etc:

#d/dz (theta1) == theta1dash:
m <- 0.3+0.31i
z <- seq(from=1,to=2+1i,len=7)
delta <- 0.001
deriv.numer <- (theta1(z=z+delta,m=m) - theta1(z=z,m=m))/delta
deriv.exact <- theta1dash(z=z+delta/2,m=m)
test(deriv.numer-deriv.exact,abs.error=1e-7)

#d/dz (theta1dash) == theta1dashdash:
deriv.numer <- (theta1dash(z=z+delta,m=m) - theta1dash(z=z,m=m))/delta
deriv.exact <- theta1dashdash(z=z+delta/2,m=m)
test(deriv.numer-deriv.exact,abs.error=1e-7)

#d/dz (theta1dashdash) == theta1dashdashdash:
deriv.numer <- (theta1dashdash(z=z+delta,m=m) - theta1dashdash(z=z,m=m))/delta
deriv.exact <- theta1dashdashdash(z=z+delta/2,m=m)
test(deriv.numer-deriv.exact,abs.error=2e-7)


# Example 13, page 668, LHS:
test(sigma.laurent(z=0.4 + 1.3i,g=c(8,4),nmax=8)-(0.278080 + 1.272785i),abs.error=6e-8)

# Example 13, page 668, RHS:
test(sigma.laurent(z=0.8 + 0.4i,g=c(7,6),nmax=8)-(0.81465765 + 0.38819473i),abs.error=1e-8)


# Check  P() against Some Maple outputs (I just made up the arguments):
test(P(1+0.3i,g=c(1+1i,2-0.33i),give.all.3=TRUE)-(0.8231651984-0.3567903513i),abs.error=1e-10)
test(P(-4-4i,g=c(0.3123+10i,0.1-0.2222i),give.all.3=TRUE)-(-1.118985985-1.038221043i),abs.error=1e-9)
test(P(10+2i,g=c(1,4+0i),give.all.3=TRUE)-(2.021264367-0.9875939553i),abs.error=1e-10)


# check e18.10.9, p650:
test(e18.10.9(parameters(g=c(1,3+0.2i))), abs.error=2e-14)
test(e18.10.9(parameters(g=c(1,3+  0i))), abs.error=1e-14)
test(e18.10.9(parameters(g=c(1,0.1+0i))), abs.error=1e-14)


# check that P'2=4P^3-g2*P-g3:
g <- c(1.44+0.1i, -0.3+0.99i)
g2 <- g[1]
g3 <- g[2]
u <- parameters(g=g)
z <- seq(from= 10-14i, to=-10+20i, len=777)
p <- P(z,g)
pd <- Pdash(z,g)
test(4*p^3-g2*p-g3-pd^2, 2e-11)


# check that (P')^2 =4(P-e1)(P-e2)(P-e3):
test(pd^2-4*(p-u$e[1])*(p-u$e[2])*(p-u$e[3]))


#now some tests of eta() and eta.series():
 z <- seq(from=1+1i,to=10+0.6i,len=99)
test(eta(z)-eta.series(z),abs.error=2e-14)
test(eta(z+1)-eta(z)*exp(pi*1i/12),1e-14)
test(eta(1i)-gamma(1/4)/2/pi^(3/4),abs.error=1e-15)
test(theta3(0,q=exp(pi*1i*z))-eta((z+1)/2)^2/eta(1+z),abs.error=4e-15)


#now test J() and lambda() for being unimodular:
 M <- matrix(c(5,4,16,13),2,2)
 z <- seq(from=1+1i,to=3+3i,len=10)
 test(J(z)-J(M %mob% z,maxiter=100),1e-7)
 test(lambda(z)-lambda(M %mob% z,maxiter=100),1e-12)

# some identities for lambda function:
 z <- seq(from= -1+0.42i,to=10+7i,len=20)
 test(lambda(z)-lambda(z+2))
 test(lambda(z+1)-lambda(z)/(lambda(z)-1))


# and one for J():
test(J(1i+-10:10)-1,abs.error=2e-15)


#standard example of divisor function:
test(divisor(140)-336)

#divisor() is multiplicative:
test(divisor(11*12)-divisor(11)*divisor(12))

#Euler's generalization of Fermat's little theorem:
test((2^totient(15)) %%15 - 1)

#totient(p)= p-1 for prime p:
test(1+totient(primes(100)) - primes(100))

#totient() is multiplicative:
test(totient(25)*totient(1:14)-totient(25)*totient(1:14))

#mobius() is multiplicative:
test(mobius(23*1:10)-mobius(23)*mobius(1:10))

#Numerical verification of Mobius inversion theorem, using f(n)=1.

mobius.invert <- function(n){
  f <- factorize(n)
  d <- unique(apply(f^t(expand.grid(lapply(1:length(f),function(...){0:1}))),2,prod))
  sum(mobius(d)*divisor(n/d,k=0))
}

jj <- c(1:10,1000:1030)
test(sapply(jj,mobius.invert)-1)
