pg.dist <-
function(variable, alpha.start=1, beta.start=0.3, epsylon = 10^(-8)){
alpha.old = c(); beta.old = c()
alpha.new = c(); beta.new = c()
alpha.old = alpha.start
beta.old = beta.start
diff = 0.1
n.iter=0
while(diff>epsylon) {
# E-step 
t = (variable+alpha.old)/(1+beta.old)
s = digamma(variable+alpha.old)-log(1+beta.old)
# M-step
beta.new = alpha.old/mean(t)
alpha.new = alpha.old-(digamma(alpha.old)+log(beta.new)-mean(s))/(trigamma(alpha.old))
diff = max(abs(alpha.old-alpha.new), abs(beta.old-beta.new ))
beta.old = beta.new 
alpha.old = alpha.new
n.iter = n.iter + 1
}
outlist = list(alpha=alpha.new, beta=beta.new, theta=1/beta.new, n.iter=n.iter)
    class(outlist) = "pg.dist"
    return(outlist)
}
