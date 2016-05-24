pl.dist <-
function(variable, p.start=0.1, epsylon = 10^(-8)){
p.old = c(); p.new = c(); t=c(); t.avg=c()
p.old = p.start
diff = 0.1
n.iter=0
while(diff>epsylon) {
# E-step 
t = ((variable+p.old+3)*(variable+1))/((variable+p.old+2)*(p.old+1))
t.avg=mean(t)
# M-step
p.new = (-(t.avg-1)+sqrt((t.avg)^2+6*t.avg+1))/(2*t.avg)
diff = abs(p.old-p.new)
p.old = p.new
n.iter = n.iter + 1
}
outlist = list(p=p.new, n.iter=n.iter)
    class(outlist) = "pl.dist"
    return(outlist)

}
