corr.nn4on<-function(p, ON.cor){
w=rnorm(100000,0,1)
y=rnorm(100000,0,1)
pvec=p
yord = numeric(length(y))
for (r in 1:length(pvec)){
if (r !=length(pvec)) {
t1 = qnorm(pvec[r])
t2 = qnorm(pvec[r+1] )
yord[(t1<y)&(y<=t2)]= r
} else {
yord[y>qnorm(pvec[r])]= r
}}
yord=yord+1
c=cor(yord[order(yord)],w[order(w)])/cor(y[order(y)],w[order(w)]) 
corrected=ON.cor/c 
return(corrected)
}