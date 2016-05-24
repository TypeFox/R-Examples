noe.compute.cgh<-function(tn,ta,tb){
# Compute vectors c, g, h for Noe's recursion
problemct<-0
# Check that ta and tb are increasing and that a<=b
problemct<-problemct+is.numeric(ta[1]>tb[1])
i<-1
while (i < (tn+1)){
       i<-i+1
       problemct<-problemct+is.numeric(ta[i-1]>ta[i])
       problemct<-problemct+is.numeric(tb[i-1]>tb[i])
       problemct<-problemct+is.numeric(ta[i]>tb[i])}
if (problemct>0){
       print("noe.compute.cgh improper ordering")
       print(problemct)
       return}
# Merge ta & tb into tc and compute tg & thia+ib-1

tc<-tg<-th<-rep(0,(2*tn+2))
tg[1]<-0
th[1]<-1
th[2]<-1
tc[1]<-0
a<-c(ta,1+tb[tn])
b<-c(tb,1+ta[tn])
ia<-ib<-1
while ((ia+ib-1) <= 2*tn){
       m<-ia+ib-1
       if (a[ia]<b[ib]) {
              tc[m+1]<-a[ia]
              tg[m+1]<-tg[m]+1
              th[m+2]<-th[m+1]
              ia<-ia+1} 
       else { tc[m+1]<-b[ib]
              tg[m+1]<-tg[m]
              th[m+2]<-th[m+1]+1
              ib<-ib+1}
tc[2*tn+2]<-max(max(a),max(b))}
ans<-matrix(0,(2*tn+2),3)
ans[,1]<-tc
ans[,2]<-tg
ans[,3]<-th
ans}
