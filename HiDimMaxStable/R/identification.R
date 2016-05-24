
select.mean<-function(data,t) {
    m<-apply(data,2,mean)
    s<-data[,m>t,drop=FALSE]
    list(selected=s/t,
         u=t(t(s)/m[m>t]),
         gamma=m[m>t]/t)
}

