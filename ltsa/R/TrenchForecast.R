`TrenchForecast` <-
function(z,r,zm,n,maxLead,UpdateAlgorithmQ=TRUE){
nz<-length(z)
if (n<0 || n>nz) 
    stop("error: invalid forecast origin")
if (length(r)< (nz+maxLead-1))
        stop("error: length(r) must be >= ", nz+maxLead-1)
stopifnot(maxLead > 0)
m<-nz-n
zc<-z-zm
zk<-zc[1:n]
zf<-vf<-matrix(numeric(maxLead*(m+1)),ncol=maxLead)
gk<-t(matrix(r[n + 1 + outer(1:maxLead,1:n,"-")],ncol=maxLead,byrow=TRUE))
GI<-TrenchInverse(toeplitz(r[1:n]))
gkGI<-crossprod(t(gk),GI)
zf[1,]<-zm+gkGI%*%zk
for (j in 1:maxLead)
    vf[1,j]<-r[1] - sum(gkGI[j,]*gk[j,])
if (m > 0){
    for (tt in 1:m){
            gk<-t(matrix(r[n+tt+1+outer(1:maxLead,1:(n+tt),"-")],ncol=maxLead,byrow=TRUE))
            zk<-c(zk,zc[n+tt])
            if (UpdateAlgorithmQ)
                GI<-ToeplitzInverseUpdate(GI, r[1:(n-1+tt)], r[n+tt])
            else
                GI<-TrenchInverse(toeplitz(r[1:(n+tt)]))
            gkGI<-crossprod(t(gk),GI)
            zf[tt+1,]<-zm + gkGI%*%zk
            for (j in 1:maxLead)
                vf[tt+1,j]<-r[1] - sum(gkGI[j,]*gk[j,])
        }
    }
dimnames(zf)<-dimnames(vf)<-list(n:(n+m), 1:maxLead)
list(Forecasts=zf,SDForecasts=sqrt(vf))
}

