csdl<-function(z, psi, L, ridge=NULL, ndx=round(L/3),
    DL=FALSE, diff.varying=FALSE){ #heat.power
        r<-z
        attr(r,"psi")<-psi
        attr(r,"L")<-L
        attr(r,"ndx")<-ndx
        attr(r,"DL")<-DL
        attr(r,"ridge")<-ridge
        attr(r,"diff.varying")<-diff.varying
        attr(r,"heat.power")<-1
        return(r)
        }
