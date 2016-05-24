findP <-
function(y, step=0.01) {

        Call <- match.call()

        m=matrix(ncol=3,nrow=floor((1-step-step)/step)+1);
        n=0;
        for(i in seq(step,1-step,by=step)){

                n=n+1;

                z=y;
                z[,1]=y[,1]^(i);
                z[,2]=y[,2]^(i);

                d=distance(z[,1],z[,2]);
                D=summary(d);

                m[n,]=c(i,(D[5]-D[2]),abs(0.5-D[3]));
        }

	maxIQR=m[m[,2]==max(m[,2]),1];
        minMed=m[m[,3]==min(m[,3]),1];
        rlist <- list(values=m,maxIQR=maxIQR,minMed=minMed);

        rlist$call <- Call
        class(rlist) <- 'findP'
        rlist
}

