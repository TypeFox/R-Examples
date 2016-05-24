vspprofile<-function(M=50, N=50, maxdepth = 1000, deltobs=20, noise = 2e-04, M1 = c(9000, -6, .001))
    {
######  vertical seismic profile for inverse problem testing
        
        ###  this is the quadratic mode: depth versus velocity
        
        depth = seq(from=0, to=maxdepth, by=1)

        ### invert to get the velocity
        vee = ( -M1[2] + sqrt(  M1[2]^2 - 4*M1[3]*(M1[1]-depth) ) )/ (2*M1[3] )  
        dobs = seq(from=deltobs, to=maxdepth, by=deltobs)
        tee = vector(length = length(dobs) )
        for(i in 1:length(dobs) )
            {
                tee[i] = sum(1/vee[depth<=dobs[i]] )
            }
       
###  time with noise: ?
        t2 = tee + noise*rnorm(length(tee) )

        dd=maxdepth/N;
        G=matrix( rep(0, length=M*N), nrow=M,ncol= N);
        ##  k=1;
        for(i in 1:M)
            {
                G[i, 1:(i) ] = rep(1,length=i)*dd;
            }

        return(list(G=G, tee=tee, t2=t2, depth=depth, vee=vee, N=N, M=M,maxdepth = maxdepth, deltobs=deltobs, noise = noise, M1 =M1) )

    }

