setDesignG<-function()
    {
        #####  create a 94 by 256 design matrix for 2D simple tomography
        m = 256;

        G = matrix(rep(0, 94*m), ncol=m, nrow=94)

        for(i in 1:16)
            {
                for(j in seq(from= (i-1)*16+1, to= i*16))
                    {
                        G[i,j] = 1.;
                    }
            }
### % Design matrix for the row scan
        for(i in 1:16)
            {
                for(j in seq(from=i, by=16, to=240+i))
                    {
                        G[i+16,j] = 1.;
                    }
            }

### % G matrix for the SW to NE diagonal scan, upper part
        for(i in 1:16)
            {
                
                for(j in 0:(i-1))
                    {
                        G[i+32,i+j*15] = sqrt(2.);
                    }
            }

### % G matrix for the SW to NE diagonal scan, lower part
        for(i in 1:15)
            {
                for(j in (0:(15-i)))
                    {
                        G[i+48,(i+1)*16+j*15] = sqrt(2.);
                    }
            }

### % G matrix for the NW to SE diagonal scan, lower part
        for(i in 1:16)
            {
                for(j in 0:(i-1))
                    {
                        G[i+63,17-i+17*j] = sqrt(2.);
                    }
            }

### % G matrix for the NW to SE diagonal scan, upper part
        for(i in 1:15)
            {
                for(j in 0:(15-i))
                    {
                        G[i+79,(i*16)+1+17*j] = sqrt(2.);
                    }
            }

        return(G)
    }
