lookback2z = function(t = 1, steps = 10, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Gyr"){
    temp = function(t, steps, H, M, L, K, units){
        if(units=="Gyr"){
            t = t/(3.16887646E-17)
        }
        zseq = seq(0,50000,len=100)
        tseq = lookback(z=zseq, H=H, M=M, L=L, K=K, units="s")
        zmid = which.min(abs(t-tseq))
        for(i in 1:steps){
            zmin = max((zmid-1),1)
            zmax = min((zmid+1),100)
            zseq = seq(zseq[zmin],zseq[zmax],len=100)
            tseq = lookback(z=zseq, H=H, M=M, L=L, K=K, units="s")
            zmid = which.min(abs(t-tseq))
        }
        return(zseq[zmid])
    }
    return(unlist(lapply(t,temp,steps=steps,H=H,M=M,L=L,K=K,units=units)))
}

