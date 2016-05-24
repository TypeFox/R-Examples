subroutine recurse(fy,xispd,tpm,nreps,epsilon,lns,nstate,nis,cis,
                   wrk,xlc,ntot,nxi,alpha,beta,gamma,xi,xisum)

implicit double precision(a-h,o-z)
logical cis
dimension xispd(nstate,nis), xlc(ntot), lns(nreps)
dimension tpm(nstate,nstate), wrk(nstate,nstate)
dimension fy(nstate,ntot), alpha(nstate,ntot), beta(nstate,ntot)
dimension gamma(nstate,ntot), xi(nstate,nstate,nxi), xisum(nstate,nstate)

# Set zero and one:
zero = 0.d0
one  = 1.d0

# Run through the replicates.
kstop = 0
do k = 1,nreps {
	n = lns(k)
        nm1 = n - 1
	kstart = 1 + kstop
        if(cis) {
            kis = 1
        } else {
            kis = k
        }

# Update the alpha's.
	call afun(fy(1,kstart),xispd(1,kis),tpm,epsilon,n,nstate,wrk,
                  xlc(kstart),alpha(1,kstart))

# Update the beta's.
	call bfun(fy(1,kstart),xispd(1,kis),tpm,epsilon,n,nstate,wrk,
                  beta(1,kstart))

# Update the gamma's.
	call gfun(alpha(1,kstart),beta(1,kstart),epsilon,n,
                  nstate,wrk,gamma(1,kstart))

# Update the xi's.
	call xfun(alpha(1,kstart),beta(1,kstart),fy(1,kstart),
                  tpm,epsilon,n,nstate,nm1,wrk,xi(1,1,kstart-k+1))

# Increment kstop.
kstop = kstop + lns(k)
}

kstop = kstop - nreps
# Sum up the xi's.
do i = 1,nstate {
	do j = 1,nstate {
		xisum(i,j) = zero
		do k = 1,kstop {
			xisum(i,j) = xisum(i,j) + xi(i,j,k)
		}
	}
}

return
end
