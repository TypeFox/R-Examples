subroutine afun(fy,xispd,tpm,epsilon,n,nstate,wrk,xlc,alpha)
implicit double precision(a-h,o-z)
dimension wrk(nstate), xispd(nstate), xlc(n)
dimension fy(nstate,n), tpm(nstate,nstate), alpha(nstate,n)

# Set some constants
one  = 1.d0
zero = 0.d0

# Set the value to give to ``log-likelihood constant'', xlc(...)
# if this is indeterminate --- i.e. less than epsilon.
# Possible choices: -1, 1, or epsilon.
dummy = epsilon

# Update the initial alpha.
tsum = zero
do j = 1,nstate {
	wrk(j) =  fy(j,1)*xispd(j)
	tsum = tsum + wrk(j)
}

if(tsum < epsilon) {
	xlc(1) = dummy
	do j = 1,nstate {
		alpha(j,1) = one/nstate
	}
}
else {
	xlc(1) = tsum
	do j = 1,nstate {
		alpha(j,1) = wrk(j)/tsum
	}
}

# Run through the remaining n-1 of the alphas (recursing!).
do kt = 2,n {
	tsum = zero
	ktm = kt - 1
	do j = 1,nstate {
		wrk(j) = zero
		do i = 1,nstate {
			wrk(j) = wrk(j) + alpha(i,ktm)*tpm(i,j)
		}
		wrk(j) = fy(j,kt)*wrk(j)
		tsum = tsum + wrk(j)
	}
	if(tsum < epsilon) {
		xlc(kt) = dummy
		do j = 1,nstate {
			alpha(j,kt) = one/nstate
		}
	}
	else {
		xlc(kt) = tsum
		do j = 1,nstate {
			alpha(j,kt) = wrk(j)/tsum
		}
	}
}

return
end
