subroutine sknotl(x,n,knot,k)
implicit double precision(a-h,o-z) 
integer     n,k,ndk,j
double precision x(n),knot(n+6),a1,a2,a3,a4


     # Allocate knots acording to the cutoffs given below


	# Cutoff constants

	 a1 = log(50e0)/log(2e0)  ; a2 = log(100e0)/log(2e0)
	 a3 = log(140e0)/log(2e0) ; a4 = log(200e0)/log(2e0)

	# Cutoff Criteria

        if(n<50)                    { ndk = n }
        else if (n>=50  & n<200)    { ndk = 2.**(a1+(a2-a1)*(n-50.)/150.)  }
        else if (n>=200 & n<800)    { ndk = 2.**(a2+(a3-a2)*(n-200.)/600.)  }
        else if (n>=800 & n<3200)   { ndk = 2.**(a3+(a4-a3)*(n-800.)/2400.)  }
        else if (n>=3200)           { ndk = 200. + (n-3200)**.2 }
		 
		 k = ndk + 6

     
     # Allocate Knots  ( note no account is taken of any weighting vector )

	do j=1,3   { knot(j) = x(1) }
	do j=1,ndk { knot(j+3) = x( 1 + (j-1)*(n-1)/(ndk-1) ) }
	do j=1,3   { knot(ndk+3+j) = x(n) }

return
end
