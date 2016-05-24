# Subroutine to list the r choose n subsets of {1,2,...,r} in an order such that
# adjacent subsets have only one element swapped.  Algorithm modeled on the pascal
# algorithm of Limin Xiang and Kazuo Ushijima (2001) "On O(1) Time Algorithms for
# Combinatorial Generation," Computer Journal, 44(4), 292-302. 
# http://comjnl.oxfordjournals.org/cgi/reprint/44/4/292

# Translated into ratfor:  12 February, 2008 R. Koenker.

subroutine combin(r,n,m,a,c,e,Last)

integer r,n,m,t,k,j,M0,Mj
integer A(n,m),c(r),e(r),Last(r)
logical odd
	
	M0 = r-n
	t = n+1
	k = 1
	j = 0
	c(1) = 0
	repeat{
		j = j + 1
		c(j) = j
		e(j) = j - 1
		if(odd(j))
			Last(j) = M0 + j
		else
			Last(j) = j + 1
		if(j == n) break
		}
	do i = 1,n
		A(i,1) = c(i)
	if(n < r) {
	repeat {
		k = k + 1
		S = c(j)
		Mj = M0 + j
		e(n+1) = n
		if(odd(j)){
			if(c(j) == Mj) {
				c(j) = c(j-1) + 1
				Last(j+1) = c(j) + 1
				}
			else
				c(j) = c(j) + 1
			}
		else { 
			if(c(j) == c(j-1) + 1)
				c(j) = Mj
			else {
				Last(j+1) = c(j)
				c(j) = c(j) - 1
				}
			}
		if(c(j) == Last(j)) {
			Last(j) = S
			e(j+1) = e(j)
			e(j) = j-1
			}
		if( (j < n) & (c(j) == Mj)) {
			t = j
			j = e(t+1)
			e(t+1) = t
			}
		else {
			if(t == j) t = t + 1
			if(t < e(n+1)) j = t
			else j = e(n+1)
			}
		do i = 1,n
			A(i,k) = c(i)
		if(j == 0) break
		}
	}
return
end


logical function odd(j) 
integer j
odd = (mod(j,2) == 1)
return
end
