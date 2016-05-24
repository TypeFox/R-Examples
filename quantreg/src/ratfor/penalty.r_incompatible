# Routines to compute the X contribution of the total variation penalty term
# Given a Delaunay triangulation structure as provided by tripack 

# Written in the nearly extinct ratfor language of the lost tribe of Murray Hill
# With luck it can be translated into the dreaded fortran understood by g77.

# Roger Koenker:  Last Modified October 29, 2002.

# The subroutines employ the Renka tripack library and strongly rely on
# the counterclockwise ordering of the elements in the adjacency list.
# the vector bnd is used to identify and neglect edges on the boundary
# In future versions one might want to have the option of making some other
# provisions for these edges.


subroutine penalty(n,m,q,x,y,bnd,tlist,tlptr,tlend,rax,jax,ned,eps,ierr)

integer n,m,q,lp,lpl,ned,ierr
integer bnd(n),tlist(q),tlptr(q),tlend(n),n4(4),p4(4),jax(m)

double precision x(n),y(n),rax(m),eps
double precision x4(4),y4(4),g4(4)

logical orient


#loop over the vertices: i is index of one end of the edge j is the other end
ned = 0
do i=1,n{
	lpl = tlend(i)
	lp  = lpl
	repeat{
		lp = tlptr(lp)
		j  = iabs(tlist(lp))
		if(j > i){
			n4(1) = i
			n4(2) = j
			call fadjs(n4,n,q,tlist,tlptr,tlend)
			if(bnd(i)*bnd(j) == 0){
				ned = ned + 1
				do k = 1,4{
					x4(k) = x(n4(k))
					y4(k) = y(n4(k))
					}
				if(orient(x4,y4)){
					call iswap(1,n4(3),1,n4(4),1)
					call dswap(1,x4(3),1,x4(4),1)
					call dswap(1,y4(3),1,y4(4),1)
					}
				call ggap(x4,y4,g4,eps,ierr)
				if(ierr == 1) return
				call srtpai(n4,1,p4,1,4)
				do k = 1,4{
					rax((ned - 1)*4 + k) = g4(p4(k))
					jax((ned - 1)*4 + k) = n4(p4(k))
					}
				if(ned*4 > m) return
				}
			}
		if(lp == lpl) break
		}
	}
return
end
logical function orient(x,y)
double precision x(4), y(4)
orient = (y(2) -y(1))*(x(3)-x(4))+(x(1)-x(2))*(y(3)-y(4)) > 0
return
end 
subroutine fadjs(n4,n,q,tlist,tlptr,tlend)
# Subroutine to find matching adjacent vertices for the triogram edges
# On input:
#       n4[1:2]  contain the indices of the edge of interest
#       tlist,tlptr,tlendd  is the (tripack) structure describing the triangulation
# On output:
#       n4[3:4]  contains indices of the two adjacent vertices 
#
# Adjacency tlist is in counter-clockwise order so we want to find the two
# vertices that are immediately above and below  n1 in the n0 tlist
#
# Roger Koenker June 4, 2002
#

integer n,q,vp,vpl,v,v0,match
integer n4(4),tlist(q),tlptr(q),tlend(n)

# Check whether edge is on the boundary

match = 0
vpl = tlend(n4(1))
vp = vpl
k = 0
repeat{
       k = k+1
       vp = tlptr(vp)
       v  = tlist(vp)
       if(k>1 & iabs(v) == n4(2)){
               n4(3) = iabs(v0)
               match = 1
               next
               }
       if(match > 0){
               n4(4) = iabs(v)
               break
               }
       v0 = v
       }
return
end

subroutine ggap(x,y,g,eps,ierr)
double precision x(4),y(4),g(4),w(2,4),h(2),D1,D2,eps
# Triogram package:  Roger Koenker  June 4, 2002
# given four (x,y) pairs for an edge compute contribution to the penalty
# ierr returns 1 if the edge is degenerate, ie either determinant is 0.
	D1 = -x(2) * y(1) + x(3) * y(1) + x(1) * y(2) -
	      x(3) * y(2) - x(1) * y(3) + x(2) * y(3)
	D2 = -x(2) * y(1) + x(4) * y(1) + x(1) * y(2) -
	     x(4) * y(2) - x(1) * y(4) + x(2) * y(4)
	if(dabs(D1) < eps | dabs(D2) < eps) {
		ierr = 1
		return
		}
        h(1) =  -(y(1) - y(2)) 
	h(2) =   (x(1) - x(2))
        w(1, 1) = (y(2) - y(3))/D1 - (y(2) - y(4))/D2
        w(2, 1) = (x(3) - x(2))/D1 - (x(4) - x(2))/D2
        w(1, 2) = (y(3) - y(1))/D1 - (y(4) - y(1))/D2
        w(2, 2) = (x(1) - x(3))/D1 - (x(1) - x(4))/D2
        w(1, 3) = (y(1) - y(2))/D1
        w(2, 3) = (x(2) - x(1))/D1
        w(1, 4) = (y(2) - y(1))/D2
        w(2, 4) = (x(1) - x(2))/D2
	do i = 1,4{
		g(i) = h(1)*w(1,i)+h(2)*w(2,i)
		}
ierr = 0
return
end
