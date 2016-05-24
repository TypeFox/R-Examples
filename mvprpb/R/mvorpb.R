mvorpb <-
function( dim.p , m.tgt , v.tgt , n.itr , it.rg ){

#	Function test.orthant
#		Evaluate an orthant probability
#	Arguments:
#		dim.p: Dimension to evaluate (Integer, Scalar)
#		m.tgt: Mean (Real vector, length: dim.p)
#		v.tgt: Covariance Matrix (Real square matrix of size dim.p)
#		n.itr: Number of intervals for numerical integration
#		it.rg: Maximum point of the numerical integration range
#
#	In this function, two functions,  f.num.cumint.wdiff and f.comb.gen.01, are called, which are
#	defined in this file.


f.num.cumint.wdiff <-  function( f.val , diff.val , int.wid ){
#	Numericl Integration using Derivitives at Grid Points
#	This function is called from test.orthant

	val.low <-   f.val[ - length(f.val) ]
	val.up  <-   f.val[ - 1 ]
	val.accm <-  c( 0 , 0.5*cumsum( ( val.low + val.up )*int.wid ) 
				- cumsum(diff(diff.val)*int.wid^2)/12
		)
	return(val.accm)
}	# enf of function f.num.cumint.wdiff

f.comb.gen.01 <-  function( q.in , candi , pre.val ){
#	Generate Combinations
#	The technique of recursive call is used.
#	This function is called from test.orthant

	if(q.in<=1){
		rep.count <-  length(candi)
		ret.val <-  rbind( matrix(rep(pre.val,rep.count),ncol=rep.count) , matrix(candi,nrow=1) )
	}else{
	if(length(candi)<=q.in){
		ret.val <-  matrix( c(pre.val ,candi) , ncol=1)
	}else{
		ret.val <- cbind(	f.comb.gen.01(q.in-1 , candi[-1] , c(pre.val,candi[1])) ,
				f.comb.gen.01(q.in   , candi[-1] , pre.val )   )
	}}   # end if

	return( ret.val )
	} # enf of function f.comb.gen.01



#
#	Step 1
#		Initial Setting

	scale.go <-  1 / sqrt( diag( v.tgt ) )
	m.wrk <-  m.tgt * scale.go 
	v.go <-  diag( scale.go ) %*% v.tgt %*% diag( scale.go )

#	Cholesky Decomposition
#			t(b.vec) %*% b.vec === v.go

	v.go <-  0.5 * ( v.go + t(v.go) )
	b.vec <-  chol( v.go )


	m.a  <-  solve(t(b.vec))%*%m.wrk
	c.vec.pre <-  t(solve( b.vec ))

	int.plane.vec <-  apply( b.vec , 1 , sum )
	int.plane.vec <-   int.plane.vec / sqrt(sum(int.plane.vec^2))
	c.vec.adj <-  t(c.vec.pre) %*% int.plane.vec
	c.vec <-  sweep( c.vec.pre , 2 , c.vec.adj , "/")


	intg.lat.p <-  ( 0:n.itr )* it.rg / n.itr
	intg.lat.d <-  diff( intg.lat.p )

	wk.mat <-  matrix( NA , nrow=length(intg.lat.p) , ncol=(2^dim.p-1) )
	wk.mat.diff <- wk.mat

	wk.orthpnt.lensq <- rep( NA ,  2^dim.p-1)
	wk.cross.mean <- rep( NA ,  2^dim.p-1)

#
#	Step 2
#		Evaluate vector w^J for each subspace

w.base <- matrix(NA , ncol=dim.p , nrow=dim.p )
w.plane.vec <- int.plane.vec
wk.orthpnt <- matrix(NA,nrow=dim.p,ncol=(2^dim.p))
wk.orthpnt[,2^dim.p-1] <- int.plane.vec
w.pos.comb <- 0
w.dim <- 0
i.w <- 0
w.stack <- rep(NA , dim.p)

while(w.dim >= 0){
	i.w <- i.w + 1
	w.u <- b.vec[ , i.w ]

	if(w.dim>0){
		w.seq <- 1:w.dim
		w.coef1 <- t( w.base[,w.seq,drop=FALSE]) %*% matrix( w.u , ncol=1)
		w.u <- w.u - w.base[,w.seq,drop=FALSE] %*% matrix(w.coef1 , ncol=1 )
      }
		w.u <- w.u / sqrt(sum(w.u^2))
		w.coef2 <- sum( w.plane.vec * w.u )
		w.pos.comb.new <- w.pos.comb + 2^(i.w - 1)
		wk.orthpnt[ ,  2^dim.p -1- w.pos.comb.new ] <- w.plane.vec - w.coef2 * w.u

	if( i.w < dim.p ){
		w.plane.vec <-  wk.orthpnt[ ,  2^dim.p -1- w.pos.comb.new ] 
		w.dim <- w.dim + 1
		w.stack[ w.dim ] <- i.w
		w.base[ , w.dim ] <- w.u
		w.pos.comb <- w.pos.comb.new

	} else {
		w.dim <- w.dim - 1
		if(w.dim > -1){
		i.w <- w.stack[ w.dim +1]
		w.pos.comb <- w.pos.comb - 2^(i.w-1)
		w.plane.vec <- wk.orthpnt[,2^dim.p-1-w.pos.comb]
		}
	} # end if
	}	# end while

wk.orthpnt <- sweep(wk.orthpnt , 2,c(t(wk.orthpnt)%*% int.plane.vec ),"/"  )

#	Step 3
#		Evaluation of Two Dimensional Nodes

	for( w.i in seq( 2 , dim.p )){
	for( w.j in seq( 1 , w.i - 1 )){
	w.vi <-  c.vec[ , w.i ]
	w.vj <-  c.vec[ , w.j ]
	w.v.dif <-  w.vi - w.vj
	w.v.dif.len.s <-  sum( w.v.dif^2 )
	w.v.dif.len  <-  sqrt( w.v.dif.len.s )

	wk.1 <-  (-1) * sum( w.vj *w.v.dif )/w.v.dif.len.s
	w.mean.posval <-  sum(w.v.dif * m.a )/w.v.dif.len
	w.vec.base <-   wk.1*w.vi + (1-wk.1)*w.vj 

	w.speed.i <-  sum((w.vi - w.vec.base)*w.v.dif)/w.v.dif.len
	w.speed.j <-  sum((w.vj - w.vec.base)*w.v.dif)/w.v.dif.len
	w.pos.i <-  intg.lat.p * w.speed.i - w.mean.posval
	w.pos.j <-  intg.lat.p * w.speed.j-  w.mean.posval

	w.prb.sec <-  (pnorm(  w.pos.i ) - pnorm(  w.pos.j  ) )

	w.prb.sec.diff <- (
			w.speed.i * dnorm( w.pos.i )
			-
			w.speed.j * dnorm( w.pos.j )
			)

	w.respos <-  2^(w.i-1) + 2^(w.j - 1)
	wk.mat[    , w.respos ] <-  w.prb.sec
	wk.mat.diff[ , w.respos ] <-  w.prb.sec.diff

	wk.orthpnt.lensq[ w.respos] <-  sum(w.vec.base^2)
  	wk.cross.mean[ w.respos] <- sum(  w.vec.base * m.a )

	}}	#	next w.i , w.j

#	Step 4
#		Evaluation of Intermediate Dimension

if( dim.p > 2){
for(w.dim in seq(3,dim.p)){

	w.pat.accm  <-  f.comb.gen.01(w.dim , seq(dim.p) , NULL)

for(w.i in seq(1 , ncol(w.pat.accm))){

	w.pat.go <-  w.pat.accm[ , w.i ]
	w.mat.x <-  c.vec[ , w.pat.go ]

	w.prb.accm <-  0
	w.prb.accm.diff <-  0

	w.highdim.pos <-  sum( 2^(w.pat.go - 1 ) )
	w.highdim.orthpnt <- wk.orthpnt[ , w.highdim.pos ]

	wk.orthpnt.lensq[ w.highdim.pos ] <- sum( w.highdim.orthpnt^2 )
	wk.cross.mean[ w.highdim.pos ] <-sum(  w.highdim.orthpnt * m.a )


for(w.j in seq(w.dim)){
#		Evaluate each integral in Step 4

	w.valno.ineq <-  w.pat.go[ - w.j ]
  	w.lowdim.pos <-  sum(2^(w.valno.ineq-1))
	w.lowdim.base <-  wk.orthpnt[,w.lowdim.pos]

	w.shiftspeed <- sqrt( wk.orthpnt.lensq[w.lowdim.pos] - wk.orthpnt.lensq[ w.highdim.pos ]  )

	if(abs(w.shiftspeed)>1e-6){
	w.mean.posval <-  (wk.cross.mean[w.lowdim.pos]-wk.cross.mean[ w.highdim.pos ])/w.shiftspeed

	w.val.lowdim <-  wk.mat[,w.lowdim.pos]
	w.phi.arg <- intg.lat.p*w.shiftspeed-w.mean.posval
	w.phi <- dnorm(w.phi.arg )
	w.prb.intg.pre   <-  w.shiftspeed * w.val.lowdim * w.phi
	w.prb.intg.diff.pre <-   w.shiftspeed * w.phi * (wk.mat.diff[ ,w.lowdim.pos ]
				- w.phi.arg * w.shiftspeed * w.val.lowdim 
				)

	w.prb.sec <-   f.num.cumint.wdiff( w.prb.intg.pre , w.prb.intg.diff.pre ,  intg.lat.d )

	if(sum(b.vec[,w.pat.go[w.j]] * w.highdim.orthpnt)>0){
	w.prb.accm <-  w.prb.accm + w.prb.sec
	w.prb.accm.diff  <-  w.prb.accm.diff  + w.prb.intg.pre
	} else {
	w.prb.accm <-  w.prb.accm -w.prb.sec
	w.prb.accm.diff  <-  w.prb.accm.diff  - w.prb.intg.pre

	}
	}	# end if
	}	# next w.j

	wk.mat[ , w.highdim.pos ] <-  w.prb.accm
	wk.mat.diff[ , w.highdim.pos ] <-  w.prb.accm.diff

}	# next w.i
}}	# next w.dim  end if

#
#	Step 5
#		Highest Dimension

	wk.respos <-  2^dim.p - 1
	wk.3 <-  wk.mat[ , wk.respos ]
	wk.mean.posval <-  sum( int.plane.vec * m.a )

	ww.phi.arg <- intg.lat.p - wk.mean.posval
	w.phi <- dnorm( ww.phi.arg )
	
	w.prb.intg.pre   <-  wk.3* w.phi
	w.prb.intg.diff.pre <-   w.phi * (wk.mat.diff[ , wk.respos ] - ww.phi.arg *  wk.3 )

	res.prb <-   f.num.cumint.wdiff( w.prb.intg.pre ,w.prb.intg.diff.pre  , intg.lat.d)
	ret.val <-  res.prb[ length(res.prb) ]
	attr(ret.val , "error-itg-rg" ) <- pnorm( - (it.rg - w.mean.posval) )
	return( ret.val  )
	}	# enf of function test.orthant

