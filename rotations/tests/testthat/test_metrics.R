rs<-rcayley(20)
Rs<-genR(rs)

context("Distance metrics")
expect_equal(rot.dist(Rs,method='intrinsic'),abs(rs))
expect_equal(rot.dist(Rs,method='projected'),sqrt(8)*sin(abs(rs)/2))
expect_equal(rot.dist(Rs,method='intrinsic'),mis.angle(Rs))

S<-genR(rvmises(1))
Rs<-genR(rs,S=S)
expect_equal(rot.dist(Rs,S,method='intrinsic'),abs(rs))
expect_equal(rot.dist(Rs,S,method='projected'),sqrt(8)*sin(abs(rs)/2))
