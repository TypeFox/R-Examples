LPiTrack <-
function(xy_mat, m = c(3, 5, 15), p = 10)
{
    if(ncol(xy_mat) != 2)
        stop("\n A matrix with two columns is required \n")
 r.list <- scale(Mod(diff(xy_mat[,1]) + 1i*(diff(xy_mat[,2]))))

 r.coefs <- LPTime(cbind(r.list[-length(r.list)], diff(r.list)),
                   exo=NULL, m = 3, p)

 dif.mat <- diff(xy_mat)
 dif.mat2 <- diff(dif.mat)
 tran2.coefs <-LPTime(dif.mat,exo = NULL ,m[1],p)
 tran3.coefs <- LPTime(dif.mat2, exo=NULL, m[1], p)

    comoment.r <-  LP.comoment(r.list[-length(r.list)], diff(r.list), m[2])
    comoment.coefs <- rbind(LP.comoment(xy_mat[,1],  xy_mat[,2], m[2]),
                            LP.comoment(dif.mat[,1], dif.mat[,2]), m[2])


moments.r <- LP.moment(r.list, m[3])
moments.difr <- LP.moment(diff(r.list), m[3])

moment.coefs <-rbind( apply(xy_mat, 2, LP.moment, m[3]),
                      apply(dif.mat, 2, LP.moment, m[3]),
                      apply(dif.mat2, 2, LP.moment), m[3])

c(r.coefs,
  tran2.coefs,
  tran3.coefs,
  comoment.r,
  comoment.coefs,
  moment.coefs,
  moments.r,
  moments.difr)
}
