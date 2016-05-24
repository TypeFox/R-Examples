box_qp_f <-
function (Q, u,b,rho,Maxiter=10^3,tol=10^-4) {
qq=nrow(Q);
obj_vals<-rep(0,Maxiter);

mode(Q)="double";mode(u)="double";mode(b)="double";
mode(Maxiter)="integer";mode(tol)="double";
mode(rho)="double";
mode(obj_vals)="double";
mode(qq)="integer";

junk<-.Fortran("box_qp_f",Q,uu=u,b,rho,Maxiter,tol,qq,grad_vec=double(qq),PACKAGE="dpglasso")

return(list(grad_vec=junk$grad_vec,u=junk$uu))

                                                        }
