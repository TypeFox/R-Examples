NPMLE.Indep <-
function(x.trunc,y.trunc,x.fix = median(x.trunc), y.fix = median(y.trunc),plotX = TRUE){

m=length(x.trunc)

x.ox=sort(unique(x.trunc))
y.oy=sort(unique(y.trunc))
nx=length(x.ox);ny=length(y.oy)

dHx1=1;Hxn=0;Ay_1=0;dAyn=1 ### initial constraint ####

###############  Starting values ####################
dHx_lyn=numeric(nx);dAy_lyn=numeric(ny)

for(i in 1:nx){
  tx=x.ox[nx-i+1]
  Rx=sum( (x.trunc<=tx)&(y.trunc>=tx) )
  dHx_lyn[nx-i+1]=sum(tx==x.trunc)/Rx
}

for(i in 1:ny){
  ty=y.oy[i]
  Ry=sum( (x.trunc<=ty)&(y.trunc>=ty) )
  dAy_lyn[i]=sum(ty==y.trunc)/Ry
}

dL_lyn=log(c(dHx_lyn[-1],dAy_lyn[-ny]))

############  Fitting Semi-Survival Independence copula  #################
l.indep=function(dL){
    dHx=c(dHx1,exp(dL[1:(nx-1)]));dAy=c(exp(dL[nx:(nx+ny-2)]),dAyn)
    Hx=c(rev(cumsum(rev(dHx)))[-1],Hxn);Ay_=c(Ay_1,cumsum(dAy)[-ny])

    prop=0
    for(i in 1:nx){
      temp_y=(y.oy>=x.ox[i])
      dHxi=dHx[i];dAyi=dAy[temp_y]
      prop=prop+sum(exp(-Ay_[temp_y])*exp(-Hx[i])*dAyi*dHxi)
    }

    l=-m*log(prop)
    for(i in 1:m){
      x_num=sum(x.ox<=x.trunc[i]);y_num=sum(y.oy<=y.trunc[i])
      Hxi=Hx[x_num];Ay_i=Ay_[y_num]
      dHxi=dHx[x_num];dAyi=dAy[y_num]
      l=l-Hxi-Ay_i+log(dHxi)+log(dAyi)
    }
    -l
}

res=nlm(l.indep,p=dL_lyn,hessian=TRUE)
dL=res$estimate;conv=res$code
dHx=c(dHx1,exp(dL[1:(nx-1)]))
dAy=c(exp(dL[nx:(nx+ny-2)]),dAyn)
Hx_=rev(cumsum(rev(dHx)));Ay=cumsum(dAy)
Fx=exp(-Hx_);Sy=exp(-Ay)
Hx=c(Hx_[-1],Hxn)
Ay_=c(Ay_1,Ay[-ny])
V_dL=solve(res$hessian)
V=diag(exp(dL))%*%V_dL%*%diag(exp(dL))
ML=-res$minimum
iter=res$iterations
Grad=res$gradient
Min_eigen=min(eigen(res$hessian)$value)

Hx_fix = Ay_fix = numeric(length(x.fix))
SE_Hx = SE_Ay = numeric(length(x.fix))
for (i in 1:length(x.fix)) {
    temp.x = c(x.ox >= x.fix[i])
    Hx_fix[i] = Hx_[max(sum(x.ox <= x.fix[i]), 1)]
    SE_Hx[i] = sqrt((temp.x[-1]) %*% V[1:(nx - 1), 1:(nx - 1)] %*% (temp.x[-1]))
}
for (i in 1:length(y.fix)) {
    temp.y = c(y.oy <= y.fix[i])
    Ay_fix[i] = Ay[max(sum(y.oy <= y.fix[i]), 1)]
    SE_Ay[i] = sqrt((temp.y[-ny]) %*% V[nx:(nx+ny-2),nx:(nx+ny-2)] %*% (temp.y[-ny]))
}
Fx_fix = exp(-Hx_fix)
Sy_fix = exp(-Ay_fix)
SE_Fx = Fx_fix * SE_Hx
SE_Sy = Sy_fix * SE_Ay

if (plotX == TRUE) {
  plot(x.ox, Fx, type = "s", xlab = "x", ylab = "Fx(x): Distribution function of X")
}

list(Hx = c(estimate=Hx_fix, Hx_se = SE_Hx), 
Ay = c(estimate=Ay_fix, Ay_se = SE_Ay),
Fx = c(estimate=Fx_fix, Fx_se = SE_Fx), 
Sy = c(estimate=Sy_fix, Sy_se = SE_Sy),
ML=ML, convergence = conv,iteration=iter,Grad=sum(Grad^2),MinEigen=Min_eigen
)

}
