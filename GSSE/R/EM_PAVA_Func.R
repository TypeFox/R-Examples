
#' @export

EM_PAVA_Func <- function(q, x, delta, timeval, p, ep=1e-5, maxiter = 400){

    # library(zoo)
    # library(Iso)

    n = length(x); # sample size

    timeval.new = c(timeval,x);
    timeval.new = sort(timeval.new);               # sorting all time points
    timeval.indices = match(timeval,timeval.new)   # indices of timeval points
    events.indices = match(x,timeval.new)          # indices of event times x

    numt.new = length(timeval.new)
    Fest = matrix(0,nrow=p,ncol=numt.new);
    Fest0 = Fest;

    for(l in 1:p) Fest0[l,] = pnorm(timeval.new,mean=mean(timeval.new))

    Fest = Fest0;  tol = 1;   count = 0;
    all.deaths = (sum(delta)==n);

    while( tol > ep & count < maxiter ){

      count = count+1;
      u = matrix(0, nrow=n, ncol = numt.new);
      v = u;

      for(i in 1:n){
        u.numerator = q[1,i] * Fest0[1,];
        u.denominator = q[1,i] * Fest0[1,] + q[2,i] * Fest0[2,];

        u[i,] = u.numerator/u.denominator;
        u[i,] = na.approx( u[i,], na.rm = FALSE, rule = 2 );

        v.numerator = q[1,i] * (1-Fest0[1,]);
        v.denominator = q[1,i] * (1-Fest0[1,]) + q[2,i] * (1-Fest0[2,]);

        v[i,] = v.numerator / v.denominator;
        v[i,] = na.approx(v[i,],na.rm=FALSE,rule=2);
      }

      s = matrix(0, nrow = p, ncol = numt.new);

      if(all.deaths == FALSE){
        w.tmp = matrix(0, nrow = n, ncol = numt.new);

        for(j in 1:numt.new){
          indicator = rep(0,n);
          indicator[ x > timeval.new[j] ] = 1;
          bottom = q[1,] * (1 - Fest0[1, events.indices]) + q[2,] * (1 - Fest0[2, events.indices]);
          top = q[1,] * (1 - Fest0[1, j]) + q[2,] * (1 - Fest0[2, j]);
          temp = (1 - delta) * top / bottom;
          temp[temp == Inf] = 0;

          if( sum(is.na(temp)) == length(temp) ){ temp = rep(0, length(temp) ) }

          if( sum(is.na(temp)) == ( length(temp) - 1 ) ){
            temp = na.locf(temp, na.rm = FALSE)
          }

          temp = na.approx(temp, na.rm = FALSE, rule = 2);

          indicator[ x <= timeval.new[j] ] = temp[ x <= timeval.new[j] ];
          w.tmp[, j] = 1 - indicator;
        }

        for(j in 1:numt.new){
          s[1,j] = sum( u[,j] * w.tmp[,j] ) + sum( v[,j] * (1-w.tmp[,j]) );
          s[2,j] = sum( (1-u[,j]) * w.tmp[,j] ) + sum( (1-v[,j]) * (1-w.tmp[,j]) );
          Fest[1,j] = sum( u[,j] * w.tmp[,j] )/s[1, j];
          Fest[2,j] = sum( (1-u[,j]) * w.tmp[,j] )/s[2, j];
        }

      } else {
        for(j in 1:numt.new){
          tem = rep(0,n);
          tem[ x <= timeval.new[j] ] = 1;
          s[1, j] = sum(u[,j] * tem) + sum( v[,j] * (1-tem) );
          s[2, j] = sum((1-u[,j]) * tem) + sum( (1-v[,j]) * (1-tem) );
          Fest[1, j] = sum(u[,j]*tem)/s[1,j];
          Fest[2, j] = sum((1-u[,j])*tem)/s[2,j];
        }
      }

      for(l in 1:p){ Fest[l,] = pava(Fest[l,], s[l,]) }

      tol = data.frame(apply(Fest-Fest0,1,abs));
      tol = apply(tol,2,mean);
      tol = sum(tol);
      Fest0 = Fest;
    }

    Fest.all = Fest;
    Fest = Fest[,timeval.indices];
    list(Fest=Fest,Fest.all=Fest.all);
}
