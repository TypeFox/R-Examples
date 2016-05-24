
#' @export

Sieve_NPMLE_Switch <- function(Y, p0G, Delta, px=NULL, Grid=NULL, Knot, degree = 3, maxiter = 400, ep = 1e-5){

    #library(splines)
    #library(zoo)

    if(length(px) == 0)
    {
      px = unique(Y[Delta==1]);
    }

    if(length(Grid) == 0)
    {
      Grid = unique(Y[Delta==1]);
    }

    n = length(Y);
    knot =  Knot;
    knots0 = quantile( (Y/max(Y))[Delta==1], probs = seq(1/(knot+1), knot/(knot+1), 1/(knot+1)) );
    Converge_Index1 = 0;  Iter = 0;  Num = 2;  knots = knots0;

    while(Iter == 0 & (Num <= knot)){

      Bt = bs( Y/max(Y), degree = degree, knots = knots, intercept = TRUE);
      nc = dim(Bt)[2];
      # ----  Initial Value  ---- #
      oldlambda.hat = Delta*pnorm(Y, mean= mean(Y))/10; oldbeta.hat = rep(0, nc);

      # ----  E-M  Iteratoion  ---- #
      iter = 1; error=1;

      while( iter <= maxiter & error > ep ){
        # -- E-Step -- #
        temp0 = (Bt%*%oldbeta.hat); Etemp0=exp(temp0);
        temp1 = p0G*exp(Delta*temp0-cumsum(oldlambda.hat*Etemp0));
        q = temp1/(temp1+(1-p0G)*exp(-cumsum(oldlambda.hat)));
        # -- M-setp -- #
        temp2 = cumsum(q[seq(n,1,-1)]); temp2 = temp2[seq(n,1,-1)];
        temp3 = seq(n,1,-1)-temp2;
        Frac = temp2*Etemp0/(temp2*Etemp0+temp3);
        grad = t(Delta*(q-Frac))%*%Bt;
        H = t(Bt)%*%(matrix((Frac^2-Frac)*Delta, n, dim(Bt)[2], byrow=FALSE)*Bt);

        if(det(H) > 0.00001 & det(H) < -0.00001 ) beta.hat = oldbeta.hat - solve(H)%*%t(grad) else {
          beta.hat = oldbeta.hat - solve(H - diag(0.0001, nc))%*%t(grad) };

        temp4 = temp2*exp(Bt %*% beta.hat) + temp3;
        lambda.hat = Delta/temp4;

        if( sum(is.na(lambda.hat)) == 0 & sum(abs(beta.hat) > 650 ) == 0 & sum(lambda.hat > 10e+04)==0 ){
          Iter = 1;

          error = sqrt( sum((oldbeta.hat - beta.hat)^2) + sum((lambda.hat - oldlambda.hat)^2) );
          if(error < ep){Converge_Index1 = 1};

          iter = iter + 1;
          oldbeta.hat = beta.hat;
          oldlambda.hat = lambda.hat;

          } else {
            Iter = 0; knots1 = c(knots0[1:(length(knots0) - Num)], mean( rev(knots0) [1:Num] ) );
            knots = knots1;  Num=Num+1;  break }

      }
    }

    lambda_2.hat = lambda.hat;

  # ----------------------------- #
    knot = Knot;
    knots0 = quantile( (Y/max(Y))[Delta==1], probs = seq(1/(knot+1), knot/(knot+1), 1/(knot+1)) );
    Converge_Index2 = 0; Iter = 0; Num = 2; knots = knots0;

    while(Iter == 0 & (Num <= knot)){

      Bt = bs( Y/max(Y), degree = degree, knots = knots, intercept = TRUE);
      Bt1 = bs( px/max(c(px,Y)), degree = degree, knots = knots, intercept = TRUE);
      nc = dim(Bt)[2];

      oldlambda.hat = Delta*pnorm(Y, mean= mean(Y))/10; oldbeta.hat = rep(0, nc);
      maxiter = 400; ep = 1e-5; iter = 1; error=1;

      while( iter <= maxiter & error > ep ){
        # -- E-Step -- #
        temp0 = -(Bt%*%oldbeta.hat); Etemp0=exp(temp0);
        temp1 = (1-p0G)*exp(Delta*temp0-cumsum(oldlambda.hat*Etemp0));
        q = (p0G*exp(-cumsum(oldlambda.hat)))/( temp1 + p0G*exp(-cumsum(oldlambda.hat)) );
        # -- M-setp -- #
        temp2 = cumsum( (1-q)[seq(n,1,-1)]); temp2 = temp2[seq(n,1,-1)];
        temp3 = seq(n,1,-1)-temp2;
        Frac = -temp2*Etemp0/(temp2*Etemp0+temp3);
        grad = t( Delta*( -(1-q) - Frac ) ) %*% Bt;
        H = t(Bt)%*%(matrix((-Frac^2+Frac)*Delta, n, dim(Bt)[2], byrow=FALSE)*Bt);

        if(det(H) > 0.00001 & det(H) < -0.00001 ) beta.hat = oldbeta.hat - solve(H)%*%t(grad) else {
          beta.hat = oldbeta.hat - solve(H - diag(0.0001, nc))%*%t(grad) };

        temp4 = temp2*exp(- Bt %*% beta.hat) + temp3;
        lambda.hat = Delta/temp4;

        if( sum(is.na(lambda.hat)) == 0 & sum(abs(beta.hat) > 650 ) == 0 & sum(lambda.hat > 10e+04)==0 ){
          Iter = 1;
          error = sqrt( sum((oldbeta.hat - beta.hat)^2) + sum((lambda.hat - oldlambda.hat)^2) );
          if(error < ep){Converge_Index2 = 1};
          iter = iter + 1;
          oldbeta.hat = beta.hat;
          oldlambda.hat = lambda.hat;
        } else {
          Iter = 0; knots1 = c(knots0[1:(length(knots0) - Num)], mean( rev(knots0) [1:Num] ) );
          knots = knots1;  Num=Num+1;  break }
      }

    }

    lambda_1.hat = lambda.hat;

    Lambda2 = cumsum(lambda_2.hat);
    Lambda1 = cumsum(lambda_1.hat);
    Lamb1 = apply( matrix(Grid, ncol=1), 1, function(x){ if(sum(Y <= x) >= 1){ max( Lambda1[Y <= x] )} else {0} } );
    Lamb2 = apply( matrix(Grid, ncol=1), 1, function(x){ if(sum(Y <= x) >= 1){ max( Lambda2[Y <= x] )} else {0} } );

   Result = list( lamb1.hat = lambda_1.hat,
                  lamb2.hat = lambda_2.hat,
                  Lamb1 = Lamb1,
                  Lamb2 = Lamb2,
                  # beta.est = beta.est,
                  Converge = c(Converge_Index1, Converge_Index2)
                );

 return(Result)

}
