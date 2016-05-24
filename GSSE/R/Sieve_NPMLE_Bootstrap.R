
#' @export

Sieve_NPMLE_Bootstrap <- function(fam_ID, Y0, Delta0, p0G0, fix_t1, fix_t2, Grid, Knot, degree =3, Bn, maxiter = 400, ep = 1e-5 )
{

    #require(zoo);
    #require(splines);

   Data = cbind(Y0, Delta0, p0G0);
   fam_ID = unique(fam_ID);

   if(length(Grid) == 0)
   {
     Grid = sort( unique(Y0[Delta0==1]) );
   }

  # --------------------------------------------------------------------------------- #

   Boot.L1 = Boot.L2 = matrix(0, ncol=length(Grid), nrow=Bn);
   Fix_t_F1 = matrix(0, ncol=length(fix_t1), nrow=Bn);
   Fix_t_F2 = matrix(0, ncol=length(fix_t2), nrow=Bn);

   Bsim = 1;

  while( Bsim <= Bn){

    B.Ind = sample( seq(1,length(fam_ID)), size=length(fam_ID), replace = TRUE);
    B.fam_ID = fam_ID[B.Ind];
    B.Data = NULL;

    for(i in 1:length(B.fam_ID) ){
      B.Data = rbind( B.Data, Data[fam_ID == B.fam_ID[i],] );
    }

    n = dim(B.Data)[1];

    OY = B.Data[,1];     ind = order(OY); Y = OY[ind];
    ODelta = B.Data[,2]; Delta = ODelta[ind];
    Op0G = B.Data[,3];   p0G = Op0G[ind];
    #Ofamid = Data[,4];   famid = Ofamid[ind];

    knot = Knot;
    knots0 = quantile( (Y/max(Y))[Delta==1], probs = seq(1/(knot+1), knot/(knot+1), 1/(knot+1)) );
    Converge_Index1 = 0; Iter = 0; Num = 2; knots = knots0;

    while(Iter == 0 & (Num <= knot)){

      Bt = bs( Y/max(Y), degree = degree, knots = knots, intercept = TRUE);
      nc = dim(Bt)[2];
      oldlambda.hat = pnorm(Y, mean= mean(Y))/10; oldbeta.hat = rep(0, nc);
      maxiter = 400; ep = 1e-5; iter = 1; error=1;

      while( iter <= maxiter & error > ep ){
        temp0 = (Bt%*%oldbeta.hat); Etemp0=exp(temp0);
        temp1 = p0G*exp(Delta*temp0-cumsum(oldlambda.hat*Etemp0));
        q = temp1/(temp1+(1-p0G)*exp(-cumsum(oldlambda.hat)));
        q[is.na(q)] = 0;

        temp2 = cumsum(q[seq(n,1,-1)]); temp2 = temp2[seq(n,1,-1)];
        temp3 = seq(n,1,-1) - temp2;
        Frac = temp2*Etemp0/(temp2*Etemp0+temp3);
        grad = t(Delta*(q-Frac)) %*% Bt;
        H = t(Bt)%*%(matrix((Frac^2-Frac)*Delta, n, dim(Bt)[2], byrow=FALSE)*Bt);

        if(det(H) > 0.00001 & det(H) < -0.00001 ) beta.hat = oldbeta.hat - solve(H)%*%t(grad) else {
           beta.hat = oldbeta.hat - solve(H - diag(0.0001, nc))%*%t(grad) };

        temp4 = temp2*exp(Bt %*% beta.hat) + temp3;
        lambda.hat = Delta/temp4;

        if( sum(is.na(lambda.hat)) == 0 & sum(abs(beta.hat) > 650 ) == 0 & sum(lambda.hat > 10e+04)==0 ){
          Iter = 1 } else {
            Iter = 0; knots1 = c(knots0[1:(length(knots0) - Num)], mean( rev(knots0) [1:Num] ) );
            knots = knots1;  Num=Num+1;  break }

        error = sqrt( sum((oldbeta.hat - beta.hat)^2) + sum((lambda.hat - oldlambda.hat)^2) );
        if(error < ep){Converge_Index1 = 1};
        iter = iter + 1;
        oldbeta.hat = beta.hat;
        oldlambda.hat = lambda.hat;
      }
    }

    F_lambda_2.hat = lambda.hat;
  # ------------------------------ #

    knot =  Knot;
    knots0 = quantile( (Y/max(Y))[Delta==1], probs = seq(1/(knot+1), knot/(knot+1), 1/(knot+1)) );
    Converge_Index2 = 0; Iter = 0; Num = 2; knots = knots0;

    while(Iter == 0 & (Num <= knot)){

      Bt = bs( Y/max(Y), degree = degree, knots = knots, intercept = TRUE);
      nc = dim(Bt)[2];

      oldlambda.hat = Delta*pnorm(Y, mean= mean(Y))/10;  oldbeta.hat = rep(0, nc);
      iter = 1; error=1;

      while( iter <= maxiter & error > ep ){
        # E-Step  ##
        temp0 = (Bt%*%oldbeta.hat); Etemp0=exp(temp0);
        temp1 = (1-p0G)*exp(Delta*temp0-cumsum(oldlambda.hat*Etemp0));
        temp2 = p0G*exp(-cumsum(oldlambda.hat));
        q = temp2/(temp1 + temp2);

        # M-setp  ##
        temp3 = cumsum( (1-q)[seq(n,1,-1)] ); temp3 = temp3[seq(n,1,-1)];
        temp4 = seq(n,1,-1) - temp3;
        Frac = temp3*Etemp0/(temp3*Etemp0+temp4);
        grad = t( Delta*((1-q) -  Frac ) ) %*% Bt;
        H = t(Bt)%*%(matrix((Frac^2-Frac)*Delta, n, dim(Bt)[2], byrow=FALSE)*Bt);

        if(det(H) > 0.00001 & det(H) < -0.00001 ) beta.hat = oldbeta.hat - solve(H)%*%t(grad) else {
          beta.hat = oldbeta.hat - solve(H - diag(0.0001, nc))%*%t(grad) };

        temp5 = temp3*exp(Bt %*% beta.hat) + temp4;
        lambda.hat = Delta/temp5;

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

    S_lambda_1.hat = lambda.hat;

    # ------------------------------------ #

    Lambda2 = cumsum(S_lambda_1.hat);
    Lambda1 = cumsum(F_lambda_2.hat);

    Lamb1 = apply( matrix(Grid, ncol=1), 1,
                   function(x){ if(sum(Y <= x) >= 1){ max( Lambda1[Y <= x] )} else{0}}
                 );

    Lamb2 = apply( matrix(Grid, ncol=1), 1,
                   function(x){ if(sum(Y <= x) >= 1){ max( Lambda2[Y <= x] )} else{0}}
                 );

    Sieve_F1.hat = 1 - exp(-Lambda1);
    Sieve_F2.hat = 1 - exp(-Lambda2);

    Sieve_lamb = c(F_lambda_2.hat, S_lambda_1.hat);

   if( (sum(is.na(Sieve_lamb)) == 0) & (sum(Sieve_lamb == Inf) < 0.5 ) )
   {
      Boot.L1[Bsim, ] = Lamb1;
      Boot.L2[Bsim, ] = Lamb2;
      Fix_t_F1[Bsim, ] = apply( matrix(fix_t1, ncol=1), 1,
                                function(x){ max( Sieve_F1.hat[Y <= x] ) }  );
      Fix_t_F2[Bsim, ] = apply( matrix(fix_t2, ncol=1), 1,
                                function(x){ max( Sieve_F2.hat[Y <= x] ) }  );
      Bsim = Bsim +1;
   }

 }

 # --------------------------------------------------------------------------------- #

 Result = list( Boot.L1 = Boot.L1,
                Boot.L2 = Boot.L2,
                SE_F1_fix_t = apply(Fix_t_F1, 2, sd),
                SE_F2_fix_t = apply(Fix_t_F2, 2, sd)   );

 return(Result)

}

