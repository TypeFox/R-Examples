if(getRversion() >= "2.15.1") utils::globalVariables(c("left","right","coefs"));

.adaptQuadratic<-function(coef, central){
  vertex = -coef[2]/(2*coef[3]);
  verticaloffset = coef[1] + coef[2]*vertex + coef[3]*vertex*vertex;
  temp = c(coef[1]-verticaloffset, coef[2], coef[3]);
  alpha = 1/(temp[1]+temp[2]*central+temp[3]*central*central);
  result = c(alpha*temp[1], alpha*temp[2], alpha*temp[3]);
  return(result);
}

## _____________________________________________________________________________________________________________

.gaussianQuadraticError<-function(params, miframe, lower, upper, side){
  
  
  if(side == "left"){ c = upper; } # impose y(upper) = 1
  else{ c = lower; } # impose y(lower) = 1
  temp = c(0,0,0);
  temp[1] = params[1];  temp[2] = params[2];  temp[3] = c;
  params = temp;
  
  penalty = 0;
  #if(exp((-1/2)*abs((0-params[3])/(params[1]+0.001))^params[2]) > 0.05) penalty = penalty + 1E10;
  #if(exp((-1/2)*abs((1-params[3])/(params[1]+0.001))^params[2]) > 0.05) penalty = penalty + 1E10;
  predicted = sapply(X=miframe[,1], FUN=function(x) exp((-1/2)*abs((x-params[3])/(params[1]+0.001))^params[2]) );
  
  predicted[predicted<0] = 0;
  observed = miframe[,2];
  squaredif = (predicted - observed)^2;
  return(sum(squaredif) + penalty);
}

## _____________________________________________________________________________________________________________

.quadraticQuadraticError<-function(params,miframe,lower,upper,side){ #upper and lower are the extremes of the X-axis interval where we are doing regression
  ## y = params[1] + params[2]*x + params[3]*x^2  
  
  # impose constraints
  a = params[1]; b = params[2]; c = -1;
  if(side=="left"){ # impose y(upper) = 1
    c = (1 - a - b*upper)/(upper*upper);
  }else{  # impose y(lower) = 1
    c = (1 - a - b*lower)/(lower*lower);
  }
  params = c(0,0,0);
  params[1] = a; params[2] = b; params[3] = c;
  
  
  penalty = 0;
  
  discriminante = params[2]*params[2]-4*params[3]*params[1];
  if(discriminante < 0){ penalty = penalty + 2E10; }
  else{
    sol1 = (-params[2]-sqrt(discriminante))/(2*params[3]);
    sol2 = (-params[2]+sqrt(discriminante))/(2*params[3]);
    if(sol1 > sol2){ temp = sol1; sol1 = sol2; sol2 = temp; } # now sol1 < sol2
    if(side == "left"){
      if(abs(sol2) < 10^(-5)) sol2 = 0;
      lower = sol2;
      if(sol2 < 0 || sol2 > upper){ penalty = penalty+1E10; }
    }
    else{ 
      if(abs(sol1) < 10^(-5)) sol1 = 0;
      upper = sol1;
      if(sol1 < lower || sol1 > 1) { penalty = penalty+1E10; } 
    }
  }
  
  # This is to guarantee that it is monotonous in the independent value regression interval
  if(lower <= -params[2]/(2*params[3]) && -params[2]/(2*params[3]) <= upper){ penalty = penalty+1E10; } 
  
  predicted = sapply(X=miframe[,1], FUN=function(x) params[1] + params[2]*x + params[3]*x^2 );
  predicted[predicted<0] = 0;
  observed = miframe[,2];
  squaredif = (predicted - observed)^2;
  return(sum(squaredif)+penalty); 
}

## _____________________________________________________________________________________________________________

.computeIntersectionX<-function(a,b,c,d,lower,upper,side,quadratic=FALSE,mostrar=FALSE){ # Finds x such that a + b*x + c*x^2 + d*x^3 == 0
  if(quadratic){ z = c(a,b,c); }
  else{  z = c(a,b,c,d); }
  roots = polyroot(z);
  res = -1;
  roots = roots[abs(Im(roots))<10^(-5) & lower <= Re(roots) & Re(roots) <= upper]; # roots that fall inside this interval
  if(length(roots) == 0){  return(); }
  else{
    roots = Re(roots);
    roots = sort(roots,decreasing=TRUE);    
    if(side=="left"){ return(roots[1]); } # left side -> return the greatest that falls inside the interval
    else{ return(roots[length(roots)]); }  # right side -> return the smallest that falls inside the interval
  }
}

## _____________________________________________________________________________________________________________

.hasCubicTurningPoints<-function(coef,lower,upper){
  derivative = c(0,0,0);
  derivative[1] = coef[2];
  derivative[2] = 2*coef[3];
  derivative[3] = 3*coef[4];
  roots = polyroot(derivative);
  roots = roots[abs(Im(roots))<10^(-5) & lower <= Re(roots) & Re(roots) <= upper]; # turning points that fall inside this interval
  
  ## Check whether zero-deritative points are really extremes or just saddle points
  if(length(roots)>0){ 
    roots = Re(roots);
    for(i in 1:length(roots)){
      xleft = roots[i]-0.0000001;
      xright = roots[i] + 0.000001;
      if((coef[2]+2*coef[3]*xleft+3*coef[4]*xleft^2)*(coef[2]+2*coef[3]*xright+3*coef[4]*xright^2) < 0){ # -1 --> different sign in derivate
         return(TRUE);
      }
    }    
  }
  return(FALSE); 
}

## _____________________________________________________________________________________________________________

.cubicQuadraticError<-function(params,miframe,lower,upper,side){
  ## y = params[1] + params[2]*x + params[3]*x^2 + params[4]*x^3
  # This is to guarantee that it is monotonous in the independent value regression interval
  
  # impose constraints
  a = params[1]; b = params[2]; c = params[3]; d = -1;
  if(side == "left"){    # impose y(upper) = 1
    d = (1 - a - b*upper-c*upper*upper)/(upper*upper*upper);
  } 
  else{ # impose y(lower) = 1
    d = (1 - a - b*lower-c*lower*lower)/(lower*lower*lower);
  }
  params = c(0,0,0,0);
  params[1] = a; params[2] = b; params[3] = c; params[4] = d;
  
  
  if(abs(params[4])<0.001) return(Inf);
  if(side == "left"){ 
    zeroat = .computeIntersectionX(params[1],params[2],params[3],params[4], 0, upper,"left"); 
    if(is.null(zeroat)) return(Inf);
    lower = zeroat;
  }
  else{
    zeroat = .computeIntersectionX(params[1],params[2],params[3],params[4], lower, 1,"right"); 
    if(is.null(zeroat)) return(Inf);
    upper = zeroat;
  }  
  penalty = 0;
  if(.hasCubicTurningPoints(params,lower,upper)){ penalty = 1E10; }

  predicted = sapply(X=miframe[,1], FUN=function(x) params[1] + params[2]*x + params[3]*x^2 + params[4]*x^3);
  #predicted[predicted<0] = 0;
  observed = miframe[,2];
  squaredif = (predicted - observed)^2;
  return(sum(squaredif)+penalty); 
}

## _____________________________________________________________________________________________________________

# This is a bit tricky. Simulate that the coefficients are local static variables inside the membership functions, 
# as the membership functions actually get only one argument x   
# The result is an already corrected left (or right) side membership function f:[0,1]->[0,1]
.makeFn<-function(params,regressiontype,coreLeftThreshold,coreRightThreshold,lower,upper,side){

    e<-new.env();    
    e$coefs = params;
    e$coreLeft = coreLeftThreshold;
    e$coreRight = coreRightThreshold;
    
    e$left=lower;
    e$right=upper;
    e$side=side    
        
    if(regressiontype == "quadratic"){  # y = a + b*x + c*x^2
      f<-function(x){ # At every call, argument x is always in [0, 1]
        xreal = left + (right - left)*x;
        result = coefs[1] + coefs[2]*xreal + coefs[3]*xreal^2;
        if(side == "left"){ result[x<0.001] = 0; result[x>0.999] = 1; }
        else{               result[x<0.001] = 1; result[x>0.999] = 0; }
        result[result<0] = 0;
        return(result);
      }
      environment(f) = e;    
      return(f);
    }else{
      if(regressiontype == "cubic"){ # y = a + b*x + c*x^2 + d*x^3
        f<-function(x){
          xreal = left + (right - left)*x;
          result = coefs[1] + coefs[2]*xreal + coefs[3]*xreal^2 + coefs[4]*xreal^3;
          if(side == "left"){ result[x<0.001] = 0; result[x>0.999] = 1; }
          else{               result[x<0.001] = 1; result[x>0.999] = 0; }
          result[result<0] = 0;
          return(result);
        }
        environment(f) = e;
        return(f);
      }else{
        if(regressiontype == "linear"){ 
          f<-function(x){ # At every call, argument x is always in [0, 1]
            xreal = left + (right - left)*x;
            result = coefs[1] + coefs[2]*xreal;
            if(side == "left"){ result[x<0.001] = 0; result[x>0.999] = 1; }
            else{               result[x<0.001] = 1; result[x>0.999] = 0; }
            result[result<0] = 0;
            return(result);
          }
          environment(f) = e;    
          return(f);
        }
        else{ 
          if(regressiontype == "inverse"){  # y = a + b*(1/x)
            f<-function(x){
              xreal = left + (right - left)*x;
              result = coefs[1] + coefs[2]*(1/xreal);
              result[result<0] = 0;
              return(result);
            }
            environment(f) = e;
            return(f);
          }
          else{
            if(regressiontype == "gaussian"){
              f<-function(x){
                xreal = left + (right - left)*x;
                result = exp((-1/2)*abs((xreal-coefs[3])/(coefs[1]+0.001))^coefs[2]);
                if(side == "left"){ result[x<0.001] = 0; result[x>0.999] = 1; }
                else{               result[x<0.001] = 1; result[x>0.999] = 0; }
                result[result<0] = 0;
                return(result);
              }
              environment(f) = e;
              return(f);
            }
            else{
              if(regressiontype == "spline"){
                f<-function(x){
                  myfunction = coefs;
                  xreal = left + (right - left)*x;
                  result = myfunction(xreal,deriv=0);
                  if(side == "left"){ result[x<0.001] = 0; result[x>0.999] = 1; }
                  else{               result[x<0.001] = 1; result[x>0.999] = 0; }
                  result[result<0] = 0;
                  return(result);                  
                }
                environment(f) = e;
                return(f);
              }
              else{
                stop("ERROR: not found regression type:",regressiontype);
              }
            }
          }
        }
      }
    }
}

## _____________________________________________________________________________________________________________

.membFnRegression<-function(miframeizq, miframeder, regressiontype, linguistic, verbose=FALSE){

    # 4 points a1 <= a2 <= a3 <= a4 which define the Fuzzy Number
    a1 = miframeizq[1,1];
    a2 = miframeizq[length(miframeizq[,1]),1];
    a3 = a2;
    a4 = miframeder[1,1];
    
    if(linguistic){
      a3 = miframeder[length(miframeder[,1]),1];
    }
        
    names(miframeizq) = c("w","y"); # w = independent variable for regression to obtain left and right side membership functions
    names(miframeder) = c("w","y"); # y = dependent variable for regression
    
    if(!linguistic){    
      condizq = abs(miframeizq$w - a2) > 0.005;
      condder = abs(miframeder$w - a2) > 0.005;
      
      miframeizq = miframeizq[condizq | condder,];
      miframeder = miframeder[condizq | condder,];
    }
    
    coreLeftThreshold = miframeizq[length(miframeizq[,1]),1];
    coreRightThreshold = miframeder[length(miframeder[,1]),1];
    
    if(!linguistic){
      miframeizq = rbind(miframeizq,c(a2,0.999));
      miframeder = rbind(miframeder,c(a2,0.999));
    }
    
    fuzzy = NA;

    ########################################################################
    ##                  QUADRATIC REGRESSION                              ##
    ########################################################################
    if(regressiontype == "quadratic"){

        ## LEFT SIDE ---------------------------------------------------

          opt1 = nls(y~I(a+b*w+(1-a-b*a2)/(a2*a2)*w^2),miframeizq,start=list(a=100,b=1));          
          coefizq = summary(opt1)$coefficients[,1];
          c = (1-coefizq[1]-coefizq[2]*a2)/(a2*a2); # the curve has been fitted assuming the imposed constraints
          coefizq = append(coefizq,c);
          
          tentative = .computeIntersectionX(coefizq[1],coefizq[2],coefizq[3],-1,quadratic=TRUE,0,a2,"left",mostrar=FALSE);      
          if(!is.null(tentative)){ 
            a1 = tentative; 
            tentativebestval = .quadraticQuadraticError(coefizq,miframeizq,0,a2,"left");
          }
          else{
            coefizq = .adaptQuadratic(coefizq, a2);
            tentativebestval = .quadraticQuadraticError(coefizq,miframeizq,0,a2,"left");
            a1 = -coefizq[2]/(2*coefizq[3]);
          }

          # Use differential evolution (approximate) to compute an alternative solution
          lowerInitial = c(coefizq[1]-100,coefizq[2]-100); upperInitial = c(coefizq[1]+100,coefizq[2]+100); 
          initialpopleft = .generateInitialPopulationQuadraticRegression(NP=50, a2, "left")
     
          opt = DEoptim(.quadraticQuadraticError, lowerInitial, upperInitial, control = DEoptim.control(trace=FALSE,itermax=200,CR=0.8,NP=50,
                      initialpop = initialpopleft), miframeizq, 0, a2, "left");
          coefizq2 = opt$optim$bestmem;
          c = (1-coefizq2[1]-coefizq2[2]*a2)/(a2*a2); # the curve has been fitted assuming the imposed constraints
          coefizq2 = append(coefizq2,c);
          izqbestval=opt$optim$bestval;
          
          # Check which of the solutions is better
          if(tentativebestval < izqbestval){   }
          else{ 
            if(izqbestval >= 1E10){ print("Left side quadratic membership function has not been found"); return(); }
            else{ 
              coefizq = coefizq2; 
              discriminante = coefizq[2]*coefizq[2]-4*coefizq[3]*coefizq[1];
              sol1 = (-coefizq[2]-sqrt(discriminante))/(2*coefizq[3]);
              sol2 = (-coefizq[2]+sqrt(discriminante))/(2*coefizq[3]);
              a1 = max(sol1,sol2);              
              #a1 = .computeIntersectionX(coefizq[1],coefizq[2],coefizq[3],-1,quadratic=TRUE,0,a2,"left",mostrar=FALSE);
              if(abs(a1) < 10^(-5)) a1 = 0;
            }
          }         

          
          ## RIGHT SIDE ---------------------------------------------------
          
          opt2 = nls(y~I(a+b*w+(1-a-b*a3)/(a3*a3)*w^2),miframeder,start=list(a=100,b=1));
          coefder = summary(opt2)$coefficients[,1];
          c = (1-coefder[1]-coefder[2]*a3)/(a3*a3); # the curve has been fitted assuming the imposed constraints
          coefder = append(coefder,c);
          tentative = .computeIntersectionX(coefder[1],coefder[2],coefder[3],-1,quadratic=TRUE,a3,1,"right",mostrar=FALSE);
          
          if(!is.null(tentative)){ 
              a4 = tentative; 
              tentativebestval = .quadraticQuadraticError(coefder,miframeder,a3,1,"right");
          }
          else{
            coefder = .adaptQuadratic(coefder, a2);
            tentativebestval = .quadraticQuadraticError(coefder,miframeder,a3,1,"right");
            a4 = -coefder[2]/(2*coefder[3]);
          }
          
          # Use differential evolution (approximate) to compute an alternative solution
          lowerInitial = c(coefder[1]-300,coefder[2]-300); upperInitial = c(coefder[1]+300,coefder[2]+300);     
          initialpopright = .generateInitialPopulationQuadraticRegression(NP=50,a3, "right");
          
          opt = DEoptim(.quadraticQuadraticError, lowerInitial, upperInitial, control = DEoptim.control(trace=FALSE,itermax=200,CR=0.8,NP=50,
                      initialpop = initialpopright), miframeder, a3, 1, "right");
          coefder2 = opt$optim$bestmem;
          c = (1-coefder2[1]-coefder2[2]*a3)/(a3*a3); # the curve has been fitted assuming the imposed constraints
          coefder2 = append(coefder2,c);
          
          derbestval=opt$optim$bestval;
           
          # Check which of the solutions is better
          if(tentativebestval < derbestval){   }
          else{ 
            if(derbestval >= 1E10){ print("Right side quadratic membership function has not been found"); return(); }
            else{ 
              coefder = coefder2;
              a4 = .computeIntersectionX(coefder[1],coefder[2],coefder[3],-1,quadratic=TRUE,a3,1,"right",mostrar=FALSE); 
            }
          } 

    }else{
      ########################################################################
      ##                      CUBIC REGRESSION                              ##
      ########################################################################
      if(regressiontype == "cubic"){
         
         ## LEFT SIDE ---------------------------------------------------

          opt1 = nls(y~I(a+b*w+c*w^2+((1-a-b*a2-c*a2*a2)/(a2*a2*a2))*w^3),miframeizq,start=list(a=100,b=1,c=1));
          coefizq = summary(opt1)$coefficients[,1];
          d = (1-coefizq[1]-coefizq[2]*a2-coefizq[3]*a2*a2)/(a2*a2*a2); # the curve has been fitted assuming the imposed constraints
          coefizq = append(coefizq,d);
          
          tentative = .computeIntersectionX(coefizq[1],coefizq[2],coefizq[3],coefizq[4],quadratic=FALSE,0,a2,"left",mostrar=FALSE);
          hasTurningPoints = FALSE;
          if(!is.null(tentative)){ 
            hasTurningPoints = .hasCubicTurningPoints(coefizq,tentative,a2);            
          }
          check = FALSE;
          if(!is.null(tentative) && !hasTurningPoints) {
            check = TRUE;
            a1 = tentative; 
            tentativebestval = .cubicQuadraticError(coefizq,miframeizq,0,a2,"left");
          }
          # use differential evolution (approximate) to check whether we can find a better solution

           lowerInitial = c(coefizq[1]-500,coefizq[2]-500,coefizq[3]-500); upperInitial = c(coefizq[1]+500,coefizq[2]+500,coefizq[3]+500);
           opt2 = DEoptim(.cubicQuadraticError, lowerInitial, upperInitial, control = DEoptim.control(trace=FALSE,itermax=200,CR=0.8), 
                          miframeizq, 0, a2,"left");
           coefizq2 = opt2$optim$bestmem;
           d = (1-coefizq2[1]-coefizq2[2]*a2-coefizq2[3]*a2*a2)/(a2*a2*a2); # the curve has been fitted assuming the imposed constraints
           coefizq2 = append(coefizq2,d);

           izqbestval=opt2$optim$bestval;
           if(check && tentativebestval < izqbestval){ }
           else{
             if(izqbestval >= 1E10){ print("Left side cubic membership function has not been found"); return(); }
             else{
               coefizq = coefizq2;
               a1 = .computeIntersectionX(coefizq[1],coefizq[2],coefizq[3],coefizq[4],0,a2,"left",mostrar=TRUE);
             }
           }
       
         ## RIGHT SIDE ----------------------------------------------------
          
          opt1 = nls(y~I(a+b*w+c*w^2+((1-a-b*a2-c*a2*a2)/(a2*a2*a2))*w^3),miframeder,start=list(a=100,b=1,c=1));
          coefder = summary(opt1)$coefficients[,1];
          d = (1-coefder[1]-coefder[2]*a2-coefder[3]*a2*a2)/(a2*a2*a2); # the curve has been fitted assuming the imposed constraints
          coefder = append(coefder,d);
          
          tentative = .computeIntersectionX(coefder[1],coefder[2],coefder[3],coefder[4],quadratic=FALSE,a3,1,"right",mostrar=FALSE);
          hasTurningPoints = FALSE;
          if(!is.null(tentative)){ 
            hasTurningPoints = .hasCubicTurningPoints(coefder,a3,tentative);
          }
          check = FALSE;
          if(!is.null(tentative) && !hasTurningPoints){ 
              a4 = tentative; 
              check = TRUE;
              tentativebestval = .cubicQuadraticError(coefder,miframeder,a2,1,"right");
          }
          
         lowerInitial = c(coefder[1]-500,coefder[2]-500,coefder[3]-500); upperInitial = c(coefder[1]+500,coefder[2]+500,coefder[3]+500);
         opt = DEoptim(.cubicQuadraticError, lowerInitial, upperInitial, control = DEoptim.control(trace=FALSE,itermax=200,CR=0.8), 
                        miframeder, a3, 1,"right");
         coefder2 = opt$optim$bestmem;  
         derbestval = opt$optim$bestval;
         d = (1-coefder2[1]-coefder2[2]*a3-coefder2[3]*a3*a3)/(a3*a3*a3); # the curve has been fitted assuming the imposed constraints
         coefder2 = append(coefder2,d);
         
         if(check && tentativebestval < derbestval){   }
         else{
           if(derbestval >= 1E10){print(paste("Right side cubic membership function has not been found:",derbestval)); return(); }
           else{
             coefder = coefder2;
             a4 = .computeIntersectionX(coefder[1],coefder[2],coefder[3],coefder[4],a3,1,"right");
           }
         }
          
      }else{
        ########################################################################
        ##                      LINEAR REGRESSION                             ##
        ########################################################################      
        if(regressiontype == "linear"){
           # constrained regression
            regizq = nls(y~I(1-b*a2+b*w),miframeizq,start=list(b=1));
            regder = nls(y~I(1-b*a3+b*w),miframeder,start=list(b=1));
            coefizq = summary(regizq)$coefficients;
            coefder = summary(regder)$coefficients;
            coefizq = append(coefizq, 1-coefizq[1]*a2, after=0);
            coefder = append(coefder, 1-coefder[1]*a3, after=0);
        }
        else{
          ########################################################################
          ##                   GAUSSIAN REGRESSION                              ##
          ########################################################################        
          if(regressiontype == "gaussian"){             
             lowerInitial = c(0.01,1); upperInitial = c(10,3);
             
             opt = DEoptim(.gaussianQuadraticError, lowerInitial, upperInitial, control = DEoptim.control(trace=FALSE,itermax=200,CR=0.8), 
                            miframeizq, 0, a2, "left");
                            
             coefizq = opt$optim$bestmem;
             # the center parameter that imposes y(a2) = 1
             coefizq = append(coefizq, a2);

             opt = DEoptim(.gaussianQuadraticError, lowerInitial, upperInitial, control = DEoptim.control(trace=FALSE,itermax=200,CR=0.8), 
                            miframeder, a3, 1, "right");
             coefder = opt$optim$bestmem;
             coefder = append(coefder, a3); # the center parameter that imposes y(a3) = 1             
             
          }
          else{
            ########################################################################
            ##                   CUBIC HYMAN SPLINE INTERPOLATION                 ##
            ########################################################################           
            if(regressiontype == "spline"){
              coefizq = splinefun(miframeizq[,1],miframeizq[,2],method="hyman");
              coefder = splinefun(miframeder[,1],miframeder[,2],method="hyman");
            }
            else{
              ########################################################################
              ##                   PIECEWISE LINEAR INTERPOLATION                   ##
              ########################################################################             
              if(regressiontype == "piecewise"){
                leftknots = miframeizq[2:length(miframeizq[,1])-1,]; # length-1 to avoid the central point (vertex) which is not a knot
                alphas = leftknots[,2];
                leftknots = leftknots[,1];
                leftknots = sort(leftknots);
                rightknots = miframeder[2:length(miframeder[,1])-1,1];
                rightknots = sort(rightknots);
    
                fuzzy = PiecewiseLinearFuzzyNumber(a1,a2,a3,a4,knot.n = length(leftknots),
                           knot.alpha = alphas, knot.left = leftknots, knot.right = rightknots);
              }
              else{
                stop('ERROR: unknown regression type:',regressiontype);
              }
            }
          }
        }
      }
    }

    #newcoefizq = .coefficientCorrection(a1,a2, coefizq, regressiontype);
    #newcoefder = .coefficientCorrection(a3,a4, coefder, regressiontype);
    
    if(regressiontype != "piecewise"){
      fizq = .makeFn(coefizq,regressiontype,coreLeftThreshold,coreRightThreshold,a1,a2,side="left");
      fder = .makeFn(coefder,regressiontype,coreLeftThreshold,coreRightThreshold,a3,a4,side="right");
              
      fuzzy = FuzzyNumber(a1,a2,a3,a4, left = fizq, right = fder);
    }

    return(list(fuzzy,miframeizq,miframeder));
}