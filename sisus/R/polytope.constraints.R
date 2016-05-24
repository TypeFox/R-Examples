polytope.constraints <-
function# Polytope constraints (incorporates both Isotope and Additional linear constraints)
### internal function for sisus
(lc
### internal variable
, n.sources
### internal variable
, n.isotopes
### internal variable
, biomass.values
### internal variable
, simplex.include.sw
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # unpack the list of additional linear constraints
  b1.constraints = lc$b1.constraints;
  A1.constraints = lc$A1.constraints;
  b2.constraints = lc$b2.constraints;
  A2.constraints = lc$A2.constraints;
  b3.constraints = lc$b3.constraints;
  A3.constraints = lc$A3.constraints;
       # if (debug.level > 0) { # DEBUG display constraints
       #   print(b1.constraints);
       #   print(A1.constraints);
       #   print(b2.constraints);
       #   print(A2.constraints);
       #   print(b3.constraints);
       #   print(A3.constraints);
       # }

  ## Define linear system of equalities and inequalities
  # These equalities define the linear system.  The first constraint that the proportions sum to one (simplex).
  # 1:  Ax=b;

  # whether or not to include simplex constraint
  if (simplex.include.sw == 1) { b1.simplex = 1;    A1.simplex = rep(1,   n.sources); } # include simplex constrait (usually)
  if (simplex.include.sw == 0) { b1.simplex = NULL; A1.simplex = rep(NULL,n.sources); } # do not include simplex

  # b1 is the rhs for the equality
  b1 = rbind(b1.constraints, b1.simplex); for (i.mixture in 1:n.isotopes) { b1 = rbind(b1, 0); }; b1 = as.vector(b1);
  A1 = rbind(A1.constraints, A1.simplex); for (i.mixture in 1:n.isotopes) { A1 = rbind(A1,as.numeric(biomass.values[,i.mixture])) ; }

  # these inequalities say the proportions are in the unit hypercube
  # 2:  Ax<=b;
  b2 = rbind(b2.constraints);
    for (i.sources in 1:n.sources) { b2 = rbind(b2, 1); }; b2 = as.vector(b2);
  A2 = rbind(A2.constraints, diag(1,n.sources)); # linear constraints and p's no more than 1

  # 3:  Ax>=b;
  b3 = rbind(b3.constraints);
    for (i.sources in 1:n.sources) { b3 = rbind(b3, 0); }; b3 = as.vector(b3);
  A3 = rbind(A3.constraints, diag(1,n.sources)); # linear constraints and p's no less than 0

  # check that b vector is positive, flipping signs in both b and A where necessary
  fix.neg = 1-as.numeric(b1<0)*2; b1 = b1*fix.neg; A1 = A1*fix.neg;
  fix.neg = 1-as.numeric(b2<0)*2; b2 = b2*fix.neg; A2 = A2*fix.neg;
  fix.neg = 1-as.numeric(b3<0)*2; b3 = b3*fix.neg; A3 = A3*fix.neg;

       # if (debug.level > 0) { # DEBUG display constraints
       #   print(b1);
       #   print(A1);
       #   print(b2);
       #   print(A2);
       #   print(b3);
       #   print(A3);
       # }

  # create a list to return with all constraints
  Ab = new.env();
  Ab$b1 = b1;
  Ab$A1 = A1;
  Ab$b2 = b2;
  Ab$A2 = A2;
  Ab$b3 = b3;
  Ab$A3 = A3;
  return( as.list(Ab) );
  ### internal variable
}
