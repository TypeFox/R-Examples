additional.linear.constraints <-
function# Additional Linear Constraints
### internal function for sisus
(constraints.cols
### internal variable
, constraints.type
### internal variable
, constraints.RHS.b
### internal variable
, constraints.sources
### internal variable
, n.sources
### internal variable
, lc.include.sw
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # create vectors and matricies of additional linear constraints
  b1.constraints = NULL; A1.constraints = NULL; # Ax =b;
  b2.constraints = NULL; A2.constraints = NULL; # Ax<=b;
  b3.constraints = NULL; A3.constraints = NULL; # Ax>=b;

  if (lc.include.sw == 1 ) {
    for (j.con in seq(1,-1+constraints.cols)) { # additional linear constraints

      # if blank column, then the type is NA, which is not available for comparisons, so set to character "NA" to skip
      if (sum((constraints.type[j.con] != ""), na.rm=TRUE) == 0) { constraints.type[j.con] = "NA"; };


      # Equal
      if (constraints.type[j.con] == "Equal") {
        #cat("Equ ",j.con,"\n");
        for (i.con in seq(1,n.sources)) {
          if (!is.na(constraints.sources[i.con,j.con])) {
            b1.constraints = rbind(b1.constraints, constraints.sources[i.con,j.con]);
            A1.constraints.Equal = rep(0,n.sources);
              A1.constraints.Equal[i.con] = 1;
            A1.constraints = rbind(A1.constraints, A1.constraints.Equal);
          }
        }
      } # Equal

      # Minimum
      if (constraints.type[j.con] == "Minimum") {
        #cat("Min ",j.con,"\n");
        for (i.con in seq(1,n.sources)) {
          if (!is.na(constraints.sources[i.con,j.con])) {
            b3.constraints = rbind(b3.constraints, constraints.sources[i.con,j.con]);
            A3.constraints.Minimum = rep(0,n.sources);
              A3.constraints.Minimum[i.con] = 1;
            A3.constraints = rbind(A3.constraints, A3.constraints.Minimum);
          }
        }
      } # Minimum

      # Maximum
      if (constraints.type[j.con] == "Maximum") {
        #cat("Max ",j.con,"\n");
        for (i.con in seq(1,n.sources)) {
          if (!is.na(constraints.sources[i.con,j.con])) {
            b2.constraints = rbind(b2.constraints, constraints.sources[i.con,j.con]);
            A2.constraints.Maximum = rep(0,n.sources);
              A2.constraints.Maximum[i.con] = 1;
            A2.constraints = rbind(A2.constraints, A2.constraints.Maximum);
          }
        }
      } # Maximum

      # EQ
      if (constraints.type[j.con] == "EQ") {
        #cat("EQ  ",j.con,"\n");
        b1.constraints = rbind(b1.constraints, constraints.RHS.b[j.con]);
        A1.constraints.EQ = rep(0,n.sources);
        for (i.con in seq(1,n.sources)) {
          if (!is.na(constraints.sources[i.con,j.con])) {
            A1.constraints.EQ[i.con] = constraints.sources[i.con,j.con];
          }
        }
        A1.constraints = rbind(A1.constraints, A1.constraints.EQ);
      } # EQ

      # GE
      if (constraints.type[j.con] == "GE") {
        #cat("GE  ",j.con,"\n");
        b3.constraints = rbind(b3.constraints, constraints.RHS.b[j.con]);
        A3.constraints.GE = rep(0,n.sources);
        for (i.con in seq(1,n.sources)) {
          if (!is.na(constraints.sources[i.con,j.con])) {
            A3.constraints.GE[i.con] = constraints.sources[i.con,j.con];
          }
        }
        A3.constraints = rbind(A3.constraints, A3.constraints.GE);
      } # GE

      # LE
      if (constraints.type[j.con] == "LE") {
        #cat("LE  ",j.con,"\n");
        b2.constraints = rbind(b2.constraints, constraints.RHS.b[j.con]);
        A2.constraints.LE = rep(0,n.sources);
        for (i.con in seq(1,n.sources)) {
          if (!is.na(constraints.sources[i.con,j.con])) {
            A2.constraints.LE[i.con] = constraints.sources[i.con,j.con];
          }
        }
        A2.constraints = rbind(A2.constraints, A2.constraints.LE);
      } # LE

    } # for additional linear constraints
  } # if

  # create a list to return with all constraints
  lc = new.env();
  lc$b1.constraints = b1.constraints;
  lc$A1.constraints = A1.constraints;
  lc$b2.constraints = b2.constraints;
  lc$A2.constraints = A2.constraints;
  lc$b3.constraints = b3.constraints;
  lc$A3.constraints = A3.constraints;
  return( as.list(lc) );
  ### internal variable
}
