write.Ab <-
function# write A matrix and b vector
### writes out the linear constraints to the process_info.txt file
(Ab
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # matrix sizes
  n.row.A1 = dim(Ab$A1)[1];  n.col.A1 = dim(Ab$A1)[2];
  n.row.A2 = dim(Ab$A2)[1];  n.col.A2 = dim(Ab$A2)[2];
  n.row.A3 = dim(Ab$A3)[1];  n.col.A3 = dim(Ab$A3)[2];

  p.o = "\n"; write.out(p.o);
  p.o = "\nLinear System values, solving for p\n"; write.out(p.o);
  # 1:  Ax=b;
  p.o = "\n  These linear equalities are the primary linear system\n"; write.out(p.o);
  p.o = "\n    Ap = b\n\n"; write.out(p.o);
  for (j.A1 in seq(1,n.row.A1)) {
    row.text = "    ";
    for (i.A1 in seq(1,n.col.A1)) {
      #row.text <- paste(row.text, Ab$A1[j.A1, i.A1], " p", i.A1);
      row.text <- paste(row.text, sprintf("%+8.7g",Ab$A1[j.A1,i.A1]), " p", i.A1, sep="");
      if (i.A1 < n.col.A1) {
        row.text <- paste(row.text, " + ", sep="");
      } else {
        row.text <- paste(row.text, " = ", sep="");
      }
    }
    row.text <- paste(row.text, sprintf("%+8.7g",Ab$b1[j.A1]), "\n", sep="");
    p.o = row.text; write.out(p.o);
  }

  p.o = "\n  These inequalities say the proportions are between 0 and 1,\n    or come from additional specified linear constraints\n"; write.out(p.o);
  # 3:  Ax>=b;
  p.o = "\n    Ap >= b\n\n"; write.out(p.o);
  for (j.A3 in seq(1,n.row.A3)) {
    row.text = "    ";
    for (i.A3 in seq(1,n.col.A3)) {
      #row.text <- paste(row.text, Ab$A3[j.A3, i.A3], " p", i.A3);
      row.text <- paste(row.text, sprintf("%+8.7g",Ab$A3[j.A3,i.A3]), " p", i.A3, sep="");
      if (i.A3 < n.col.A3) {
        row.text <- paste(row.text, " + ", sep="");
      } else {
        row.text <- paste(row.text, " >= ", sep="");
      }
    }
    row.text <- paste(row.text, sprintf("%+8.7g",Ab$b3[j.A3]), "\n", sep="");
    p.o = row.text; write.out(p.o);
  }
  # 2:  Ax<=b;
  p.o = "\n    Ap <= b\n\n"; write.out(p.o);
  for (j.A2 in seq(1,n.row.A2)) {
    row.text = "    ";
    for (i.A2 in seq(1,n.col.A2)) {
      #row.text <- paste(row.text, Ab$A2[j.A2, i.A2], " p", i.A2);
      row.text <- paste(row.text, sprintf("%+8.7g",Ab$A2[j.A2,i.A2]), " p", i.A2, sep="");
      if (i.A2 < n.col.A2) {
        row.text <- paste(row.text, " + ", sep="");
      } else {
        row.text <- paste(row.text, " <= ", sep="");
      }
    }
    row.text <- paste(row.text, sprintf("%+8.7g",Ab$b2[j.A2]), "\n", sep="");
    p.o = row.text; write.out(p.o);
  }

  ### internal variable
}
