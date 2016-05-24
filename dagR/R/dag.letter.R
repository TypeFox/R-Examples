dag.letter <-
function(dag, letter, x, y, alt.symb = TRUE)
{ # function to draw the letters in the DAG;
  # alt.symb is new since v1.1.2;

 if(is.null(dag$symbols)) dag$symbols<-rep(NA, length(dag$x));

 if(letter==1)
 { if( (alt.symb==FALSE) || (is.na(dag$symbols[1])) ) {
     text(x, y, "X");
   } else text(x, y, bquote(.(dag$symbols[1])));
 } else
 if(letter==length(dag$x))
 { if( (alt.symb==FALSE) || (is.na(dag$symbols[letter])) ) {
     text(x, y, "Y");
   } else text(x, y, bquote(.(dag$symbols[letter])));
 } else
 {
  i_c<-0; # covariable counter for subscripts
  i_u<-0; # unknown covs counter for...
  for(i1 in 2:letter)
  {
    if(dag$names[i1]=="unknown" || dag$cov.types[i1]==2)
    { i_u<-i_u+1;
    } else
    { i_c<-i_c+1;
    }
  }
  
  if(dag$names[letter]=="unknown" || dag$cov.types[letter]==2) 
  { if( (alt.symb==FALSE) || (is.na(dag$symbols[letter])) )
    { text(x, y, bquote(U[.(i_u)]));
    } else
    { text(x, y, bquote(.(dag$symbols[letter])));
    }
  } else
  { if(is.in(as.numeric(letter), dag$adj)==TRUE)
    # the above use of as.numeric() was required due to some strange behaviour:
    #  without it, in some cases, this block would be skipped despite letter
    #  being in dag$adj; for instance with demo.dag2, if C1,3,5,7, or 9 was
    #  adjusted, they would not appear so in the write.paths() output, whereas
    #  both the manual use of dag.letter() and the coding for the nodes in dag.draw
    #  would produce the correct barred/underlined output. however, other DAGs,
    #  e.g. demo.dag1, were not affected!
    { if( (alt.symb==FALSE) || (is.na(dag$symbols[letter])) ) 
      { text(x, y, bquote(underline(bar(C))[.(i_c)]));
      } else  
      { text(x, y, bquote(underline(bar(.(dag$symbols[letter])))));
      }
    } else 
    { if( (alt.symb==FALSE) || (is.na(dag$symbols[letter])) ) 
      { text(x, y, bquote(C[.(i_c)]));
      } else
      { text(x, y, bquote(.(dag$symbols[letter])));
      }
    }
  }
 }
}

