dag.legend <-
function(dag, lx=-0.15, ly=-0.075, alt.symb = TRUE)
{ # write legend
  # alt.symb is new in version 1.1.2;

  if(is.null(dag$symbols)) {
    dag$symbols<-rep(NA, length(dag$x));
    writeLines("WARNING! DAG's from version 1.1.2 on should contain a $symbols array. This has been added only temporarily by dag.legend.");
  }

  i1<-0;
  i_c<-0; # covariable counter for subscripts
  i_u<-0; # unknown covs counter for...
  nodes<-length(dag$x);
  while (i1 < nodes)
  {
    i1<-i1+1;
    text(lx, ly-0.05*i1, i1);
    if(dag$names[i1]=="unknown" || dag$cov.types[i1]==2)
    { i_u<-i_u+1;
      if( (alt.symb==FALSE) || (is.na(dag$symbols[i1])) ) {
        text(lx+0.1, ly-0.05*i1, bquote(U[.(i_u)]));
      } else text(lx+0.1, ly-0.05*i1, bquote(.(dag$symbols[i1])));
    } else
    { if(i1==1)
      { if( (alt.symb==FALSE) || (is.na(dag$symbols[1])) ) {
          text(lx+0.1, ly-0.05*i1, "X");
        } else text(lx+0.1, ly-0.05*i1, bquote(.(dag$symbols[1])));
      } else
      { if(i1==nodes)
        { if( (alt.symb==FALSE) || (is.na(dag$symbols[nodes])) ) {
            text(lx+0.1, ly-0.05*i1, "Y");
          } else text(lx+0.1, ly-0.05*i1, bquote(.(dag$symbols[nodes])));
        } else
        { i_c<-i_c+1;
          if( (alt.symb==FALSE) || (is.na(dag$symbols[i1])) ) {
            text(lx+0.1, ly-0.05*i1, bquote(C[.(i_c)]));
          } else text(lx+0.1, ly-0.05*i1, bquote(.(dag$symbols[i1])));
        }
      }
    }
    text(lx+0.15, ly-0.05*i1, dag$names[i1], pos=4);
  }
}

