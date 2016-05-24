dag.draw <-
function(dag, legend=TRUE, paths=TRUE, numbering=FALSE, p=FALSE, alt.symb=TRUE, noxy=0,...)
{ # the is.numeric() condition is "new" and avoid error when paths have not been evaluated;
  # the y.low definition for legend==F && paths==T is "new" to ensure display of paths;

  # the noxy option is new to version 1.1.2;  

  # the alternative symbols options is new to version 1.1.2;
  #  it requires $symbols of the dag not to be NULL, so this is checked and
  #  a warning may be issued;
  # it would have been nice to write the custom symbols with something like
  #  boxed.labels or shadowtext, but these do not allow bar/underlining etc.

  if(is.null(dag$symbols)) {
    dag$symbols<-rep(NA, length(dag$x));
    dag$warning<-"WARNING! DAG's from version 1.1.2 on should contain a $symbols array. This has been added by dag.draw";
  }

  y.low<-(-0.25);
  if( legend==TRUE && (paths==TRUE && is.numeric(dag$pathsN)) )
  { y.low<-y.low-max(c(0.05*(length(dag$x)-4), 0.15+0.065*(dag$pathsN-5)));
  } else  
  { if(legend==TRUE)
    { y.low<-y.low-0.05*(length(dag$x)-4);
    } else
    if( paths==TRUE && is.numeric(dag$pathsN) )
    { y.low<-y.low-(0.15+0.065*(dag$pathsN-5));
    }
  }

  y.high<-max(c(1, max(dag$y)+0.25));
  x.left<-min(c(-0.25, min(dag$x)-0.25));
  x.right<-max(c(1.25, max(dag$x)+0.25));
  plot(x=c(x.left, x.right),
       y=c(y.low, y.high), type="n", axes=FALSE, xlab='', ylab='');

  if( (alt.symb==FALSE) || (is.na(dag$symbols[1])) ) {
    text(dag$x[1], dag$y[1], "X");
  } else text(dag$x[1], dag$y[1], bquote(.(dag$symbols[1])));

  if( (alt.symb==FALSE) || (is.na(dag$symbols[length(dag$x)])) ) {
    text(dag$x[length(dag$x)], dag$y[length(dag$y)], "Y");
  } else text(dag$x[length(dag$x)], dag$y[length(dag$y)],
              bquote(.(dag$symbols[length(dag$x)])));

  if( (noxy==0) || (noxy==2) )
  { garrows(dag$x[1], dag$y[1], dag$x[length(dag$x)], dag$y[length(dag$y)], dag$xgap, dag$ygap, dag$len);
    if(noxy==0) text(0.5*(dag$x[1]+dag$x[length(dag$x)]),
                     0.5*(dag$y[1]+dag$y[length(dag$y)])+dag$ygap/2, "?");
  }

# write covariate symbols
  i1<-1;
  i_c<-0; # covariable counter for subscripts
  i_u<-0; # unknown covs counter for...
  nodes<-length(dag$x);
  while (i1 < nodes-1)
  {
    i1<-i1+1;
    if(dag$names[i1]=="unknown" || dag$cov.types[i1]==2)
    { i_u<-i_u+1;
      if( (alt.symb==FALSE) || (is.na(dag$symbols[i1])) ) {
        text(dag$x[i1], dag$y[i1], bquote(U[.(i_u)]));
      } else text(dag$x[i1], dag$y[i1], bquote(.(dag$symbols[i1])));

    } else
    { i_c<-i_c+1;
      if(is.in(i1, dag$adj)==TRUE)
      { if( (alt.symb==FALSE) || (is.na(dag$symbols[i1])) ) {
          text(dag$x[i1], dag$y[i1], bquote(underline(bar(C))[.(i_c)]));
        } else text(dag$x[i1], dag$y[i1], bquote(underline(bar(.(dag$symbols[i1])))));
      } else
      { if( (alt.symb==FALSE) || (is.na(dag$symbols[i1])) ) {
          text(dag$x[i1], dag$y[i1], bquote(C[.(i_c)]));
        } else text(dag$x[i1], dag$y[i1], bquote(.(dag$symbols[i1])));
      }
    }
  }

# draw legend
  if(legend==TRUE) dag.legend(dag, alt.symb = alt.symb);

# write paths
  if(paths==TRUE) write.paths(dag, alt.symb = alt.symb);
  
# drawing arcs or associations
  i2<-1;
  while (i2 <= length(dag$arc[,1]))
  {
    if(dag$arc.type[i2]==0)
    { garrows(dag$x[dag$arc[i2,1]], dag$y[dag$arc[i2,1]], dag$x[dag$arc[i2,2]], dag$y[dag$arc[i2,2]], dag$xgap, dag$ygap, dag$len);
    } else if(dag$arc.type[i2]==1)
    {
      if(is.na(dag$curve.x[i2])==TRUE)
      { x0<-dag$x[dag$arc[i2,1]]; y0<-dag$y[dag$arc[i2,1]];
        x1<-dag$x[dag$arc[i2,2]]; y1<-dag$y[dag$arc[i2,2]];
        dag$curve.x[i2]<- -(y1-y0)/15+(x0+x1)/2;
        dag$curve.y[i2]<-  (x1-x0)/15+(y0+y1)/2;
      }
      smoothArc(A=c(dag$x[dag$arc[i2,1]], dag$y[dag$arc[i2,1]]),
                B=c(dag$x[dag$arc[i2,2]], dag$y[dag$arc[i2,2]]),
                C=c(dag$curve.x[i2], dag$curve.y[i2]),
                gap=max(dag$xgap, dag$ygap), p=p); 
    }
    if(numbering==TRUE)
    { if(dag$arc.type[i2]!=1)
      { text(0.5*(dag$x[dag$arc[i2,1]]+dag$x[dag$arc[i2,2]]),
             0.5*(dag$y[dag$arc[i2,1]]+dag$y[dag$arc[i2,2]]),
             i2);
      } else
      { text(0.5*(dag$x[dag$arc[i2,1]]+dag$curve.x[i2]),
             0.5*(dag$y[dag$arc[i2,1]]+dag$curve.y[i2]),
             i2);
      }
    }
   i2<-i2+1;
  }
  return(dag);
}

