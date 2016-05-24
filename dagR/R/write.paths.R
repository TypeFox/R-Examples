write.paths <-
function(dag, px=0.5, py=-0.060, alt.symb = TRUE)
{ # writes the paths in the DAG graph;
  # alt.symb is new in version 1.1.2;
 if(is.null(dag$paths)==FALSE)
 {
  #if(dag$pathsN==1){ dag$paths<-matrix(dag$paths, nrow=1); }

  if(is.null(dag$symbols)) {
    dag$symbols<-rep(NA, length(dag$x));
    writeLines("WARNING! DAG's from version 1.1.2 on should contain a $symbols array. This has been added only temporarily by write.paths().");
  }

  i1<-0;                         # path number
  max.l<-0;                      # maximal path length
  while (i1 < dag$pathsN)
  {
    i1<-i1+1;
    cur.node<-1;                 # current node
    dag.letter(dag, 1, x=px, y=py-i1*0.065, alt.symb = alt.symb);    
    i2<-1; # arc counter
    arc.typ<-'';
    while(is.na(dag$paths[i1, i2])==FALSE)
    {
      if( (dag$arc[dag$paths[i1, i2], 1] != cur.node) )
      { cur.node<-dag$arc[dag$paths[i1, i2], 1];
        arc.typ<-'<';
      } else
      { cur.node<-dag$arc[dag$paths[i1, i2] ,2];
        arc.typ<-'>';
      }
      if(dag$arc.type[dag$paths[i1, i2]]==1)
      { arc.typ<-'-';
      }
      text(px+i2*0.1-0.05, py-i1*0.065, arc.typ);
      dag.letter(dag, cur.node, x=px+i2*0.1, y=py-i1*0.065, alt.symb = alt.symb);
      i2<-i2+1;
      if(i2>max.l){max.l<-i2;}
    }
  }
  for(i1 in 1:dag$pathsN)
  { text(px+max.l*0.1-0.05, py-i1*0.065, dag$path.status[i1], pos=4);
  }
 }
}

