dag.move <-
function(dag)
{ # first mouse click identifies node or smoothArc to be moved;
  # right click stops the moving;
  dag<-dag.draw(dag, p=TRUE);

  nodes<-length(dag$x);
  arcs<-length(dag$curve.x);
  
  i.move<-NULL;
  while(is.null(i.move)==TRUE || is.na(i.move)==TRUE)
  { i.move<-identify(x=c(dag$x, dag$curve.x),
                     y=c(dag$y, dag$curve.y), n=1, plot=FALSE);
  }
 
  new.xy<-locator(n=1);
  while(is.null(new.xy)==FALSE)
  {
    if(i.move>nodes)
    { dag$curve.x[i.move-nodes]<-new.xy$x;
      dag$curve.y[i.move-nodes]<-new.xy$y;
    } else
    { dag$x[i.move]<-new.xy$x;
      dag$y[i.move]<-new.xy$y;
    }
    dag<-dag.draw(dag, p=TRUE);
    new.xy<-locator(n=1);
  }
  return(dag);
}

