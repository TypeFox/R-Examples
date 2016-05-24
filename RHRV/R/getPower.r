getPower=function(wx,nodes,wf)
{   
  # sum contributions
  power=0;
  l=length(nodes)

  for (j in 2:l)
  {

      # index for wavelets coefficients returned by waveslim::modwt
      depth=nodes[[j]][1]
      numberNode=nodes[[j]][2]
      index=sum(2^(1:depth-1))+numberNode;
      power=align(wx[[index]]^2,depth,numberNode,wf)+power;

  }


  return(power);
}

