sym.interval.pc <-
function(sym.data, method = c("vertex","centers"),maxit,plot,scale,center)  {
  idn <- all(sym.data$sym.var.types == sym.data$sym.var.types[1])
  
  if (idn == FALSE) 
    stop("All variables have to be of the same type")
  method <- match.arg(method)
  
  if ((sym.data$sym.var.types[1] != "$C") && (sym.data$sym.var.types[1] != "$I")) 
    stop("Variables have to be continuos or Interval")
  else if (sym.data$sym.var.types[1] == "$C") 
    res <- principal.curve(sym.data$data, plot.true = plot,maxit = maxit)
  else if (sym.data$sym.var.types[1] == "$I") {
    
    vertex<-vertex.interval(sym.data)
    individuals<-scale(as.matrix(vertex$vertex),scale = scale , center = center)
    
    
    if (method == "centers") {
      
      centers<-centers.interval(sym.data)
      
      res <- principal.curve(as.matrix(centers), plot.true = plot,maxit = maxit)
      
      n<-dim(individuals)
      
      projection.matrix<-matrix(data = NA,nrow = n[1] , ncol = n[2])
      
      distance.vector<-rep(NA,n[1])
      
      lambda<-rep(NA,n[1])
      
      orthogonal.projection<-rep(NA,n[1])
      
      for(i in 1:n[1])
      {
        neig<-neighbors.vertex(as.matrix(individuals[i,]),res$s,2)
        v<- -neig$neighbors[1,] + neig$neighbors[2,]
        vp<- -neig$neighbors[1,]  + individuals[i,]
        proy<-sum(v*vp)/(norm.vect(v)^2)*v
        proy.point<-neig$neighbors[1,] + proy
        projection.matrix[i,]<-proy.point 
        orthogonal.projection[i]<-sum((vp-proy)*v)
        distance.vector[i]<-norm.vect(vp-proy)
        
        lambda1<-res$lambda[neig$order[1:2]]
        if(lambda1[1] <= lambda1[2])
        {
          lambda[i]<- -lambda1[1] + norm.vect(proy)
        }
        else{
          lambda[i]<- lambda1[1] - norm.vect(proy)
        }
        
      }
      
      res.var.ind<-variance.princ.curve(data = individuals,curve = projection.matrix)
      res.var.mid<-variance.princ.curve(data = as.matrix(centers),curve = res$s)
      res.var<-list(res.var.ind = res.var.ind , res.var.mid = res.var.mid)
      colnames(projection.matrix)<-sym.data$sym.var.names
      res.limits<-sym.interval.pc.limits(sym.data = sym.data,prin.curve = projection.matrix,
                                         num.vertex = vertex$num.vertex,
                                         lambda = lambda,res.var$res.var.mid$var.order)
      
      num.vars<-sym.data$M
      variables<-rep("X",num.vars)
      for(i in 1:num.vars){
        variables[res.var.ind$var.order[i]]<-paste0("prin_surface_" , as.character(i))
      }
      
      colnames(projection.matrix)<-variables
      projection.matrix<-projection.matrix[,res.var.ind$var.order]
      correl<-cor(x = projection.matrix,y = vertex$vertex)
    }
    else if (method == "vertex") {
      res <- principal.curve(individuals, plot.true = plot, maxit = maxit) 
      
      res.var<-variance.princ.curve(data = individuals,curve = res$s)
      
      res.limits<-sym.interval.pc.limits(sym.data = sym.data,prin.curve = res$s,
                                         num.vertex = vertex$num.vertex,
                                         lambda = res$lambda,res.var$var.order)
      num.vars<-sym.data$M
      variables<-rep("X",num.vars)
      for(i in 1:num.vars){
        variables[res.var$var.order[i]]<-paste0("prin_surface_" , as.character(i))
      }
      
      colnames(res$s)<-variables
      res$s<-res$s[,res.var$var.order]
      correl<-cor(x = res$s,y = vertex$vertex)
    }
    
    
    return(list(prin.curve = res , sym.prin.curve = res.limits , var.curve = res.var,cor.ps = correl))
    
  } 
  
  return(TRUE)
}
