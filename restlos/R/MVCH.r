MVCH <- 
function(data, ps=.75, pf=.2, k=1000, a.poi=2, del.poi=1)
{
  #####################################
  # iterative Fast-pface#
  #####################################
  # function arguments
  # -------------------
  # ps:proportion of data selected per step; default 0.75
  # pf: proportion of final subset; default 0.2
  # k:maximum number of iterations; default 1,000
  
  fmvch.pface_sim<-function(data, h, in.subs=50, d.poi, alpha=.95, nc=1, a.poi, d.mat.ov ){
    #############################
    #Fast-MVCH - pface#
    #############################
    # function arguments
    # -------------------
    # data:input data set (as matrix)
    # h: size of the robust subsample to be identified; default: (n+d+1)/2
    # in.subs:number of initial subsamples (of size h) to be drawn; default: 250
    
   
    
    ggst.step<-function(testpts,calpts){
      
	  dist.list <- disthull(p = matrix(data[testpts,], ncol= ncol(data)), D = matrix(data[calpts,], ncol= ncol(data)), ind.test = testpts, ind.cal = calpts)
      
	  distis <- dist.list$dists
      
      add.pts <- order(distis, na.last=T)[1:max(1,min(a.poi, h-length(calpts),length(distis)))]
	  

	  calpts <- c(calpts, testpts[add.pts])
      
	  n.sub <- calpts
      
    }
    
    
    f.pee <- function(n.sub, nc=1, d.poi=1){
	
      d.poi <- max(1,d.poi)
    
		repeat{
			f <- length(n.sub) - h
			repeat{
			  if(file.exists(paste(getwd(),"/qhull_out.txt",sep=""))){
				file.remove(paste(getwd(),"/qhull_out.txt",sep=""))
			  }
			  hull <- convhulln(data[n.sub,], "QJ Fa Pp TO 'qhull_out.txt'")
			  if(is.list(hull)) break
			}
			c.hull <- hull$hull
			d.fac <- nrow(c.hull)
			
			res <- scan(paste(getwd(),"/qhull_out.txt",sep=""),strip.white=T,quiet=T)
			dims <- res[1]
			
			fac.area <- res[-1]
			file.remove(paste(getwd(),"/qhull_out.txt",sep=""))
			
			del.facet <- order(fac.area, decreasing = T)[1:min(d.fac, nc)]
			
			pot.de.pts <- unique(as.numeric(c.hull[del.facet,]))
			
			del.poi <- unique(sample(pot.de.pts, min(f, d.poi, length(pot.de.pts))))
			
			n.sub <- n.sub[-del.poi]
			f <- length(n.sub) - h
			if(f <= 0) break
      }
      return(n.sub)
    }



		
	disthull <- function(D, p, ind.test, ind.cal){
		# D ... matrix of subset coordinates (n1 x d)
		# p ... matrix of points to be tested (n2 x d)
		#
		
		if(ncol(D) != ncol(p)) stop("matrices must have same dimension")

		######################################
		# distance to face function
		######################################

		dist.face <- function(normal, p.mat, poi){
			# normal .... vector of normal (length d)
			# p.mat ... matrix of normal points (dim. d x d)
			# poi ... vector of point to be tested (length d)
			
			if(!is.matrix(p.mat)) stop("p.mat must be matrix")
			
			if(nrow(normal) != ncol(p.mat) & ncol(p.mat) != length(poi)) stop("all arguments must have same dimension")
			
			
			d <- ncol(p.mat)
			
			A <- cbind(normal, t(p.mat))
			#b <- p[1,]
			A <- rbind(A, c(0,rep(1,d)))
			b <- c(poi , 1)
			b <- matrix(b, ncol=1)

			res <- Solve(A=A, B=b)

			# check
			chk <- all(res[-1] > 0) & res[1] > 0
			
			if(chk){
				return(res[1])
				}
			else{
				return(NA)
				}

		}

		con.hull <- convhulln(D, options="Tv n TO 'qhull.out'")		

		path <- paste(getwd(), '/qhull.out', sep="")
		norms <- read.table(path, fill=T)
		norms <- norms[-c(1:2),]
		norms <- norms[,-ncol(norms)]

		con.list <- vector("list", nrow(con.hull))

		for(i in 1:nrow(con.hull)){

			con.list[[i]] <- list(poi.mat = D[con.hull[i,],], normal = norms[i,])

		}
		
		if(nrow(p) > 1){
			dist.mat <- t(sapply(con.list, function(y) apply(p, 1, function(x) dist.face(normal = matrix(as.numeric(y$normal), ncol=1), p.mat = y$poi.mat, poi = x))))
			}
		else{
			dist.mat <- matrix(sapply(con.list, function(y) apply(p, 1, function(x) dist.face(normal = matrix(as.numeric(y$normal), ncol=1), p.mat = y$poi.mat, poi = x))), ncol=1)
	
		}
		
		
		dist.vec <- numeric(nrow(p))

		for(i in 1:ncol(dist.mat)){ 
			
			x <- dist.mat[,i]
				
			if(all(is.na(x))){
				dist.vec[i] <- min(d.mat.ov[ind.test[i], ind.cal])
			}
			else{
				dist.vec[i] <- min(x, na.rm = T)
			}

		}

		return(list(dists = dist.vec, dist.mat = dist.mat, c.hull = con.hull))
	}

        
    n<-nrow(data)
    len <- 1:n
    d <- ncol(data)
    
	 
    i=1
    k=1
    tris <- delaunayn(data, options="Fa Pp Fn TO 'qhull_out.txt'")
    
    tmp <- as.list(readLines(paste(getwd(),"/qhull_out.txt",sep="")))
    
    rdela_area <- as.numeric(unlist(tmp[(nrow(tris)+3):length(tmp)]))
    
    file.remove(paste(getwd(),"/qhull_out.txt",sep=""))
    
    n.tris <- nrow(tris)
	
	in.subs <- min(floor(sqrt(n.tris)), in.subs)
	
	int.subs <- vector("list", in.subs)
	
    if(n.tris > in.subs){
      probs <- 1/(2^(1:n.tris))
      probs <- probs[order(rdela_area)]
      x <- sample(1:n.tris, in.subs, prob=probs)
      for(i in 1:in.subs){
        int.subs[[i]] <- tris[x[i],]
      }
    }else{
      in.subs=n.tris
      for(i in 1:in.subs){
        int.subs[[i]] <- tris[i,]
      }
    }
    
    fin.subs <- vector("list", in.subs)
    vol.subs <- vector("list", in.subs)
    

    for(i in 1:in.subs){
      j <- 1
      n.sub <- int.subs[[i]]
	  
      repeat{
        j<-j+1
        rem.pts <- setdiff(len,n.sub)
        
        n.sub <- ggst.step(calpts=n.sub, testpts=rem.pts)
        
        if(length(n.sub) >= h | j > 100) {
          
          rem.pts <- setdiff(len,n.sub)
		  
		  if(length(rem.pts) > 0){
		  
			repeat{
				if(file.exists(paste(getwd(),"/qhull_out.txt",sep=""))){
				  file.remove(paste(getwd(),"/qhull_out.txt",sep=""))
				}
				hull <- convhulln(data[n.sub,],"QJ FA n Tv Fa Fn Pp TO 'qhull_out.txt'")
				if(is.list(hull)) break
			  }
			  
			  c.hull <- hull$hull
			  
			  res <- scan(paste(getwd(),"/qhull_out.txt",sep=""),strip.white=T,quiet=T)
			  dims <- res[1:2]
			  
			  #Normals
			  
			  nrmls <- matrix(res[3:(2+prod(dims))], dims[2], byrow=T)
			  aN <- nrmls[, d+1]
			  nrmls <- nrmls[ ,-(d+1)] * -1
			  
			  file.remove(paste(getwd(), "/qhull_out.txt", sep=""))
			  
			  if(!is.matrix(data[rem.pts,])==T){
				cx <- 1
			  }
			  else{ 
				cx <- nrow(data[rem.pts,])
				}
			  
			  nt <- nrow(nrmls)
			  			  
			  test.mat <- data[rem.pts,] %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE)
			  #test.mat<-matrix(1, cx, nt)
			  # -1 : outside a facet
			  # 1 : inside a facet
			  
			  tol <- sqrt(.Machine$double.eps)
			  test.mat[abs(test.mat) <= tol] <- 0
			  
			  test.mat <- sign(test.mat)
			  
			  ind.insider <- which(apply(test.mat,1,function(x) all(x %in% c(0,1))))
			  if(any(ind.insider)){
				n.sub <- unique(c(n.sub, rem.pts[ind.insider]))
				rem.pts <- rem.pts[-ind.insider]
				test.mat <- test.mat[-ind.insider,]
			  }
			}  
          if(length(n.sub) > h){
			
			d.poi.2 <- max(d.poi, ceiling((length(n.sub) - h)/2))
            n.sub <- f.pee(n.sub, d.poi = d.poi.2, nc = ceiling(d.poi.2/ncol(data)))
          }
          
          repeat{
            ch.x <- convhulln(data[n.sub,],"Pp QJ FA TO 'qhull_out.txt'")
            if(is.list(ch.x)) break
          }
          
          
          
          break
        }
        
      }
      fin.subs[[i]] <- n.sub
      
      vol.subs[[i]] <- ch.x$vol
    }
    
    
    mins <- which.min(unlist(vol.subs))
    subs <- fin.subs[[mins]]
    
    n.sub <- subs
    vols <- min(unlist(vol.subs))
    
	
	
    return(list(set=n.sub, vol=vols))
    
  }
  
  
  n <- nrow(data)
  lfsub <- ceiling(n*pf)
  lssub <- floor(ps*n)
  i<-1
  subi<-1:n
  
  d.mat.ov.a <- as.matrix(dist(data)) 
  
  
  repeat{
	
	d.mat.ov.aa <- d.mat.ov.a[subi, subi]
	colnames(d.mat.ov.aa) <- rownames(d.mat.ov.aa) <- 1:length(subi)
	
    ergi <- fmvch.pface_sim(data=data[subi,], h=lssub, a.poi=a.poi, d.poi=del.poi, d.mat.ov = d.mat.ov.aa)
    
	
    subi <- subi[ergi$set]
    lssub <- ceiling(length(subi)*ps)
    i<-i+1
    if(i == 2){
		set.1 <- subi
		}
    if(lssub <= lfsub | lssub<=ncol(data) | i>k) break
    
  }
  tmp.mode <- colMeans(data[subi,])
  
  return(list(mode=tmp.mode, set=subi, vol=ergi$vol, set.1 = set.1))
  
}
