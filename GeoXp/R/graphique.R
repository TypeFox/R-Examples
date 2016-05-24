graphique<-
function (var1, var2, var3, obs, num, graph = "", couleurs = "",
    symbol = 16, labvar = "", nbcol = 10, alpha1, W,
    Xpoly, Ypoly, F, G, opt1 = 1, opt2 = 1, quantiles = 0,
    labmod = "", direct=NULL, inertie, label = 0, kernel, obsq, locmoran = FALSE,
    bin = NULL,cex.lab=1, buble=FALSE,cbuble=NULL,legmap=NULL,legends=list(FALSE,FALSE),
    xlim,ylim)
{
    dev.set(num)
  #  par(mar = c(4.5, 4.5, 0.5, 0.5))

####################################################
# Initialisation des bubbles
####################################################

if(buble && (length(cbuble)!=0))
{
  if(min(cbuble)==0)
   {cbuble[which(cbuble==0)]<-min(cbuble[which(cbuble!=0)])/2}

  if(length(legmap)>0)
   {
    if(as.numeric(legmap[3])==0)
     {legmap[3]<-min(cbuble)}
   }
}
else
   {cbuble=rep(0.7,length(var1))
    cbuble[which(obs==TRUE)]=1 
   }
 
####################################################
# Graphique par graphique
####################################################
        
    if (graph == "Histogram") 
    {
        if(is.null(bin)) bin<-"count"
        mnv <- min(var1)
        mxv <- max(var1)
        h <- (mxv - mnv)/(nbcol)
        cpt <- NULL
        absc <- mnv
        
       # analysis of the first bar 
        x1 <- mnv + h 
        
        cnt = length(which((var1 >= absc[1]) & (var1 <= x1)))
        
        cpt <- c(cpt, cnt)
        absc <- c(absc, x1)
        
        # analysis of the others bars            
        for (i in 2:nbcol) 
        {
          x1 <- mnv + h * i
          cnt <- length(which((var1 > absc[i]) & (var1 <= x1))) 
                                   
          cpt <- c(cpt, cnt)
          absc <- c(absc, x1)
        }
       
       sum.cpt<-sum(cpt) 


        cpt <- c(cpt, 0)
        if(bin[1]=="percent")
         {cpt<-cpt/sum.cpt
          if(labvar[2]=="") labvar[2]="percent"
         }
         else
          {if(bin[1]=="density")
           {cpt<-cpt/sum.cpt/h
           if(labvar[2]=="") labvar[2]="density"
           }
          }
        
        if(labvar[2]=="") labvar[2]="counts"
        plot(absc, cpt, "n", xlim = c(mnv - h, mxv + h), ylim = c(0,max(cpt)),
        xlab = labvar[1], ylab = labvar[2],frame.plot=FALSE)
        
        for (i in 1:nbcol) 
         {
          rect(absc[i], 0, absc[i + 1], cpt[i], col = couleurs)
         }
         
        if (length(var1[obs]) != 0) 
        {
            vrob <- var1[obs]
            cpt1 <- NULL
           
           # analysis of the first bar
            i <- 1
            cnt <- length(which((vrob >= absc[1]) & (vrob <= absc[2]))) 
            cpt1 <- c(cpt1, cnt)
            
            # analysis of the others bars
            for (i in 2:nbcol) 
            {
             cnt <- length(which((vrob > absc[i]) & (vrob <= absc[i + 1])))
             cpt1 <- c(cpt1, cnt)
            }
            
            cpt1 <- c(cpt1, 0)
            
        if(bin[1]=="percent")
         {cpt1<-cpt1/sum.cpt
         }
         else
          {if(bin[1]=="density")
           {cpt1<-cpt1/sum.cpt/h
           }
          }
          
            for (i in 1:nbcol) 
               {
                  rect(absc[i], 0, absc[i + 1], cpt1[i], col = "yellow")
                  rect(absc[i], 0, absc[i + 1], cpt1[i], density = 4,
                  angle = 45, lwd = 1, col = "red")
               }
        }        
    }
    
    
    
    if (graph == "histo.nb") 
     {
      if(is.null(bin)) bin<-"count"
      nblist<-unlist(var2)
        
      mnv <- min(nblist)
      mxv <- max(nblist)
        
      h <- (mxv - mnv)/(nbcol)
      cpt <- NULL
      absc <- mnv
      
      # analysis of 1st bar  
      x1 <- mnv + h  
      cnt<-length(which((nblist >= absc[1]) & (nblist <= x1)))
      cpt <- c(cpt, cnt)
      absc <- c(absc, x1)  
      
        for (i in 2:nbcol) 
        {
          x1 <- mnv + h * i        
          cnt<-length(which((nblist> absc[i]) & (nblist <= x1))) 
          cpt <- c(cpt, cnt)
          absc <- c(absc, x1)
        }
        
      sum.cpt<-sum(cpt)
      cpt <- c(cpt, 0)
      
       if(bin[1]=="percent")
        {cpt<-cpt/sum.cpt
         if(labvar[2]=="") labvar[2]="percent"
        }
       else
        {if(bin[1]=="density")
          {cpt<-cpt/sum.cpt/h
           if(labvar[2]=="") labvar[2]="density"
          }
        }
          
      plot(absc, cpt, "n", xlim = c(mnv - h, mxv + h), ylim = c(0,max(cpt)),
      xlab = labvar[1], ylab = labvar[2])
      
       for (i in 1:nbcol) 
        {
         rect(absc[i], 0, absc[i + 1], cpt[i], col=couleurs[1])
        }

       if (sum(as.integer(obs)) != 0) 
       {         
        ident.select<-which(obs==TRUE,arr.ind=TRUE)
        n.ident.select<-dim(ident.select)[1]
            
        vrob<-NULL
          
          for (j in 1:n.ident.select)
           {vrob <- c(vrob,var2[[ident.select[j,1]]][which(var1[[ident.select[j,1]]]==ident.select[j,2])])}
             
            cpt1 <- NULL
            
            cnt <- length(which((vrob >= absc[1]) & (vrob <= absc[2])))
            cpt1 <- c(cpt1, cnt)
            
          
            for (i in 2:nbcol) 
            {
              cnt<-length(which((vrob > absc[i]) & (vrob <= absc[i + 1]))) 
              cpt1 <- c(cpt1, cnt)
            }
            
            cpt1 <- c(cpt1, 0)
            
           if(bin[1]=="percent")
            {cpt1<-cpt1/sum.cpt
            }
           else
            {if(bin[1]=="density")
             {cpt1<-cpt1/sum.cpt/h
             }
            }
           
            for (i in 1:nbcol) 
            {
             rect(absc[i], 0, absc[i + 1], cpt1[i], density = 4,
             angle = 45, lwd = 1, col = "red")

            }
        #    lines(density(vrob),col='red1',lwd=1.5,lty=2)
        }
           #     lines(density(nblist),col=couleurs[2],lwd=2,lty=2)
                
    }
    
     
    
    if ((graph == "Barplot") || (graph == "Cluster")) 
    {
    
      if(is.null(bin)) bin<-"count"
      r <- table(var1)
      nomsr <- names(r)
      r.total<-r
      
      if(bin[1]=="percent") r<-r/sum(r)
       
      if (labmod[1] == "") 
       {
        labmod <- nomsr
       }
        
        
      if (length(couleurs) != length(levels(as.factor(var1)))) 
        {
         couleurs=rep(couleurs[1],length(levels(as.factor(var1))))
        }
        
        g <- barplot(r, xlim = c(0, length(r)), names.arg = labmod,
        width = 0.8, col = couleurs, lwd = 1, xlab = labvar[1],
        ylab = labvar[2])

         
        if (length(var1[obs]) != 0) 
        {
         t.obs <- table(var1[obs])
         nomst <- names(t.obs)
         if(bin[1]=="percent") t.obs<-t.obs/sum(r.total)
               
         for (i in 1:length(t.obs)) 
          {
           quit <- FALSE
           j <- 1
          
             while ((j <= length(r)) && (!quit)) 
             {
              if (nomst[i] == nomsr[j]) 
                  {rect(g[j] - 0.4, 0, g[j] + 0.4, t.obs[[i]],
                   col = "yellow")                
                   rect(g[j] - 0.4, 0, g[j] + 0.4, t.obs[[i]],
                   col = "red", density = 4, angle = 45,lwd = 1)
                   quit <- TRUE
                  }
              j <- j + 1
             }
          }
        }
    }
    
    
    if ((graph == "bar.nb") ) 
    {
      r <- table(card(var1))
      nomsr <- names(r)
     
      if (labmod[1] == "") 
      {
            labmod <- nomsr
      }
               
       g <- barplot(r, xlim = c(0, length(r)), names.arg = labmod,
       width = 0.8, col =couleurs, xlab = labvar[1], ylab = labvar[2])

        if (sum(as.integer(obs)) != 0) 
        {           
         ident.select<-which(obs==TRUE,arr.ind=TRUE)
         ind.b<-unique(ident.select[,1])
         t <- table(card(var1)[ind.b])
         
         nomst <- names(t)
            
           for (i in 1:length(t)) 
           {
             quit <- FALSE
             j <- 1
             
              while ((j <= length(r)) && (!quit)) 
              {
               if (nomst[i] == nomsr[j]) 
                  {
                   rect(g[j] - 0.4, 0, g[j] + 0.4, t[[i]],
                   col = "red", density = 4, angle = 45,lwd = 1)
                        
                   quit <- TRUE
                  }
               j <- j + 1
              }
           }
        }
    }    
    


    if (graph == "Scatterplot" || graph=="Acp1")
    {
        if (length(quantiles)!=0) couleur.quant <- heat.colors(length(quantiles))

        if(graph=="Acp1")
        {
         labvar[1] <- paste("component ", labvar[1], " : ", round(inertie[direct[1]], 0), "%", sep = "")
         labvar[2] <- paste("component ", labvar[2], " : ", round(inertie[direct[2]], 0), "%", sep = "")
        }
        
        layout(matrix(c(1, 3, 0, 2), 2, 2, byrow = TRUE), c(1,4), c(5, 1) )
        par(mar = c(2, 1, 2, 2))
        boxplot(var2, axes = FALSE)
        title(ylab = labvar[2], line = 0)
        par(mar = c(1, 2, 1, 2))
        boxplot(var1, horizontal = TRUE, axes = FALSE)
        title(xlab = labvar[1], line = 0)
        par(mar = c(2, 2, 2, 2))
        plot(var1, var2, "n", xlab = "", ylab = "")
        points(var1[!obs], var2[!obs], col = couleurs, pch = 16,cex=0.8)

        if(graph=="Acp1")
        {
        segments(min(var1), 0, max(var1), 0, col = "black")
        segments(0, min(var2), 0, max(var2), col = "black")
        }
        
        if (opt1 & graph == "Scatterplot" )
        {
            xg <- seq(min(var1), max(var1), length = 100)
            reg <- lm(var2 ~ var1)
            lines(xg, reg$coefficients[1] + reg$coefficients[2] * xg)
        }
        
        if(quantiles)
        {
              etendue <- diff(range(var1))
               fitmax  <- rqss(var2 ~ qss(var1, constraint= "N",lambda=etendue), tau = alpha1, control=sfn.control(warn.mesg=FALSE))    
            plot.rqss(fitmax,add=TRUE,rug=FALSE,titles="",col='blue') 
           # Dat.nls <- nlrq(var2 ~ SSlogis(var1, Asym, mid, scal),tau=alpha1)
           # lines(sort(var1), predict(Dat.nls, newdata=list(var1=sort(var1))),col='blue')
        }
        
        if (length(var1[obs]) != 0) 
        {
         points(var1[obs], var2[obs], col = "red", pch = symbol,cex = 1.2)

          if(graph=="Acp1")
           {
            if(is.logical(label))
            {
             if (label)
             {
              msg <- paste(round(labmod, 0), "%", sep = "")
              text(var1[obs] + 0.02, var2[obs] + 0.02, msg[obs], cex = cex.lab, font=3,adj=c(0.75,-0.75))
             }
           }
         }
        }
       if(dev.cur()!=3)
       {
       layout(1)
       par(mar=c(5.1,4.1,4.1,2.1))
       }
    }
    
    if (graph == "Boxplot") 
    {
        r <- boxplot(var1, boxwex = 0.5, xlab = labvar[1], ylab = labvar[2],
        col=couleurs)
        
        mat <- r$stats
        out <- r$out
        quit1 <- FALSE
        quit2 <- FALSE
        quit3 <- FALSE
        quit4 <- FALSE
        if (length(var1[obs]) != 0) 
        {
        
        n1 <- length(which((var1[obs] >= mat[1, 1]) & (var1[obs] < mat[2, 1])))
        n2 <- length(which((var1[obs] >= mat[2, 1]) & (var1[obs] < mat[3, 1])))
        n3 <- length(which((var1[obs] >= mat[3, 1]) & (var1[obs] < mat[4, 1])))
        n4 <- length(which((var1[obs] >= mat[4, 1]) & (var1[obs] < mat[5, 1]))) 
        
        ifelse(length(out)==0,n5<-0,n5<-length(which((var1[obs] >= mat[5, 1]) & (var1[obs] <= max(out)))))
        ifelse(length(out)==0,n6<-0,n6<-length(which((var1[obs] < mat[1, 1]) & (var1[obs] >= min(out))))) 
              
        if(n1>0 && (!quit1)) 
        {
         rect(0.92, mat[1, 1], 1.08, mat[2, 1], col = "red",
         density = 4, angle = 45, lwd = 1)
                   
         quit1 <- !quit1
         }
               
        if (n2 >0 && (!quit2) ) 
        {
         rect(0.875, mat[2, 1], 1.125, mat[3, 1],
         col = "red", density = 4, angle = 45, lwd = 1)
                                        
         quit2 <- !quit2
         }
                
        if (n3>0 && (!quit3)) 
         {
          rect(0.875, mat[3, 1], 1.125, mat[4, 1],
          col = "red", density = 4, angle = 45, lwd = 1)
                    
          quit3 <- !quit3
          }
          
        if (n4>0 && (!quit4)) 
        {
         rect(0.92, mat[4, 1], 1.08, mat[5, 1], col = "red",
         density = 4, angle = 45, lwd = 1)
                  
         quit4 <- !quit4
         }
     
        if(n5>0)
         {points(rep(1,n5), var1[obs][which((var1[obs] >= mat[5, 1]) & (var1[obs] <= max(out)))],
          col = "red", pch = symbol)}
          
        if(n6>0)
         {points(rep(1,n6), var1[obs][which((var1[obs] < mat[1, 1]) & (var1[obs] >= min(out)))],
          col = "red", pch = symbol)}
                     
        }
    }
    
        if (graph == "bp2") {
         if (length(var1[obs]) != 0) {
         boxplot(c(var1,var1[obs])~c(rep(1,length(var1)),rep(2,length(var1[obs]))), varwidth=TRUE, xlab = labvar,names=c("Full sample","Selected sample"))           
        }
    }
    
    
    if (graph == "Polyboxplot") 
    {
     x <- factor(var2)
     
     if (labmod[1] == "") 
     {
      labmod <- levels(x)
     }
               
      if (length(couleurs) != length(levels(as.factor(var2)))) 
      {
       couleurs=rep(couleurs[1],length(levels(as.factor(var2))))
      }

      if (length(symbol) != length(levels(as.factor(var2)))) 
      {
       symbol=rep(symbol[1],length(levels(as.factor(var2))))
      }
              
     r <- boxplot(var1 ~ x, boxwex = 0.8, xlab = labvar[1],
     ylab = labvar[2], col = couleurs, pch=symbol, names = labmod,varwidth = bin)
     
     r <- boxplot(var1 ~ x, boxwex = 0.8, xlab = labvar[1],
     ylab = labvar[2], col = couleurs, pch=symbol, plot = FALSE,varwidth = bin)
     

     mat <- r$stats
     out <- r$out
     quit1 <- FALSE
     quit2 <- FALSE
     quit3 <- FALSE
     quit4 <- FALSE
     
     if (length(var1[obs]) != 0) 
     {
      quit <- FALSE
  #    k <- 1 
      for (k in 1:length(r$n)) 
      {   
       for (j in 1:length(var1[obs])) 
       {
        if (as.character(var2[obs][j]) == as.character(r$names[k])) 
        {
         if ((var1[obs][j] >= mat[1, k]) && (var1[obs][j] < mat[2, k])) 
         {
           rect(k - 0.2, mat[1, k], k + 0.2, mat[2,k], col = "red", density = 4, angle = 45,
           lwd = 1)
         }
        
         if ((var1[obs][j] >= mat[2, k]) && (var1[obs][j] < mat[3, k])) 
         {
          rect(k - 0.4, mat[2, k], k + 0.4, mat[3,k], col = "red", density = 4, angle = 45,
          lwd = 1)
         }
    
         if ((var1[obs][j] >= mat[3, k]) && (var1[obs][j] < mat[4, k])) 
         {
          rect(k - 0.4, mat[3, k], k + 0.4, mat[4,k], col = "red", density = 4, angle = 45,
          lwd = 1)
         }
        
         if ((var1[obs][j] >= mat[4, k]) && (var1[obs][j] <= mat[5, k])) 
         {
          rect(k - 0.2, mat[4, k], k + 0.2, mat[5,k], col = "red", density = 4, angle = 45,
          lwd = 1)
         }

         if (length(as.vector(out)) > 0) 
         {
          if ((var1[obs][j] >= mat[5, k]) && (var1[obs][j] <= max(out)) && 
          (as.character(var2[obs][j]) == as.character(r$names[k]))) 
          {
           points(k, var1[obs][j], col = "red", pch = symbol)
           }
          
          if ((var1[obs][j] < mat[1, k]) && (var1[obs][j] >= min(out))) 
          {
           points(k, var1[obs][j], col = "red", pch = symbol)
          }
         }
        }
       }
   #    k <- k + 1
       }
      }
    }
    
    if (graph == "Densityplot1") 
    {
        ypartsup <- 0
        h <- (alpha1/100) * (max(var1) - min(var1))/2
        dens <- bkde(var1, kernel = kernel, bandwidth = h, gridsize = 100)
        ysup <- max(dens$y) + 0.1 * max(dens$y)
      
        if (length(var1[obs]) != 0) 
        {
          h2 <- (alpha1/100) * (max(var1[obs]) - min(var1[obs]))/2
          denspart <- bkde(var1[obs], kernel = kernel, bandwidth = h2,gridsize = 100)
          ypartsup <- max(denspart$y) + 0.1 * max(denspart$y)
        }
       
        plot(dens, type = "l", col = couleurs[1], ylim = c(0, max(ysup,ypartsup)),
        xlab = labvar[1], ylab = labvar[2])
       
        if (length(var1[obs]) != 0) 
        {
          points(denspart, type = "l", col = "red", lty=2)
        }
    }
    
    if (graph == "Densityplot2") 
    {
        h <- (alpha1/100) * (max(var1) - min(var1))/2
        aire = 0
        
        dens <- bkde(var1, kernel = kernel, bandwidth = h, gridsize = 100)
        
        plot(dens, type = "l", col = couleurs[1], ylim = c(0, max(dens$y) +
        0.1 * max(dens$y)), xlab = labvar[1], ylab = labvar[2])

        interv.k<-NULL

        if(length(Xpoly)>0)
        {
         for (k in 1:length(Xpoly))
         { i <- 1
             while ((Xpoly[[k]][1] > dens$x[i]) && (i < 99)) 
             {
              i <- i + 1
             }
          i1 <- i
          i <- 1
            
            while ((Xpoly[[k]][2] > dens$x[i]) && (i < 99)) 
            {
             i <- i + 1
            }
      
           interv.k<-setdiff(union(interv.k,i1:i),intersect(interv.k,i1:i) )

        }   
    
    i<-1
    abs.poly<-c(dens$x[interv.k[1]],dens$x[interv.k[1]])
    ord.poly<-c(0,dens$y[interv.k[1]])
    
        while(i <= (length(interv.k)-1)) 
        {
         if(interv.k[i+1]==(interv.k[i]+1))
          {abs.poly<-c(abs.poly,dens$x[interv.k[i+1]])
           ord.poly<-c(ord.poly,dens$y[interv.k[i+1]])            
        
           aire = aire + (dens$x[interv.k[i] + 1] - dens$x[interv.k[i]]) * dens$y[interv.k[i]] +
           (dens$y[interv.k[i] + 1] - dens$y[interv.k[i]]) * (dens$x[interv.k[i] + 1] -
            dens$x[interv.k[i]])/2         
          
          i<-i+1
          }        
         else
         {polygon(c(abs.poly,dens$x[interv.k[i]]),c(ord.poly,0), col = "red",density=10, angle=45)
          
          abs.poly<-c(dens$x[interv.k[i+1]],dens$x[interv.k[i+1]])
          ord.poly<-c(0,dens$y[interv.k[i+1]])       
          
            aire = aire + (dens$x[interv.k[i+1]] - dens$x[interv.k[i+1]-1]) * dens$y[interv.k[i+1]-1] +
            (dens$y[interv.k[i+1]] - dens$y[interv.k[i+1]-1]) * (dens$x[interv.k[i+1]] -
             dens$x[interv.k[i+1]-1])/2
         
          i<-i+1   
         }               
        }
    
    polygon(c(abs.poly,dens$x[interv.k[i]]),c(ord.poly,0), col = "red",density=10, angle=45)
             
    #    for (i in 1:(length(interv.k)-1)) 
    #    {              
    #     polygon(c(dens$x[interv.k[i]], dens$x[interv.k[i]], dens$x[interv.k[i] + 1], dens$x[interv.k[i] + 1]),
    #     c(0, dens$y[interv.k[i]], dens$y[interv.k[i] + 1], 0), col = "red",border=0)
         
    #     aire = aire + (dens$x[interv.k[i] + 1] - dens$x[interv.k[i]]) * dens$y[interv.k[i]] +
    #     (dens$y[interv.k[i] + 1] - dens$y[interv.k[i]]) * (dens$x[interv.k[i] + 1] -
    #      dens$x[interv.k[i]])/2
                  
    #     if(interv.k[i+1]!=interv.k[i]+1)
    #       {polygon(c(dens$x[interv.k[i+1]-1], dens$x[interv.k[i+1]-1], dens$x[interv.k[i+1]], dens$x[interv.k[i+1]]),
    #        c(0, dens$y[interv.k[i+1]-1], dens$y[interv.k[i+1]], 0), col = "red",border=0)
                
    #        aire = aire + (dens$x[interv.k[i+1]] - dens$x[interv.k[i+1]-1]) * dens$y[interv.k[i+1]-1] +
    #        (dens$y[interv.k[i+1]] - dens$y[interv.k[i+1]-1]) * (dens$x[interv.k[i+1]] -
    #         dens$x[interv.k[i+1]-1])/2
    #      }     
    #    }
      }  
    }
    
    if ((graph == "Angleplot")||(graph == "pairwise"))
    {

      ifelse(quantiles,color.quant <- colors()[grep("red", colors())],color.quant<-couleurs[1])
       
      diag(var1) <- NA
      diag(var2) <- NA
      vect1 <- as.vector(var1[which(!is.na(var1))])
      vect2 <- as.vector(var2[which(!is.na(var2))])
      vect <- cbind(vect1, vect2)
      res <- sort(vect[, 1], index.return = TRUE)
      x <- res$x
      y <- vect[res$ix, 2]
      z <- seq(1, length(x), by = (length(x)/1000))
      z <- round(z)
      

        if(quantiles) 
        {
          fitmax  <- rqss(vect2 ~ qss(vect1, constraint= "N", lambda=diff(range(vect1))), tau = alpha1, control=sfn.control(warn.mesg=FALSE))
        }

        if(quantiles) 
        {
            plot(vect1, vect2, "n", xlab = labvar[1], ylab = labvar[2],
            xlim = c(0, max(vect1)),axes=FALSE)
            if (graph == "Angleplot")
            {axis(1, c(0,pi/4,pi/2,3*pi/4,pi), c("0",expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi)))}
            else
            {axis(1)}
            axis(2)                             
            points(vect1[which(vect2> predict(fitmax, newdata=list(vect1=vect1)))],
            vect2[which(vect2 > predict(fitmax, newdata=list(vect1=vect1)))],
            col = couleurs, pch = 16, cex = 0.8)
            
            #lines(sort(vect1), predict(Dat.nls, newdata=list(variable.x=sort(vect1))),col='blue')
            plot.rqss(fitmax,add=TRUE,rug=FALSE,titles="")  
        }
        else 
        {
          plot(vect1, vect2, "n", xlab=labvar[1], ylab = labvar[2],
          xlim = c(0, max(vect1)),axes=FALSE)
            if (graph == "Angleplot")
            {axis(1, c(0,pi/4,pi/2,3*pi/4,pi), c("0",expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi)))}
            else
            {axis(1)}
            axis(2)
          points(vect1, vect2, col = couleurs, pch = 16, cex = 0.8)
        }
       
       if (length(var1[obs]) != 0) 
         {points(var1[obs], var2[obs], col = "red", pch = symbol,cex = 1)}

       if((labvar[1]!="Rank")&(length(direct)>0))
       {abline(v=direct)}
    }


    if (graph == "pairwise2")
    {

      diag(var1) <- NA
      diag(var2) <- NA
      
      vect1 <- as.vector(var1[which(!is.na(var1))])
      vect2 <- as.vector(var2[which(!is.na(var2))])
      vect <- cbind(vect1, vect2)
      res <- sort(vect[, 1], index.return = TRUE)
      x <- res$x
      y <- vect[res$ix, 2]
      z <- seq(1, length(x), by = (length(x)/1000))
      z <- round(z)


      xg <- seq(min(vect[z, 1]), max(vect[z, 1]), length = 100)

          plot(vect1, vect2, "n", xlab=labvar[1], ylab = labvar[2],
          xlim = c(0, max(vect1)),axes=FALSE)
          if(length(bin)!=0)
          {abline(v=bin)}
          axis(1)
          axis(2)
          points(vect1, vect2, col = couleurs, pch = 16, cex = 0.8)


        if(length(var1[obs]) != 0)
         {points(var1[obs], var2[obs], col = "red", pch = symbol,cex = 1)}

    }
    
    
    if (graph == "Variocloud") 
    {
        ifelse(quantiles,couleur <- heat.colors(length(quantiles)),couleur<-couleurs[1])

        msg <- paste("Studied variable : ", labvar[1], sep = "")
        
        
        diag(var1) <- NA
        diag(var2) <- NA
        diag(var3) <- NA
        
        vect1 <- as.vector(var1[which(!is.na(var1))])
        vect2 <- as.vector(var2[which(!is.na(var2))])
        vect3 <- as.vector(var3[which(!is.na(var2))])
        
        vect1 <- vect1[which(!(vect1 == -10))]
        vect2 <- vect2[which(!(vect2 == -10))]
        vect3 <- vect3[which(!(vect3 == -10))]
        vect <- cbind(vect1, vect2, vect3)
        
        if(length(xlim)!=2) xlim = c(0, max(vect1))
           
        res <- sort(vect[, 1], index.return = TRUE)
        x <- res$x
        y <- vect[res$ix, 2]
        y2 <- vect[res$ix, 3]
        z <- seq(1, length(x), by = (length(x)/1000))
        z <- round(z)
    
        vect1 <- vect[,1]
        vect2 <- vect[,2]
       
        etendue <- max(vect[z, 1]) - min(vect[z, 1])
                
        if (quantiles) 
        {
            fitmax  <- rqss(vect2 ~ qss(vect1, constraint= "N",lambda=etendue), tau = alpha1, control=sfn.control(warn.mesg=FALSE))                 
        }
       
        n <- length(x)
        d <- max(x[2:length(x)] - x[1:(length(x) - 1)])
        xg <- seq(min(vect[z, 1]), max(vect[z, 1]), length = 100)
        
        ifelse(length(bin) == 0, xg2 <- seq(min(vect1), max(vect1), length = 1/(max(c(30/n,0.08, d/etendue)))),xg2 <- sort(bin) )
        
        xg3 <- c(xg2[1],(xg2[1:(length(xg2) - 1)] + xg2[2:length(xg2)])/2)

        if(quantiles) 
        {
         plot(vect1, vect2, "n", xlab = "Distance", ylab = "semivariance",
         xlim = xlim, ylim=ylim)
        
           points(vect1[which(vect2> predict(fitmax, newdata=list(vect1=vect1)))],
            vect2[which(vect2 > predict(fitmax, newdata=list(vect1=vect1)))],
            col = couleurs, pch = 16, cex = 0.8)
                  
          plot.rqss(fitmax,add=TRUE,rug=FALSE, titles="", col = couleur[1]) 
        }
        else 
        {
            plot(vect1, vect2, "n", xlab = "Distance", ylab = "semivariance", 
            xlim = xlim, ylim=ylim)
            
            points(vect1, vect2, col = couleurs, pch = 16, cex = 0.8)
        }
        
        if (opt1 == 2) 
        {
           dis = matrix(0, ncol = (length(xg2) ))
               
            dis[1]<-sum(y[which(x <= xg3[2])])/length(which(x <= xg3[2]))
            
            for (i in 2:(length(xg3)-1)) 
             {
               dis[i] <- sum(y[which((x <= xg3[i+1]) & (xg3[i-1] < x ))])/length(which((x <= xg3[i+1]) & (xg3[i-1] < x )))
             }
        
             dis[length(xg3)]<-sum(y[which((x <= xg3[length(xg3)]) & (xg3[length(xg3)-1] < x ))])/length(which((x <= xg3[length(xg3)]) & (xg3[length(xg3)-1] < x )))
           
            points(xg3, dis, col = "red", pch = 1)
            lines(xg3, dis, col = "red")
      
        }

       
        if (opt2 == 2) 
        {
            dis2 = matrix(0, ncol = (length(xg2)))
            Nh = matrix(0, ncol = (length(xg2)))
                 
            dis2[1]<-sum(y2[which(x <= xg3[2])])
            Nh[1] <- length(which(x <= xg3[2]))
            
            for (i in 2:(length(xg3)-1)) 
             {
               dis2[i] <- sum(y2[which((x <= xg3[i+1]) & (xg3[i-1] < x ))])
               Nh[i] <- length(which((x <= xg3[i+1]) & (xg3[i-1] < x )))
             }
        
             dis2[length(xg3)]<-sum(y2[which((x <= xg3[length(xg3)]) & (xg3[length(xg3)-1] < x ))])
             Nh[length(xg3)] <- length(which((x <= xg3[length(xg3)]) & (xg3[length(xg3)-1] < x )))


             dis2 <- (dis2/Nh)^4/(0.457 + 0.494/Nh)/2
             
        
             points(xg3, dis2, col = "red", pch = 2)
             lines(xg3, dis2, col = "red", lty = 2)
        

        }     
   
        
        if (length(var1[obs]) != 0) 
        {
          points(var1[obs], var2[obs], col = "red", pch = symbol,
          cex = 1)
        }
    }
    
    if (graph == "Lorentz") 
    {
        G <- G[cumsum(as.data.frame(table(F))$Freq)]
        F <- F[cumsum(as.data.frame(table(F))$Freq)]
        var1u = unique(var1)
        var1u <- c(0, var1u)
        plot(F, G, pch=symbol,  col = couleurs,xlab=labvar[1],ylab=labvar[2])
        segments(0, 0, 1, 1, col = "blue")
        lines(F, G, col = couleurs)
        if (length(var1[obs]) != 0) 
        {
            imax = which.min(F <= Xpoly)
            imax = imax + 1
            vsort <- sort(var1u)
            vsort <- c(0, vsort)
            segments(F[imax - 1], 0, F[imax - 1], G[imax - 1],
            col = "red")
            
            segments(F[imax - 1], G[imax - 1], 0, G[imax - 1],
                col = "red")
            if (F[2] > G[2]) 
            {
                msg <- paste("f = ", abs(round(F[imax - 1], 2)))
                text(0.15, 0.98, msg, cex = 0.8)
                msg <- paste("g = ", abs(round(G[imax - 1], 2)))
                text(0.15, 0.94, msg, cex = 0.8)
                msg <- paste("expectile = ", abs(round(var1u[which(var1u ==
                  vsort[imax])], 2)))
                text(0.15, 0.9, msg, cex = 0.8)
            }
            else 
            {
                msg <- paste("f = ", abs(round(F[imax - 1], 2)))
                text(0.8, 0.12, msg, cex = 0.8)
                msg <- paste("g = ", abs(round(G[imax - 1], 2)))
                text(0.8, 0.08, msg, cex = 0.8)
                msg <- paste("expectile = ", abs(round(var1u[which(var1u ==
                  vsort[imax])], 2)))
                text(0.8, 0.04, msg, cex = 0.8)
            }
        }
    }
    
    if (graph == "VLorentz") 
    {
        G <- G[cumsum(as.data.frame(table(F))$Freq)]
        F <- F[cumsum(as.data.frame(table(F))$Freq)]
        var1u = unique(var1)
        var1u <- c(0, var1u)
        var1u <- sort(var1u)
        plot(F, G, pch=symbol, col = couleurs,xlab=labvar[1],ylab=labvar[2])
        segments(0, 0, 1, 1, col = "blue")
        lines(F, G, col = couleurs)
        if ((length(var1[obs]) != 0) && (length(var1[obs])!=length(var1))) 
        {
            imax = which.min(var1u <= Xpoly)
            segments(F[imax - 1], 0, F[imax - 1], G[imax - 1],
                col = "red")
            segments(F[imax - 1], G[imax - 1], 0, G[imax - 1],
                col = "red")
            if (F[2] > G[2]) 
            {
                msg <- paste("f = ", round(F[imax - 1], 2))
                text(0.15, 0.98, msg, cex = 0.8)
                msg <- paste("g = ", round(G[imax - 1], 2))
                text(0.15, 0.94, msg, cex = 0.8)
                msg <- paste("expectile = ", Xpoly)
                text(0.15, 0.9, msg, cex = 0.8)
            }
            else 
            {
                msg <- paste("f = ", round(F[imax - 1], 2))
                text(0.8, 0.12, msg, cex = 0.8)
                msg <- paste("g = ", round(G[imax - 1], 2))
                text(0.8, 0.08, msg, cex = 0.8)
                msg <- paste("expectile = ", Xpoly)
                text(0.8, 0.04, msg, cex = 0.8)
            }
        }
        
        if(length(var1[obs])==length(var1))
        {segments(F[length(F)], 0, F[length(F)], G[length(F)],
          col = "red")
         segments(F[length(F)], G[length(F)], 0, G[length(F)],
         col = "red")
        
         msg <- paste("f = ", round(F[length(F)], 2))
         text(0.8, 0.12, msg, cex = 0.8)
         msg <- paste("g = ", round(G[length(F)], 2))
         text(0.8, 0.08, msg, cex = 0.8)
         msg <- paste("expectile = ", Xpoly)
         text(0.8, 0.04, msg, cex = 0.8)
                
        }        
    }

    
    if (graph == "Moran") 
    {
       plot(var1, var2, "n", xlab = labvar[1], ylab = labvar[2])
        #,xlim=c(-12,12),ylim=c(-4.5,4.5)
      segments(min(var1), mean(var2), max(var1), mean(var2), col = "black", lty=2)
      segments(mean(var1), min(var2), mean(var1), max(var2), col = "black", lty=2)
      points(var1[!obs], var2[!obs], col = couleurs[obsq[!obs]], pch = symbol[obsq[!obs]],
      cex=cbuble[!obs])
     

        if(!is.null(bin))
        {if(bin) abline(lm(var2 ~ var1)) }

        
        if (length(var1[obs]) != 0) 
        {
         points(var1[obs], var2[obs], col = "red", pch = symbol[obsq[obs]], cex = cbuble[obs])
          
         if (locmoran) 
         {
          ilocal <- var3[obs]
          msg <- paste(round(ilocal, 2))
          text(var1[obs] + 0.02, var2[obs] + 0.02, msg, cex = cex.lab,font=3,adj=c(0.75,-0.75));  
         }
        }
  
       if(legends[[1]])
       {
        legend(legends[[3]]$x,legends[[3]]$y, c(legmap[4],legmap[5],legmap[6]),
        title=legmap[7],pch = 16,cex=cex.lab, 
        pt.cex=c(as.numeric(legmap[1]),as.numeric(legmap[2]),as.numeric(legmap[3])))
       }  
        
    }
      
    if (graph == "Quadrant") 
    {
        plot(var1, var2, "n", xlab = labvar[1], ylab = labvar[2])
        segments(min(var1), mean(var2), max(var1), mean(var2), col = "black",lty=2)
        segments(mean(var1), min(var2), mean(var1), max(var2), col = "black",lty=2)
        points(var1[!obs], var2[!obs], col = couleurs[obsq[!obs]], pch = symbol[obsq[!obs]],
        cex=cbuble[!obs])

        if(!is.null(bin))
         {if(bin) abline(lm(var2 ~ var1))}
        
        if (length(var1[obs]) != 0) 
        {
         points(var1[obs], var2[obs], col = "red", pch = symbol[obsq[obs]],
         cex=cbuble[obs])
        }
        
       if(legends[[1]])
       {
        legend(legends[[3]]$x,legends[[3]]$y, c(legmap[4],legmap[5],legmap[6]),
        title=legmap[7],pch = 16,cex=cex.lab, 
        pt.cex=c(as.numeric(legmap[1]),as.numeric(legmap[2]),as.numeric(legmap[3])))
       }  
        
    }


    if (graph == "Neighbourplot") 
    {
        plot(var1, var1, "n", xlab = labvar[1], ylab = labvar[2])
        if (opt1) lines(c(0, max(var1)), c(0, max(var1)))
        
        indt <- which(obs == TRUE, arr.ind = TRUE)
        indf <- which(obs == FALSE, arr.ind = TRUE)
        if (length(indf) > 0) 
        {
            for (i in 1:nrow(indf)) 
            {
                if (W[indf[i, 1], indf[i, 2]] != 0) 
                {
                  points(var1[indf[i, 1]], var1[indf[i, 2]],
                    col = couleurs[1], pch = 16, cex=0.8)
                }
            }
        }
        if (length(var1[obs]) != 0) 
        {
           for (i in 1:nrow(indt))
            {
              if (W[indt[i, 1], indt[i, 2]] != 0) 
              {
               points(var1[indt[i, 1]], var1[indt[i, 2]],
               col = "red", pch = symbol[1], cex = 1)
              }
            }
        }
    }
    

   if (graph == "Acp2")
    {
        centre <- rep(0, length(var1))
        msg1 <- paste("component ", labvar[1],
        " : ", round(inertie[direct[1]], 0), "%", sep = "")
        msg2 <- paste("component ", labvar[2],
        " : ", round(inertie[direct[2]], 0), "%", sep = "")

        if ((max(abs(var1)) > 1) || (max(abs(var2)) > 1))
        {
            plot(var1, var2, "n", xlab = msg1, ylab = msg2)
        }
        else
        {
            plot(-1:1, -1:1, "n", xlab = msg1, ylab = msg2)
            segments(-1, 0, 1, 0, col = "black")
            segments(0, -1, 0, 1, col = "black")
        }
        segments(centre, centre, var1, var2, col = couleurs)
        msg <- paste(legmap, " : ", round(labmod, 0), "%", sep = "")
        text(var1 +0.02, var2 + runif(length(var2),-0.02,0.02), msg,
        cex=cex.lab,font=3,adj=c(0.75,-0.75))
        x = seq(-1, 1, 0.02)
        y = sqrt(1 - x^2)
        lines(x, y)
        lines(x, -y)
    }


}
