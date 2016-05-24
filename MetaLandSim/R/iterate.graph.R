iterate.graph <-
  function (iter, mapsize, dist_m, areaM, areaSD, Npatch,
            disp, span, par1="none", par2=NULL, par3=NULL, 
            par4=NULL, par5=NULL, method="percentage",
            parm, nsew="none", succ="none", param_df, kern, conn,
            colnz, ext, beta1=NULL, b=1, c1=NULL, c2=NULL,
            z=NULL, R=NULL, graph=TRUE)
  {
    ma <- matrix(nrow = span, ncol = iter)
    md <- matrix(nrow = span, ncol = iter)
    np <- matrix(nrow = span, ncol = iter)
    occ <- matrix(nrow =span, ncol = iter)
    trn <- matrix(nrow =span, ncol = iter)
    for(i in 1:iter)
    {
      rland1 <- rland.graph(mapsize, dist_m, areaM, areaSD, Npatch,
                            disp, plotG=FALSE)
      span1 <- span.graph(rland1, span, par1, par2, par3, par4, par5)
      span2 <- length(span1)
      sim <- simulate_graph(rland1, span1, simulate.start=TRUE, method, parm,
                            nsew, succ, param_df, kern, conn, colnz, ext, beta1,
                            b, c1, c2, z, R)
      marea <- list.stats(sim, "mean_area", plotG=FALSE)
      ma[1:length(marea), i] <- marea
      mdistance <- list.stats(sim, "mean_nneigh", plotG=FALSE)
      md[1:length(mdistance), i] <- mdistance
      numberpatches <- list.stats(sim, "n_patches", plotG=FALSE)
      np[1:length(numberpatches), i] <- numberpatches
      occupation1 <- list.stats(sim, "occupation", plotG=FALSE)
      occ[1:length(occupation1), i] <- occupation1
      turnover1 <- list.stats(sim, "turnover", plotG=FALSE)
      trn[1:length(turnover1), i] <- turnover1
      
      cat("Completed iteration",i," of ",iter,"\n")
      
    }
    ma <- as.data.frame (ma)
    ma[, iter+1] <- rowMeans(ma[, 1:iter])
    ma[, iter+2] <- apply(as.matrix(ma[, 1:iter]), 1, sd)
    ma <- na.omit(ma)
    md <- as.data.frame (md)
    md[, iter+1] <- rowMeans(md[, 1:iter])
    md[,iter+2] <- apply(as.matrix(md[, 1:iter]), 1, sd)
    md <- na.omit(md)
    np <- as.data.frame (np)
    np[, iter+1] <- rowMeans(np[, 1:iter])
    np[,iter+2] <- apply(as.matrix(np[,1:iter]),1, sd)
    np <- na.omit(np)
    occ <- as.data.frame (occ)
    occ[, iter+1] <- rowMeans(occ[, 1:iter])
    occ[, iter+2] <- apply(as.matrix(occ[, 1:iter]), 1, sd)
    occ <- na.omit(occ)
    trn <- as.data.frame (trn)
    trn[, iter+1] <- rowMeans(trn[, 1:iter])
    trn[, iter+2] <- apply(as.matrix(trn[, 1:iter]), 1, sd)
    trn <- na.omit(trn)
    if(graph == TRUE)
    {
      time_vector <- 1:span
      g_area <- gvisLineChart(data.frame(time_step=time_vector, mean_area=ma[, iter+1]),
                              xvar="time_step", yvar="mean_area",
                              options=list(title="Mean area (Ha)", width=600, height=300,
                                           curveType="function", legend="none",
                                           titleTextStyle="{colour:'black', fontName:'Courier', fontSize:16}",
                                           vAxis="{title: 'hectares'}", hAxis="{title: 'time steps'}",
                                           series="[{color: '#006400'}]", backgroundColor="#D1EEEE")) 
        g_dist <- gvisLineChart(data.frame(time_step=time_vector, mean_distance=md[, iter+1]),
                                xvar="time_step", yvar="mean_distance",
                                options=list(title="Mean distance to nearest habitat patch (m)", width=600, height=300,
                                curveType="function", legend="none",
                                titleTextStyle="{colour:'black', fontName:'Courier', fontSize:16}",
                                vAxis="{title: 'meters'}", hAxis="{title: 'time steps'}",
                                series="[{color:'#0000FF'}]", backgroundColor="#D1EEEE"))
      g_npatches <- gvisLineChart(data.frame(time_step=time_vector, npatches=np[, iter+1]),
                                  xvar="time_step", yvar="npatches",
                                  options=list(title="Number of patches", width=600, height=300,
                                               curveType="function", legend="none", 
                                               titleTextStyle="{colour:'black', fontName:'Courier', fontSize:16}",
                                               vAxis="{title: 'number of patches'}", hAxis="{title: 'time steps'}",
                                               series="[{color:'#8B0000'}]", backgroundColor="#D1EEEE"))
      g_occ <- gvisLineChart(data.frame(time_step=time_vector, occ=occ[, iter+1]),
                             xvar="time_step", yvar="occ",
                             options=list(title="Species patch occupancy (%)", width=600, height=300,
                                          curveType="function", legend="none",
                                          titleTextStyle="{colour:'black', fontName:'Courier', fontSize:16}",
                                          vAxis="{title: '% of patch occupancy'}", hAxis="{title: 'time steps'}",
                                          series="[{color: '#FF4500'}]", backgroundColor="#D1EEEE"))
      g_trn <- gvisLineChart(data.frame(time_step=time_vector, trn=trn[, iter+1]),
                             xvar="time_step", yvar="trn",
                             options=list(title="Occupancy turnover (%)", width=600,height=300,
                                          curveType="function", legend="none",
                                          titleTextStyle="{colour:'black', fontName:'Courier', fontSize:16}",
                                          vAxis="{title: '% of turnover'}", hAxis="{title: 'time steps'}",
                                          series="[{color:'#8B4789'}]", backgroundColor="#D1EEEE"))
      ln1 <- gvisMerge(g_area, g_dist, horizontal=TRUE)
      ln2 <- gvisMerge(g_npatches, g_occ, horizontal=TRUE)
      ln.final <- gvisMerge(gvisMerge(ln1, ln2, horizontal=FALSE), g_trn, horizontal=FALSE)
      paste("<div><span>Metapopulation persistence in a dynamic landscape (parameter 1 = ", par1,"; see help for further description)</span><br />", sep="") -> ln.final$html$caption
      paste("\n<!-- htmlFooter -->\n<span> \n ",R.Version()$version.string,"&#8226; <a href=\"http://code.google.com/p/google-motion-charts-with-r/\">googleVis-",packageVersion("googleVis"),"</a>\n &#8226; MetaLandSim-",packageVersion("MetaLandSim"),"\n &#8226; <a href=\"https://developers.google.com/terms/\">Google Terms of Use</a> &#8226; <a href=\"https://google-developers.appspot.com/chart/interactive/docs/gallery/linechart.html#Data_Policy\">Data Policy</a>\n</span></div>\n</body>\n</html>\n",sep="") -> ln.final$html$footer
      plot(ln.final)
    }
    iter_vector <- rep("iter",length=iter)
    number_vector <- seq(1:iter)
    names_vector <- paste(iter_vector,number_vector,sep="")
    names(ma)[1:iter] <- names_vector
    names(ma)[iter+1] <- "mean"
    names(ma)[iter+2] <- "SD"
    names(md)[1:iter] <- names_vector
    names(md)[iter+1] <- "mean"
    names(md)[iter+2] <- "SD"
    names(np)[1:iter] <- names_vector
    names(np)[iter+1] <- "mean"
    names(np)[iter+2] <- "SD"
    names(occ)[1:iter] <- names_vector
    names(occ)[iter+1] <- "mean"
    names(occ)[iter+2] <- "SD"
    names(trn)[1:iter] <- names_vector
    names(trn)[iter+1] <- "mean"
    names(trn)[iter+2] <- "SD"
    output <- list(mean_area=ma, mean_distance=md, number_patches=np,
                   occupancy=occ, turnover=trn)
    return(output)
  }
