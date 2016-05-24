AggregateFireFracs <-
   function(vegChangesLocal, fireAreaFracsLocal, vt2pvtlutLocal, pvts2aggregate) {
      years = vegChangesLocal[[2]]
      nYrs = length(years)
      vts = vegChangesLocal[[3]]
      nVTs = length(vts)
      
      ndxsOfVTs2aggregate = c()
      nPVTs = length(pvts2aggregate)
      pvts2aggregateNdx = 1
      while (pvts2aggregateNdx<=nPVTs) {
         tgtPVT = pvts2aggregate[[pvts2aggregateNdx]]
         found = FALSE
         vtNdx = 0
         while (vtNdx<nVTs && !found) {
            vtNdx = vtNdx + 1
            vt = vt2pvtlutLocal$VT[vtNdx]
            pvt = levels(vt2pvtlutLocal$PVT)[vt2pvtlutLocal$PVT[vtNdx]]
            cat(c("tgtPVT, vtNdx, vt, pvt = ", tgtPVT, vtNdx, vt, pvt, "\n"))
            
            if (pvt==tgtPVT) {
               ndxsOfVTs2aggregate = c(ndxsOfVTs2aggregate, vtNdx)
               pvts2aggregateNdx = pvts2aggregateNdx + 1
               found = TRUE
               cat(c("found vt ", vt, " at vtNdx ", vtNdx, " for tgtPVT ", tgtPVT, "\n"))
            }
         }
         stopifnot(found)
      }
      stopifnot(length(ndxsOfVTs2aggregate)==nPVTs)
      
      aggVTfracs = array(0, nYrs)
      aggFireAreaFracs = array(0, nYrs)
      vtFracs = vegChangesLocal[[4]]
      stopifnot(dim(vtFracs)[1]==nVTs)
      stopifnot(dim(vtFracs)[2]==nYrs)
      for (aggNdx in 1:nPVTs) {
         vtNdx = ndxsOfVTs2aggregate[aggNdx]
         vt = vt2pvtlutLocal$VT[vtNdx]
         cat(c("vtNdx, vt = ", vtNdx, vt, "\n"))
         aggVTfracs = aggVTfracs + vtFracs[vtNdx,]
         aggFireAreaFracs = aggFireAreaFracs + fireAreaFracsLocal[,vtNdx]*vtFracs[vtNdx,]
      }
      
      aggFireFracs = array(0, nYrs)
      for (yrNdx in 1:nYrs) {
         if (aggVTfracs[yrNdx]>0) aggFireFracs[yrNdx] = aggFireAreaFracs[yrNdx]/aggVTfracs[yrNdx]
         else aggFireFracs[yrNdx] = NA
      }
      
      return(aggFireFracs)
   }