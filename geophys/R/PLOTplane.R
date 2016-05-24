PLOTplane <-
function(Rp, planecol="brown")
      {
        if(missing(planecol)) {  planecol="brown"  }
        points(Rp[,1], Rp[,2], pch=21, col='red', bg='yellow')

        segments(Rp[1,1],Rp[1,2], Rp[2,1],Rp[2,2], col=planecol)
        segments(Rp[2,1],Rp[2,2], Rp[3,1],Rp[3,2], col=planecol)
        segments(Rp[3,1],Rp[3,2], Rp[1,1],Rp[1,2], col=planecol)
      }

