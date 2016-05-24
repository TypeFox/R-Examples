PFoutput <-
function(PF, stas=NULL, sol=NULL, format=0 )
  {
###  save pickfile information to disk in a variety of formats

    if(missing(format)) format = 0
    ###   format = 0,1,2,3,4
    ###   format <= 0 means save in all formats

    #### sol = solution to earthquake location problem (lat,lon, z, t0)
    ###   stas = station location list

    

    if(any( 0 %in% format )  )
      {
        twpx = PF
        twpx$phase[twpx$phase=="Y"] = "P"
        PF =  INITpickfile(stas=stas, src=NULL, WPX=twpx)
      }

    
    if(!is.null(sol))
      {
        upf = UPdateEQLOC(PF,  sol, stas=stas)
      }
    else
      {
        upf = PF

      }

    dout = c(upf$LOC$yr, upf$LOC$jd, upf$LOC$hr, upf$LOC$mi,  upf$LOC$sec)
    fout1 = PCfiledatetime(dout, 0)
    foutp = paste(fout1,"uwpf", sep="." )
    upf$filename =   foutp
    output=upf$filename
    
   foutcvs = paste(fout1,"wpx", sep="." )
    fout2 = paste(fout1,"RDATA", sep="." )

    if(any( c(1,0) %in% format ) )   save(file=fout2, twpx)
    if(any( c(2,0) %in% format ) )   RSEIS::writeUWpickfile(upf, output=output)
    if(any( c(3,0) %in% format ) )
      {
        write.csv(twpx, file =foutcvs, row.names = FALSE)
      }
  }
