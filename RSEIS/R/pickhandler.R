pickhandler<-function(i1=1, ppick=0, kzap="Y", err=NA, res=0, ycol=rgb(0,0,1) , pol=0, flg=0, onoff=1, NPX=1, WPX=WPX, NH=NH)
  {
#######   used in swig  for handling picks
    if(missing(i1)) i1=1
    if(missing(kzap)) kzap="Y"
    if(missing(err)) err=NA
    if(missing(res)) res=NA
    if(missing(ycol)) ycol=rgb(0,0,1)
    if(missing(pol)) pol= 0
    if(missing(flg)) flg= 0
    if(missing(onoff)) onoff=1
    if(missing(NPX)) NPX=1

    if(missing(NH))
      {
        print("Missing RSEIS (NH) list in pickhandler: need datetime")
        return(NULL)
      }
    
    if(missing(WPX))
      {
        WPX = cleanWPX()
      }
    
    ## print(i1)
    
    asec = NH$info$sec[i1]+NH$info$msec[i1]/1000+NH$info$t1[i1]-NH$info$off[i1]+ppick
    pic1 = recdate(NH$info$jd[i1], NH$info$hr[i1], NH$info$mi[i1], asec)
    
    WPX$tag[NPX]=paste(sep=".",NH$STNS[i1],  NH$COMPS[i1])
    WPX$name[NPX]=NH$STNS[i1]
    WPX$comp[NPX]=NH$COMPS[i1]
    WPX$c3[NPX]=NH$OCOMPS[i1]
    WPX$phase[NPX]=kzap
    
    WPX$err[NPX]=err
    WPX$pol[NPX]=pol
    WPX$flg[NPX]=flg
    WPX$res[NPX]= res
    WPX$dur[NPX]= res

    WPX$yr[NPX]=NH$info$yr[i1]
    WPX$mo[NPX]= NH$info$mo[i1]
    WPX$dom[NPX]=NH$info$dom[i1]
    WPX$jd[NPX]=pic1$jd
    WPX$hr[NPX]=pic1$hr
    WPX$mi[NPX]=pic1$mi
    WPX$sec[NPX]=pic1$sec
    WPX$col[NPX]=ycol
    WPX$onoff[NPX] = onoff

    return(WPX)
    
  }
