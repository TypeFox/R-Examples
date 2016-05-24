`getGEOmap` <-
function(fn)
{
  ####   fn = path to 2 files with root names =fn, the suffexes are .strks and .pnts
   ###   get a leesmap from an ascii file.
  ###  the ascii file is gernated from a leesmap using  driveleesmap input.leesmap

## proposed map code scheme
## UW map codes for plotting

## a = major coasts  islands  lakes
## b = intermediate coasts islands lakes
## c = minor coasts islands lakes
## d =  very minor coasts islands lakes
## e =  major rivers
## f =  intermediate rivers
## g =  minor rivers
## h =  very minor rivers
## i = political borders
## j = major faults
## k = minor faults
## l = geology formation
## m = major high ways
## n = minor roads
## p = plates
## o = other
 
  fstks = paste(sep='', fn, ".strks")
  fpnts = paste(sep='', fn, ".pnts")

  stks =  scan(file=fstks, fill=TRUE, list(nam='', num=0, index=0, col=0, style=0, code=''), quiet=TRUE)
  pnts=scan(file=fpnts, list(lat=0, lon=0), quiet=TRUE)

  stks$code[stks$code==''] = "o"

  return(list(STROKES=stks, POINTS=pnts))

}

