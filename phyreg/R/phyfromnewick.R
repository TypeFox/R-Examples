phyfromnewick <-
function(file="",str,datatype="Branch_lengths") {  # reads a newick file or string, and converts to internal format
                                                                                                                      # doesn't deal with names of higher nodes (see Wikipedia)
                                                                                                                      # calls its own internal copy of uw to cut out single-daughter higher nodes
                                                                                                                      # essentially ignores all quote marks. Could do better!
###### Its own copy of uw()  -- NO LONGER. I can now have it at top level and NOT have build / check complain about documentation

######### end of definition of uw -- beginning of body of phyfromnewick

if (!missing(file)) str<-scan(file,sep="\n",what="character")

qchr<-'"'
obchr<-"("
cbchr<-")"
comchr<-","
colchr<-":"
scchr<-";"
sqchr<-"'"

symbset<-c(obchr,cbchr,comchr,colchr,scchr)


emptystr<-""
Heightsstr<-"Heights"

str<-strsplit(str,split=NULL)[[1]]
str<-str[str!=" " & str!=qchr & str!=sqchr]  ### new This renders the outer if the the outer loop unnecessary as there none of thsoe characters

fop<-paste0("Failure message(s) from phylofromnewick:\n")

cs<-c("");
nob<-0;nco<-0;ncb<-0;curexob<-0;
jj<-1
for (ii in 1:length(str)) {
 chr<-str[[ii]];
 if (chr=="(") { nob<-nob+1;curexob<-curexob+1} else
    if (chr==")") { ncb<-ncb+1; curexob<-curexob-1 } else 
    if (chr==",")  {nco<-nco+1} 
 if (curexob<0) {
    fop<-paste(fop, "Bracketing is wrong in input string,")
     fop<-paste(fop, "after ", nob, " open brackets, ", ncb, " close brackets, and ", nco, "commas.")
     cat(fop); stop("leaving phylofromnewick")
  }
 if (chr!=" " & chr!="\n") {cs[jj]<-chr; jj<-jj+1 }
}
 if (curexob!=0) {
    fop<-paste0(fop, "Bracketing is wrong at end of input string,")
    fop<-paste0(fop, "after ", nob, " open brackets, ",  ncb, " close brackets, and ",  nco, "commas.")
    cat(fop);stop("leaving phylofromnewick")
}

parent<-rep(-1,nob)
final<-rep(-1,nob)
finalparent<-rep(-1,nob)
finalparentname<-paste0("Higher_Node_", as.character(seq(1:(nob)))) 
specparent<-rep(-1,1+nco)
specname<-paste0("Species_", as.character(seq(1:(1+nco))))  # 
finalname<-rep("specdata",nob) ### new
specdata<-rep("nodata",1+nco)
finaldata<-rep("nodata", nob)

getstuff<-function(i) {  # note uses variables cs, symbset, from environment
                                      # i is the first position where the stuff, if any, actually starts.
                                      # it can be called for i==end+1 and behaves well
  if (i>=(2+length(cs))) stop("getstuff alarmed to be asked to read past the end of the string")
  if (i==(1+length(cs))) return(list(stuff="nodata",nexti=i))
  k<-0
  while (!(cs[[i+k]] %in% symbset) & !((i+k)==length(cs))) k<-k+1;
  if (!(cs[[i+k]] %in% symbset) & ((i+k)==length(cs))) k<-k+1
  if (k>0L) stuff<-paste0(cs[i:(i+k-1)],collapse="") else stuff<-"nodata"
  return(list(stuff=stuff,nexti=(i+k)))
 }
 
 ## After , we have illegal ; and legal [ new species , A ) : & no action ( ]
 ## After ( we have illegal ; ) and legal [ new higher node and new species  , : A new higher node only ( ]
 ## After ) we have illegal ( and legal [ data for higher node : name of higher node A complete higher node ) , ; ]
 ## After : we have illegal ( : and legal [ data in dataslot A no data in dataslot ) , ; ]
 ## Of course, there is no after ;
 ## A is always absorbed up to the next syntactic symbol by getstuff
  
  
curhinode<-0; nexthinode<-1;curhinodefinal<-1;nextspecno<-1;
ii<-1;
repeat {                                   #  /* Do level 1 */

chr<-cs[[ii]]; if(chr==scchr) break
if (ii==length(cs)) nextchr<-scchr else nextchr<-cs[[ii+1]]

#print(paste("Start",ii,chr,nextchr))

if (chr==obchr) {                                   #      /* Do level 2 */
  section<-obchr
  parent[nexthinode]<-curhinode
  curhinode<-nexthinode
  nexthinode<-nexthinode+1
  if (nextchr!=obchr) {  ## Start level 3
    g<-getstuff(ii+1)
    specparent[nextspecno]<-curhinode;
    if (g$stuff!="nodata") specname[nextspecno]<-g$stuff
    nextspecno<-nextspecno+1;
    inspecies<-1;
    ii<-g$nexti} else          ## End level 3
    ii<-ii+1                       ## Start and end level 3
}    else                                                 #  /* End level 2 */
if (chr==cbchr) {                                      #   /* Do level 2 */
  section<-cbchr
  final[curhinode]<-curhinodefinal;
  finalparent[curhinodefinal]<-parent[curhinode];
  curhinode<-parent[curhinode];
  if (!(nextchr %in% symbset)) {
      g<-getstuff(ii+1)
      if (g$stuff!="nodata") finalparentname[curhinodefinal]<-g$stuff
      ii<-g$nexti-1; nextchr<-cs[[ii+1]]
    }
    if (nextchr==colchr) {
       g<-getstuff(ii+2)
       if (g$stuff!="nodata") finaldata[curhinodefinal]<-g$stuff  # -1 of long duration removed from index -- fingers crossed!
       ii<-g$nexti-1  # So we can do an unconditional ii<-ii+1 after the bracket
      }
      curhinodefinal<-curhinodefinal+1;
      ii<-ii+1
     inspecies<-0;
 }     else                                                    #/* End level 2 */
   if (chr==comchr) {                          #  /* Do level 2 */
  section<-comchr
  if (nextchr!=obchr) {                                    # Start level 3
     g<-getstuff(ii+1)
    specparent[nextspecno]<-curhinode;
    if (g$stuff!="nodata") specname[nextspecno]<-g$stuff
    nextspecno<-nextspecno+1;
    inspecies<-1;
     ii<-g$nexti
   }   else                                                 # End level 3
   ii<-ii+1                                         # Start and end of level 3
 }     else                                            #  /* End level 2 */
if (chr==colchr) {                                    #      /* Do level 2 */
  section<-colchr
  g<-getstuff(ii+1)
  name<-paste0(g$stuff,collapse="")
  if (g$stuff!="nodata") specdata[nextspecno-1]<-g$stuff
  ii<-g$nexti
 }       # now there *is* no "else"           #  /* End level 2 */

if (ii>length(cs)) break  ## replaces the SAS do until. Note ii is changed within the loop
                                         ## but in fact a string that ends properly (;) should not end here
                                         
#print(paste("End",ii,chr,nextchr))
#if (exists("g")) print(paste(g$stuff,g$nexti))

}  #  /* End level 1 */

prephy<-c(specparent, head(finalparent,-1))

prephy2 <-final[prephy];

phy<-1+nco+prephy2;

## Remember to return phy!

outspecies<-cbind(Name=specname , Number=seq(1,length(specdata)), AddData=specdata)

#  Remember to return outspecies


outhigher<-cbind(Name=finalparentname,Number=c((length(specdata)+1):length(phy), max(phy)), AddData=finaldata);

# Remember to return outhigher

  noheights<-FALSE
  if (datatype=="Heights")   {
         warnlevel<-options("warn"=-1)
         outheights<-as.numeric(specdata,finaldata)  
         options(warn=warnlevel$warn)
         if (any(is.na(outheights))) noheights<-TRUE
         }  else {
         warnlevel<-options("warn"=-1)
         tempdata<-as.numeric(c(specdata,head(finaldata,-1)));
         options(warn=warnlevel$warn)
         if (any(is.na(tempdata))) noheights<-TRUE else
           {   tempdata<-c(tempdata,0);
               hts<-rep(0,length(phy)+1);
               for (ii in length(phy):1) hts[ii]<-tempdata[ii]+hts[phy[ii]]
               hts<-max(hts)-hts;
               outheights<-hts;
            }
          }
  
  myuw<-uw(phy,rep(1,length(specdata)),length(specdata))
  if (!noheights) newheights<-outheights[myuw$on]
  
  if (noheights) return(list(phy=myuw$txp, firstphy=phy, sd=outspecies,oh=outhigher)) else
      return(list(phy=myuw$txp, hts=newheights,firstphy=phy, sd=outspecies,oh=outhigher))

}
