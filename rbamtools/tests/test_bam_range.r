
coords<-getRefCoords(reader,"chr1")
rg<-bamRange(reader,coords)


if(size(rg)!=2216)
    stop("[test_bam_range.r] size of range must be 2216!")
