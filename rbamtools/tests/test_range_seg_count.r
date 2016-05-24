

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Count alignments in specified genomic segments
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
coords <- c(0, 0, 2e4)
segments <- seq(14000, 20000, 20)
segcount<-rangeSegCount(reader, coords, segments)
segcount
dfr<-as.data.frame(segcount)

if(nrow(dfr)!=301)
    stop("[test_range_seg_count.r] data.frame has wrong number of lines!")

if(sum(dfr$count)!=2112)
    stop("[test_range_seg_count.r] Wrong total number of aligns!")
