

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# test_bam_align: General
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
align <- getNextAlign(reader)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Test for accessing new align FLAG (0x800)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

f <- flag(align)
suppAlign(align)
suppAlign(align) <- TRUE
if(!suppAlign(align))
    stop("[test_bam_align.r] Setting suppAlign failed")

suppAlign(align) <- FALSE
if(flag(align)!=f)
    stop("[test_bam_align.r]
            Setting and re-setting suppAlign changed other flags")
