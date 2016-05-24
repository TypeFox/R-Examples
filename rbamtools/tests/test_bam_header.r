
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# test_bam_header: General
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
header<-getHeader(reader)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Section: headerReadGroup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
## Create headerReadGroup by simply parsing text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

hrg <- c("@RG\tID:rg1\tCN:seqCenter1\tDT:01.01.2011\tSM:sm1",
                    "@RG\tID:rg2\tCN:seqCenter2\tDT:01.01.2012\tSM:sm2")
object <- new("headerReadGroup", hrg)

if(!identical(object@ID,paste("rg", 1:2, sep="")))
    stop("[test_bam_header.r] Setting ID slot in headerReadGroup failed")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Create headerReadGroup by user interface
# and check for consistency of returned header text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

hrg <- new("headerReadGroup")
hrg <- addReadGroup(hrg, list(ID="rg1", CN="sct1", FO="fo1"))
hrg <- addReadGroup(hrg, list(ID="rg2", CN="sct2", FO="fo2", LB="lb2"))
hrg <- addReadGroup(hrg, list(ID="rg3", CN="sct3", LB="lb3"))
hrg2 <- new("headerReadGroup", getHeaderText(hrg))

if(!identical(hrg, hrg2))
    stop("[test_bam_header.r] Adding headerReadGroup's failed")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Section: headerProgram
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

prog <- new("headerProgram")
setVal(prog, "ID", "pg.001")
setVal(prog, "PN", "Program name")
setVal(prog, "CL", "Command line")
setVal(prog, "DS", "Description")

if(!identical(getVal(prog,"ID"), "pg.001"))
    stop("[test_bam_header.r] Setting ID in headerProgram failed")

if(!identical(getVal(prog,"DS"), "Description"))
    stop("[test_bam_header.r] Setting Description in header Program failed")

# Set new Description segment
header<-getHeader(reader)
htxt<-getHeaderText(header)
prog <-headerProgram(htxt)

setVal(prog, "DS", "Description")
headerProgram(htxt) <- prog

## Convert to binary header representation and retrieve Description segment
bh <- bamHeader(htxt)
if(!identical(getVal(headerProgram(getHeaderText(bh)),"DS"),"Description"))
    stop("[test_bam_header.r] Writing Program description to bamHeader failed")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# END OF FILE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
