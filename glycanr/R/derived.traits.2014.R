# Derived traits for Glycan peaks in IgG for UPLC
# based on paper from 2014.
#
# @references
# Jennifer E. Huffman et al. 
# "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
# \url{http://dx.doi.org/10.1074/mcp.M113.037465}
igg.uplc.derived.traits.2014 <- function(data, print.exp.names=FALSE) {
    if(print.exp.names){
        return(paste0("GP", 1:24))
    }
    
    # derived glycans
    dteh(data$IGP24 <- with(data, (GP16 + GP18 + GP23)/(GP16 + GP18 + GP23 + GP8 + GP9 + GP14)) * 100)
    dteh(data$IGP25 <- with(data, (GP19 + GP24)/(GP19 + GP24 + GP10 + GP11 + GP15)) * 100)
    dteh(data$IGP26 <- with(data, (GP16 + GP18 + GP23)/(GP16 + GP18 + GP23 + GP4 + GP8 + GP9 + GP14)) * 100)
    dteh(data$IGP27 <- with(data, (GP19 + GP24)/(GP19 + GP24 + GP6 + GP10 + GP11 + GP15)) * 100)
    dteh(data$IGP28 <- with(data, GP16/(GP16 + GP8 + GP9)) * 100)
    dteh(data$IGP29 <- with(data, GP18/(GP18 + GP14 + GP23)) * 100)
    dteh(data$IGP30 <- with(data, GP23/(GP23 + GP14 + GP18)) * 100)
    dteh(data$IGP31 <- with(data, GP19/(GP19 + GP15 + GP24)) * 100)
    dteh(data$IGP32 <- with(data, GP24/(GP24 + GP15 + GP19)) * 100)
    dteh(data$IGP33 <- with(data, (GP16 + GP18 + GP19)/(GP23 + GP24)))
    dteh(data$IGP34 <- with(data, (GP16 + GP18)/GP23))
    dteh(data$IGP35 <- with(data, GP19/GP24))
    dteh(data$IGP36 <- with(data, (GP19 + GP24)/(GP16 + GP18 + GP23)))
    dteh(data$IGP37 <- with(data, GP19/(GP16 + GP18)))
    dteh(data$IGP38 <- with(data, GP19/(GP16 + GP18 + GP19)) * 100)
    dteh(data$IGP39 <- with(data, GP24/GP23))
    dteh(data$IGP40 <- with(data, GP24/(GP23 + GP24)) * 100)
    
    # neutral glycans
    dteh(GPn <- with(data, GP1 + GP2 + GP4 + GP6 + GP7 + GP8 + GP9 + GP10 + GP11 + GP12 + GP13 + GP14 + GP15))
    dteh(data$IGP41 <- with(data, GP1/GPn) * 100)
    dteh(data$IGP42 <- with(data, GP2/GPn) * 100)
    dteh(data$IGP43 <- with(data, GP4/GPn) * 100)
    dteh(data$IGP44 <- with(data, GP5/GPn) * 100)
    dteh(data$IGP45 <- with(data, GP6/GPn) * 100)
    dteh(data$IGP46 <- with(data, GP7/GPn) * 100)
    dteh(data$IGP47 <- with(data, GP8/GPn) * 100)
    dteh(data$IGP48 <- with(data, GP9/GPn) * 100)
    dteh(data$IGP49 <- with(data, GP10/GPn) * 100)
    dteh(data$IGP50 <- with(data, GP11/GPn) * 100)
    dteh(data$IGP51 <- with(data, GP12/GPn) * 100)
    dteh(data$IGP52 <- with(data, GP13/GPn) * 100)
    dteh(data$IGP53 <- with(data, GP14/GPn) * 100)
    dteh(data$IGP54 <- with(data, GP15/GPn) * 100)
    
    # neutral derived glycans
    dteh(data$IGP55 <- with(data, (IGP41 + IGP42 + IGP43 + IGP45)))
    dteh(data$IGP56 <- with(data, (IGP46 + IGP47 + IGP48 + IGP49 + IGP50)))
    dteh(data$IGP57 <- with(data, (IGP51 + IGP52 + IGP53 + IGP54)))
    dteh(data$IGP58 <- with(data, (IGP41 + IGP43 + IGP45 + IGP47 + IGP48 + IGP49 + IGP50 + IGP53 + IGP54)))
    dteh(data$IGP59 <- with(data, (IGP41 + IGP43 + IGP45)/IGP55) * 100)
    dteh(data$IGP60 <- with(data, (IGP47 + IGP48 + IGP49 + IGP50)/IGP56) * 100)
    dteh(data$IGP61 <- with(data, (IGP53 + IGP54)/IGP57) * 100)
    dteh(data$IGP62 <- with(data, (IGP41 + IGP43 + IGP47 + IGP48 + IGP53)))
    dteh(data$IGP63 <- with(data, (IGP41 + IGP43)/IGP55) * 100)
    dteh(data$IGP64 <- with(data, (IGP47 + IGP48)/IGP56) * 100)
    dteh(data$IGP65 <- with(data, IGP53/IGP57) * 100)
    dteh(data$IGP66 <- with(data, (IGP45 + IGP49 + IGP50 + IGP54)))
    dteh(data$IGP67 <- with(data, IGP45/IGP55) * 100)
    dteh(data$IGP68 <- with(data, (IGP49 + IGP50)/IGP56) * 100)
    dteh(data$IGP69 <- with(data, IGP54/IGP57) * 100)
    dteh(data$IGP70 <- with(data, IGP66/IGP62) * 100)
    dteh(data$IGP71 <- with(data, IGP66/IGP58) * 100)
    dteh(data$IGP72 <- with(data, IGP62/(IGP52 + IGP66)))
    dteh(data$IGP73 <- with(data, IGP52/(IGP62 + IGP66)) * 1000)
    dteh(data$IGP74 <- with(data, IGP54/IGP53))
    dteh(data$IGP75 <- with(data, IGP54/(IGP53 + IGP54)) * 100)
    dteh(data$IGP76 <- with(data, IGP53/(IGP52 + IGP54)))
    dteh(data$IGP77 <- with(data, IGP52/(IGP53 + IGP54)) * 1000)
    
    return(data)
}


# Derived traits for Glycan peaks in IgG for LCMS
# based on paper from 2014.
#
# @references
# Jennifer E. Huffman et al. 
# "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
# \url{http://dx.doi.org/10.1074/mcp.M113.037465}
igg.lcms.derived.traits.2014 <- function(data, print.exp.names=FALSE) {
    if(print.exp.names){
        return(c("IgG1_G0F","IgG1_G1F","IgG1_G2F","IgG1_G0FN","IgG1_G1FN",
               "IgG1_G2FN","IgG1_G1FS1","IgG1_G2FS1","IgG1_G1FNS1","IgG1_G2FNS1",
               "IgG1_G0","IgG1_G1","IgG1_G2","IgG1_G0N","IgG1_G1N",
               "IgG1_G2N","IgG1_G1S1","IgG1_G2S1","IgG1_G1NS1","IgG1_G2NS1",
               "IgG2_G0F","IgG2_G1F","IgG2_G2F","IgG2_G0FN","IgG2_G1FN",
               "IgG2_G2FN","IgG2_G1FS1","IgG2_G2FS1","IgG2_G1FNS1","IgG2_G2FNS1",
               "IgG2_G0","IgG2_G1","IgG2_G2","IgG2_G0N","IgG2_G1N",
               "IgG2_G2N","IgG2_G1S1","IgG2_G2S1","IgG2_G1NS1","IgG2_G2NS1",
               "IgG4_G0F","IgG4_G1F","IgG4_G2F","IgG4_G0FN","IgG4_G1FN",
               "IgG4_G2FN","IgG4_G1FS1","IgG4_G2FS1","IgG4_G1FNS1","IgG4_G2FNS1",
               "IgG4_G0","IgG4_G1","IgG4_G2","IgG4_G0N","IgG4_G1N",
               "IgG4_G2N","IgG4_G1S1","IgG4_G2S1","IgG4_G1NS1","IgG4_G2NS1"))
    }
    
    # =======================================
    # IgG1 derived traits
    # =======================================

    dteh(data$LC_IGP21 <- with(data, IgG1_G0F+IgG1_G1F+IgG1_G2F+IgG1_G0FN+IgG1_G1FN+IgG1_G2FN+IgG1_G1FS1+IgG1_G2FS1+IgG1_G1FNS1+IgG1_G2FNS1))
    dteh(data$LC_IGP22 <- with(data, IgG1_G0FN+IgG1_G1FN+IgG1_G2FN+IgG1_G1FNS1+IgG1_G2FNS1+IgG1_G0N+IgG1_G1N+IgG1_G2N+IgG1_G1NS1+IgG1_G2NS1))
    dteh(data$LC_IGP23 <- with(data, (IgG1_G1F+IgG1_G1FN+IgG1_G1FS1+IgG1_G1FNS1+IgG1_G1+IgG1_G1N+IgG1_G1S1+IgG1_G1NS1)*0.5+(IgG1_G2F+IgG1_G2FN+IgG1_G2FS1+IgG1_G2FNS1+IgG1_G2+IgG1_G2N+IgG1_G2S1+IgG1_G2NS1)))
    dteh(data$LC_IGP24 <- with(data, IgG1_G1FS1+IgG1_G2FS1+IgG1_G1FNS1+IgG1_G2FNS1+IgG1_G1S1+IgG1_G2S1+IgG1_G1NS1+IgG1_G2NS1))
    dteh(data$LC_IGP25 <- with(data, LC_IGP24/(2*LC_IGP23)*100))
    dteh(data$LC_IGP26 <- with(data, (IgG1_G1S1+IgG1_G2S1)/(IgG1_G1+IgG1_G1S1+IgG1_G2+IgG1_G2S1)*100))
    dteh(data$LC_IGP27 <- with(data, (IgG1_G1S1+IgG1_G2S1)/(IgG1_G0+IgG1_G1+IgG1_G1S1+IgG1_G2+IgG1_G2S1)*100))
    dteh(data$LC_IGP28 <- with(data, IgG1_G1S1/(IgG1_G1+IgG1_G1S1)*100))
    dteh(data$LC_IGP29 <- with(data, IgG1_G2S1/(IgG1_G2+IgG1_G2S1)*100))
    dteh(data$LC_IGP30 <- with(data, (IgG1_G1NS1+IgG1_G2NS1)/(IgG1_G1N+IgG1_G1NS1+IgG1_G2N+IgG1_G2NS1)*100))
    dteh(data$LC_IGP31 <- with(data, (IgG1_G1NS1+IgG1_G2NS1)/(IgG1_G0N+IgG1_G1N+IgG1_G1NS1+IgG1_G2N+IgG1_G2NS1)*100))
    dteh(data$LC_IGP32 <- with(data, IgG1_G1NS1/(IgG1_G1N+IgG1_G1NS1)*100))
    dteh(data$LC_IGP33 <- with(data, IgG1_G2NS1/(IgG1_G2N+IgG1_G2NS1)*100))
    dteh(data$LC_IGP34 <- with(data, (IgG1_G1FS1+IgG1_G2FS1)/(IgG1_G1F+IgG1_G1FS1+IgG1_G2F+IgG1_G2FS1)*100))
    dteh(data$LC_IGP35 <- with(data, (IgG1_G1FS1+IgG1_G2FS1)/(IgG1_G0F+IgG1_G1F+IgG1_G1FS1+IgG1_G2F+IgG1_G2FS1)*100))
    dteh(data$LC_IGP36 <- with(data, IgG1_G1FS1/(IgG1_G1F+IgG1_G1FS1)*100))
    dteh(data$LC_IGP37 <- with(data, IgG1_G2FS1/(IgG1_G2F+IgG1_G2FS1)*100))
    dteh(data$LC_IGP38 <- with(data, (IgG1_G1FNS1+IgG1_G2FNS1)/(IgG1_G1FN+IgG1_G1FNS1+IgG1_G2FN+IgG1_G2FNS1)*100))
    dteh(data$LC_IGP39 <- with(data, (IgG1_G1FNS1+IgG1_G2FNS1)/(IgG1_G0FN+IgG1_G1FN+IgG1_G1FNS1+IgG1_G2FN+IgG1_G2FNS1)*100))
    dteh(data$LC_IGP40 <- with(data, IgG1_G1FNS1/(IgG1_G1FN+IgG1_G1FNS1)*100))
    dteh(data$LC_IGP41 <- with(data, IgG1_G2FNS1/(IgG1_G2FN+IgG1_G2FNS1)*100))
    dteh(data$LC_IGP42 <- with(data, (IgG1_G1NS1+IgG1_G2NS1)/(IgG1_G1S1+IgG1_G2S1)))
    dteh(data$LC_IGP43 <- with(data, (IgG1_G1FNS1+IgG1_G2FNS1)/(IgG1_G1FS1+IgG1_G2FS1)))
    dteh(data$LC_IGP44 <- with(data, (IgG1_G1NS1+IgG1_G2NS1)/(IgG1_G1S1+IgG1_G1NS1+IgG1_G2S1+IgG1_G2NS1)))
    dteh(data$LC_IGP45 <- with(data, (IgG1_G1FNS1+IgG1_G2FNS1)/(IgG1_G1FS1+IgG1_G1FNS1+IgG1_G2FS1+IgG1_G2FNS1)))

    # neutral
    dteh(IgG1.neutral <- with(data, IgG1_G0F+IgG1_G1F+IgG1_G2F+IgG1_G0FN+IgG1_G1FN+IgG1_G2FN+IgG1_G0+IgG1_G1+IgG1_G2+IgG1_G0N+IgG1_G1N+IgG1_G2N),
                                 mess="Not all neutral glycans for IgG1 where found in the data frame")

    dteh(data$LC_IGP46 <- with(data, IgG1_G0F/IgG1.neutral) * 100)
    dteh(data$LC_IGP47 <- with(data, IgG1_G1F/IgG1.neutral) * 100)
    dteh(data$LC_IGP48 <- with(data, IgG1_G2F/IgG1.neutral) * 100)
    dteh(data$LC_IGP49 <- with(data, IgG1_G0FN/IgG1.neutral) * 100)
    dteh(data$LC_IGP50 <- with(data, IgG1_G1FN/IgG1.neutral) * 100)
    dteh(data$LC_IGP51 <- with(data, IgG1_G2FN/IgG1.neutral) * 100)
    dteh(data$LC_IGP52 <- with(data, IgG1_G0/IgG1.neutral) * 100)
    dteh(data$LC_IGP53 <- with(data, IgG1_G1/IgG1.neutral) * 100)
    dteh(data$LC_IGP54 <- with(data, IgG1_G2/IgG1.neutral) * 100)
    dteh(data$LC_IGP55 <- with(data, IgG1_G0N/IgG1.neutral) * 100)
    dteh(data$LC_IGP56 <- with(data, IgG1_G1N/IgG1.neutral) * 100)
    dteh(data$LC_IGP57 <- with(data, IgG1_G2N/IgG1.neutral) * 100)

    # neutral - derived traits
    dteh(data$LC_IGP58 <- with(data, (LC_IGP52+LC_IGP46+LC_IGP49+LC_IGP55)))
    dteh(data$LC_IGP59 <- with(data, (LC_IGP53+LC_IGP47+LC_IGP50+LC_IGP56)))
    dteh(data$LC_IGP60 <- with(data, (LC_IGP54+LC_IGP48+LC_IGP51+LC_IGP57)))
    dteh(data$LC_IGP61 <- with(data, (LC_IGP46+LC_IGP49+LC_IGP47+LC_IGP50+LC_IGP48+LC_IGP51)))
    dteh(data$LC_IGP62 <- with(data, (LC_IGP46+LC_IGP49)/LC_IGP58*100))
    dteh(data$LC_IGP63 <- with(data, (LC_IGP47+LC_IGP50)/LC_IGP59*100))
    dteh(data$LC_IGP64 <- with(data, (LC_IGP48+LC_IGP51)/LC_IGP60*100))
    dteh(data$LC_IGP65 <- with(data, (LC_IGP46+LC_IGP47+LC_IGP48)))
    dteh(data$LC_IGP66 <- with(data, LC_IGP46/LC_IGP58*100))
    dteh(data$LC_IGP67 <- with(data, LC_IGP47/LC_IGP59*100))
    dteh(data$LC_IGP68 <- with(data, LC_IGP48/LC_IGP60*100))
    dteh(data$LC_IGP69 <- with(data, (LC_IGP49+LC_IGP50+LC_IGP51)))
    dteh(data$LC_IGP70 <- with(data, LC_IGP49/LC_IGP58*100))
    dteh(data$LC_IGP71 <- with(data, LC_IGP50/LC_IGP59*100))
    dteh(data$LC_IGP72 <- with(data, LC_IGP51/LC_IGP60*100))
    dteh(data$LC_IGP73 <- with(data, (LC_IGP55+LC_IGP56+LC_IGP57+LC_IGP49+LC_IGP50+LC_IGP51)))
    dteh(data$LC_IGP74 <- with(data, (LC_IGP55+LC_IGP49)/LC_IGP58*100))
    dteh(data$LC_IGP75 <- with(data, (LC_IGP56+LC_IGP50)/LC_IGP59*100))
    dteh(data$LC_IGP76 <- with(data, (LC_IGP57+LC_IGP51)/LC_IGP60*100))
    dteh(data$LC_IGP77 <- with(data, (LC_IGP55+LC_IGP56+LC_IGP57)))
    dteh(data$LC_IGP78 <- with(data, LC_IGP55/LC_IGP58*100))
    dteh(data$LC_IGP79 <- with(data, LC_IGP56/LC_IGP59*100))
    dteh(data$LC_IGP80 <- with(data, LC_IGP57/LC_IGP60*100))
    dteh(data$LC_IGP81 <- with(data, LC_IGP65/LC_IGP77))
    dteh(data$LC_IGP82 <- with(data, LC_IGP69/LC_IGP65))
    dteh(data$LC_IGP83 <- with(data, LC_IGP69/LC_IGP61*100))
    dteh(data$LC_IGP84 <- with(data, LC_IGP69/LC_IGP73*100))
    dteh(data$LC_IGP85 <- with(data, LC_IGP65/LC_IGP73))
    dteh(data$LC_IGP86 <- with(data, LC_IGP77/LC_IGP61*1000))

    # =======================================
    # IgG2 derived traits
    # =======================================

    dteh(data$LC_IGP107 <- with(data, IgG2_G0F+IgG2_G1F+IgG2_G2F+IgG2_G0FN+IgG2_G1FN+IgG2_G2FN+IgG2_G1FS1+IgG2_G2FS1+IgG2_G1FNS1+IgG2_G2FNS1))
    dteh(data$LC_IGP108 <- with(data, IgG2_G0FN+IgG2_G1FN+IgG2_G2FN+IgG2_G1FNS1+IgG2_G2FNS1+IgG2_G0N+IgG2_G1N+IgG2_G2N+IgG2_G1NS1+IgG2_G2NS1))
    dteh(data$LC_IGP109 <- with(data, (IgG2_G1F+IgG2_G1FN+IgG2_G1FS1+IgG2_G1FNS1+IgG2_G1+IgG2_G1N+IgG2_G1S1+IgG2_G1NS1)*0.5+(IgG2_G2F+IgG2_G2FN+IgG2_G2FS1+IgG2_G2FNS1+IgG2_G2+IgG2_G2N+IgG2_G2S1+IgG2_G2NS1)))
    dteh(data$LC_IGP110 <- with(data, IgG2_G1FS1+IgG2_G2FS1+IgG2_G1FNS1+IgG2_G2FNS1+IgG2_G1S1+IgG2_G2S1+IgG2_G1NS1+IgG2_G2NS1))
    dteh(data$LC_IGP111 <- with(data, LC_IGP110/(2*LC_IGP109)*100))
    dteh(data$LC_IGP112 <- with(data, (IgG2_G1S1+IgG2_G2S1)/(IgG2_G1+IgG2_G1S1+IgG2_G2+IgG2_G2S1)*100))
    dteh(data$LC_IGP113 <- with(data, (IgG2_G1S1+IgG2_G2S1)/(IgG2_G0+IgG2_G1+IgG2_G1S1+IgG2_G2+IgG2_G2S1)*100))
    dteh(data$LC_IGP114 <- with(data, IgG2_G1S1/(IgG2_G1+IgG2_G1S1)*100))
    dteh(data$LC_IGP115 <- with(data, IgG2_G2S1/(IgG2_G2+IgG2_G2S1)*100))
    dteh(data$LC_IGP116 <- with(data, (IgG2_G1NS1+IgG2_G2NS1)/(IgG2_G1N+IgG2_G1NS1+IgG2_G2N+IgG2_G2NS1)*100))
    dteh(data$LC_IGP117 <- with(data, (IgG2_G1NS1+IgG2_G2NS1)/(IgG2_G0N+IgG2_G1N+IgG2_G1NS1+IgG2_G2N+IgG2_G2NS1)*100))
    dteh(data$LC_IGP118 <- with(data, IgG2_G1NS1/(IgG2_G1N+IgG2_G1NS1)*100))
    dteh(data$LC_IGP119 <- with(data, IgG2_G2NS1/(IgG2_G2N+IgG2_G2NS1)*100))
    dteh(data$LC_IGP120 <- with(data, (IgG2_G1FS1+IgG2_G2FS1)/(IgG2_G1F+IgG2_G1FS1+IgG2_G2F+IgG2_G2FS1)*100))
    dteh(data$LC_IGP121 <- with(data, (IgG2_G1FS1+IgG2_G2FS1)/(IgG2_G0F+IgG2_G1F+IgG2_G1FS1+IgG2_G2F+IgG2_G2FS1)*100))
    dteh(data$LC_IGP122 <- with(data, IgG2_G1FS1/(IgG2_G1F+IgG2_G1FS1)*100))
    dteh(data$LC_IGP123 <- with(data, IgG2_G2FS1/(IgG2_G2F+IgG2_G2FS1)*100))
    dteh(data$LC_IGP124 <- with(data, (IgG2_G1FNS1+IgG2_G2FNS1)/(IgG2_G1FN+IgG2_G1FNS1+IgG2_G2FN+IgG2_G2FNS1)*100))
    dteh(data$LC_IGP125 <- with(data, (IgG2_G1FNS1+IgG2_G2FNS1)/(IgG2_G0FN+IgG2_G1FN+IgG2_G1FNS1+IgG2_G2FN+IgG2_G2FNS1)*100))
    dteh(data$LC_IGP126 <- with(data, IgG2_G1FNS1/(IgG2_G1FN+IgG2_G1FNS1)*100))
    dteh(data$LC_IGP127 <- with(data, IgG2_G2FNS1/(IgG2_G2FN+IgG2_G2FNS1)*100))
    dteh(data$LC_IGP128 <- with(data, (IgG2_G1NS1+IgG2_G2NS1)/(IgG2_G1S1+IgG2_G2S1)))
    dteh(data$LC_IGP129 <- with(data, (IgG2_G1FNS1+IgG2_G2FNS1)/(IgG2_G1FS1+IgG2_G2FS1)))
    dteh(data$LC_IGP130 <- with(data, (IgG2_G1NS1+IgG2_G2NS1)/(IgG2_G1S1+IgG2_G1NS1+IgG2_G2S1+IgG2_G2NS1)))
    dteh(data$LC_IGP131 <- with(data, (IgG2_G1FNS1+IgG2_G2FNS1)/(IgG2_G1FS1+IgG2_G1FNS1+IgG2_G2FS1+IgG2_G2FNS1)))

    # neutral
    dteh(IgG2.neutral <- with(data, IgG2_G0F+IgG2_G1F+IgG2_G2F+IgG2_G0FN+IgG2_G1FN+IgG2_G2FN+IgG2_G0+IgG2_G1+IgG2_G2+IgG2_G0N+IgG2_G1N+IgG2_G2N),
                                 mess="Not all neutral glycans for IgG2 where found in the data frame")

    dteh(data$LC_IGP132 <- with(data, IgG2_G0F/IgG2.neutral) * 100)
    dteh(data$LC_IGP133 <- with(data, IgG2_G1F/IgG2.neutral) * 100)
    dteh(data$LC_IGP134 <- with(data, IgG2_G2F/IgG2.neutral) * 100)
    dteh(data$LC_IGP135 <- with(data, IgG2_G0FN/IgG2.neutral) * 100)
    dteh(data$LC_IGP136 <- with(data, IgG2_G1FN/IgG2.neutral) * 100)
    dteh(data$LC_IGP137 <- with(data, IgG2_G2FN/IgG2.neutral) * 100)
    dteh(data$LC_IGP138 <- with(data, IgG2_G0/IgG2.neutral) * 100)
    dteh(data$LC_IGP139 <- with(data, IgG2_G1/IgG2.neutral) * 100)
    dteh(data$LC_IGP140 <- with(data, IgG2_G2/IgG2.neutral) * 100)
    dteh(data$LC_IGP141 <- with(data, IgG2_G0N/IgG2.neutral) * 100)
    dteh(data$LC_IGP142 <- with(data, IgG2_G1N/IgG2.neutral) * 100)
    dteh(data$LC_IGP143 <- with(data, IgG2_G2N/IgG2.neutral) * 100)

    # neutral - derived traits
    dteh(data$LC_IGP144 <- with(data, (LC_IGP138+LC_IGP132+LC_IGP135+LC_IGP141)))
    dteh(data$LC_IGP145 <- with(data, (LC_IGP139+LC_IGP133+LC_IGP136+LC_IGP142)))
    dteh(data$LC_IGP146 <- with(data, (LC_IGP140+LC_IGP134+LC_IGP137+LC_IGP143)))
    dteh(data$LC_IGP147 <- with(data, (LC_IGP132+LC_IGP135+LC_IGP133+LC_IGP136+LC_IGP134+LC_IGP137)))
    dteh(data$LC_IGP148 <- with(data, (LC_IGP132+LC_IGP135)/LC_IGP144*100))
    dteh(data$LC_IGP149 <- with(data, (LC_IGP133+LC_IGP136)/LC_IGP145*100))
    dteh(data$LC_IGP150 <- with(data, (LC_IGP134+LC_IGP137)/LC_IGP146*100))
    dteh(data$LC_IGP151 <- with(data, (LC_IGP132+LC_IGP133+LC_IGP134)))
    dteh(data$LC_IGP152 <- with(data, LC_IGP132/LC_IGP144*100))
    dteh(data$LC_IGP153 <- with(data, LC_IGP133/LC_IGP145*100))
    dteh(data$LC_IGP154 <- with(data, LC_IGP134/LC_IGP146*100))
    dteh(data$LC_IGP155 <- with(data, (LC_IGP135+LC_IGP136+LC_IGP137)))
    dteh(data$LC_IGP156 <- with(data, LC_IGP135/LC_IGP144*100))
    dteh(data$LC_IGP157 <- with(data, LC_IGP136/LC_IGP145*100))
    dteh(data$LC_IGP158 <- with(data, LC_IGP137/LC_IGP146*100))
    dteh(data$LC_IGP159 <- with(data, (LC_IGP141+LC_IGP142+LC_IGP143+LC_IGP135+LC_IGP136+LC_IGP137)))
    dteh(data$LC_IGP160 <- with(data, (LC_IGP141+LC_IGP135)/LC_IGP144*100))
    dteh(data$LC_IGP161 <- with(data, (LC_IGP142+LC_IGP136)/LC_IGP145*100))
    dteh(data$LC_IGP162 <- with(data, (LC_IGP143+LC_IGP137)/LC_IGP146*100))
    dteh(data$LC_IGP163 <- with(data, (LC_IGP141+LC_IGP142+LC_IGP143)))
    dteh(data$LC_IGP164 <- with(data, LC_IGP141/LC_IGP144*100))
    dteh(data$LC_IGP165 <- with(data, LC_IGP142/LC_IGP145*100))
    dteh(data$LC_IGP166 <- with(data, LC_IGP143/LC_IGP146*100))
    dteh(data$LC_IGP167 <- with(data, LC_IGP151/LC_IGP163))
    dteh(data$LC_IGP168 <- with(data, LC_IGP155/LC_IGP151))
    dteh(data$LC_IGP169 <- with(data, LC_IGP155/LC_IGP147*100))
    dteh(data$LC_IGP170 <- with(data, LC_IGP155/LC_IGP159*100))
    dteh(data$LC_IGP171 <- with(data, LC_IGP151/LC_IGP159))
    dteh(data$LC_IGP172 <- with(data, LC_IGP163/LC_IGP147*1000))

    # =======================================
    # IgG4 derived traits
    # =======================================

    dteh(data$LC_IGP183 <- with(data, (IgG4_G0FN+IgG4_G1FN+IgG4_G2FN+IgG4_G1FNS1+IgG4_G2FNS1) ))
    dteh(data$LC_IGP184 <- with(data, (IgG4_G1F+IgG4_G1FN+IgG4_G1FS1+IgG4_G1FNS1)*0.5+(IgG4_G2F+IgG4_G2FN+IgG4_G2FS1+IgG4_G2FNS1)))
    dteh(data$LC_IGP185 <- with(data, (IgG4_G1FS1+IgG4_G2FS1+IgG4_G1FNS1+IgG4_G2FNS1)))
    dteh(data$LC_IGP186 <- with(data, LC_IGP185/(2*LC_IGP184)*100))
    dteh(data$LC_IGP187 <- with(data, (IgG4_G1FS1+IgG4_G2FS1)/(IgG4_G1F+IgG4_G1FS1+IgG4_G2F+IgG4_G2FS1)*100))
    dteh(data$LC_IGP188 <- with(data, (IgG4_G1FS1+IgG4_G2FS1)/(IgG4_G0F+IgG4_G1F+IgG4_G1FS1+IgG4_G2F+IgG4_G2FS1)*100))
    dteh(data$LC_IGP189 <- with(data, IgG4_G1FS1/(IgG4_G1F+IgG4_G1FS1)*100))
    dteh(data$LC_IGP190 <- with(data, IgG4_G2FS1/(IgG4_G2F+IgG4_G2FS1)*100))
    dteh(data$LC_IGP191 <- with(data, (IgG4_G1FNS1+IgG4_G2FNS1)/(IgG4_G1FN+IgG4_G1FNS1+IgG4_G2FN+IgG4_G2FNS1)*100))
    dteh(data$LC_IGP192 <- with(data, (IgG4_G1FNS1+IgG4_G2FNS1)/(IgG4_G0FN+IgG4_G1FN+IgG4_G1FNS1+IgG4_G2FN+IgG4_G2FNS1)*100))
    dteh(data$LC_IGP193 <- with(data, IgG4_G1FNS1/(IgG4_G1FN+IgG4_G1FNS1)*100))
    dteh(data$LC_IGP194 <- with(data, IgG4_G2FNS1/(IgG4_G2FN+IgG4_G2FNS1)*100))
    dteh(data$LC_IGP195 <- with(data, (IgG4_G1FNS1+IgG4_G2FNS1)/(IgG4_G1FS1+IgG4_G2FS1)))
    dteh(data$LC_IGP196 <- with(data, (IgG4_G1FNS1+IgG4_G2FNS1)/(IgG4_G1FS1+IgG4_G1FNS1+IgG4_G2FS1+IgG4_G2FNS1)))
  
    # neutral
    dteh(IgG4.neutral <- with(data, IgG4_G0F+IgG4_G1F+IgG4_G2F+IgG4_G0FN+IgG4_G1FN+IgG4_G2FN),
                                 mess="Not all neutral glycans for IgG4 where found in the data frame")

    dteh(data$LC_IGP197 <- with(data, IgG4_G0F/IgG4.neutral) * 100)
    dteh(data$LC_IGP198 <- with(data, IgG4_G1F/IgG4.neutral) * 100)
    dteh(data$LC_IGP199 <- with(data, IgG4_G2F/IgG4.neutral) * 100)
    dteh(data$LC_IGP200 <- with(data, IgG4_G0FN/IgG4.neutral) * 100)
    dteh(data$LC_IGP201 <- with(data, IgG4_G1FN/IgG4.neutral) * 100)
    dteh(data$LC_IGP202 <- with(data, IgG4_G2FN/IgG4.neutral) * 100)

    # neutral - derived traits
    dteh(data$LC_IGP203 <- with(data, LC_IGP197 + LC_IGP200))
    dteh(data$LC_IGP204 <- with(data, LC_IGP198 + LC_IGP201))
    dteh(data$LC_IGP205 <- with(data, LC_IGP199 + LC_IGP202))

    return(data)
}

# Translate names between computer readable and human readable
# for derived traits of IgG with LCMS
# based on paper from 2014.
#
# Translates names between computer readable and human readable
# for derived traits of IgG with LCMS
#
# @author Ivo Ugrina
# @param orignames vector; type string
# @param to type of translation. If \code{inverse} is used everything will be
#   translated. For \code{computer} names will be translated to computer
#   readable, and for \code{human} names will be translated to human readable.
# @return Returns a character vector with original and translated names
# @references
# Jennifer E. Huffman et al. 
# "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
# \url{http://dx.doi.org/10.1074/mcp.M113.037465}
ildt.translate.2014 <- function(orignames, to="inverse") {

    allnames <- "computer	human
LC_IGP1	IgG1_G0F
LC_IGP2	IgG1_G1F
LC_IGP3	IgG1_G2F
LC_IGP4	IgG1_G0FN
LC_IGP5	IgG1_G1FN
LC_IGP6	IgG1_G2FN
LC_IGP7	IgG1_G1FS1
LC_IGP8	IgG1_G2FS1
LC_IGP9	IgG1_G1FNS1
LC_IGP10	IgG1_G2FNS1
LC_IGP11	IgG1_G0
LC_IGP12	IgG1_G1
LC_IGP13	IgG1_G2
LC_IGP14	IgG1_G0N
LC_IGP15	IgG1_G1N
LC_IGP16	IgG1_G2N
LC_IGP17	IgG1_G1S1
LC_IGP18	IgG1_G2S1
LC_IGP19	IgG1_G1NS1
LC_IGP20	IgG1_G2NS1
LC_IGP21	IgG1 Fucosylation
LC_IGP22	IgG1 Bisecting_GlcNAc
LC_IGP23	IgG1 Galactosylation
LC_IGP24	IgG1 Sialylation
LC_IGP25	IgG1 SA per Gal
LC_IGP26	IgG1 GS1/(G+GS1)
LC_IGP27	IgG1 GS1/(G0+G+GS1)
LC_IGP28	IgG1 G1S1/(G1+G1S1)
LC_IGP29	IgG1 G2S1/(G2+G2S1)
LC_IGP30	IgG1 BGS1/(BG+BGS1)
LC_IGP31	IgG1 BGS1/(BG0+BG+BGS1)
LC_IGP32	IgG1 BG1S1/(BG1+BG1S1)
LC_IGP33	IgG1 BG2S1/(BG2+BG2S1)
LC_IGP34	IgG1 FGS1/(FG+FGS1)
LC_IGP35	IgG1 FGS1/(F+FG+FGS1)
LC_IGP36	IgG1 FG1S1/(FG1+FG1S1)
LC_IGP37	IgG1 FG2S1/(FG2+FG2S1)
LC_IGP38	IgG1 FBGS1/(FBG+FBGS1)
LC_IGP39	IgG1 FBGS1/(FB+FBG+FBGS1)
LC_IGP40	IgG1 FBG1S1/(FBG1+FBG1S1)
LC_IGP41	IgG1 FBG2S1/(FBG2+FBG2S1)
LC_IGP42	IgG1 BS1/S1
LC_IGP43	IgG1 FBS1/FS1
LC_IGP44	IgG1 BS1/(S1+BS1)
LC_IGP45	IgG1 FBS1/(FS1+FBS1)
LC_IGP46	IgG1_G0Fn
LC_IGP47	IgG1_G1Fn
LC_IGP48	IgG1_G2Fn
LC_IGP49	IgG1_G0FNn
LC_IGP50	IgG1_G1FNn
LC_IGP51	IgG1_G2FNn
LC_IGP52	IgG1_G0n
LC_IGP53	IgG1_G1n
LC_IGP54	IgG1_G2n
LC_IGP55	IgG1_G0Nn
LC_IGP56	IgG1_G1Nn
LC_IGP57	IgG1_G2Nn
LC_IGP58	IgG1G0n
LC_IGP59	IgG1G1n
LC_IGP60	IgG1G2n
LC_IGP61	IgG1 Fn total
LC_IGP62	IgG1  FG0n total/G0n
LC_IGP63	IgG1  FG1n total/G1n
LC_IGP64	IgG1  FG2n total/G2n
LC_IGP65	IgG1 Fn
LC_IGP66	IgG1 FG0n/G0n
LC_IGP67	IgG1 FG1n/G1n
LC_IGP68	IgG1 FG2n/G2n
LC_IGP69	IgG1 FBn
LC_IGP70	IgG1 FBG0n/G0n
LC_IGP71	IgG1 FBG1n/G1n
LC_IGP72	IgG1 FBG2n/G2n
LC_IGP73	IgG1 Bn total
LC_IGP74	IgG1 BG0n total/G0n
LC_IGP75	IgG1 BG1n total/G1n
LC_IGP76	IgG1 BG2n total/G2n
LC_IGP77	IgG1 Bn
LC_IGP78	IgG1 BG0n/G0n
LC_IGP79	IgG1 BG1n/G1n
LC_IGP80	IgG1 BG2n/G2n
LC_IGP81	IgG1 Fn/Bn
LC_IGP82	IgG1 FBn/Fn
LC_IGP83	IgG1 FBn/Fn total
LC_IGP84	IgG1 FBn/Bn total
LC_IGP85	IgG1 Fn/Bn total
LC_IGP86	IgG1 Bn/Fn total 
LC_IGP87	IgG2_G0F
LC_IGP88	IgG2_G1F
LC_IGP89	IgG2_G2F
LC_IGP90	IgG2_G0FN
LC_IGP91	IgG2_G1FN
LC_IGP92	IgG2_G2FN
LC_IGP93	IgG2_G1FS1
LC_IGP94	IgG2_G2FS1
LC_IGP95	IgG2_G1FNS1
LC_IGP96	IgG2_G2FNS1
LC_IGP97	IgG2_G0
LC_IGP98	IgG2_G1
LC_IGP99	IgG2_G2
LC_IGP100	IgG2_G0N
LC_IGP101	IgG2_G1N
LC_IGP102	IgG2_G2N
LC_IGP103	IgG2_G1S1
LC_IGP104	IgG2_G2S1
LC_IGP105	IgG2_G1NS1
LC_IGP106	IgG2_G2NS1
LC_IGP107	IgG2 Fucosylation
LC_IGP108	IgG2 Bisecting GlcNAc
LC_IGP109	IgG2 Galactosylation
LC_IGP110	IgG2 Sialylation
LC_IGP111	IgG2 SA per Gal
LC_IGP112	IgG2 GS1/(G+GS1)
LC_IGP113	IgG2 GS1/(G0+G+GS1)
LC_IGP114	IgG2 G1S1/(G1+G1S1)
LC_IGP115	IgG2 G2S1/(G2+G2S1)
LC_IGP116	IgG2 BGS1/(BG+BGS1)
LC_IGP117	IgG2 BGS1/(BG0+BG+BGS1)
LC_IGP118	IgG2 BG1S1/(BG1+BG1S1)
LC_IGP119	IgG2 BG2S1/(BG2+BG2S1)
LC_IGP120	IgG2 FGS1/(FG+FGS1)
LC_IGP121	IgG2 FGS1/(F+FG+FGS1)
LC_IGP122	IgG2 FG1S1/(FG1+FG1S1)
LC_IGP123	IgG2 FG2S1/(FG2+FG2S1)
LC_IGP124	IgG2 FBGS1/(FBG+FBGS1)
LC_IGP125	IgG2 FBGS1/(FB+FBG+FBGS1)
LC_IGP126	IgG2 FBG1S1/(FBG1+FBG1S1)
LC_IGP127	IgG2 FBG2S1/(FBG2+FBG2S1)
LC_IGP128	IgG2 BS1/S1
LC_IGP129	IgG2 FBS1/FS1
LC_IGP130	IgG2 BS1/(S1+BS1)
LC_IGP131	IgG2 FBS1/(FS1+FBS1)
LC_IGP132	IgG2_G0Fn
LC_IGP133	IgG2_G1Fn
LC_IGP134	IgG2_G2Fn
LC_IGP135	IgG2_G0FNn
LC_IGP136	IgG2_G1FNn
LC_IGP137	IgG2_G2FNn
LC_IGP138	IgG2_G0n
LC_IGP139	IgG2_G1n
LC_IGP140	IgG2_G2n
LC_IGP141	IgG2_G0Nn
LC_IGP142	IgG2_G1Nn
LC_IGP143	IgG2_G2Nn
LC_IGP144	IgG2G0n
LC_IGP145	IgG2G1n
LC_IGP146	IgG2G2n
LC_IGP147	IgG2 Fn total
LC_IGP148	IgG2  FG0n total/G0n
LC_IGP149	IgG2  FG1n total/G1n
LC_IGP150	IgG2  FG2n total/G2n
LC_IGP151	IgG2 Fn
LC_IGP152	IgG2 FG0n/G0n
LC_IGP153	IgG2 FG1n/G1n
LC_IGP154	IgG2 FG2n/G2n
LC_IGP155	IgG2 FBn
LC_IGP156	IgG2 FBG0n/G0n
LC_IGP157	IgG2 FBG1n/G1n
LC_IGP158	IgG2 FBG2n/G2n
LC_IGP159	IgG2 Bn total
LC_IGP160	IgG2 BG0n total/G0n
LC_IGP161	IgG2 BG1n total/G1n
LC_IGP162	IgG2 BG2n total/G2n
LC_IGP163	IgG2 Bn
LC_IGP164	IgG2 BG0n/G0n
LC_IGP165	IgG2 BG1n/G1n
LC_IGP166	IgG2 BG2n/G2n
LC_IGP167	IgG2 Fn/Bn
LC_IGP168	IgG2 FBn/Fn
LC_IGP169	IgG2 FBn/Fn total
LC_IGP170	IgG2 FBn/Bn total
LC_IGP171	IgG2 Fn/Bn total
LC_IGP172	IgG2 Bn/Fn total
LC_IGP173	IgG4_G0F
LC_IGP174	IgG4_G1F
LC_IGP175	IgG4_G2F
LC_IGP176	IgG4_G0FN
LC_IGP177	IgG4_G1FN
LC_IGP178	IgG4_G2FN
LC_IGP179	IgG4_G1FS1
LC_IGP180	IgG4_G2FS1
LC_IGP181	IgG4_G1FNS1
LC_IGP182	IgG4_G2FNS1
LC_IGP183	IgG4 Bisecting GlcNAc
LC_IGP184	IgG4 Galactosylation
LC_IGP185	IgG4 Sialylation
LC_IGP186	IgG4 SA per Gal
LC_IGP187	IgG4 FGS1/(FG+FGS1)
LC_IGP188	IgG4 FGS1/(F+FG+FGS1)
LC_IGP189	IgG4 FG1S1/(FG1+FG1S1)
LC_IGP190	IgG4 FG2S1/(FG2+FG2S1)
LC_IGP191	IgG4 FBGS1/(FBG+FBGS1)
LC_IGP192	IgG4 FBGS1/(FB+FBG+FBGS1)
LC_IGP193	IgG4 FBG1S1/(FBG1+FBG1S1)
LC_IGP194	IgG4 FBG2S1/(FBG2+FBG2S1)
LC_IGP195	IgG4 FBS1/FS1
LC_IGP196	IgG4 FBS1/(FS1+FBS1)
LC_IGP197	IgG4_G0Fn
LC_IGP198	IgG4_G1Fn
LC_IGP199	IgG4_G2Fn
LC_IGP200	IgG4_G0FNn
LC_IGP201	IgG4_G1FNn
LC_IGP202	IgG4_G2FNn
LC_IGP203	IgG4G0n
LC_IGP204	IgG4G1n
LC_IGP205	IgG4G2n"

    con <- textConnection(allnames)
    allnames <- read.delim(con, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    close(con)

    cn <- match(orignames, allnames$computer)
    hn <- match(orignames, allnames$human)
    if(to=="inverse"){
        newnames <- ifelse(!is.na(cn),
                           allnames$human[cn],
                           allnames$computer[hn])
    }
    if(to=="computer"){
        newnames <- ifelse(!is.na(cn),
                           allnames$computer[cn],
                           allnames$computer[hn])
    }
    if(to=="human"){
        newnames <- ifelse(!is.na(hn),
                           allnames$human[hn],
                           allnames$human[cn])
    }

    return(newnames)
}


# Translate names between computer readable and human readable
# for derived traits of IgG with UPLC
# based on an article from 2014.
#
# Translates names between computer readable and human readable
# for derived traits of IgG with UPLC 
#
# @author Ivo Ugrina
# @param orignames vector; type string
# @param to type of translation. If \code{inverse} is used everything will be
#   translated. For \code{computer} names will be translated to computer
#   readable, and for \code{human} names will be translated to human readable.
# @return Returns a character vector with original and translated names
# @references
# Jennifer E. Huffman et al. 
# "Comparative Performance of Four Methods for High-throughput Glycosylation Analysis of Immunoglobulin G in Genetic and Epidemiological Research*"
# \url{http://dx.doi.org/10.1074/mcp.M113.037465}
iudt.translate.2014 <- function(orignames, to="inverse") {

    allnames <- "computer	human
IGP1	GP1
IGP2	GP2
IGP3	GP4
IGP4	GP5
IGP5	GP6
IGP6	GP7
IGP7	GP8
IGP8	GP9
IGP9	GP10
IGP10	GP11
IGP11	GP12
IGP12	GP13
IGP13	GP14
IGP14	GP15
IGP15	GP16
IGP16	GP17
IGP17	GP18
IGP18	GP19
IGP19	GP20
IGP20	GP21
IGP21	GP22
IGP22	GP23
IGP23	GP24
IGP24	FGS/(FG+FGS)
IGP25	FBGS/(FBG+FBGS)
IGP26	FGS/(F+FG+FGS)
IGP27	FBGS/(FB+FBG+FBGS)
IGP28	FG1S1/(FG1+FG1S1)
IGP29	FG2S1/(FG2+FG2S1+FG2S2)
IGP30	FG2S2/(FG2+FG2S1+FG2S2)
IGP31	FBG2S1/(FBG2+FBG2S1+FBG2S2)
IGP32	FBG2S2/(FBG2+FBG2S1+FBG2S2)
IGP33	FtotalS1/FtotalS2
IGP34	FS1/FS2
IGP35	FBS1/FBS2
IGP36	FBStotal/FStotal
IGP37	FBS1/FS1
IGP38	FBS1/(FS1+FBS1)
IGP39	FBS2/FS2
IGP40	FBS2/(FS2+FBS2)
IGP41	GP1n
IGP42	GP2n
IGP43	GP4n
IGP44	GP5n
IGP45	GP6n
IGP46	GP7n
IGP47	GP8n
IGP48	GP9n
IGP49	GP10n
IGP50	GP11n
IGP51	GP12n
IGP52	GP13n
IGP53	GP14n
IGP54	GP15n
IGP55	G0n
IGP56	G1n
IGP57	G2n
IGP58	Fn total
IGP59	FG0n total/G0n
IGP60	FG1n total/G1n
IGP61	FG2n total/G2n
IGP62	Fn
IGP63	FG0n/G0n
IGP64	FG1n/G1n
IGP65	FG2n/G2n
IGP66	FBn
IGP67	FBG0n/G0n
IGP68	FBG1n/G1n
IGP69	FBG2n/G2n
IGP70	FBn/Fn
IGP71	FBn/Fn total
IGP72	Fn/(Bn+ FBn)
IGP73	Bn/(Fn+ FBn)
IGP74	FBG2n/FG2n
IGP75	FBG2n/(FG2n+ FBG2n)
IGP76	FG2n/(BG2n+ FBG2n)
IGP77	BG2n/(FG2n+ FBG2n)"

    con <- textConnection(allnames)
    allnames <- read.delim(con, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    close(con)

    cn <- match(orignames, allnames$computer)
    hn <- match(orignames, allnames$human)
    if(to=="inverse"){
        newnames <- ifelse(!is.na(cn),
                           allnames$human[cn],
                           allnames$computer[hn])
    }
    if(to=="computer"){
        newnames <- ifelse(!is.na(cn),
                           allnames$computer[cn],
                           allnames$computer[hn])
    }
    if(to=="human"){
        newnames <- ifelse(!is.na(hn),
                           allnames$human[hn],
                           allnames$human[cn])
    }

    return(newnames)
}

