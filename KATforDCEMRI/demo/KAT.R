##    KATforDCEMRI: a Kinetic Analysis Tool for DCE-MRI
##    Copyright 2014 Genentech, Inc.
##
##    For questions or comments, please contact
##    Gregory Z. Ferl, Ph.D.
##    Genentech, Inc.
##    Development Sciences
##    1 DNA Way, Mail stop 463A
##    South San Francisco, CA, United States of America
##    E-mail: ferl.gregory@gene.com

runme <- function(){
    data(dcemri.data, package="KATforDCEMRI")

    dir.create("KATforDCEMRI_benchmark_test")
    setwd("KATforDCEMRI_benchmark_test")

    attach(dcemri.data)

    ## SHRINK THE ROI MASK
    maskROI[,,] <- 0
    maskROI[32:42,32:42,] <- 1

    runtime1 <- system.time(KAT.checkData(file.name="KAT", vector.times=vectorTimes, map.CC=mapCC, mask.ROI=maskROI, vector.AIF=vectorAIF))
    runtime2 <- system.time(KAT(file = "KAT.RData", results_file="KAT_benchmark_test-full", range.map=1.05, cutoff.map=0.95, AIF.shift="NONE", tlag.Tofts.on=FALSE, export.matlab=FALSE))

    ## runtime3 <- system.time(KAT.checkData(file.name="KATtrunc", vector.times=vectorTimes[1:44], map.CC=mapCC[,,,1:44], mask.ROI=maskROI, vector.AIF=vectorAIF[1:44]))
    ## runtime4 <- system.time(KAT(file = "KATtrunc.RData", results_file="KAT_benchmark_test-truncated", range.map=1.05, cutoff.map=0.95))
    detach(dcemri.data)

    ## runtime <- format(runtime1[3] + runtime2[3] + runtime3[3] + runtime4[3], digits=3)
    runtime <- format(runtime1[3] + runtime2[3], digits=3)

    ## KAT.plot(F1="KAT_benchmark_test-full_slice1.RData", F2="KAT_benchmark_test-full_slice2.RData", F3="KAT_benchmark_test-truncated_slice1.RData", F4="KAT_benchmark_test-truncated_slice2.RData", export.matlab=FALSE)
    KAT.plot(F1="KAT_benchmark_test-full_slice1.RData", F2="KAT_benchmark_test-full_slice2.RData", F3="KAT_benchmark_test-full_slice3.RData", F4="KAT_benchmark_test-full_slice4.RData", export.matlab=FALSE)

    load("KAT_benchmark_test-full_slice1.RData")
    Ktrans_s1A <- as.numeric(format(dcemri.data$paramestmedianxT$Ktrans.median, digits=3))
    kep_s1A <- as.numeric(format(dcemri.data$paramestmedianxT$kep.median, digits=3))
    vb_s1A <- as.numeric(format(dcemri.data$paramestmedianxT$vb.median, digits=3))
    cvKtrans_s1A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransxT), na.rm=TRUE)/Ktrans_s1A, digits=3))
    cvkep_s1A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepxT), na.rm=TRUE)/kep_s1A, digits=3))
    cvvb_s1A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapvbxT), na.rm=TRUE)/vb_s1A, digits=3))

    Ktrans_Ts1A <- as.numeric(format(dcemri.data$paramestmedianT$Ktrans.median, digits=3))
    kep_Ts1A <- as.numeric(format(dcemri.data$paramestmedianT$kep.median, digits=3))
    cvKtrans_Ts1A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransT), na.rm=TRUE)/Ktrans_Ts1A, digits=3))
    cvkep_Ts1A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepT), na.rm=TRUE)/kep_s1A, digits=3))

    load("KAT_benchmark_test-full_slice2.RData")
    Ktrans_s2A <- as.numeric(format(dcemri.data$paramestmedianxT$Ktrans.median, digits=3))
    kep_s2A <- as.numeric(format(dcemri.data$paramestmedianxT$kep.median, digits=3))
    vb_s2A <- as.numeric(format(dcemri.data$paramestmedianxT$vb.median, digits=3))
    cvKtrans_s2A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransxT), na.rm=TRUE)/Ktrans_s2A, digits=3))
    cvkep_s2A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepxT), na.rm=TRUE)/kep_s2A, digits=3))
    cvvb_s2A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapvbxT), na.rm=TRUE)/vb_s2A, digits=3))

    Ktrans_Ts2A <- as.numeric(format(dcemri.data$paramestmedianT$Ktrans.median, digits=3))
    kep_Ts2A <- as.numeric(format(dcemri.data$paramestmedianT$kep.median, digits=3))
    cvKtrans_Ts2A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransT), na.rm=TRUE)/Ktrans_Ts2A, digits=3))
    cvkep_Ts2A <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepT), na.rm=TRUE)/kep_s2A, digits=3))

    ## load("KAT_benchmark_test-truncated_slice1.RData")
    load("KAT_benchmark_test-full_slice3.RData")
    Ktrans_s1B <- as.numeric(format(dcemri.data$paramestmedianxT$Ktrans.median, digits=3))
    kep_s1B <- as.numeric(format(dcemri.data$paramestmedianxT$kep.median, digits=3))
    vb_s1B <- as.numeric(format(dcemri.data$paramestmedianxT$vb.median, digits=3))
    cvKtrans_s1B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransxT), na.rm=TRUE)/Ktrans_s1B, digits=3))
    cvkep_s1B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepxT), na.rm=TRUE)/kep_s1B, digits=3))
    cvvb_s1B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapvbxT), na.rm=TRUE)/vb_s1B, digits=3))

    Ktrans_Ts1B <- as.numeric(format(dcemri.data$paramestmedianT$Ktrans.median, digits=3))
    kep_Ts1B <- as.numeric(format(dcemri.data$paramestmedianT$kep.median, digits=3))
    cvKtrans_Ts1B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransT), na.rm=TRUE)/Ktrans_Ts1B, digits=3))
    cvkep_Ts1B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepT), na.rm=TRUE)/kep_s1B, digits=3))

    ## load("KAT_benchmark_test-truncated_slice2.RData")
    load("KAT_benchmark_test-full_slice4.RData")
    Ktrans_s2B <- as.numeric(format(dcemri.data$paramestmedianxT$Ktrans.median, digits=3))
    kep_s2B <- as.numeric(format(dcemri.data$paramestmedianxT$kep.median, digits=3))
    vb_s2B <- as.numeric(format(dcemri.data$paramestmedianxT$vb.median, digits=3))
    cvKtrans_s2B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransxT), na.rm=TRUE)/Ktrans_s2B, digits=3))
    cvkep_s2B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepxT), na.rm=TRUE)/kep_s2B, digits=3))
    cvvb_s2B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapvbxT), na.rm=TRUE)/vb_s2B, digits=3))

    Ktrans_Ts2B <- as.numeric(format(dcemri.data$paramestmedianT$Ktrans.median, digits=3))
    kep_Ts2B <- as.numeric(format(dcemri.data$paramestmedianT$kep.median, digits=3))
    cvKtrans_Ts2B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapKtransT), na.rm=TRUE)/Ktrans_Ts2B, digits=3))
    cvkep_Ts2B <- as.numeric(format(100*sd(as.vector(dcemri.data$mapkepT), na.rm=TRUE)/kep_s2B, digits=3))

    pdf(file="KAT_demo-page1.pdf", height=11, width=8.5)

    plot(0:35, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
    text(5, 35, paste("KATforDCEMRI version ", dcemri.data$KATversion, " BENCHMARK TEST", sep=""), pos=4, font=2, col="red")
    text(5, 33, paste("date:", date(), "\n"), pos=4)
    text(5, 32, paste("Total Processing Time:", runtime, "seconds"), pos=4)
    text(5, 29, paste("sysname/release:", Sys.info()[[1]], Sys.info()[[2]], "\n"), pos=4)
    text(5, 28, paste("nodename:", Sys.info()[[4]], "\n"), pos=4)
    text(5, 27, paste("user:", Sys.info()[[7]], "\n"), pos=4)
    text(5, 25, "Estimated (inter-voxel %CV) versus True Parameter values (extended Tofts)", pos=4, font=2)
    text(5, 23, "LOW NOISE DATA", pos=4)
    text(5, 22, paste("Slice 1: Ktrans =", Ktrans_s1A, "1/min (", cvKtrans_s1A,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 21, paste("Slice 1: kep =", kep_s1A, "1/min (", cvkep_s1A, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 20, paste("Slice 1: vb =", vb_s1A, " (", cvvb_s1A, "%) [true value = 0]", sep=""), pos=4)
    text(5, 18, paste("Slice 2: Ktrans =", Ktrans_s2A, "1/min (", cvKtrans_s2A,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 17, paste("Slice 2: kep =", kep_s2A, "1/min (", cvkep_s2A, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 16, paste("Slice 2: vb =", vb_s2A, " (", cvvb_s2A, "%) [true value = 0.05]", sep=""), pos=4)
    text(5, 14, "HIGH NOISE DATA", pos=4)
    text(5, 13, paste("Slice 3: Ktrans =", Ktrans_s1B, "1/min (", cvKtrans_s1B,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 12, paste("Slice 3: kep =", kep_s1B, "1/min (", cvkep_s1B, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 11, paste("Slice 3: vb =", vb_s1B, " (", cvvb_s1B, "%) [true value = 0]", sep=""), pos=4)
    text(5, 9, paste("Slice 4: Ktrans =", Ktrans_s2B, "1/min (", cvKtrans_s2B,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 8, paste("Slice 4: kep =", kep_s2B, "1/min (", cvkep_s2B, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 7, paste("Slice 4: vb =", vb_s2B, " (", cvvb_s2B, "%) [true value = 0.05]", sep=""), pos=4)
    dev.off()


    pdf(file="KAT_demo-page2.pdf", height=11, width=8.5)

    plot(0:35, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")

    text(5, 25, "Estimated (inter-voxel %CV) versus True Parameter values (Tofts Model)", pos=4, font=2)
    text(5, 23, "LOW NOISE DATA", pos=4)
    text(5, 22, paste("Slice 1: Ktrans =", Ktrans_Ts1A, "1/min (", cvKtrans_Ts1A,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 21, paste("Slice 1: kep =", kep_Ts1A, "1/min (", cvkep_Ts1A, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 20, "Slice 1: vb = 0 (Fixed) [true value = 0]", pos=4)
    text(5, 18, paste("Slice 2: Ktrans =", Ktrans_Ts2A, "1/min (", cvKtrans_Ts2A,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 17, paste("Slice 2: kep =", kep_Ts2A, "1/min (", cvkep_Ts2A, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 16, "Slice 2: vb = 0 (Fixed) [true value = 0.05]", pos=4)
    text(5, 14, "HIGH NOISE DATA", pos=4)
    text(5, 13, paste("Slice 3: Ktrans =", Ktrans_Ts1B, "1/min (", cvKtrans_Ts1B,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 12, paste("Slice 3: kep =", kep_Ts1B, "1/min (", cvkep_Ts1B, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 11, "Slice 3: vb = 0 (Fixed) [true value = 0]", pos=4)
    text(5, 9, paste("Slice 4: Ktrans =", Ktrans_Ts2B, "1/min (", cvKtrans_Ts2B,"%) [true value = 0.22 1/min]", sep=""), pos=4)
    text(5, 8, paste("Slice 4: kep =", kep_Ts2B, "1/min (", cvkep_Ts2B, "%) [true value = 1.1 1/min]", sep=""), pos=4)
    text(5, 7, "Slice 4: vb = 0 (Fixed) [true value = 0.05]", pos=4)
    dev.off()
}

runme()
