# ===========================================================================
# File: "benchmarkAll.R"
#                        Created: 2012-11-01 18:17:51
#              Last modification: 2012-11-13 13:39:30
# Author: Bernard Desgraupes
# e-mail: <bdesgraupes@users.sourceforge.net>
# Unit test file for the R package clusterCrit.
# ===========================================================================


pkg <- "clusterCrit"

if (require("rbenchmark", quietly = TRUE)) {
    library(package=pkg, character.only = TRUE)

    dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsInternal_400_4.Rdata")
    load(file=dataPath, envir=.GlobalEnv)

    # Internal indices
    names <- c("all",getCriteriaNames(TRUE))
    df <- benchmark(
	    intCriteria(traj_400_4, part_400_4[[4]], c("all")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Ball_Hall")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Banfeld_Raftery")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("C_index")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Calinski_Harabasz")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Davies_Bouldin")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Det_Ratio")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Dunn")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Gamma")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("G_plus")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI11")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI12")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI13")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI21")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI22")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI23")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI31")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI32")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI33")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI41")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI42")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI43")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI51")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI52")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("GDI53")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Ksq_DetW")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Log_Det_Ratio")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Log_SS_Ratio")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("McClain_Rao")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("PBM")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Point_Biserial")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Ray_Turi")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Ratkowsky_Lance")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Scott_Symons")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("SD_Scat")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("SD_Dis")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("S_Dbw")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Silhouette")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Tau")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Trace_W")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Trace_WiB")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Wemmert_Gancarski")),
	    intCriteria(traj_400_4, part_400_4[[4]], c("Xie_Beni")),
	    columns = c("test","elapsed", "replications", "relative", "user.self", "sys.self"),
	    replications=100, order="elapsed"
	)
	
# 	df <- cbind(names,df)
	print(df)

    # External indices
    dataPath <- file.path(path.package(package="clusterCrit"),"unitTests","data","testsExternal100.Rdata")
    load(file=dataPath, envir=.GlobalEnv)
	
    names <- c("all",getCriteriaNames(FALSE))
    df <- benchmark(
	    extCriteria(clus_p2, clus_p3, c("all")),
	    extCriteria(clus_p2, clus_p3, c("Czekanowski_Dice")),
	    extCriteria(clus_p2, clus_p3, c("Folkes_Mallows")), 
	    extCriteria(clus_p2, clus_p3, c("Hubert")),
	    extCriteria(clus_p2, clus_p3, c("Jaccard")), 
	    extCriteria(clus_p2, clus_p3, c("Kulczynski")), 
	    extCriteria(clus_p2, clus_p3, c("McNemar")),
	    extCriteria(clus_p2, clus_p3, c("Phi")), 
	    extCriteria(clus_p2, clus_p3, c("Precision")), 
	    extCriteria(clus_p2, clus_p3, c("Rand")), 
	    extCriteria(clus_p2, clus_p3, c("Recall")), 
	    extCriteria(clus_p2, clus_p3, c("Rogers_Tanimoto")),
	    extCriteria(clus_p2, clus_p3, c("Russel_Rao")),
	    extCriteria(clus_p2, clus_p3, c("Sokal_Sneath1")),
	    extCriteria(clus_p2, clus_p3, c("Sokal_Sneath2")),
	    columns = c("elapsed", "replications", "relative", "user.self", "sys.self"),
	    replications=100
	)
	
	df <- cbind(names,df)
	print(df)

} else {
    cat("R package 'rbenchmark' cannot be loaded.\n")
}




# # #     names <- c("all",getCriteriaNames(TRUE))
# # #     len <- length(names)
# # #     
# # #     tests <- list()
# # #     for (i in 1:len) {
# # #         tests[[i]] <- expression( intCriteria(traj_400_4, part_400_4[[4]], c(names[i])) )
# # #     }
# # # 
# # #     df <- do.call(benchmark,
# # # 	    c(tests, list(replications=10,
# # #                       columns=c("elapsed", "replications", "relative", "user.self", "sys.self"))))
# # # 	
# # #     df <- cbind(names,df)
# # #     print(df)
