#
#  N.D. Crins, C. Roever, A.D. Goralczyk, T. Friede.
#  Interleukin-2 receptor antagonists for pediatric liver transplant recipients:
#  A systematic review and meta-analysis of controlled studies.
#  Pediatric Transplantation, 18(8):839-850, 2014.
#

CrinsEtAl2014 <- data.frame("publication"=c("Spada (2006)", "Ganschow (2005)",
                                            "Gibelli (2004)", "Heffron (2003)",
                                            "Schuller (2005)", "Gras (2008)"),
                            "year"=c(2006, 2005, 2004, 2003, 2005, 2008),
                            "randomized"=factor(c("yes","no")[c(1,2,2,1,2,2)],
                                                levels=c("yes","no")),
                            "control.type"=factor(c("concurrent","historical",
                                                    "historical","concurrent",
                                                    "concurrent","historical"),
                                                  levels=c("concurrent","historical")),
                            "comparison"=factor(c("IL-2RA only","delayed CNI","no/low steroids")[c(3,1,1,2,1,3)],
                                                levels=c("IL-2RA only","delayed CNI","no/low steroids")),
                            "IL2RA"=factor(c("basiliximab","daclizumab")[c(1,1,1,2,2,1)]),
                            "CNI"=factor(c("cyclosporine A","tacrolimus")[c(2,1,1,2,2,2)]),
                            "MMF"=factor(c("yes","no")[c(2,2,2,1,1,2)],
                                         levels=c("yes","no")),
                            "followup"=c(12,36,6,24,6,36),
                            "exp.AR.events"=c(4,9,16,14,3,0),
                            "exp.SRR.events"=c(NA,4,NA,2,NA,1),
                            "exp.total"=c(36,54,28,61,18,50),
                            "cont.AR.events"=c(11,29,19,15,8,3),
                            "cont.SRR.events"=c(NA,6,NA,4,NA,4),
                            "cont.total"=c(36,54,28,20,12,34),
                            stringsAsFactors=FALSE)[c(4,3,5,2,1,6),]

rownames(CrinsEtAl2014) <- as.character(1:6)
