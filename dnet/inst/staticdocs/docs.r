library(staticdocs)
library(grid)
list(
    readme = "",
    
    index = list(
        sd_section("Network analysis and visualisation functions",
            "These functions are used for network analysis and visualisation, including: identification of gene-active networks, and network-based sample classifications and visualisations on 2D sample landscape.",
            c(
                'dNetPipeline',
                'dNetReorder',
                'dNetConfidence',
                'dPvalAggregate',
                'dNetInduce',
                'dBUMfit',
                'dBUMscore',
                'dSVDsignif',
                'dFDRscore',
                'dNetFind',
                'dCommSignif',
                'dContrast',
                'visNet',
                'visNetMul',
                'visNetAnimate',
                'visNetReorder',
                'visNetArc',
                'visNetCircle',
                'visBoxplotAdv',
                'dRDataLoader',
                'dFunArgs'
            )
        ),
        sd_section("Random Walk with Restart (RWR)",
            "These functions are used for RWR-based analysis, including: calculation of network affinity and estimation of contact strength between samples.",
            c(
                'dRWR',
                'dRWRcontact',
                'dRWRpipeline',
                'dCheckParallel'
            )
        ),
        sd_section("Enrichment analysis and visualisation functions",
            "These functions are used for enrichment analysis and visualisation, including: enrichment analysis (Fisher's exact test, Hypergeometric test or Binomial test) that can account for the hierarchy of the ontology; and gene set enrichment analysis (GSEA).",
            c(
                'dEnricher',
                'dEnricherView',
                'dGSEA',
                'dGSEAview',
                'dGSEAwrite',
                'visGSEA'
            )
        ),
        sd_section("Built-in ontologies, and supporting functions for analysis and visualisation",
            "These ontologies each are represented as a direct acyclic graph (DAG). DAG is stored as an object of the class 'igraph', which can be processed by various analysis and visualisation functions.",
            c(
                "ig.GOMF",
                "ig.GOBP",
                "ig.GOCC",
                "ig.HPPA",
                "ig.HPMI",
                "ig.HPCM",
                "ig.HPMA",
                "ig.DO",
                "ig.MP",
                "dDAGinduce",
                'dDAGreverse',
                'dDAGroot',
                'dDAGtip',
                'dDAGlevel',
                'dDAGannotate',
                'dDAGancestor',
                'dDAGtermSim',
                'dDAGgeneSim',
                "visDAG"
            )
        ),
        sd_section("Built-in databases in human",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in human (Hs; tax_id=9606).",
            c(
                "org.Hs.eg",
                "org.Hs.egDO", 
                "org.Hs.egGOBP",
                "org.Hs.egGOCC",
                "org.Hs.egGOMF",
                "org.Hs.egHPMI", 
                "org.Hs.egHPCM",
                "org.Hs.egHPMA",
                "org.Hs.egHPPA",
                "org.Hs.egMP",
                "org.Hs.egPS",
                "org.Hs.egSF",  
                "org.Hs.egDGIdb",              
                "org.Hs.string"
            )
        ),
        sd_section("Built-in genesets in human",
            "These genesets are derived from the molecular signatures database (Msigdb) in human (Hs).",
            c(
                "org.Hs.egMsigdbH",
                "org.Hs.egMsigdbC1",
                "org.Hs.egMsigdbC2CGP", 
                "org.Hs.egMsigdbC2CP",
                "org.Hs.egMsigdbC2KEGG",
                "org.Hs.egMsigdbC2REACTOME",
                "org.Hs.egMsigdbC2BIOCARTA", 
                "org.Hs.egMsigdbC3TFT",
                "org.Hs.egMsigdbC3MIR",
                "org.Hs.egMsigdbC4CGN",
                "org.Hs.egMsigdbC4CM",
                "org.Hs.egMsigdbC5BP",
                "org.Hs.egMsigdbC5MF", 
                "org.Hs.egMsigdbC5CC",
                "org.Hs.egMsigdbC6",
                "org.Hs.egMsigdbC7"
            )
        ),
        sd_section("Built-in databases in mouse",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in mouse (Mm; tax_id=10090).",
            c(
                "org.Mm.eg",
                "org.Mm.egDO", 
                "org.Mm.egGOBP",
                "org.Mm.egGOCC",
                "org.Mm.egGOMF",
                "org.Mm.egHPMI", 
                "org.Mm.egHPCM",
                "org.Mm.egHPMA",
                "org.Mm.egHPPA",
                "org.Mm.egMP",
                "org.Mm.egPS",
                "org.Mm.egSF",
                "org.Mm.string"
            )
        ),
        sd_section("Built-in databases in arabidopsis",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in arabidopsis (At; tax_id=3702).",
            c(
                "org.At.eg",
                "org.At.egGOBP",
                "org.At.egGOCC",
                "org.At.egGOMF",
                "org.At.egPS",
                "org.At.egSF",
                "org.At.string"
            )
        ),
        sd_section("Built-in databases in c.elegans",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in c.elegans (Ce; tax_id=6239).",
            c(
                "org.Ce.eg",
                "org.Ce.egGOBP",
                "org.Ce.egGOCC",
                "org.Ce.egGOMF",
                "org.Ce.egPS",
                "org.Ce.egSF",
                "org.Ce.string"
            )
        ),
        sd_section("Built-in databases in fruitfly",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in fruitfly (Dm; tax_id=7227).",
            c(
                "org.Dm.eg",
                "org.Dm.egGOBP",
                "org.Dm.egGOCC",
                "org.Dm.egGOMF",
                "org.Dm.egPS",
                "org.Dm.egSF",
                "org.Dm.string"
            )
        ),
        sd_section("Built-in databases in zebrafish",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in zebrafish (Da; tax_id=7955).",
            c(
                "org.Da.eg",
                "org.Da.egGOBP",
                "org.Da.egGOCC",
                "org.Da.egGOMF",
                "org.Da.egPS",
                "org.Da.egSF",
                "org.Da.string"
            )
        ),
        sd_section("Built-in databases in rat",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in rat (Rn; tax_id=10116).",
            c(
                "org.Rn.eg",
                "org.Rn.egGOBP",
                "org.Rn.egGOCC",
                "org.Rn.egGOMF",
                "org.Rn.egPS",
                "org.Rn.egSF",
                "org.Rn.string"
            )
        ),
        sd_section("Built-in databases in chicken",
            "These databases are used for analysis (i.e. enrichment, evolution and network analysis) in chicken (Gg; tax_id=9031).",
            c(
                "org.Gg.eg",
                "org.Gg.egGOBP",
                "org.Gg.egGOCC",
                "org.Gg.egGOMF",
                "org.Gg.egPS",
                "org.Gg.egSF",
                "org.Gg.string"
            )
        ),
        sd_section("Built-in datasets",
            "These datasets are used for the demos.",
            c(
                "Hiratani_TableS1",
                "CLL",
                "TCGA_mutations"
            )
        )
    ),
    
    if(0){
    icons = list(  
        eCal = sd_icon({
          textGrob("Common", rot = 45, gp = gpar(cex = 1))
        }),
        visRunES = sd_icon({
          textGrob("Hot", rot = 45, gp = gpar(cex = 1.2))
        }),
        eView = sd_icon(inherit = "eCal"),
        visNet = sd_icon(inherit = "visRunES")
    )
    }
)
