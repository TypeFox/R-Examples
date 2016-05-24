library(staticdocs)
list(
    readme = "",
    
    index = list(
        sd_section("Functions for analysis and visualisations",
            "These analysis and visualisation functions are used to process ontologies (and annotations), to do enrichment analysis, to calculate semantic similarity between annotated domains based on ontology term semantic similarity, and to perform random walk with restart upon domain-domain (semantic) networks. Most of analyses are supported by high-performance parallel computing.",
            c(
                'dcDAGannotate',
                'dcRDataLoader',
                'dcConverter',
                'dcEnrichment',
                'visEnrichment',
                'dcDAGdomainSim',
                'dcRWRpipeline'
            )
        ),
        sd_section("Functions for built-in (customised) data building",
            "These functions are used to build objects of S4 classes 'InfoDataFrame', 'Anno' and 'Onto' from user-customised input data (domains, ontologies and annotations).",
            c(
                'dcBuildInfoDataFrame',
                'dcBuildAnno',
                'dcBuildOnto'
            )
        ),
        sd_section("Definitions for S4 classes and methods",
            "These documentations are to help understand S4 classes and methods defined in the package.",
            c(
                'InfoDataFrame-class',
                'InfoDataFrame-method',
                'AnnoData-class',
                'Anno-class',
                'Anno-method',
                'AdjData-class',
                'Onto-class',
                'Onto-method',
                'Eoutput-class',
                'Eoutput-method',
                'Dnetwork-class',
                'Dnetwork-method',
                'Cnetwork-class',
                'Cnetwork-method',
                'Coutput-class',
                'Coutput-method'
            )
        ),
        sd_section("Ontologies mainly including open biomedical ontology (obo)",
            "These ontologies each are represented as a direct acyclic graph (DAG). DAG is stored as an object of class 'Onto'.",
            c(
                "onto.GOBP",
                "onto.GOMF",
                "onto.GOCC",
                "onto.DO",
                "onto.HPPA",
                "onto.HPMI",
                "onto.HPON",
                "onto.MP",
                "onto.EC",
                "onto.KW",
                "onto.UP"
            )
        ),
        sd_section("SCOP domain superfamilies and their annotations by ontologies",
            "These R objects are about SCOP domain superfamilies (sf) and their annotations by various ontologies, derived from the dcGO database.",
            c(
                "SCOP.sf",
                "SCOP.sf2GOBP",
                "SCOP.sf2GOMF",
                "SCOP.sf2GOCC",
                "SCOP.sf2DO",
                "SCOP.sf2HPPA",
                "SCOP.sf2HPMI",
                "SCOP.sf2HPON", 
                "SCOP.sf2MP",
                "SCOP.sf2EC",
                "SCOP.sf2KW",
                "SCOP.sf2UP"
            )
        ),
        sd_section("SCOP domain families and their annotations by ontologies",
            "These R objects are about SCOP domain families (fa) and their annotations by various ontologies, derived from the dcGO database.",
            c(
                "SCOP.fa",
                "SCOP.fa2GOBP",
                "SCOP.fa2GOMF",
                "SCOP.fa2GOCC",
                "SCOP.fa2DO",
                "SCOP.fa2HPPA",
                "SCOP.fa2HPMI",
                "SCOP.fa2HPON", 
                "SCOP.fa2MP",
                "SCOP.fa2EC",
                "SCOP.fa2KW",
                "SCOP.fa2UP"
            )
        ),
        sd_section("Pfam domains and their annotations by ontologies",
            "These R objects are about Pfam domains (Pfam) and their annotations by Gene Ontology (GO).",
            c(
                "Pfam",
                "Pfam2GOBP",
                "Pfam2GOMF",
                "Pfam2GOCC"
            )
        ),
        sd_section("InterPro domains and their annotations by ontologies",
            "These R objects are about InterPro domains (InterPro) and their annotations by Gene Ontology (GO).",
            c(
                "InterPro",
                "InterPro2GOBP",
                "InterPro2GOMF",
                "InterPro2GOCC"
            )
        ),
        sd_section("Rfam RNA families and their annotations by ontologies",
            "These R objects are about Rfam RNA families (Rfam) and their annotations by Gene Ontology (GO).",
            c(
                "Rfam",
                "Rfam2GOBP",
                "Rfam2GOMF",
                "Rfam2GOCC"
            )
        ),
        sd_section("Complete domains (domain-ome) in Eukaryotic tree of life (eTOL)",
            "These databases and functions are used for domain-centric genome analysis in Eukaryotes. Note: these domains are defined as SCOP domain superfamilies.",
            c(
                "Ancestral_domainome",
                "eTOL",
                "dcAncestralML",
                "dcAncestralMP",
                "dcSubtreeClade",
                "dcSubtreeTips",
                "dcTreeConnectivity",
                "dcDuplicated"
            )
        ),
        sd_section("Functions for domain-centric ontology creation and ontology term predictions",
            "These functions are used for creating domain-centric ontology, which in turn is used for predicting domain-centric ontology terms from input protein domain architectures.",
            c(
                "dcSplitArch",
                "dcFunArgs",
                "dcSparseMatrix",
                "dcList2Matrix",
                "dcSupraBetter",
                "dcAlgo",
                "dcAlgoPropagate",
                "dcAlgoPredict",
                "dcAlgoPredictMain",
                "dcAlgoPredictGenome",
                "dcAlgoPredictPR",
                "dcRWRpredict",
                "dcNaivePredict"
            )
        ),
        sd_section("Databases used for ontology term predictions",
            "These databases contain ontology annotations (along with hypergeometric scores) for domains and domain combinations, which are then used for domain-centric ontology term predictions.",
            c(
                "Feature2GOBP.sf",
                "Feature2GOMF.sf",
                "Feature2GOCC.sf",
                "Feature2HPPA.sf",
                "Feature2GOBP.pfam",
                "Feature2GOMF.pfam",
                "Feature2GOCC.pfam",
                "Feature2HPPA.pfam",
                "Feature2GOBP.interpro",
                "Feature2GOMF.interpro",
                "Feature2GOCC.interpro",
                "Feature2HPPA.interpro"
            )
        )

    )

)
