#' Database of TSD information for marine turtles
#' @title Database of TSD information for marine turtles
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @author Maria Sousa Martins \email{maria.esmartins@@gmail.com}
#' @docType data
#' @name STSRE_TSD
#' @encoding UTF-8
#' @description Database of TSD information for turtles\cr
#' The columns are:\cr
#' * Species: Name de the species in binominal nommenclature\cr
#' * Country: From which country the eggs come from\cr
#' * Area: Name of the beach or region the eggs come from\cr
#' * RMU: For marine turtles, name of the RMU for this population; see Wallace, B.P., DiMatteo, A.D., Hurley, B.J., Finkbeiner, E.M., Bolten, A.B., Chaloupka, M.Y., Hutchinson, B.J., Abreu-Grobois, F.A., Amorocho, D., Bjorndal, K.A., Bourjea, J., Bowen, B.W., Duenas, R.B., Casale, P., Choudhury, B.C., Costa, A., Dutton, P.H., Fallabrino, A., Girard, A., Girondot, M., Godfrey, M.H., Hamann, M., Lopez-Mendilaharsu, M., Marcovaldi, M.A., Mortimer, J.A., Musick, J.A., Nel, R., Seminoff, J.A., Troeng, S., Witherington, B., Mast, R.B., 2010. Regional management units for marine turtles: a novel framework for prioritizing conservation and research across multiple scales. Plos One 5, e15465.\cr
#' * Incubation.temperature: Nominal incubation temperature\cr
#' * Fluctuation: How much the temperature could fluctuate around nominal temperature\cr
#' * Precision: What is the precision of the measure of temperature\cr
#' * Correction.factor: Difference between the incubator temperature and the eggs temperature\cr
#' * IP.min: Shorter incubation period\cr
#' * IP.max: Longer incubation period\cr
#' * IP.amplitude: IP.max-IP.min if available\cr
#' * IP.mean: Mean incubation periods\cr
#' * IP.mode: Mode for incubation periods\cr
#' * IP.SE: Standard error for incubation periods\cr
#' * IP.SD: Standard deviation for incubation periods\cr
#' * IP.pm: If incubation period are indicated with a plus-minus term with no more precision\cr
#' * Total: Total number of eggs incubated\cr
#' * Hatched: Number of hatchlings\cr
#' * Intersexes: Number of individuals intersexes or ambiguous for sex phenotype\cr
#' * Males: Number of individuals indentified as males\cr
#' * Females: Number of individuals indentified as females\cr
#' * Sexed: Number of sexed individuals\cr
#' * Clutch: Identity of clutch\cr
#' * Reference: Bibliographic reference\cr
#' * Note: Diverse information for this incubation\cr
#' * Version: Date of the version of this database
#' @references Binckley C.A., Spotila J.R., Wilson K.S. and Paladino F.V. (1998) Sex determination and sex ratios of Pacific Leatherback Turtles, Dermochelys coriacea. Copeia, 1998(2), 291-300.
#' @references Crastz F. (1982) Embryological stages of the marine turtle Lepidochelys olivacea (Eschscholtz). Rev. Biol. Trop., 30, 113-120.
#' @references De Souza, R.R., Vogt, R.C., 1994. Incubation temperature influences sex and hatchling size in the neotropical turtle Podocnemis unifilis. Journal of Herpetology 28, 453-464.
#' @references Girondot M. 1999 Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1(3), 479-486.
#' @references Gutzke, W. H. N., & Paukstis, G. L. (1983). Influence of the hydric environment on sexual differentiation in turtles. J. Exp. Zool., 226, 467-469. [Chrysemys picta]
#' @references Hewavisenthi, S. (1999) Influence of incubation environment on the development of the flatback turtle (Natator depressus). Ph.D. Thesis, Central Queensland University.
#' @references Hewavisenthi, S. & Parmenter, C.J. (2000) Hydric environment and sex determination in the flatback turtle (Natator depressus Garman) (Chelonia : Cheloniidae). Australian Journal of Zoology, 48, 653-659.
#' @references Hulin, V., Delmas, V., Girondot, M., Godfrey, M. H. and Guillon, J.-M. 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia, 160, 493-506.
#' @references Limpus C.J., Reed P.C. and Miller J.D. (1985) Temperature dependent sex determination in Queensland sea turtles: Intraspecific variation in Caretta caretta. In Grigg G., Shine R. and Ehmann H. (eds) Biology of Australian frogs and reptiles. New South Wales, Australia: Royal Zoological Society, pp 343-351.
#' @references Lopez Correa, J.Y., 2010. Diferenciacion gonadica en crias de Lepidochelys olivacea (Eschscholtz, 1829) (Testudinata: Cheloniidae), Instituto Politicnico Nacional, Centro Interdisciplinaria de Ciencias Marinas, La Paz, BCS, p. 108.
#' @references Marcovaldi M.A., Godfrey M.H. and Mrosovsky N. (1997) Estimating sex ratios of loggerhead turtles in Brazil from pivotal incubation temperatures. Canadian Journal of Zoology-Revue Canadienne De Zoologie, 75, 755-770.
#' @references McCoy C.J., Vogt R.C. and Censky E.J. (1983) Temperature-controlled sex determination in the sea turtle Lepidochelys olivacea. J. Herpetol., 17(4), 404-406.
#' @references Merchant-Larios, H, Diaz-Hernandez, V, Marmolejo-Valencia A. 2010. Gonadal morphogenesis and gene expression in reptiles with temperature-dependent sex determination. Sexual Development 4:50-61.
#' @references Merchant-Larios, H., Villalpando-Fierro, I. & Centeno-Urruiza, B. 1989. Gonadal morphogenesis under controlled temperature in the sea turtle Lepidochelys olivacea. Herpetol. Monographs, 3, 43-61.
#' @references Merchant-Larios, H., Ruiz-Ramirez, S., Moreno-Mendoza, N. & Marmolejo-Valencia, A. 1997. Correlation among thermosensitive period, estradiol response, and gonad differentiation in the sea turtle Lepidochelys olivacea. General and Comparative Endocrinology, 107, 373-385.
#' @references Michel-Morfinu J.E., Gomez Munoz V.M. and Navarro Rodriguez C. (2001) Morphometric model for sex assessment in hatchling olive ridley sea turtles. Chelonian Conservation and Biology, 40(1), 53-58.
#' @references Mrosovsky N., Kamel S., Rees A.F. and Margaritoulis D. (2002) Pivotal temperature for loggerhead turtles (Caretta caretta) from Kyparissia Bay, Greece. Canadian Journal of Zoology-Revue Canadienne De Zoologie, 80(12), 2118-2124.
#' @references Mrosovsky N., Dutton P.H. and Whitmore C.P. (1984) Sex ratios of two species of sea turtle nesting in Suriname. Can. J. Zool., 62, 2227-2239.
#' @references Mrosovsky N., Bass A., Corliss L.A., Richardson J.I. and Richardson T.H. (1992) Pivotal and beach temperature for hawksbill turtles nesting in Antigua. Can. J. Zool., 70, 1920-1925.
#' @references Rimblot F., Fretey J., Mrosovsky N., Lescure J. and Pieau C. (1985) Sexual differentiation as a function of the incubation temperature of eggs in the sea-turtle Dermochelys coriacea (Vandelli, 1761). Amphibia-Reptilia, 85(6), 83-92.
#' @references Rimblot-Baly F., Lescure J., Fretey J. and Pieau C. (1986-1987) Sensibilite a la temperature de la differenciation sexuelle chez la tortue Luth, Dermochelys coriacea (Vandelli, 1761); application des donnees de l'incubation artificielle a l'etude de la sex-ratio dans la nature. Annales des Sciences Naturelles, Zoologie, 8, 277-290.
#' @references Ruiz Garcia N.A. (2014) Efectos de la temperatura sobre el desarrollo embrionario y el desempero de crias de la tortuga golfina, Lepidochelys olivacea. Tesis Nivel Maestria.
#' @references Valenzuela, N., 2001. Constant, shift, and natural temperature effects on sex determination in Podocnemis expansa. Ecology 82, 3010-3024.
#' @references Tokunaga, S., Iwakiri, Y., Nakajima, Y., (1999). Temperature-dependent sex determination of a sea turtle, Caretta caretta, from Miyazaki,Japan. Bull. Kitakyushu Mus. Nat. Hist. 18, 147-156.
#' @references Wibbels T., Rostal D.C. and Byles R. (1998) High pivotal temperature in the sex determination of the olive ridley sea turtle, Lepidochelys olivacea, from Playa Nancite, Costa Rica. Copeia, 1998(4), 1086-1088.
#' @references Yntema C.L. and Mrosovsky N. (1980) Sexual differentiation in hatchling loggerheads (Caretta caretta) incubated at different controlled temperatures. Herpetologica, 36(1), 33-36.
#' @references Yntema C.L. and Mrosovsky N. (1982) Critical periods and pivotal temperatures for sexual differentiation in loggerhead sea turtles. Canadian Journal of Zoology-Revue Canadienne De Zoologie, 60(5), 1012-1016.
#' @keywords datasets
#' @family Functions for temperature-dependent sex determination
#' @usage STSRE_TSD
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(STSRE_TSD)
#' STSRE_TSD$Version[1]
#' totalIncubation_Lo <- subset(STSRE_TSD, Species=="Lepidochelys olivacea" & Sexed!=0)
#' tot_Lo <- with(totalIncubation_Lo, tsd(males=Males, females=Females, 
#'  temperatures=Incubation.temperature, par=c(P=29, S=-0.01), xlim=c(25, 35)))
#'  predict(tot_Lo)
#' }
#' @format A dataframe with raw data.
NULL
