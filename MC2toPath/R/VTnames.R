VTnames <-
function (baseCalibration) 
{ if (baseCalibration=="CONUS") {
	c(
	"cold barren",
	"tundra aka alpine",
	"taiga-tundra",
	"boreal needleleaf forest",
	"boreal woodland",
	"subalpine forest",
	"maritime needleleaf forest",
	"temperate needleleaf forest",
	"temperate deciduous broadleaf forest",
	"cool mixed forest",
	"temperate warm mixed forest",
	"temperate needleleaf woodland",
	"temperate deciduous broadleaf woodland",
	"temperate cool mixed woodland",
	"temperate warm mixed woodland",
	"C3 shrubland",
	"C3 grassland",
	"temperate desert and semidesert",
	"subtropical needleleaf forest",
	"subtropical deciduous broadleaf forest",
	"warm evergreen broadleaf forest",
	"subtropical mixed forest",
	"subtropical needleleaf woodland",
	"subtropical deciduous broadleaf woodland",
	"subtropical evergreen broadleaf woodland",
	"subtropical mixed woodland",
	"C4 shrubland",
	"C4 grassland",
	"subtropical desert and semidesert",
	"tropical evergreen broadleaf forest",
	"tropical deciduous woodland",
	"tropical savanna",
	"unused33",
	"unused34",
	"tropical desert",
	"moist temperate needleleaf forest",
	"unused37",
	"subalpine meadow",
	"water and wetlands",
	"natural barren",
	"developed",
	"larch forest",
   "unused43",
   "unused44",
   "unused45",
   "unused46",
   "unused47",
   "unused48",
   "dry temperate needleleaf forest")
	}
  else if (baseCalibration=="WCR")
  {
    c(
      "cold barren",
      "tundra aka alpine",
      "taiga-tundra",
      "boreal needleleaf forest",
      "boreal woodland",
      "subalpine forest",
      "maritime needleleaf forest",
      "temperate needleleaf forest",
      "temperate deciduous broadleaf forest",
      "cool mixed forest",
      "temperate warm mixed forest",
      "temperate needleleaf woodland",
      "temperate deciduous broadleaf woodland",
      "temperate cool mixed woodland",
      "temperate warm mixed woodland",
      "C3 shrubland",
      "C3 grassland",
      "temperate desert and semidesert",
      "subtropical needleleaf forest",
      "subtropical deciduous broadleaf forest",
      "warm evergreen broadleaf forest",
      "subtropical mixed forest",
      "subtropical needleleaf woodland",
      "subtropical deciduous broadleaf woodland",
      "subtropical evergreen broadleaf woodland",
      "subtropical mixed woodland",
      "C4 shrubland",
      "C4 grassland",
      "subtropical desert and semidesert",
      "tropical evergreen broadleaf forest",
      "tropical deciduous woodland",
      "tropical savanna",
      "unused33",
      "unused34",
      "tropical desert",
      "cool needleleaf forest",
      "unused37",
      "subalpine meadow",
      "water and wetlands",
      "natural barren",
      "developed",
      "larch forest",
      "Sitka spruce zone",
      "western hemlock zone",
      "Pacific silver fir zone",
      "mountain hemlock zone",
      "subalpine fir zone",
      "subalpine parkland zone")
  }
  else if (baseCalibration=="WWETAC")
  {
    c(
      "ice", #1
      "tundra aka alpine", #2
      "taiga-tundra", #3
      "boreal EN forest", #4
      "boreal mixed woodland", #5
      "subalpine", #6
      "maritime EN forest", #7
      "temperate EN forest", #8
      "temperate deciduous broadleaf forest", #9
      "temperate cool mixed forest", #10
      "temperate warm mixed forest", #11
      "temperate EN woodland", #12
      "temperate deciduous broadleaf woodland", #13
      "temperate cool mixed woodland", #14
      "temperate warm mixed woodland", #15
      "temperate shrubland", #16
      "C3 grassland", #17
      "temperate desert", #18
      "subtropical EN forest", #19
      "subtropical deciduous broadleaf forest", #20
      "subtropical evergreen broadleaf forest", #21
      "subtropical mixed forest", #22
      "subtropical EN woodland", #23
      "subtropical deciduous broadleaf woodland", #24
      "subtropical evergreen broadleaf woodland", #25
      "subtropical mixed woodland", #26
      "subtropical shrubland", #27
      "C4 grassland", #28
      "subtropical desert", #29
      "tropical evergreen broadleaf forest", #30
      "tropical deciduous woodland", #31
      "tropical savanna", #32
      "tropical shrubland", #33
      "tropical grassland", #34
      "tropical desert", #35
      "cool needleleaf forest", #36
      "coniferous xeromorphic woodland", #37
      "subalpine meadow", #38
      "water", #39
      "natural barren", #40
      "developed") #41
  }
  else stopifnot(FALSE)
}
