is.integer(ID)
!duplicated(ID)
ID > 0 & ID < 1754

is.numeric(Latitude)
Latitude < 0

is.numeric(Longitude)
Longitude < 180 & Longitude > -180

is.null(Longitude) == is.null(Latitude)

is.character(Adm1)
is.character(Adm2)
is.character(Adm3)
is.character(Country)
is.integer(Altitude)

is.null(Adm1) == is.null(Longitude)  # 
is.null(Adm2) == is.null(Longitude)
is.null(Adm3) == is.null(Longitude)
is.null(Country) == is.null(Longitude)
is.null(Altitude) == is.null(Longitude)

is.numeric(pH)
# sapply(pH, is.withinRange,0,14)
pH >= 0  # pH bigger than
pH <= 14  # pH lesser than



is.numeric(Conductivity)
Conductivity >= 0

is.numeric(CaCO3)
CaCO3 >= 0

is.numeric(Sand)
is_within_range(Sand, 0, 100)

is.numeric(Lime)
is_within_range(Lime, 0, 100)

is.numeric(Clay)
is_within_range(Clay, 0, 100)

is.character(Soil_texture)
is_one_of(Soil_texture, c("Loam", "Clay Loam", "Sandy Loam", "Sandy Clay Loam", "Clay", "Sand", "Loamy Sand", "Sandy Clay", "Silt Loam", "Silty Clay Loam", "Silty Clay", "Organic Soil"))

is.numeric(P)
P >= 0 
