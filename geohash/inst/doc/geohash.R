## ---- eval=FALSE---------------------------------------------------------
#  # Encode!
#  gh_encode(42.60498046875, -5.60302734375, 5)
#  # [1] "ezs42"

## ---- eval=FALSE---------------------------------------------------------
#  # Decode!
#  gh_decode("ezs42")
#  # lat      lng      lat_error  lng_error
#  # 42.60498 42.60498 0.02197266 0.02197266

## ---- eval=FALSE---------------------------------------------------------
#  # Get all neighbours!
#  gh_neighbours("ezs42")
#  #  north northeast  east southeast south southwest  west northwest
#  #  ezs48     ezs49 ezs43     ezs41 ezs40     ezefp ezefr     ezefx

## ---- eval=FALSE---------------------------------------------------------
#  # Get the southeast neighbours!
#  southeast("ezs42")
#  # [1] "ezs41"

