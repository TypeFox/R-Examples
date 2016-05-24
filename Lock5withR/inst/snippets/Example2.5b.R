tally(Response ~ Gender, data = OneTrueLove)
tally( ~ Response | Gender, data = OneTrueLove)
tally(Gender ~ Response, data = OneTrueLove)
tally( ~ Gender | Response, data = OneTrueLove)

