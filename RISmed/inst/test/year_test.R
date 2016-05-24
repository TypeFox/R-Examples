library(RISmed)

query <- EUtilsSummary("Duke[Author] AND Kakefuda[Author] AND 1981[Date - Publication] AND 67[Volume] AND 449[Pagination]")

query@count
query@querytranslation

query_result <- EUtilsGet(query)

Author(query_result)
Year(query_result)