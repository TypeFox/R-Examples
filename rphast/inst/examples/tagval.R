tags <- c("tag1 \"val 1a\"; tag2 \"val 2a\" \"val2a.1\" 123; tag3 \"val3a\"",
          "tag1 \"val 1b\"; tag2 \"val 2b\"; tag4 \"val4b\"",
          "tag3 \"val3a\" 1; tag4 2;")
tagval(tags, "tag1")
tagval(tags, "tag2")
tagval(tags, "tag3")
tagval(tags, "tag4")
tagval(tags, "notag")
rm(tags)
