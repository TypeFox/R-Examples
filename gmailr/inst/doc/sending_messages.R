## ---- include=FALSE------------------------------------------------------
library(gmailr)

## ----sending_messages_simple---------------------------------------------
mime() %>%
  to("james.f.hester@gmail.com") %>%
  from("me@somewhere.com") %>%
  text_body("Gmailr is a very handy package!") -> text_msg

## ----sending_messages_simple_print---------------------------------------
strwrap(as.character(text_msg))

## ----sending_messages_html-----------------------------------------------
mime() %>%
  to("james.f.hester@gmail.com") %>%
  from("me@somewhere.com") %>%
  html_body("<b>Gmailr</b> is a <i>very</i> handy package!") -> html_msg

## ----sending_messages_attachments_2--------------------------------------
write.csv(file = "iris.csv", iris)

html_msg %>%
  subject("Here are some flowers") %>%
  attach_file("iris.csv") -> file_attachment

## ----sending_messages_attachments_1--------------------------------------
html_msg %>% attach_part(part = charToRaw("attach me!"), name = "please") -> simple_attachment

## ----sending_messages_create_draft, eval=FALSE---------------------------
#  create_draft(file_attachment)

## ----sending_messages_insert_message, eval=FALSE-------------------------
#  insert_message(file_attachment)

## ----sending_messages_file_attachment, eval=FALSE------------------------
#  insert_message(file_attachment)

## ----sending_messages_send_draft, eval=FALSE-----------------------------
#  my_drafts <- drafts()
#  
#  send_draft(id(my_drafts, "draft_id")[1])

## ----sending_messages_send_message, eval=FALSE---------------------------
#  send_message(file_attachment)

## ----sending_messages_clenup, include=FALSE------------------------------
unlink("iris.csv")

