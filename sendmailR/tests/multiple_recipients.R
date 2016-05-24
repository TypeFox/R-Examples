library("sendmailR")

sendmail(from="from@example.org",
         to="to1@example.org",
         subject="foo",
         msg="bar",
         control=list(transport="debug"))

sendmail(from="from@example.org",
         to=c("to1@example.org", "to2@example.org"),
         subject="foo",
         msg="bar",
         control=list(transport="debug"))

sendmail(from="from@example.org",
         to=c("to1@example.org", "to2@example.org"),
         cc="cc1@example.org",
         subject="foo",
         msg="bar",
         control=list(transport="debug"))

sendmail(from="from@example.org",
         to=c("to1@example.org", "to2@example.org"),
         cc="cc1@example.org",
         subject="foo",
         msg="bar",
         control=list(transport="debug"))
         
sendmail(from="from@example.org",
         to="to1@example.org",
         cc=c("cc1@example.org", "cc2@example.org"),
         subject="foo",
         msg="bar",
         control=list(transport="debug"))

sendmail(from="from@example.org",
         to="to1@example.org",
         bcc="bcc1@example.org",
         subject="foo",
         msg="bar",
         control=list(transport="debug"))

sendmail(from="from@example.org",
         to="to1@example.org",
         bcc=c("bcc1@example.org", "bcc2@example.org"),
         subject="foo",
         msg="bar",
         control=list(transport="debug"))
