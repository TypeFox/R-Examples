convert_from_PCS2003_4l_to_ISCO88_3d <-
function(PCS2003_4l, data){
data$ISCO88_3d[PCS2003_4l %in% c("111a")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("111b")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("111c")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("111d")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("111e")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("111f")]<-613
data$ISCO88_3d[PCS2003_4l %in% c("121a")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("121b")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("121c")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("121d")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("121e")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("121f")]<-613
data$ISCO88_3d[PCS2003_4l %in% c("122a")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("122b")]<-614
data$ISCO88_3d[PCS2003_4l %in% c("122c")]<-615
data$ISCO88_3d[PCS2003_4l %in% c("131a")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("131b")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("131c")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("131d")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("131e")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("131f")]<-613
data$ISCO88_3d[PCS2003_4l %in% c("211a")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("211b")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("211c")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("211d")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("211e")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("211f")]<-714
data$ISCO88_3d[PCS2003_4l %in% c("211g")]<-722
data$ISCO88_3d[PCS2003_4l %in% c("211h")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("211i")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("211j")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("212a")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("212b")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("212c")]<-720
data$ISCO88_3d[PCS2003_4l %in% c("212d")]<-720
data$ISCO88_3d[PCS2003_4l %in% c("213a")]<-743
data$ISCO88_3d[PCS2003_4l %in% c("214a")]<-742
data$ISCO88_3d[PCS2003_4l %in% c("214b")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("214c")]<-734
data$ISCO88_3d[PCS2003_4l %in% c("214d")]<-711
data$ISCO88_3d[PCS2003_4l %in% c("214e")]<-730
data$ISCO88_3d[PCS2003_4l %in% c("214f")]<-740
data$ISCO88_3d[PCS2003_4l %in% c("215a")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("215b")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("215c")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("215d")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("216a")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("216b")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("216c")]<-744
data$ISCO88_3d[PCS2003_4l %in% c("217a")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("217b")]<-933
data$ISCO88_3d[PCS2003_4l %in% c("217c")]<-514
data$ISCO88_3d[PCS2003_4l %in% c("217d")]<-743
data$ISCO88_3d[PCS2003_4l %in% c("217e")]<-740
data$ISCO88_3d[PCS2003_4l %in% c("218a")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("219a")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("221a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("221b")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("222a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("222b")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223b")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223c")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223d")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223e")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223f")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223g")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("223h")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("224a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("224b")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("224c")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("224d")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("225a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("226a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("226b")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("226c")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("227a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("227b")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("227c")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("227d")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("231a")]<-121
data$ISCO88_3d[PCS2003_4l %in% c("232a")]<-121
data$ISCO88_3d[PCS2003_4l %in% c("233a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("233b")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("233c")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("233d")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("311a")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("311b")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("311c")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("311d")]<-244
data$ISCO88_3d[PCS2003_4l %in% c("311e")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("311f")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("312a")]<-242
data$ISCO88_3d[PCS2003_4l %in% c("312b")]<-242
data$ISCO88_3d[PCS2003_4l %in% c("312c")]<-241
data$ISCO88_3d[PCS2003_4l %in% c("312d")]<-241
data$ISCO88_3d[PCS2003_4l %in% c("312e")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("312f")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("312g")]<-242
data$ISCO88_3d[PCS2003_4l %in% c("313a")]<-411
data$ISCO88_3d[PCS2003_4l %in% c("331a")]<-111
data$ISCO88_3d[PCS2003_4l %in% c("332a")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("332b")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("333a")]<-242
data$ISCO88_3d[PCS2003_4l %in% c("333b")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("333c")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("333d")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("333e")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("333f")]<-247
data$ISCO88_3d[PCS2003_4l %in% c("334a")]<-10
data$ISCO88_3d[PCS2003_4l %in% c("335a")]<-111
data$ISCO88_3d[PCS2003_4l %in% c("341a")]<-232
data$ISCO88_3d[PCS2003_4l %in% c("341b")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("342a")]<-231
data$ISCO88_3d[PCS2003_4l %in% c("342b")]<-231
data$ISCO88_3d[PCS2003_4l %in% c("342c")]<-231
data$ISCO88_3d[PCS2003_4l %in% c("342d")]<-231
data$ISCO88_3d[PCS2003_4l %in% c("342e")]<-200
data$ISCO88_3d[PCS2003_4l %in% c("342f")]<-200
data$ISCO88_3d[PCS2003_4l %in% c("342g")]<-200
data$ISCO88_3d[PCS2003_4l %in% c("342h")]<-200
data$ISCO88_3d[PCS2003_4l %in% c("343a")]<-244
data$ISCO88_3d[PCS2003_4l %in% c("344a")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("344b")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("344c")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("344d")]<-222
data$ISCO88_3d[PCS2003_4l %in% c("351a")]<-243
data$ISCO88_3d[PCS2003_4l %in% c("352a")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("352b")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("353a")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("353b")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("353c")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("354a")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("354b")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("354c")]<-245
data$ISCO88_3d[PCS2003_4l %in% c("354d")]<-347
data$ISCO88_3d[PCS2003_4l %in% c("354e")]<-347
data$ISCO88_3d[PCS2003_4l %in% c("354f")]<-347
data$ISCO88_3d[PCS2003_4l %in% c("354g")]<-231
data$ISCO88_3d[PCS2003_4l %in% c("371a")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("372a")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("372b")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("372c")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("372d")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("372e")]<-242
data$ISCO88_3d[PCS2003_4l %in% c("372f")]<-243
data$ISCO88_3d[PCS2003_4l %in% c("373a")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("373b")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("373c")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("373d")]<-123
data$ISCO88_3d[PCS2003_4l %in% c("374a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("374b")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("374c")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("374d")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("375a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("375b")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376b")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376c")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376d")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376e")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376f")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("376g")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("377a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("380a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("381a")]<-221
data$ISCO88_3d[PCS2003_4l %in% c("381b")]<-221
data$ISCO88_3d[PCS2003_4l %in% c("381c")]<-221
data$ISCO88_3d[PCS2003_4l %in% c("382a")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("382b")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("382c")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("382d")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("383a")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("383b")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("383c")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("384a")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("384b")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("384c")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("385a")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("385b")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("385c")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("386a")]<-211
data$ISCO88_3d[PCS2003_4l %in% c("386b")]<-211
data$ISCO88_3d[PCS2003_4l %in% c("386c")]<-211
data$ISCO88_3d[PCS2003_4l %in% c("386d")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("386e")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("387a")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("387b")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("387c")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("387d")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("387e")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("387f")]<-221
data$ISCO88_3d[PCS2003_4l %in% c("388a")]<-213
data$ISCO88_3d[PCS2003_4l %in% c("388b")]<-213
data$ISCO88_3d[PCS2003_4l %in% c("388c")]<-213
data$ISCO88_3d[PCS2003_4l %in% c("388d")]<-213
data$ISCO88_3d[PCS2003_4l %in% c("388e")]<-213
data$ISCO88_3d[PCS2003_4l %in% c("389a")]<-122
data$ISCO88_3d[PCS2003_4l %in% c("389b")]<-314
data$ISCO88_3d[PCS2003_4l %in% c("389c")]<-314
data$ISCO88_3d[PCS2003_4l %in% c("421a")]<-233
data$ISCO88_3d[PCS2003_4l %in% c("421b")]<-233
data$ISCO88_3d[PCS2003_4l %in% c("422a")]<-232
data$ISCO88_3d[PCS2003_4l %in% c("422b")]<-232
data$ISCO88_3d[PCS2003_4l %in% c("422c")]<-232
data$ISCO88_3d[PCS2003_4l %in% c("422d")]<-235
data$ISCO88_3d[PCS2003_4l %in% c("422e")]<-235
data$ISCO88_3d[PCS2003_4l %in% c("423a")]<-334
data$ISCO88_3d[PCS2003_4l %in% c("423b")]<-334
data$ISCO88_3d[PCS2003_4l %in% c("424a")]<-334
data$ISCO88_3d[PCS2003_4l %in% c("425a")]<-243
data$ISCO88_3d[PCS2003_4l %in% c("431a")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("431b")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("431c")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("431d")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("431e")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("431f")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("431g")]<-323
data$ISCO88_3d[PCS2003_4l %in% c("432a")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("432b")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("432c")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("432d")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("433a")]<-321
data$ISCO88_3d[PCS2003_4l %in% c("433b")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("433c")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("433d")]<-322
data$ISCO88_3d[PCS2003_4l %in% c("434a")]<-346
data$ISCO88_3d[PCS2003_4l %in% c("434b")]<-346
data$ISCO88_3d[PCS2003_4l %in% c("434c")]<-346
data$ISCO88_3d[PCS2003_4l %in% c("434d")]<-333
data$ISCO88_3d[PCS2003_4l %in% c("434e")]<-333
data$ISCO88_3d[PCS2003_4l %in% c("434f")]<-333
data$ISCO88_3d[PCS2003_4l %in% c("434g")]<-333
data$ISCO88_3d[PCS2003_4l %in% c("435a")]<-346
data$ISCO88_3d[PCS2003_4l %in% c("435b")]<-346
data$ISCO88_3d[PCS2003_4l %in% c("441a")]<-246
data$ISCO88_3d[PCS2003_4l %in% c("441b")]<-348
data$ISCO88_3d[PCS2003_4l %in% c("451a")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451b")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451c")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451d")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451e")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451f")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451g")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("451h")]<-344
data$ISCO88_3d[PCS2003_4l %in% c("452a")]<-345
data$ISCO88_3d[PCS2003_4l %in% c("452b")]<-10
data$ISCO88_3d[PCS2003_4l %in% c("461a")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("461b")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("461c")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("461d")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("461e")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("461f")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("462a")]<-131
data$ISCO88_3d[PCS2003_4l %in% c("462b")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("462c")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("462d")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("462e")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("463a")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("463b")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("463c")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("463d")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("463e")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("464a")]<-342
data$ISCO88_3d[PCS2003_4l %in% c("464b")]<-244
data$ISCO88_3d[PCS2003_4l %in% c("465a")]<-347
data$ISCO88_3d[PCS2003_4l %in% c("465b")]<-313
data$ISCO88_3d[PCS2003_4l %in% c("465c")]<-313
data$ISCO88_3d[PCS2003_4l %in% c("466a")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("466b")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("466c")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("467a")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("467b")]<-343
data$ISCO88_3d[PCS2003_4l %in% c("467c")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("467d")]<-341
data$ISCO88_3d[PCS2003_4l %in% c("468a")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("468b")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("471a")]<-321
data$ISCO88_3d[PCS2003_4l %in% c("471b")]<-321
data$ISCO88_3d[PCS2003_4l %in% c("472a")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("472b")]<-214
data$ISCO88_3d[PCS2003_4l %in% c("472c")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("472d")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("473a")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("473b")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("473c")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("474a")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("474b")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("474c")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("475a")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("475b")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("476a")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("476b")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("477a")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("477b")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("477c")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("477d")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("478a")]<-312
data$ISCO88_3d[PCS2003_4l %in% c("478b")]<-312
data$ISCO88_3d[PCS2003_4l %in% c("478c")]<-312
data$ISCO88_3d[PCS2003_4l %in% c("478d")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("479a")]<-321
data$ISCO88_3d[PCS2003_4l %in% c("479b")]<-311
data$ISCO88_3d[PCS2003_4l %in% c("480a")]<-613
data$ISCO88_3d[PCS2003_4l %in% c("480b")]<-615
data$ISCO88_3d[PCS2003_4l %in% c("481a")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("481b")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("482a")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("483a")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("484a")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("484b")]<-812
data$ISCO88_3d[PCS2003_4l %in% c("485a")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("485b")]<-825
data$ISCO88_3d[PCS2003_4l %in% c("486a")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("486b")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("486c")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("486d")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("486e")]<-700
data$ISCO88_3d[PCS2003_4l %in% c("487a")]<-413
data$ISCO88_3d[PCS2003_4l %in% c("487b")]<-933
data$ISCO88_3d[PCS2003_4l %in% c("488a")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("488b")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("521a")]<-414
data$ISCO88_3d[PCS2003_4l %in% c("521b")]<-414
data$ISCO88_3d[PCS2003_4l %in% c("522a")]<-412
data$ISCO88_3d[PCS2003_4l %in% c("523a")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("523b")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("523c")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("523d")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("524a")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("524b")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("524c")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("524d")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("525a")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("525b")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("525c")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("525d")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("526a")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("526b")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("526c")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("526d")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("526e")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("531a")]<-516
data$ISCO88_3d[PCS2003_4l %in% c("531b")]<-516
data$ISCO88_3d[PCS2003_4l %in% c("531c")]<-516
data$ISCO88_3d[PCS2003_4l %in% c("532a")]<-10
data$ISCO88_3d[PCS2003_4l %in% c("532b")]<-10
data$ISCO88_3d[PCS2003_4l %in% c("532c")]<-10
data$ISCO88_3d[PCS2003_4l %in% c("533a")]<-10
data$ISCO88_3d[PCS2003_4l %in% c("533b")]<-915
data$ISCO88_3d[PCS2003_4l %in% c("533c")]<-915
data$ISCO88_3d[PCS2003_4l %in% c("534a")]<-915
data$ISCO88_3d[PCS2003_4l %in% c("534b")]<-915
data$ISCO88_3d[PCS2003_4l %in% c("541a")]<-422
data$ISCO88_3d[PCS2003_4l %in% c("541b")]<-422
data$ISCO88_3d[PCS2003_4l %in% c("541c")]<-422
data$ISCO88_3d[PCS2003_4l %in% c("541d")]<-422
data$ISCO88_3d[PCS2003_4l %in% c("542a")]<-411
data$ISCO88_3d[PCS2003_4l %in% c("542b")]<-411
data$ISCO88_3d[PCS2003_4l %in% c("543a")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("543b")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("543c")]<-410
data$ISCO88_3d[PCS2003_4l %in% c("543d")]<-419
data$ISCO88_3d[PCS2003_4l %in% c("543e")]<-419
data$ISCO88_3d[PCS2003_4l %in% c("543f")]<-419
data$ISCO88_3d[PCS2003_4l %in% c("543g")]<-419
data$ISCO88_3d[PCS2003_4l %in% c("543h")]<-419
data$ISCO88_3d[PCS2003_4l %in% c("544a")]<-411
data$ISCO88_3d[PCS2003_4l %in% c("545a")]<-421
data$ISCO88_3d[PCS2003_4l %in% c("545b")]<-421
data$ISCO88_3d[PCS2003_4l %in% c("545c")]<-421
data$ISCO88_3d[PCS2003_4l %in% c("545d")]<-421
data$ISCO88_3d[PCS2003_4l %in% c("546a")]<-511
data$ISCO88_3d[PCS2003_4l %in% c("546b")]<-413
data$ISCO88_3d[PCS2003_4l %in% c("546c")]<-413
data$ISCO88_3d[PCS2003_4l %in% c("546d")]<-511
data$ISCO88_3d[PCS2003_4l %in% c("546e")]<-511
data$ISCO88_3d[PCS2003_4l %in% c("551a")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("552a")]<-421
data$ISCO88_3d[PCS2003_4l %in% c("553a")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("553b")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("553c")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554a")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554b")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554c")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554d")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554e")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554f")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554g")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554h")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("554j")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("555a")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("556a")]<-522
data$ISCO88_3d[PCS2003_4l %in% c("561a")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("561b")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("561c")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("561d")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("561e")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("561f")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("562a")]<-514
data$ISCO88_3d[PCS2003_4l %in% c("562b")]<-514
data$ISCO88_3d[PCS2003_4l %in% c("563a")]<-513
data$ISCO88_3d[PCS2003_4l %in% c("563b")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("563c")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("564a")]<-915
data$ISCO88_3d[PCS2003_4l %in% c("564b")]<-910
data$ISCO88_3d[PCS2003_4l %in% c("621a")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("621b")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("621c")]<-833
data$ISCO88_3d[PCS2003_4l %in% c("621d")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("621e")]<-710
data$ISCO88_3d[PCS2003_4l %in% c("621f")]<-710
data$ISCO88_3d[PCS2003_4l %in% c("621g")]<-811
data$ISCO88_3d[PCS2003_4l %in% c("622a")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("622b")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("622c")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("622d")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("622e")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("622f")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("622g")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("623a")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("623b")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("623c")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("623d")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("623e")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("623f")]<-722
data$ISCO88_3d[PCS2003_4l %in% c("623g")]<-722
data$ISCO88_3d[PCS2003_4l %in% c("624a")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("624b")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("624c")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("624d")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("624e")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("624f")]<-812
data$ISCO88_3d[PCS2003_4l %in% c("624g")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("625a")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("625b")]<-740
data$ISCO88_3d[PCS2003_4l %in% c("625c")]<-815
data$ISCO88_3d[PCS2003_4l %in% c("625d")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("625e")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("625f")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("625g")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("625h")]<-816
data$ISCO88_3d[PCS2003_4l %in% c("626a")]<-812
data$ISCO88_3d[PCS2003_4l %in% c("626b")]<-812
data$ISCO88_3d[PCS2003_4l %in% c("626c")]<-814
data$ISCO88_3d[PCS2003_4l %in% c("627a")]<-826
data$ISCO88_3d[PCS2003_4l %in% c("627b")]<-826
data$ISCO88_3d[PCS2003_4l %in% c("627c")]<-826
data$ISCO88_3d[PCS2003_4l %in% c("627d")]<-814
data$ISCO88_3d[PCS2003_4l %in% c("627e")]<-734
data$ISCO88_3d[PCS2003_4l %in% c("627f")]<-825
data$ISCO88_3d[PCS2003_4l %in% c("628a")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("628b")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("628c")]<-722
data$ISCO88_3d[PCS2003_4l %in% c("628d")]<-722
data$ISCO88_3d[PCS2003_4l %in% c("628e")]<-916
data$ISCO88_3d[PCS2003_4l %in% c("628f")]<-700
data$ISCO88_3d[PCS2003_4l %in% c("628g")]<-700
data$ISCO88_3d[PCS2003_4l %in% c("631a")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("632a")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("632b")]<-711
data$ISCO88_3d[PCS2003_4l %in% c("632c")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("632d")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("632e")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("632f")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("632g")]<-714
data$ISCO88_3d[PCS2003_4l %in% c("632h")]<-714
data$ISCO88_3d[PCS2003_4l %in% c("632j")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("632k")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("633a")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("633b")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("633c")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("633d")]<-724
data$ISCO88_3d[PCS2003_4l %in% c("634a")]<-721
data$ISCO88_3d[PCS2003_4l %in% c("634b")]<-722
data$ISCO88_3d[PCS2003_4l %in% c("634c")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("634d")]<-731
data$ISCO88_3d[PCS2003_4l %in% c("635a")]<-743
data$ISCO88_3d[PCS2003_4l %in% c("636a")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("636b")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("636c")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("636d")]<-512
data$ISCO88_3d[PCS2003_4l %in% c("637a")]<-732
data$ISCO88_3d[PCS2003_4l %in% c("637b")]<-700
data$ISCO88_3d[PCS2003_4l %in% c("637c")]<-743
data$ISCO88_3d[PCS2003_4l %in% c("637d")]<-700
data$ISCO88_3d[PCS2003_4l %in% c("641a")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("641b")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("642a")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("642b")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("643a")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("644a")]<-832
data$ISCO88_3d[PCS2003_4l %in% c("651a")]<-833
data$ISCO88_3d[PCS2003_4l %in% c("651b")]<-833
data$ISCO88_3d[PCS2003_4l %in% c("652a")]<-833
data$ISCO88_3d[PCS2003_4l %in% c("652b")]<-833
data$ISCO88_3d[PCS2003_4l %in% c("653a")]<-413
data$ISCO88_3d[PCS2003_4l %in% c("654a")]<-831
data$ISCO88_3d[PCS2003_4l %in% c("654b")]<-831
data$ISCO88_3d[PCS2003_4l %in% c("654c")]<-831
data$ISCO88_3d[PCS2003_4l %in% c("655a")]<-831
data$ISCO88_3d[PCS2003_4l %in% c("656a")]<-834
data$ISCO88_3d[PCS2003_4l %in% c("656b")]<-834
data$ISCO88_3d[PCS2003_4l %in% c("656c")]<-834
data$ISCO88_3d[PCS2003_4l %in% c("671a")]<-931
data$ISCO88_3d[PCS2003_4l %in% c("671b")]<-931
data$ISCO88_3d[PCS2003_4l %in% c("671c")]<-931
data$ISCO88_3d[PCS2003_4l %in% c("671d")]<-931
data$ISCO88_3d[PCS2003_4l %in% c("672a")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("673a")]<-821
data$ISCO88_3d[PCS2003_4l %in% c("673b")]<-821
data$ISCO88_3d[PCS2003_4l %in% c("673c")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("674a")]<-815
data$ISCO88_3d[PCS2003_4l %in% c("674b")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("674c")]<-827
data$ISCO88_3d[PCS2003_4l %in% c("674d")]<-812
data$ISCO88_3d[PCS2003_4l %in% c("674e")]<-820
data$ISCO88_3d[PCS2003_4l %in% c("675a")]<-826
data$ISCO88_3d[PCS2003_4l %in% c("675b")]<-828
data$ISCO88_3d[PCS2003_4l %in% c("675c")]<-800
data$ISCO88_3d[PCS2003_4l %in% c("676a")]<-933
data$ISCO88_3d[PCS2003_4l %in% c("676b")]<-933
data$ISCO88_3d[PCS2003_4l %in% c("676c")]<-932
data$ISCO88_3d[PCS2003_4l %in% c("676d")]<-933
data$ISCO88_3d[PCS2003_4l %in% c("676e")]<-800
data$ISCO88_3d[PCS2003_4l %in% c("681a")]<-712
data$ISCO88_3d[PCS2003_4l %in% c("681b")]<-713
data$ISCO88_3d[PCS2003_4l %in% c("682a")]<-723
data$ISCO88_3d[PCS2003_4l %in% c("683a")]<-741
data$ISCO88_3d[PCS2003_4l %in% c("684a")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("684b")]<-913
data$ISCO88_3d[PCS2003_4l %in% c("685a")]<-700
data$ISCO88_3d[PCS2003_4l %in% c("691a")]<-833
data$ISCO88_3d[PCS2003_4l %in% c("691b")]<-612
data$ISCO88_3d[PCS2003_4l %in% c("691c")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("691d")]<-611
data$ISCO88_3d[PCS2003_4l %in% c("691e")]<-613
data$ISCO88_3d[PCS2003_4l %in% c("691f")]<-614
data$ISCO88_3d[PCS2003_4l %in% c("692a")]<-615
return(data)
}
