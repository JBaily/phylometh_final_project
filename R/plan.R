plan <- drake_plan(
  
  #Getting a tree and data
  #Tree is constructed with 16S sequences from my samples and relative 
  #abundance data is also from my samples. 
  
  #The tree was constructed and data aligned using the SILVA SINA program, the 
  #reference for which is in "data/arb-silva/SILVA_SINA_documentation.txt", 
  #for inquiring minds. 

  tree = read.tree("data/project.tree"),
  tree_test = plot.phylo(tree, cex = 0.4),
  tree_image = plot_tree(file=file_out("results/initial_tree.pdf"), tree, 
                         type = "phylogram", label = TRUE, cex = 0.4),
  data.test.1 = read.csv(file="data/phylometh_project_data_test.csv", sep = ",",
                           stringsAsFactors=FALSE, row.names = 1, header= TRUE),
  data.test.2 = convert.test.data(data.test.1), 
  data.test.sd = apply(data.test.2, 2, sd),
  print(data.test.sd),
  data.test.sd.ordered = sort.int(data.test.sd, decreasing=TRUE),
  print(data.test.sd.ordered),
  which.values.comp = c(names(data.test.sd.ordered)[1], 
                        names(data.test.sd.ordered)[2]),
  best.values.comp = c(data.test.sd.ordered[1], data.test.sd.ordered[2]),
  print(which.values.comp),
  print(best.values.comp),
  
  to.keep = c(which.values.comp[1],which.values.comp[2]), 
  data.test.3 = data.test.2[,to.keep],
  data.test.4 = assign.states(data.test.3), 
  
  #Cleaning the data -- making sure that taxa match and there is no missing data
  cleaned.test = CleanData(tree, data.test.4$state.data), 
  
  #Making ancestral reconstructions with phangorn
  cleaned.test.phyDat = phangorn::phyDat(cleaned.test$data, type="USER", 
                                         levels = c("11","12","13","21","22",
                                                    "23","31","32","33")), 
  anc.p.test = phangorn::ancestral.pars(cleaned.test$phy, 
                                   cleaned.test.phyDat),
  
  plot.test.p = print_plotAnc(cleaned.test$phy, anc.p.test, 
                                type="test_pars_anc_recon",cex = 0.4),
  anc.p.test.results = do.call(rbind, anc.p.test),
  print.results(anc.p.test.results, "anc_p_test"),
  
  anc.ml.test = phangorn::ancestral.pml(phangorn::pml(cleaned.test$phy, 
                                                 cleaned.test.phyDat), 
                                        type="ml"),
  plot.test.ml = print_plotAnc(cleaned.test$phy, anc.ml.test, 
                               type="test_ml_anc_recon", cex=0.4),
  anc.ml.test.results = do.call(rbind, anc.ml.test),
  print.results(anc.ml.test.results, "anc_ml_test"),
  
  
  #######################################
  ## Now with hand-picked time points. ## 
  #######################################
  to.keep.pick = c("Value.4", "Value.6"), 
  data.pick.3 = data.test.2[,to.keep.pick],
  data.pick.4 = assign.states(data.pick.3), 
  
  cleaned.pick = CleanData(tree, data.pick.4$state.data), 
  
  cleaned.pick.phyDat = phangorn::phyDat(cleaned.pick$data, type="USER", 
                                         levels = c("11","12","13","21","22",
                                                    "23","31","32","33")), 
  anc.p.pick = phangorn::ancestral.pars(cleaned.pick$phy, 
                                        cleaned.pick.phyDat),
  plot.pick.p = print_plotAnc(cleaned.pick$phy, anc.p.pick, 
                              type="pick_pars_anc_recon",cex = 0.4),
  anc.p.pick.results = do.call(rbind, anc.p.pick),
  print.results(anc.p.pick.results, "anc_p_pick"),
  
  anc.ml.pick = phangorn::ancestral.pml(phangorn::pml(cleaned.pick$phy, 
                                                      cleaned.pick.phyDat), 
                                        type="ml"),
  plot.pick.ml = print_plotAnc(cleaned.pick$phy, anc.ml.pick, 
                               type="pick_ml_anc_recon", cex=0.4),
  
  anc.ml.pick.results = do.call(rbind, anc.ml.pick),
  print.results(anc.ml.pick.results, "anc_ml_pick"),
  
  ############################################################################
  
  ###########################################################
  ### Running ancestral reconstructions with control data ###
  ###########################################################
  
  data.control.1 = read.csv(file="data/phylometh_project_data_comp.csv", sep = ",",
                         stringsAsFactors=FALSE, header= TRUE),
  print(data.control.1),
  data.control = assign.control(data.control.1), 
  
  #Cleaning the data -- making sure that taxa match and there is no missing data
  cleaned.control = CleanData(tree, data.control), 
  
  #Making ancestral reconstructions with phangorn
  cleaned.control.phyDat = phangorn::phyDat(cleaned.control$data, type="USER", 
                                         levels = c("GEN","NF","AO","SR","ME",
                                                    "MY")), 
  
  anc.p.control = phangorn::ancestral.pars(cleaned.control$phy, 
                                        cleaned.control.phyDat),
  plot.control.p = print_plotAnc(cleaned.control$phy, anc.p.control, 
                              type="control_pars_anc_recon",cex = 0.4),
  anc.p.control.results = do.call(rbind, anc.p.control),
  print.results(anc.p.control.results, "anc_p_control"),
  
  
  anc.ml.control = phangorn::ancestral.pml(phangorn::pml(cleaned.control$phy, 
                                                      cleaned.control.phyDat), 
                                        type="ml"),
  plot.control.ml = print_plotAnc(cleaned.control$phy, anc.ml.control, 
                               type="control_ml_anc_recon", cex=0.4),
  anc.ml.control.results = do.call(rbind, anc.ml.test),
  print.results(anc.ml.control.results, "anc_ml_control"),
  
  ############################################################################
  
  ######################################
  ### Comparing the different models ###
  ######################################
  
  #I will be using the Maximum likelihood trees for the rest of this analysis, 
  #though I believe this would work for the parsimony trees as well. 
  
  #Assigning state values--by default, the number of possible states in the 
  #test model often be greater than the number of metabolic classes of interest.
  #This means that some of the states must be assigned to the same metabolic 
  #class. Furthermore, the states must be assigned to the "best" representative
  #metabolic class that they match to. I do not have the slightest clue as to 
  #the most efficient way to make R do this to any degree of accuracy or 
  #meaningfulness, so I have done it manually by comparing the two ancestral 
  #reconstructions and seeing which color groups match the best. Let it be 
  #known, however, that this is a part of the method that can absolutely be 
  #improved upon in the future, by me or by someone else. 
  
  #GEN - 11, 22, 23, 31. NF - 13, AO - 32, SR - 33, ME - 21, & MY - 12. 
  
  #########################################
  ### Control vs. Machine-picked (test) ###
  #########################################
  
  control.2.mat = convert.control(anc.ml.control,6,
                                  c("Node","GEN","NF","AO","SR","ME","MY")),
  test.2.mat = convert.test.2(anc.ml.test,6,
                              c("Node","GEN","NF","AO","SR","ME","MY")),
  
  comp.2.OTU = convert_to_OTU_table(control.2.mat, test.2.mat),
  comp.2.samp = create_sample_matrix(comp.2.OTU, "control", "test.2"),
  
  comp.2.OTU.t = otu_table(comp.2.OTU, taxa_are_rows = TRUE),
  comp.2.samp.t = sample_data(as.data.frame(comp.2.samp)),
  
  ps.object.2 = phyloseq(comp.2.OTU.t, comp.2.samp.t),
  
  ord.nmds.bray = ordinate(ps.object.2, method = "NMDS", distance = "bray", trymax = 1000),
  nmds.bay.plot.2 = plot_ordination(ps.object.2, ord.nmds.bray, shape = "method", 
                                    color = "bacteria", title="Bray NMDS"),
  plot(nmds.bay.plot.2),
  
  ### Looking at the NMDS plot, it becomes quite apparent that nodes rarely 
  ### match up across the "method" types, so it is unlikely that the test
  ### reconstruction is a good approximation of the control one, but we will 
  ### statistically test that with a PERMANOVA with the adonis() function
  ### to be sure. 
  
  distance.bray.2 = phyloseq::distance(ps.object.2, method = "bray"),
  sample.dataframe = data.frame(sample_data(ps.object.2)),
  
  permanova.2 = adonis(distance.bray.2 ~ method, data = sample.dataframe),
  print(permanova.2),
  print.results(permanova.2, "permanova_two_sites"),
  
  ### P-value < 0.001, which means that there is a significant difference 
  ### between the two reconstructions. Thus, this particular method of using 
  ### relative abundance as a proxy for general metabolic classes is not viable.

  
  ############################################################################
  ############################################################################
  
  
  #############################   Extras   ####################################
  
  #############################################################################
  ### Comparing the all two point models (control, "test", and hand-picked) ###
  #############################################################################
  
  pick.2.mat = convert.pick.2(anc.ml.pick,6,
                              c("Node","GEN","NF","AO","SR","ME","MY")),
  
  comp.2.OTU.all = cbind(control.2.mat, test.2.mat, pick.2.mat),
  comp.2.samp.all = create_sample_matrix_all(comp.2.OTU.all, "control", "test.2", "pick.2"),
  
  comp.2.OTU.t.all = otu_table(comp.2.OTU.all, taxa_are_rows = TRUE),
  comp.2.samp.t.all = sample_data(as.data.frame(comp.2.samp.all)),
  
  ps.object.2.all = phyloseq(comp.2.OTU.t.all, comp.2.samp.t.all),
  
  ord.nmds.bray.all = ordinate(ps.object.2.all, method = "NMDS", 
                           distance = "bray", trymax = 1000),
  nmds.bay.plot.2.all = plot_ordination(ps.object.2.all, ord.nmds.bray.all, 
                                    shape = "method", color = "bacteria", 
                                    title="Bray NMDS"),
  plot(nmds.bay.plot.2.all),
  
  comp.ct = subset_samples(ps.object.2.all, method=="control" | method=="test.2"), #comparing control and test - same as earlier
  comp.cp = subset_samples(ps.object.2.all, method=="control" | method=="pick.2"), #comparing control and pick 
  comp.pt = subset_samples(ps.object.2.all, method=="test.2" | method=="pick.2"), #comparing pick and test, just for fun. 
  
  distance.bray.2.all = phyloseq::distance(ps.object.2.all, method = "bray"),
  sample.dataframe.all = data.frame(sample_data(ps.object.2.all)),
  
  distance.bray.2.ct = phyloseq::distance(comp.ct, method = "bray"),
  sample.dataframe.ct = data.frame(sample_data(comp.ct)),
  
  distance.bray.2.cp = phyloseq::distance(comp.cp, method = "bray"),
  sample.dataframe.cp = data.frame(sample_data(comp.cp)),
  
  distance.bray.2.pt = phyloseq::distance(comp.pt, method = "bray"),
  sample.dataframe.pt = data.frame(sample_data(comp.pt)),
  
  permanova.2.all = adonis(distance.bray.2.all ~ method, data = sample.dataframe.all),
  permanova.2.ct = adonis(distance.bray.2.ct ~ method, data = sample.dataframe.ct),
  permanova.2.cp = adonis(distance.bray.2.cp ~ method, data = sample.dataframe.cp),
  permanova.2.pt = adonis(distance.bray.2.pt ~ method, data = sample.dataframe.pt),
  
  print(permanova.2.all),
  print(permanova.2.ct),
  print(permanova.2.cp),
  print(permanova.2.pt),
  
  print.results(permanova.2.all, "permanova_two_sites_all"),
  print.results(permanova.2.ct, "permanova_two_sites_ct"),
  print.results(permanova.2.cp, "permanova_two_sites_cp"),
  print.results(permanova.2.pt, "permanova_two_sites_pt"),
  
  ### In a surprising turn of events, the p-value for the hand-picked time 
  ### points and the control is actually much higher than that of the computer
  ### picked time points (0.093 vs. 0.001). While I would absolutely not call 
  ### the hand-picked time points a particularly good reconstruction as 
  ### compared to the control model, I suppose I would call it serviceable? The
  ### PERMANOVA test, at least, does not think the two models are significantly
  ### different, at least. This does seem to speak to the complexity of chosing
  ### time-points, however, as selecting for the time points with the most 
  ### variability does not translate to meaningful data, but just eyeballing it
  ### does. Given that I noticed the patterns by eyeballing, that does check
  ### out. HOWEVER, I nearly completely discarded the hand-picked time points, 
  ### as the ancestral reconstruction, did not "look as good" as the machine
  ### -picked one, which serves as a caution against not confirming your 
  ### suspicions with actual statistical analysis. 
  
  
  
  #######################################################
  ### Trying the procedure with three time-point data ###
  #######################################################
  
  ### Computer-selected time points###
  
  print(data.test.sd),
  print(data.test.sd.ordered),
  which.values.comp.3 = c(names(data.test.sd.ordered)[1], 
                        names(data.test.sd.ordered)[2],
                        names(data.test.sd.ordered)[3]),
  best.values.comp.3 = c(data.test.sd.ordered[1], data.test.sd.ordered[2], 
                         data.test.sd.ordered[3]),
  print(which.values.comp.3),
  print(best.values.comp.3),
  
  to.keep.3 = c(which.values.comp.3[1],which.values.comp.3[2],
              which.values.comp.3[3]), 
  data.test.3.3 = data.test.2[,to.keep.3],
  print(data.test.3.3),
  data.test.4.3 = assign.states.3(data.test.3.3), 
  
  all.states = c(111,112,113,121,122,123,131,132,133,211,212,213,
                 221,222,223,231,232,233,311,312,313,321,322,323,
                 331,332,333),
  missing.states = data.test.4.3$missing,
  present.states = all.states[! all.states %in% missing.states],
  print(present.states),
  print(missing.states),
  
  cleaned.test.3 = CleanData(tree, data.test.4.3$state.data), 
  
  cleaned.test.3.phyDat = phangorn::phyDat(cleaned.test.3$data, type="USER", 
                                         levels = as.character(present.states)), 

  anc.ml.test.3 = phangorn::ancestral.pml(phangorn::pml(cleaned.test.3$phy, 
                                                      cleaned.test.3.phyDat), 
                                        type="ml"),
  plot.test.ml.3 = print_plotAnc(cleaned.test.3$phy, anc.ml.test.3, 
                               type="test_ml_anc_recon.three_states", cex=0.4),
  anc.ml.test.results.3 = do.call(rbind, anc.ml.test.3),
  print.results(anc.ml.test.results.3, "anc_ml_test.three_states"),
  
  test.3.mat = convert.test.3(anc.ml.test.3,6,
                              c("Node","GEN","NF","AO","SR","ME","MY")),
  
  comp.3.OTU = convert_to_OTU_table(control.2.mat, test.3.mat),
  comp.3.samp = create_sample_matrix(comp.3.OTU, "control", "test.triple"),
  
  comp.3.OTU.t = otu_table(comp.3.OTU, taxa_are_rows = TRUE),
  comp.3.samp.t = sample_data(as.data.frame(comp.3.samp)),
  
  ps.object.3 = phyloseq(comp.3.OTU.t, comp.3.samp.t),
  
  ord.nmds.bray.3 = ordinate(ps.object.3, method = "NMDS", distance = "bray", 
                             trymax = 500),
  nmds.bay.plot.3 = plot_ordination(ps.object.3, ord.nmds.bray.3, shape = "method", 
                                    color = "bacteria", title="Bray NMDS"),
  plot(nmds.bay.plot.3),

  ### You can already see from the NMDS plot that there are a lot more overlaps
  ### between nodes from the control and the triple method. They also tend to 
  ### approximate each other much better this time around as well. 
  
  distance.bray.3 = phyloseq::distance(ps.object.3, method = "bray"),
  sample.dataframe.3 = data.frame(sample_data(ps.object.3)),
  
  permanova.3 = adonis(distance.bray.3 ~ method, data = sample.dataframe.3),
  print(permanova.3),
  print.results(permanova.3, "permanova_three_sites"),
  
  ### With a p-value of 0.089, this model is slightly worse than the two site
  ### hand-picked method, which I find really interesting. The pattern of this
  ### NMDS prompted me to go back and look at the "all" NMDS and I was found
  ### that the pick and control nodes also tended to coincide about as 
  ### frequently as the three site computer model did, which I had not picked 
  ### up on the first time. It's a bit mind-boggling that just a bit of human
  ### consideration can bump up the quality of this model so much--have an 8
  ### state model be more accurate (presumably) than a 27 state model--and I do 
  ### not even know if I picked the optimal points! All of this seems to point
  ### to instituting a brute-force method across two, three, and maybe even 
  ### four site models (I suspect model quality will act a bit like a horseshoe
  ### here--more resolution of the occasionally erratic patterns might confuse
  ### the method) to find the optimal model. The only problem with this is that 
  ### I still have no idea how to program R to make the match judgements that I 
  ### do when comparing the various test reconstructions to the control. I do 
  ### not think that is particularly possible within in the span of a day, 
  ### however, but it is definitely something to think about in the future--
  ### outlining the specific thought process and decision tree that I go through 
  ### to make those assignments, that is. I also think increasing the number of
  ### sequences in the tree might help a tad, too, but that's a secondary 
  ### concern, I believe. I'd also want to consolidate the functions that only 
  ### differ in the size of the input matrix that they handle, like 
  ### convert.test.2/convert.test.3 and the "create_sample_matrix", for
  ### instance. And the ones where I consolidate the values in the test models
  ### to align with the control model, but I think that one will be much easier,
  ### as it's just a matter of adjusting the inputs of the functions, not 
  ### broadly abstracting the function as with the previous group of functions. 
  
  ### Overall, this was pretty fun and I've learned so much, both from this 
  ### project and absolutely from the class as a whole. This has to be
  ### one of the most informative classes that I have ever taken, so thank you 
  ### Brian and I hope you have a wonderful summer! 
)
