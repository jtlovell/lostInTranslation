## Title:  Parse GENESPACE output
## Author: JTLovell
## Date:   30-May 2024

################################################################################
# 1. Combine orthogroup information stored in the GENESPACE output
library(data.table)
library(ape)
library(data.table)
cbArab <- fread("orthologs/analysis/soybean_v_arabidopsis/results/combBed.txt")
cbArab[,contrast := "soybean_v_arabidopsis"]

cbRice <- fread("orthologs/analysis/maize_v_rice/results/combBed.txt")
cbRice[,contrast := "maize_v_rice"]

cbMous <- fread("orthologs/analysis/human_v_mouse/results/combBed.txt")
cbMous[,contrast := "human_v_mouse"]

cb <- rbind(cbArab, cbRice, cbMous)

################################################################################
# 2. count the number of genes and genomes in orthogroups
# -- choose genes of interest for gene trees
cbMous[,`:=`(
  ngenesSynOG = .N, ngenomesSynOG = uniqueN(genome),
  ngenesMod = sum(genome == "human"), 
  ngenesTst = sum(genome == "mouse")), 
  by = c("contrast", "og")]
cbMous[,`:=`(
  ngenesHOG = .N, ngenomesHOG = uniqueN(genome)), 
  by = c("contrast", "globHOG")]
cbMous[, isRep := ngenesMod == ngenesTst]

cbRice[,`:=`(
  ngenesSynOG = .N, ngenomesSynOG = uniqueN(genome),
  ngenesMod = sum(genome == "rice"), 
  ngenesTst = sum(genome == "maize")), 
  by = c("contrast", "globHOG")]
cbRice[,isOtl := ngenesTst >3  & ngenesMod >3 & genome == "rice" & ngenesTst < 7 & ngenesMod < 7]
cbRice[,`:=`(
  ngenesSynOG = .N, ngenomesSynOG = uniqueN(genome),
  ngenesMod = sum(genome == "rice"), 
  ngenesTst = sum(genome == "maize")), 
  by = c("contrast", "og")]
cbRice[,`:=`(
  ngenesHOG = .N, ngenomesHOG = uniqueN(genome)), 
  by = c("contrast", "globHOG")]
cbRice[, isRep := ngenesMod == 1 & ngenesTst == 2]
choose <- "1000" # pyrroline-5-carboxylate reductase
cbRice[,isOtl := ngenesTst >= 5 & ngenesMod >= 5]
choose <- "1746" # pyrroline-5-carboxylate reductase

cbArab[,`:=`(
  ngenesSynOG = .N, ngenomesSynOG = uniqueN(genome),
  ngenesMod = sum(genome == "arabidopsis"), 
  ngenesTst = sum(genome == "soybean")), 
  by = c("contrast", "og")]
cbArab[,`:=`(
  ngenesHOG = .N, ngenomesHOG = uniqueN(genome)), 
  by = c("contrast", "globHOG")]
cbArab[, isRep := ngenesMod == 4 & ngenesTst == 4]
choose <- "OG0000960" # pyrroline-5-carboxylate reductase
cbArab[,isOtl := ngenesTst >= 5 & ngenesMod >= 5]
(sum(cbArab$ngenesTst >= 4 & cbArab$ngenesMod >= 4) / sum(cbArab$ngenesSynOG >= 0))

################################################################################
# 3. Make gene trees for representative orthogroups 
pepArab <- readAAStringSet("orthologs/analysis/soybean_v_arabidopsis/peptide/arabidopsis.fa")
pepSoyb <- readAAStringSet("orthologs/analysis/soybean_v_arabidopsis/peptide/soybean.fa")
cbrep <- subset(cbArab, og == "3263")
gsMod <- cbrep$id[cbrep$genome == "arabidopsis"]
gsTst <- cbrep$id[cbrep$genome == "soybean"]
peprep <- c(pepArab[gsMod], pepSoyb[gsTst])

writeXStringSet(peprep, file = "orthologs/analysis/asrep.unal.fa")
system("mafft orthologs/analysis/asrep.unal.fa > orthologs/analysis/asrep.al.fa")
system("fasttree orthologs/analysis/asrep.al.fa > orthologs/analysis/asrep.al.tre")

cbotl <- subset(cbArab, og == "4393")
gsMod <- cbotl$id[cbotl$genome == "arabidopsis"]
gsTst <- cbotl$id[cbotl$genome == "soybean"]
peprep <- c(pepArab[gsMod], pepSoyb[gsTst])

writeXStringSet(peprep, file = "orthologs/analysis/asotl.unal.fa")
system("mafft orthologs/analysis/asotl.unal.fa > orthologs/analysis/asotl.al.fa")
system("fasttree orthologs/analysis/asotl.al.fa > orthologs/analysis/asotl.al.tre")

pdf("orthologs/atTrees.pdf", height = 3, width = 6)
par(mfrow = c(1, 2))
tre <- read.tree("orthologs/analysis/asrep.al.tre")
plot(ladderize(tre), cex = .5, type = "u")
tre2 <- read.tree("orthologs/analysis/asotl.al.tre")
plot(ladderize(tre2), cex = .5, type = "u")
dev.off()

pepRice <- readAAStringSet("orthologs/analysis/maize_v_rice//peptide/rice.fa")
pepMaiz <- readAAStringSet("orthologs/analysis/maize_v_rice/peptide/maize.fa")
cbrep <- subset(cbRice, og == "1000")
gsMod <- cbrep$id[cbrep$genome == "rice"]
gsTst <- cbrep$id[cbrep$genome == "maize"]
peprep <- c(pepRice[gsMod], pepMaiz[gsTst])

writeXStringSet(peprep, file = "orthologs/analysis/rmrep.unal.fa")
system("mafft orthologs/analysis/rmrep.unal.fa > orthologs/analysis/rmrep.al.fa")
system("fasttree orthologs/analysis/rmrep.al.fa > orthologs/analysis/rmrep.al.tre")

cbotl <- subset(cbRice,  globHOG == "N0.HOG0000932 OG0000639")
gsMod <- cbotl$id[cbotl$genome == "rice"]
gsTst <- cbotl$id[cbotl$genome == "maize"]
peprep <- c(pepRice[gsMod], pepMaiz[gsTst])

writeXStringSet(peprep, file = "orthologs/analysis/asotl.unal.fa")
system("mafft orthologs/analysis/asotl.unal.fa > orthologs/analysis/asotl.al.fa")
system("fasttree orthologs/analysis/asotl.al.fa > orthologs/analysis/asotl.al.tre")

pdf("orthologs/zrTrees.pdf", height = 3, width = 6)
par(mfrow = c(1, 2))
tre <- read.tree("orthologs/analysis/rmrep.al.tre")
plot(ladderize(tre), cex = .5, type = "u")
tre2 <- read.tree("orthologs/analysis/asotl.al.tre")
plot(ladderize(tre2), cex = .5, type = "u")
dev.off()

cbRice[,`:=`(ngenesSynOG = .N, ngenomesSynOG = uniqueN(genome)), by = c("contrast", "og")]
cbRice[,`:=`(ngenesHOG = .N, ngenomesHOG = uniqueN(genome)), by = c("contrast", "globHOG")]

cbArab[,`:=`(ngenesSynOG = .N, ngenomesSynOG = uniqueN(genome)), by = c("contrast", "og")]
cbArab[,`:=`(ngenesHOG = .N, ngenomesHOG = uniqueN(genome)), by = c("contrast", "globHOG")]

cb[,pavType := ifelse(ngenes == 2 & ngenomes == 2, "1to1",
                      ifelse(ngenomes == 1, "PAV", "CNV"))]
subset(cb, globHOG == "N0.HOG0000522 OG0000305" & contrast == "human_v_mouse")

cbHog <- subset(cb, !duplicated(paste(contrast, og)))
cbCnt <- cbHog[,list(n = .N), by = c("contrast", "pavType")]
cbq <- cbHog[,list(q90 = quantile(ngenes, .95)), by = "contrast"]

cbCnt[,pavType := factor(pavType, levels = c("1to1", "PAV", "CNV"))]

pdf("orthologs/cnvCnt.pdf", height = 2, width = 6)
cbCnt[,xlab := sprintf("%s %s",contrast, sum(n)), by = "contrast"]
ggplot(cbCnt, aes(fill = pavType, x = xlab, y = n, group = pavType))+
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(~contrast, scale = "free", nrow = 3)+
  scale_fill_manual(values = c("salmon", "khaki", "dodgerblue4"))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())+
  coord_flip()
dev.off()

cbHog[,nGeneCut := ifelse(ngenes > 10, 10, ngenes)]
ggplot(cbHog, 
       aes(fill = contrast, x = nGeneCut, group = contrast))+
  geom_histogram(binwidth = 1)+
  facet_grid(contrast ~ .)

with(subset(cbHog, contrast == "human_v_mouse"), sum(ngenes > 7)/sum(ngenes > 0))
# HLA-DRB is among the 98.4th percentile of largest gene families
pepHum <- readAAStringSet("orthologs/analysis/human_v_mouse/peptide/human.fa")
pepMou <- readAAStringSet("orthologs/analysis/human_v_mouse/peptide/mouse.fa")
cbHLA <- subset(cb, globHOG == "N0.HOG0000337 OG0000214")
gsHum <- cbHLA$id[cbHLA$genome == "human"]
gsMou <- cbHLA$id[cbHLA$genome == "mouse"]
pepHLA <- c(pepHum[gsHum], pepMou[gsMou])

writeXStringSet(pepHLA, file = "orthologs/analysis/hla.unal.fa")
system("mafft orthologs/analysis/hla.unal.fa > orthologs/analysis/hla.al.fa")
system("fasttree orthologs/analysis/hla.al.fa > orthologs/analysis/hla.al.tre")


pdf("orthologs/hlatre.pdf", height = 4, width = 4)
tr <- read.tree("orthologs/analysis/hla.al.tre")
plot(ladderize(tr), type = "u")
dev.off()
treu <- read.tree("orthologs/analysis/human_v_mouse/orthofinder/Results_May30/Resolved_Gene_Trees/OG0000666_tree.txt")
plot(ladderize(tr), type = "u")


################################################################################
# 4. make the timetree phylogeny
pdf("orthologs/timetre.pdf", height = 2, width = 2)
tr <- read.tree("/Users/lovell/Downloads/sp.nwk")
plot(tr, cex = .5)
dev.off()


cbi <- subset(cb, contrast == "soybean_v_arabidopsis")
cbi[,ngm := sum(genome == "soybean"), by = "globHOG"]
cbi[,nat := sum(genome == "arabidopsis"), by = "globHOG"]
cbg <- subset(cbi, ngm == nat & nat == 4)

quantile(cbi$ngenes[!duplicated(cbi$globHOG)], .984)
with(subset(cbHog, contrast == "human_v_mouse"), sum(ngenes > 7)/sum(ngenes > 0))
# HLA-DRB is among the 98.4th percentile of largest gene families
pepHum <- readAAStringSet("orthologs/analysis/human_v_mouse/peptide/human.fa")
pepMou <- readAAStringSet("orthologs/analysis/human_v_mouse/peptide/mouse.fa")
cbHLA <- subset(cb, globHOG == "N0.HOG0000337 OG0000214")
gsHum <- cbHLA$id[cbHLA$genome == "human"]
gsMou <- cbHLA$id[cbHLA$genome == "mouse"]
pepHLA <- c(pepHum[gsHum], pepMou[gsMou])


tre <- read.tree("orthologs/analysis/soybean_v_arabidopsis/orthofinder/Results_May30/Resolved_Gene_Trees/OG0000468_tree.txt")
pdf("orthologs/atgmtreOG0000468.pdf", height = 4, width = 4)
plot(ladderize(tre), type = "u")
dev.off()

################################################################################
# 5. Compare same/different annotation methods
cb1 <- fread("orthologs/analysis/withinRice/results/combBed.txt")
cb2 <- fread("orthologs/analysis/withinSoybean/results/combBed.txt")
cb1[,contrast := "diff"]
cb2[,contrast := "same"]
cbc <- rbind(cb1, cb2)
tp <- cbc[,list(type = ifelse(uniqueN(genome) == 2 & .N == 2, "1:1",
                              ifelse(uniqueN(genome) == 1, "PAV", "CNV"))), by = c("contrast","globHOG")]
tp <- tp[,list(n = .N), by = c("contrast", "type")]
tp[,type := factor(type, levels = rev(c("CNV", "PAV", "1:1")))]
pdf("orthologs/sameDiffBar.pdf", height = 3, width = 3)
ggplot(tp, aes(x = contrast, y = n, fill = type))+
  geom_bar(position = "stack", stat = "identity")+
  facet_wrap(~contrast, scale = "free", nrow = 2)+
  coord_flip()
dev.off()
