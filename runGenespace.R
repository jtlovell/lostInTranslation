## Title:  Run GENESPACE
## Author: JTLovell
## Date:   30-May 2024

library(GENESPACE)
path2mcscanx <- "/programs/MCScanX/"

################################################################################
# 1. Human v. mouse
basedir <- "orthologs/analysis"
genomeRepo <- "orthologs/genomeRepo/"
gids <- c("human_Hg38", "mouse_GRCm39")
gnames <- c("human", "mouse")
wd <- file.path(basedir, "human_v_mouse")
if(!dir.exists(wd))
  dir.create(wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids,
  genomeIDs = gnames, 
  faString = "peptide",
  presets = "ncbi",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx, 
  nCores = 4)
out <- run_genespace(gpar, overwrite = FALSE)

################################################################################
# 2. Rice v. Maize
basedir <- "orthologs/analysis"
genomeRepo <- "orthologs/genomeRepo/"
gids <- c("maize_B73", "rice_IRGSP")
gnames <- c("maize", "rice")
wd <- file.path(basedir, "maize_v_rice")
if(!dir.exists(wd))
  dir.create(wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids,
  genomeIDs = gnames, 
  faString = "peptide",
  presets = "ncbi",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  ploidy = c(2, 1),
  path2mcscanx = path2mcscanx, 
  nCores = 4)
out <- run_genespace(gpar, overwrite = FALSE)

################################################################################
# 3. arabidopsis v. soybean
basedir <- "orthologs/analysis"
genomeRepo <- "orthologs/genomeRepo/"
gids <- c("soybean_WM82", "arabidopsis_TAIR10")
gnames <- c("soybean", "arabidopsis")
wd <- file.path(basedir, "soybean_v_arabidopsis")
if(!dir.exists(wd))
  dir.create(wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids,
  genomeIDs = gnames,
  presets = "phytozome",
  faString = "peptide",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx, 
  ploidy = c(4, 4),
  nCores = 4)
out <- run_genespace(gpar, overwrite = FALSE)

################################################################################
# 4. Human Hg38 v. CHM13
basedir <- "orthologs/analysis"
genomeRepo <- "orthologs/genomeRepo/"
gids <- c("human_Hg38", "human_CHM13")
gnames <- c("hg38", "chm13")
wd <- file.path(basedir, "hg38_v_chm13")
if(!dir.exists(wd))
  dir.create(wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids,
  genomeIDs = gnames, 
  faString = "peptide",
  presets = "ncbi",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx, 
  nCores = 4)
out <- run_genespace(gpar, overwrite = FALSE)


################################################################################
# 5. Rice kitaake vs. niponbarre
basedir <- "orthologs/analysis"
genomeRepo <- "orthologs/genomeRepo"
gids <- c("rice_IRGSP", "rice_kitaake")
gnames <- c("niponbarre", "kitaake")
wd <- file.path(basedir, "withinRice")
if(!dir.exists(wd))
  dir.create(wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids[1],
  genomeIDs = gnames[1], 
  faString = "peptide",
  presets = "ncbi",
  genespaceWd = wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids[2],
  genomeIDs = gnames[2],
  presets = "phytozome",
  faString = "prot",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx, 
  nCores = 4)
out <- run_genespace(gpar, overwrite = FALSE)

################################################################################
# 5. maize B73 vs. MO17
basedir <- "orthologs/analysis"
genomeRepo <- "orthologs/genomeRepo"
gids <- c("soybean_WM82_phytozome", "soybean_Fiskeby_phytozome")
gnames <- c("Wm82", "Fiskeby")
wd <- file.path(basedir, "withinSoybean")
if(!dir.exists(wd))
  dir.create(wd)
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = gids,
  genomeIDs = gnames, 
  presets = "phytozome",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx, 
  nCores = 4)
out <- run_genespace(gpar, overwrite = FALSE)


