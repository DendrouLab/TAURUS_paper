stopifnot(
  require(AUCell),
  require(GSEABase),
  require(escape),
  require(GSVA),
  require(SingleCellExperiment),
  require(ggplot2),
  require(optparse),
  require(viridis),
  require(hrbrthemes),
  require(tidyverse),
  require(scales),
  require(ggpubr),
  require(ggalt),
  require(pheatmap),
  require(dplyr),
  require(randomcoloR),
  require(gridExtra),
  require(GGally)
)

rundir <- "/Users/tom/Desktop/sig_revised"

option_list <- list(
  make_option(
    c("--outdir"),
    default="output/",
    help="where should the output files be saved"
  ),
  make_option(
    c("--pseudobulk_matrix"),
    default=file.path(rundir,
                      "pseudocounts.tsv.gz"),
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
  ),
  make_option(
      c("--figwidth"),
      default=10,
      help="figure width in inches"),
  make_option(
      c("--figheight"),
      default=8,
      help="figure width in inches")

)

opt <- parse_args(OptionParser(option_list=option_list))

###read in expression matrix, metadata, and scores
exprMatrix <- read.table(opt$pseudobulk_matrix, header = TRUE, row.names = 1)
exprMatrix <- as.matrix(exprMatrix)
dim(exprMatrix)

metadata <- read.table("metadata.txt", header = TRUE, row.names = 1, sep = '\t')

a <- read.csv("MM.csv", header = FALSE)
a <- a$V1

b <- read.csv("M4.csv", header = FALSE)
b <- b$V1

c <- read.csv("M5.csv", header = FALSE)
c <- c$V1

d <- read.csv("M6.csv", header = FALSE)
d <- d$V1

e <- read.csv("Smilie.csv", header = FALSE)
e <- e$V1

###convert into GeneSets
MM <- GeneSetCollection(GeneSet(a, setName = "MM"))
M4Sets <- GeneSetCollection(GeneSet(b, setName="M4"))
M5Sets <- GeneSetCollection(GeneSet(c, setName="M5"))
M6Sets <- GeneSetCollection(GeneSet(d, setName="M6"))
Smilie <- GeneSetCollection(GeneSet(e, setName="SmilieScore"))

###generate scores
MMScore <- enrichIt(obj = exprMatrix, gene.sets = MM, groups = 1000, cores = 1)
SmilieScore <- enrichIt(obj = exprMatrix, gene.sets = Smilie, groups = 1000, cores = 1)
M4Scores <- enrichIt(obj = exprMatrix, gene.sets = M4Sets, groups = 1000, cores = 1)
M5Scores <- enrichIt(obj = exprMatrix, gene.sets = M5Sets, groups = 1000, cores = 1)
M6Scores <- enrichIt(obj = exprMatrix, gene.sets = M6Sets, groups = 1000, cores = 1)

Scores <- cbind(MMScore,SmilieScore,M4Scores,M5Scores,M6Scores)

###write out scores, so that they can be merged with adata.obs
write.csv(Scores, "Scores.csv")

###merge scores into metadata
metadata <- cbind(metadata, Scores)