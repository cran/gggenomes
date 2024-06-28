## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = TRUE, fig.align = "center")

## ----message=FALSE, fig.height=1.2, fig.width=7-------------------------------
library(gggenomes)

# a minimal seq track
s0 <- tibble::tibble(
  seq_id = c("a", "b"),
  length = c(600, 550)
)

# a minimal gene track
g0 <- tibble::tibble(
  seq_id = c("a", "a", "b"),
  start = c(50, 350, 80),
  end = c(250, 500, 450)
)

# a simple link track
l0 <- tibble::tibble(
  seq_id = c("a", "a"),
  start = c(50, 400),
  end = c(250, 480),
  seq_id2 = c("b", "b"),
  start2 = c(80, 350),
  end2 = c(300, 430)
)

p <- gggenomes(genes=g0, seqs=s0, links=l0)
p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +   # label each sequence 
  geom_gene() +        # draw genes as arrow
  geom_link()          # draw some connections between syntenic regions


## -----------------------------------------------------------------------------
# Let's use some of the bundled example data here
data(package="gggenomes")

p <- gggenomes(
  genes=emale_genes,  # a gene track, added as first feat track
  seqs=emale_seqs,    # a seq track
  feats=list(emale_tirs, emale_ngaros),  # multiple feat tracks
  links=emale_ava     # a link track
)

# inspect the tracks of the plot
p %>% track_info

# plot all tracks
p +
#  geom_link() +  # the first link track
  geom_gene() +  # the first feat track filtered for geneish feats: CDS, mRNA, ..
  geom_feat() +  # the first feat track not named "genes", here emale_tirs
  # use an extra feat track by name
  geom_feat(data=feats(emale_ngaros), color="plum3")

## ----fig.height=1.5, fig.width=7----------------------------------------------
# inspect seqs track with layout vars - note y,x,xend
p %>% pull_seqs
# inspect genes track with layout vars - note y,x,xend, but also other
# columns such as strand, feat_id or type, that are added automatically
p %>% pull_genes

## ----fig.height=1.5, fig.width=7----------------------------------------------
# some genes
g0 <- tibble::tibble(
  seq_id = c("a", "a", "b"),
  start = c(50, 350, 80),
  end = c(250, 500, 450)
)

p <- gggenomes(g0) 
p +
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +   # label each sequence 
  geom_gene()          # draw genes as arrow

## ----fig.height=1.5, fig.width=7, message=FALSE-------------------------------
# note: ex() is just a helper to get stable paths to gggenomes example data
s0 <- read_seqs(ex("emales/emales.fna"))
g0 <- read_feats(ex("emales/emales.gff"))

gggenomes(g0, s0) +
  geom_seq() + geom_gene()

# for lazy people
gggenomes(ex("emales/emales.gff")) + geom_gene() 

# and really fancy: multiple remote files, all at once
gbk_phages <- c(
  PSSP7 = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/858/745/GCF_000858745.1_ViralProj15134/GCF_000858745.1_ViralProj15134_genomic.gff.gz",
  PSSP3 = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/904/555/GCF_000904555.1_ViralProj195517/GCF_000904555.1_ViralProj195517_genomic.gff.gz")
gggenomes(gbk_phages) + geom_gene() +
  geom_seq_label()

## ----fig.height=1.5, fig.width=7----------------------------------------------
# seq track: one entry per sequence
s0 <- tibble::tibble(
  bin_id = c("A", "A", "B"),
  seq_id = c("a1","a2","b1"),
  length = c(2e5, 3e5, 5e5)
)

p <-  gggenomes(seqs=s0) 
p +
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label()    # label each sequence 
  #geom_bin_label()  # label each bin

## ----fig.height=1.5, fig.width=7----------------------------------------------
# zoom in on a longer sequence - note the scale on the x-axis
s0 <- tibble::tibble(
  seq_id = "a1",
  length = 10000,
  start = 1000,
  end = 3000
)
gggenomes(seqs=s0) + geom_seq() + geom_seq_label()


## ----fig.height=1.5, fig.width=7----------------------------------------------
# some genes
g0 <- tibble::tibble(
  seq_id = c("a"),
  start = c(1, 800),
  end = c(500, 1200),
  # NOTE: introns need to be a list-column!
  introns = list(c(50,200), c(50,200,250,300))
)

gggenomes(g0) +
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +   # label each sequence 
  geom_gene()          # draw genes as arrow

## -----------------------------------------------------------------------------
# some links
l0 <- tibble::tibble(
  seq_id = c("a", "a", "a"),
  start = c(200, 801, 1600),
  end = c(550, 1300, 1800),
  seq_id2 = c("b", "b", "b"),
  start2 = c(1100, 1, 1800),
  end2 = c(1450, 500, 1600)
)

# corresponding sequences
s1 <- tibble::tibble(
  seq_id = c("a", "b"), 
  length = c(2000, 2000),
  start = c(1, 1),
  end = c(2000, 2000)
)


gggenomes(seqs=s1, links=l0)  +
  geom_seq() +                  # draws contigs/chromosome lines  
  geom_seq_label()             # labels each sequence
 # geom_link(offset = 0.05)      # draws links between contigs

