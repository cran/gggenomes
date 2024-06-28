## ----eval=FALSE---------------------------------------------------------------
#  library(gggenomes)
#  
#  # parse sequence length and some metadata from fasta file
#  emale_seqs <- read_fai("emales.fna") %>%
#    tidyr::extract(seq_desc, into = c("emale_type", "is_typespecies"), "=(\\S+) \\S+=(\\S+)",
#      remove=F, convert=T) %>%
#    dplyr::arrange(emale_type, length)
#  
#  # plot the genomes - first six only to keep it simple for this example
#  emale_seqs_6 <- emale_seqs[1:6,]
#  p1 <- gggenomes(emale_seqs_6) +
#   geom_seq() + geom_bin_label()
#  p1

## ----eval=FALSE---------------------------------------------------------------
#  emale_genes <- read_gff("emales.gff") %>%
#    dplyr::rename(feature_id=ID) %>%                       # we'll need this later
#    dplyr::mutate(gc_cont=as.numeric(gc_cont))             # per gene GC-content
#  
#  p2 <- gggenomes(emale_seqs_6, emale_genes) +
#    geom_seq() + geom_bin_label() +
#    geom_gene(aes(fill=gc_cont)) +
#    scale_fill_distiller(palette="Spectral")
#  p2

## ----eval=FALSE---------------------------------------------------------------
#  # prefilter hits by minimum length and maximum divergence
#  emale_tirs_paf <- read_paf("emales-tirs.paf") %>%
#    dplyr::filter(seq_id1 == seq_id2 & start1 < start2 & map_length > 99 & de < 0.1)
#  emale_tirs <- bind_rows(
#    dplyr::select(emale_tirs_paf, seq_id=seq_id1, start=start1, end=end1, de),
#    dplyr::select(emale_tirs_paf, seq_id=seq_id2, start=start2, end=end2, de))
#  
#  p3 <- gggenomes(emale_seqs_6, emale_genes, emale_tirs) +
#    geom_seq() + geom_bin_label() +
#    geom_feature(size=5) +
#    geom_gene(aes(fill=gc_cont)) +
#    scale_fill_distiller(palette="Spectral")
#  p3

## ----eval=FALSE---------------------------------------------------------------
#  emale_links <- read_paf("emales.paf")
#  
#  p4 <- gggenomes(emale_seqs_6, emale_genes, emale_tirs, emale_links) +
#    geom_seq() + geom_bin_label() +
#    geom_feature(size=5, data=use_features(features)) +
#    geom_gene(aes(fill=gc_cont)) +
#    geom_link() +
#    scale_fill_distiller(palette="Spectral")
#  
#  p4 <- p4 %>% flip_bins(3:5)
#  p4

## ----eval=FALSE---------------------------------------------------------------
#  emale_gc <- thacklr::read_bed("emales-gc.tsv") %>%
#    dplyr::rename(seq_id=contig_id)
#  
#  p5 <- p4 %>% add_features(emale_gc)
#  p5 <- p5 + geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),
#      group=seq_id, linetype="GC-content"), use_features(emale_gc),
#                         fill="blue", alpha=.5)
#  p5

## ----eval=FALSE---------------------------------------------------------------
#  emale_cogs <- read_tsv("emales-cogs.tsv", col_names = c("feature_id", "cluster_id", "cluster_n"))
#  emale_cogs %<>% dplyr::mutate(
#    cluster_label = paste0(cluster_id, " (", cluster_n, ")"),
#    cluster_label = fct_lump_min(cluster_label, 5, other_level = "rare"),
#    cluster_label = fct_lump_min(cluster_label, 15, other_level = "medium"),
#    cluster_label = fct_relevel(cluster_label, "rare", after=Inf))
#  emale_cogs
#  
#  
#  p6 <- gggenomes(emale_seqs_6, emale_genes, emale_tirs, emale_links) %>%
#    add_features(emale_gc) %>%
#    add_clusters(genes, emale_cogs) %>%
#    flip_bins(3:5) +
#    geom_seq() + geom_bin_label() +
#    geom_feature(size=5, data=use_features(features)) +
#    geom_gene(aes(fill=cluster_label)) +
#    geom_link() +
#    geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),
#      group=seq_id, linetype="GC-content"), use_features(emale_gc),
#                         fill="blue", alpha=.5) +
#    scale_fill_brewer("Conserved genes", palette="Set3")
#  
#  p6

## ----eval=FALSE---------------------------------------------------------------
#  emale_blast <- read_blast("emales_mavirus-blastp.tsv")
#  emale_blast %<>%
#    dplyr::filter(evalue < 1e-3) %>%
#    dplyr::select(feature_id=qaccver, start=qstart, end=qend, saccver) %>%
#    dplyr::left_join(read_tsv("mavirus.tsv", col_names = c("saccver", "blast_hit", "blast_desc")))
#  
#  # manual annotations by MFG
#  emale_transposons <- read_gff("emales-manual.gff", types = c("mobile_element"))
#  
#  
#  p7 <- gggenomes(emale_seqs_6, emale_genes, emale_tirs, emale_links) %>%
#    add_features(emale_gc) %>%
#    add_clusters(genes, emale_cogs) %>%
#    add_features(emale_transposons) %>%
#    add_subfeatures(genes, emale_blast, transform="aa2nuc") %>%
#    flip_bins(3:5) +
#    geom_feature(aes(color="integrated transposon"),
#      use_features(emale_transposons), size=7) +
#    geom_seq() + geom_bin_label() +
#    geom_link(offset = c(0.3, 0.2), color="white", alpha=.3) +
#    geom_feature(aes(color="terminal inverted repeat"), use_features(features),
#      size=4) +
#    geom_gene(aes(fill=cluster_label)) +
#    geom_feature(aes(color=blast_desc), use_features(emale_blast), size=2,
#      position="pile") +
#    geom_ribbon(aes(x=(x+xend)/2, ymax=y+.24, ymin=y+.38-(.4*score),
#      group=seq_id, linetype="GC-content"), use_features(emale_gc),
#                         fill="blue", alpha=.5) +
#    scale_fill_brewer("Conserved genes", palette="Set3") +
#    scale_color_viridis_d("Blast hits & Features", direction = -1) +
#    scale_linetype("Graphs") +
#    ggtitle(expression(paste("Endogenous mavirus-like elements of ",
#    italic("C. burkhardae"))))
#  
#  p7

