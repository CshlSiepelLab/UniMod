#### functions for processing data ####
# process bigwig file and convert it into single base resolution
process_bw <- function(bw, strand) {
  strand(bw) <- strand
  bw$score <- abs(bw$score)
  bw <- bw[bw$score > 0]
  bw <- keepStandardChromosomes(bw, pruning.mode = "coarse")
  bw <- BRGenomics::makeGRangesBRG(bw)
  return(bw)
}

# get the sum of read counts for regions within a gene
summarise_bw <-
  function(bw, grng, col_name) {
    rc <- bw %>%
      plyranges::find_overlaps_directed(grng) %>%
      plyranges::group_by(gene_id) %>%
      plyranges::summarise(score = sum(score * width))
    colnames(rc) <- c("gene_id", col_name)
    return(rc)
  }
