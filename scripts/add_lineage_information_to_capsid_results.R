library(readr)
library(dplyr)

blast <- read_tsv(snakemake@input[['diamond']], 
                  col_names = c("qseqid", "qstart", "qend", "qlen", "qstrand",
                                "sseqid", "sstart", "send", "slen", "pident",
                                "evalue", "cigar", "qseq_translated", "full_qseq",
                                "full_qseq_mate"))
lca <- read_tsv(snakemake@input[['lca']])

blast_lca <- left_join(blast, lca, by = c("sseqid"="rep"))

write_tsv(blast_lca, snakemake@output[['lca']])
