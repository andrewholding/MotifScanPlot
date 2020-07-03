#Motif Scan plots by A N Holding

plot_motifscan<-function(gene_chr,gene_start,gene_end,gene_padding) {
    gene_length =abs(gene_start - gene_end)
    gene_dl<-paste0(gene_chr,":",gene_start-gene_padding,",",gene_end+gene_padding)
    
    xml_data<-system(paste0('curl -s "http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=',gene_dl,'"'),intern = TRUE)
    
    
    
    require(XML)
    xml_obj<-xmlParse(xml_data)
    xml_list<-xmlToList(xml_obj)
    
    #Motif - https://bioconductor.org/packages/release/bioc/vignettes/universalmotif/inst/doc/SequenceSearches.pdf
    library(universalmotif)
    library(Biostrings)
    
    #HomerMotifs 
    #ERE download from http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults.html
    motif.ere<-read_homer("motif91.motif")
    motif.hif1a<-read_homer("motif152.motif")
    motif.hif1b<-read_homer("motif153.motif")
    motif.hif2a<-read_homer("motif154.motif")
    #motif<-"NAGGTCACNNTGACC"
    #motif.ere <- create_motif("NAGGTCACNNTGACC", nsites = 15)
    
    dna <- gsub("\n","",xml_list$SEQUENCE$DNA$text)
    dnaSet<-DNAStringSet(dna)
    
    #Threshold to 1 to get every site
    motif_scan <- scan_sequences(c(motif.ere,motif.hif1a,motif.hif1b,motif.hif2a), dnaSet, RC = TRUE,  threshold =1, threshold.type = "pvalue")
    motif_df <- data.frame(motif_scan)
    
    #gene_pos
    
    
    library(ggpubr)
    p<-ggplot(data=motif_df, aes(x=start, y=2^score,group=motif,color=motif)) 
    p<-p+  geom_line()
    
    pos_start=gene_padding
    pos_finish=gene_padding+gene_length
    
    p + geom_vline(aes(xintercept = pos_start)) +
        geom_vline(aes(xintercept = pos_finish)) +
        ggtitle(gene_dl)
}



#TFF1 chr21:42,362,282-42,366,535
gene_start = 42362282
gene_end = 42366535
gene_chr = "chr21"
gene_padding=10000
plot_motifscan(gene_chr,gene_start,gene_end,gene_padding)

#SCN5A chr3:38,548,066-38,649,628
gene_start = 38548066
gene_end = 38649628
gene_chr = "chr3"
gene_padding=50000
plot_motifscan(gene_chr,gene_start,gene_end,gene_padding)

#SCN5A_prom chr3:38,649,526-38,649,672
gene_start = 38649526
gene_end = 38649672
gene_chr = "chr3"
gene_padding=5000
plot_motifscan(gene_chr,gene_start,gene_end,gene_padding)