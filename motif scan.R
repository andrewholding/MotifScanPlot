#Motif Scan plots by A N Holding
library(ggpubr)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(org.Hs.eg.db)


plot_MotifScan<-function(gene_chr,gene_start,gene_end,gene_padding,title="",strand="") {
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
    
    p<-ggplot(data=motif_df, aes(x=start, y=2^score,group=motif,color=motif)) 
    p<-p+  geom_line()
    
    pos_start=gene_padding
    pos_finish=gene_padding+gene_length
    
    p + geom_vline(aes(xintercept = pos_start)) +
        geom_vline(aes(xintercept = pos_finish)) +
        ggtitle(paste0(title, " (" ,strand, ") ", gene_dl))
}


sym2eg<-function(ids){
    list_symbol2eg <- as.character(org.Hs.egALIAS2EG[mappedkeys(org.Hs.egALIAS2EG)])
    ids <- as.character(ids)
    outlist <- list_symbol2eg[ids]
    names(outlist) <- ids
    outlist[is.na(outlist)] <- paste("unknown.", ids[is.na(outlist)], sep = "")
    outlist <- gsub("unknown.unknown.", "", outlist)
    return(outlist)
}

getGeneLocation<-function(geneSymbol) {
    cls <- columns(TxDb.Hsapiens.UCSC.hg38.knownGene)
    kts <- keytypes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    locations<-select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys=sym2eg(geneSymbol), columns=cls[c(16,20,17,21)], keytype=kts[5])
    return(locations)
}

#plot_MotifScan("chr21",42362282,42366535,10000,"TFF1")
#plot_MotifScan("chr3",38548066,38649628,50000,"SCN5A")
#plot_MotifScan("chr3",38649526,38649672,5000,"SCN5A Promoter")


plot_MotifScanByName<-function(geneSymbol, padding) {
    
    geneLocation<-getGeneLocation(geneSymbol)   
    
    plot_MotifScan(geneLocation$TXCHROM[1],geneLocation$TXSTART[1],geneLocation$TXEND[1],padding,geneSymbol,geneLocation$TXSTRAND[1])
}


pdf("output.pdf")
plot_MotifScanByName("TFF1",10000)
plot_MotifScanByName("SCN5A",50000)
plot_MotifScanByName("ATP1A1",50000)

plot_MotifScanByName("ATP1A2",50000)
plot_MotifScanByName("ATP1A3",50000)
plot_MotifScanByName("ATP1A4",50000)

plot_MotifScanByName("ATP1B1",50000)
plot_MotifScanByName("ATP1B2",50000)

plot_MotifScanByName("SLC5A1",50000)
plot_MotifScanByName("SLC5A2",50000)

plot_MotifScanByName("SLC38A1",50000)
plot_MotifScanByName("SLC38A3",50000)
plot_MotifScanByName("SLC38A5",50000)
plot_MotifScanByName("SLC38A7",50000)

plot_MotifScanByName("SLC7A5",50000)
plot_MotifScanByName("SLC7A8",50000)

plot_MotifScanByName("SLC1A2",50000)

plot_MotifScanByName("SLC16A10",50000)

plot_MotifScanByName("SLC6A19",50000)
plot_MotifScanByName("SLC6A14",50000)

plot_MotifScanByName("SLC7A6",50000)

plot_MotifScanByName("SLC8A1",50000)

plot_MotifScanByName("SLC9A1",50000)

plot_MotifScanByName("SLC4A7",50000)
plot_MotifScanByName("SLC12A2",50000)

plot_MotifScanByName("SCNN1A",50000)
plot_MotifScanByName("SCNN1B",50000)
plot_MotifScanByName("SCNN1G",50000)

plot_MotifScanByName("ASIC1",50000)
plot_MotifScanByName("ASIC2",50000)
plot_MotifScanByName("ASIC3",50000)

plot_MotifScanByName("SCN8A",50000)
plot_MotifScanByName("SCN10A",50000)
plot_MotifScanByName("SCN9A",50000)
plot_MotifScanByName("SCN2A",50000)
plot_MotifScanByName("SCN4A",50000)
plot_MotifScanByName("SCN5A",50000)

dev.off()