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
    motif.are<-read_homer("motif6.motif")
    motif.ere<-read_homer("motif91.motif")
    motif.gre<-read_homer("motif145.motif")
    
    #motif.hif1a<-read_homer("motif152.motif")
    #motif.hif1b<-read_homer("motif153.motif")
    #motif.hif2a<-read_homer("motif154.motif")
    #motif<-"NAGGTCACNNTGACC"
    #motif.ere <- create_motif("NAGGTCACNNTGACC", nsites = 15)
    
    dna <- gsub("\n","",xml_list$SEQUENCE$DNA$text)
    dnaSet<-DNAStringSet(dna)
    
    #Threshold to 1 to get every site
    motif_scan <- scan_sequences(c(motif.are,motif.gre,motif.ere), dnaSet, RC = TRUE,  threshold =1, threshold.type = "pvalue")
    motif_df <- data.frame(motif_scan)
    
    #gene_pos
    
    p<-ggplot(data=motif_df, aes(x=start, y=2^score,group=motif,color=motif)) 
    p<-p+  geom_line()
    
    pos_start=gene_padding
    pos_finish=gene_padding+gene_length
    
    return(  p + geom_vline(aes(xintercept = pos_start)) +
                 geom_vline(aes(xintercept = pos_finish)) +
                 ggtitle(paste0(title, " (" ,strand, ") ", gene_dl)))
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

plot_MotifScanByName<-function(geneSymbol, padding) {
    
    geneLocation<-getGeneLocation(geneSymbol)   
    
    p<-plot_MotifScan(geneLocation$TXCHROM[1],min(geneLocation$TXSTART),max(geneLocation$TXEND),padding,geneSymbol,geneLocation$TXSTRAND[1])
    return(p)
}


plot_MotifScanByName("TFF1",10000)
plot_MotifScanByName("TMPRSS2",10000)
plot_MotifScanByName("ACE2",10000)

geneLocation<-getGeneLocation("TMPRSS2")   
p<-plot_MotifScanByName("TMPRSS2",25000)
p+geom_vline(xintercept=geneLocation$TXSTART-min(geneLocation$TXSTART)+25000,color="black")+geom_vline(xintercept=geneLocation$TXEND-min(geneLocation$TXSTART)+25000,color="red")

