## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Autor: Juan Francisco Alba Valle & Claudia Muñoz Mesa
## Fecha: Octubre 2019

peak.col<- "callpeakcol_peaks.narrowPeak"
peak.treat<- "callpeakmTCP4_peaks.narrowPeak"

library("clusterProfiler")

library("org.At.tair.db")

library("topGO")

library("pathview")

library("ChIPseeker")

library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


## Leer fichero de picos
col.peaks <- readPeakFile(peakfile = peak.col, header=FALSE)
col.peaks

treat.peaks <- readPeakFile(peakfile = peak.treat, header=FALSE)
treat.peaks

## Definir la región que se considera promotor entorno al TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)
promoter

## Anotación de los picos
col.peak.annot<- annotatePeak(peak = col.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)
col.peak.annot

treat.peak.annot<- annotatePeak(peak = treat.peaks, 
                              tssRegion=c(-1000, 1000),
                              TxDb=txdb)
treat.peak.annot

## Crear un plot con el porcentaje de la naturaleza de los picos dentro del DNA
## % de exones, promotores etc.

plotAnnoPie(col.peak.annot)
plotDistToTSS(col.peak.annot,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

plotAnnoPie(treat.peak.annot)
plotDistToTSS(treat.peak.annot,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Convertir la anotación a data frame
control.annotation <- as.data.frame(col.peak.annot)
head(control.annotation)

target.genes.control <- unique(control.annotation$geneId[control.annotation$annotation == "Promoter"])
length(target.genes.control)

write(x = target.genes.control ,file = "target_genes_control.txt")


treat.annotation <- as.data.frame(treat.peak.annot)
head(treat.annotation)

target.genes.treat <- unique(treat.annotation$geneId[treat.annotation$annotation == "Promoter"])
length(target.genes.treat)


write(x = target.genes.treat ,file = "target_genes_treat.txt")


#### Cáculo de los genes diferenciales entre el control y el tratamiento
target.genes<-setdiff(target.genes.treat,target.genes.control)

length(target.genes)

write(x = target.genes ,file = "target.genes.txt")


# Para realizar enriquecimientos funcionales de rutas kegg usamos clusterProfiler.
# Con HOMER encontraremos los Transcription Factors Binding Sites.
# Los paquetes han sido instalados arriba


# Necesitaremos el universo del cual realizaremos los enriquecimientos

ath.genes<-as.data.frame(genes(txdb))
ath.genes<- ath.genes$gene_id

# ENRIQUECIMIENTO GO

## Sirve para saber el ketipe
idType(OrgDb = org.At.tair.db)

## la variable ont, indica el tipo de análisis de enriquecimiento que se realiza,
## "BP" biological process
## "MF" molecular function
## "CC" celular compartiment

e.bp<- enrichGO(gene = target.genes,
             OrgDb = org.At.tair.db,
             keyType = "TAIR",
             ont = "BP",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             universe = ath.genes)

head(e.bp)

e.bp.table <- as.data.frame(e.bp)

head(e.bp.table)

barplot(e.bp, showCategory=20)


jpeg(filename = "plot_e_bp.jpeg", res = )
plot.e.bp<- plotGOgraph(x = e.bp,
                         firstSigNodes = 10,
                         useInfo = "all",
                         sigForAll = TRUE,
                         useFullNames = TRUE)
dev.off()

e.mf<- enrichGO(gene = target.genes,
                OrgDb = org.At.tair.db,
                keyType = "TAIR",
                ont = "MF",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                universe = ath.genes)

head(e.mf)

e.mf.table <- as.data.frame(e.mf)

head(e.mf.table)

barplot(e.mf, showCategory=20)
 

jpeg(filename = "plot_e_mf.jpeg")
plot.e.mf<- plotGOgraph(x = e.mf,
                         firstSigNodes = 10,
                         useInfo = "all",
                         sigForAll = TRUE,
                         useFullNames = TRUE)
dev.off()

e.cc<- enrichGO(gene = target.genes,
                OrgDb = org.At.tair.db,
                keyType = "TAIR",
                ont = "CC",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                universe = ath.genes)

head(e.cc)

e.cc.table <- as.data.frame(e.cc)

head(e.cc.table)

barplot(e.cc, showCategory=20)
   

jpeg(filename = "plot_e_cc.jpeg")
plot.e.cc<- plotGOgraph(x = e.cc,
                         firstSigNodes = 10,
                         useInfo = "all",
                         sigForAll = TRUE,
                         useFullNames = TRUE)
dev.off()


# ENRIQUECIMIENTO KEGG

e.kegg<- enrichKEGG(gene = target.genes, 
                    organism = "ath", 
                    universe = ath.genes)

e.kegg.table <- as.data.frame(e.kegg)


head(e.kegg.table)

Kegg.ID<-e.kegg.table$ID

# Preparación de los genes

# Path.id de cada función (en la función anterior)

for(i in 1:length(Kegg.ID))
{
  kegg.pathview<- pathview(gene.data = ath.genes, 
                           pathway.id = Kegg.ID[i] ,
                           species = "ath",
                           gene.idtype = "TAIR")
  
}