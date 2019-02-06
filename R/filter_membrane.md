---
title: "Filter Membrane Proteins"
output: 
  html_document: 
    keep_md: yes
---





```r
require(Seurat)
require(dplyr)
require(org.Hs.eg.db)
```

# Relies on having done the standard seurat workflow before `FindAllMarkers`


```r
data("pbmc_small")
print(head(pbmc_small@var.genes))
```

```
## [1] "LTB"      "EAF2"     "CD19"     "KIAA0125" "CYB561A3" "IGLL5"
```

```r
columns(org.Hs.eg.db)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [25] "UNIGENE"      "UNIPROT"
```

```r
# "GO:0009986" means cell surface 
# "GO:0005886" means plasma membrane
# "GO:0044459" means plasma membrane part

surface_annotated <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = pbmc_small@var.genes,
    columns = "GO",
    keytype = "ALIAS") %>% 
    filter(GO %in% c("GO:0009986", "GO:0005886", "GO:0044459"))
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
head(surface_annotated)
```

```
##     ALIAS         GO EVIDENCE ONTOLOGY
## 1     LTB GO:0005886      TAS       CC
## 2    CD19 GO:0005886      TAS       CC
## 3 PIK3IP1 GO:0005886      IEA       CC
## 4    CCR7 GO:0005886      IDA       CC
## 5    CCR7 GO:0005886      TAS       CC
## 6    CCR7 GO:0009986      IDA       CC
```

```r
pbmc.markers <- FindAllMarkers(
    object = pbmc_small, genes.use = unique(surface_annotated$ALIAS), 
    only.pos = TRUE,
    min.pct = 0.25, 
    thresh.use = 1, 
    print.bar = FALSE)


pbmc.markers
```

```
##                p_val avg_logFC pct.1 pct.2    p_val_adj cluster    gene
## LTB     1.537519e-11  2.503489 1.000 0.333 3.536293e-09       0     LTB
## PIK3IP1 1.465349e-05  3.556491 0.379 0.020 3.370302e-03       0 PIK3IP1
## PTGDR   8.352422e-05  3.770585 0.417 0.044 1.921057e-02       2   PTGDR
## S1PR4   1.954469e-03  2.652548 0.583 0.176 4.495278e-01       2   S1PR4
## CA2     1.524201e-13  4.797574 0.800 0.014 3.505663e-11       3     CA2
```

