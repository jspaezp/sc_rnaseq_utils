---
title: "Tree Based Variable Selection"
output: 
  html_document: 
    keep_md: yes
---





```r
require(Seurat)
require(dplyr)
require(ranger)
require(partykit)

data("pbmc_small")
```




```r
#' Gets a data frame from a seurat object
#' 
#' Provided a seurat object, returs a data frame of the count values, being 
#' the columns each 'gene' and the rows each UMI/cell.
#'
#' It returns only the genes annotated as variable and the identity column.
#'
#' @param seurat A seurat object
#'
#' @return
#' @export
#'
#' @examples
#' > as.data.frame.Seurat(pbmc_small)[1:3,1:3]
#'                    LTB EAF2 CD19
#' ATGCCAGAACGACT 6.062788    0    0
#' CATGGCCTGTGCAT 6.714813    0    0
#' GAACCTGATGAACC 7.143118    0    0
#' @importFrom Seurat FetchData GetIdent
as.data.frame.Seurat <- function(seurat) {
    tmp <- Seurat::FetchData(seurat, vars.all = seurat@var.genes) 
    tmp <-  as.data.frame(tmp, stringsAsFactors = FALSE) 
       
    
    tmp$ident <- Seurat::GetIdent(seurat, uniq = FALSE, cells.use = rownames(tmp))
    return(tmp)
}


pbmc_df <- as.data.frame.Seurat(pbmc_small)

# Some minor renaming is required for most methods, since some genes have - in 
# their names, and we require R-valid names.

original_names <- colnames(pbmc_df)
newnames <- make.names(original_names)
colnames(pbmc_df) <- newnames
```


## Most Basic Form of Decision Tree based variable selection



```r
tree <- ctree(ident ~ ., pbmc_df)
tree
```

```
## 
## Model formula:
## ident ~ LTB + EAF2 + CD19 + KIAA0125 + CYB561A3 + IGLL5 + PIK3IP1 + 
##     KHDRBS1 + CCR7 + ACSM3 + SRSF7 + S1PR4 + LYAR + SATB1 + IL17RA + 
##     POP7 + ZNF330 + COPS6 + PPBP + PF4 + HIST1H2AC + TALDO1 + 
##     CA2 + ACRBP + TSC22D1 + VDAC3 + GNLY + PTGDR + ARHGDIA + 
##     PCMT1 + S100B
## 
## Fitted party:
## [1] root
## |   [2] PPBP <= 4.75309
## |   |   [3] GNLY <= 4.6764
## |   |   |   [4] LTB <= 3.79076: 1 (n = 28, err = 14.3%)
## |   |   |   [5] LTB > 3.79076: 0 (n = 29, err = 13.8%)
## |   |   [6] GNLY > 4.6764: 2 (n = 13, err = 15.4%)
## |   [7] PPBP > 4.75309: 3 (n = 10, err = 0.0%)
## 
## Number of inner nodes:    3
## Number of terminal nodes: 4
```

```r
varimp(tree)
```

```
##      PPBP      GNLY       LTB 
## 3.5174292 1.0646640 0.6986878
```

```r
plot(tree)
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

which could be plotted later by selecting the genes that worked best.
Also note how only imprtances are given to the genes that were used in the end.


```r
ggplot2::qplot(x = pbmc_df$PPBP, y = pbmc_df$GNLY, colour = pbmc_df$ident)
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggplot2::qplot(x = pbmc_df$PPBP, fill = pbmc_df$ident, geom = "density")
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```r
ggplot2::qplot(x = pbmc_df$GNLY, fill = pbmc_df$ident, geom = "density")
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-4-3.png)<!-- -->


## Most Basic Form of Random forest based variable selection


```r
forest <- ranger::ranger(
  ident ~ ., data =  pbmc_df,
  importance = "impurity_corrected")

forest
```

```
## Ranger result
## 
## Call:
##  ranger::ranger(ident ~ ., data = pbmc_df, importance = "impurity_corrected") 
## 
## Type:                             Classification 
## Number of trees:                  500 
## Sample size:                      80 
## Number of independent variables:  31 
## Mtry:                             5 
## Target node size:                 1 
## Variable importance mode:         impurity_corrected 
## Splitrule:                        gini 
## OOB prediction error:             18.75 %
```

```r
df_importances <- as.data.frame(
  ranger::importance_pvalues(forest)[
      ranger::importance_pvalues(forest)[,'pvalue']<0.05,])
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.

## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```r
head(df_importances[order(df_importances[, "importance"]),] )
```

```
##         importance pvalue
## PIK3IP1  0.3457385      0
## S100B    0.3569180      0
## ACRBP    0.3860325      0
## PTGDR    0.4733949      0
## S1PR4    0.5242100      0
## TSC22D1  0.5660943      0
```

```r
dim(df_importances)
```

```
## [1] 17  2
```


## Adding_weights for under-represented clusters



```r
# This add weihts to the clusters in proportion to the inverse of their abundance
# TODO find a better way to calculate the weights ...
weights_vars <- ceiling(as.numeric(100/table(pbmc_df$ident)))


forest <- ranger::ranger(
  ident ~ ., data =  pbmc_df,
  importance = "impurity_corrected", 
  class.weights = weights_vars)

forest
```

```
## Ranger result
## 
## Call:
##  ranger::ranger(ident ~ ., data = pbmc_df, importance = "impurity_corrected",      class.weights = weights_vars) 
## 
## Type:                             Classification 
## Number of trees:                  500 
## Sample size:                      80 
## Number of independent variables:  31 
## Mtry:                             5 
## Target node size:                 1 
## Variable importance mode:         impurity_corrected 
## Splitrule:                        gini 
## OOB prediction error:             12.50 %
```

```r
df_importances <- as.data.frame(
  ranger::importance_pvalues(forest)[
      ranger::importance_pvalues(forest)[,'pvalue']<0.05,])
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.

## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```r
head(df_importances[order(df_importances[, "importance"]),] )
```

```
##         importance pvalue
## S100B     9.423453      0
## PIK3IP1  10.398379      0
## KHDRBS1  11.380939      0
## VDAC3    12.042864      0
## TSC22D1  13.374409      0
## S1PR4    13.656005      0
```

```r
dim(df_importances)
```

```
## [1] 19  2
```

Note that due to the nature of the random forest, it cannot be plotted ...
but importances can


```r
qplot(x = df_importances$importance, y = rownames(df_importances))
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

.... since factors are weird in R we need to reorder the levels for the plot to
make sense ...


```r
genes <- forcats::fct_reorder(rownames(df_importances), df_importances$importance)
qplot(x = df_importances$importance, y = genes)
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-8-1.png)<!-- -->



## Packing all in a function


```r
find_cluster_classification_forest <- function(cluster_id,
                                               class_df,
                                               markers, 
                                               report_n = 4) {
    print(paste0("Cluster: ", cluster_id))
    foo2 <- markers %>% 
        unique() %>%
      {class_df[,colnames(class_df) %in% c(., 'ident')]}
    
    in_cluster <- foo2$ident == cluster_id
    out_cluster <- !in_cluster
    
    num_outcluster <- sum(out_cluster)
    num_incluster <- sum(in_cluster)

    foo2$ident <- foo2$ident == cluster_id
    
    # This section alleviates slightly the problem when the cluster is extremely
    # small compared to the non-cluster points.
    # TODO  add a warning whenever this happens...
    if (100 < (num_outcluster/num_incluster)) {
      foo2 <- foo2[c(which(foo2$ident), 
                     sample(which(!foo2$ident), num_incluster*100)), ]
    }
    
    foo2$ident <- factor(make.names(foo2$ident))
        
    forest <- ranger::ranger(
        ident ~ ., 
        data = foo2,
        importance = "impurity_corrected",
        num.trees = 2000)

    df_importances <- as.data.frame(
      ranger::importance_pvalues(forest))
    
    df_importances$marker <- rownames(df_importances)
    
    markers <- dplyr::top_n(df_importances, report_n, importance)$marker  

    return(
      list(markers = markers, 
           tree = forest, 
           cluster = cluster_id, 
           df_importances = df_importances))
}


tmp <- find_cluster_classification_forest(
    cluster_id = "1",
    class_df = pbmc_df,
    markers = pbmc_small@var.genes,
    report_n = 5 )
```

```
## [1] "Cluster: 1"
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```r
lapply(tmp, summary)
```

```
## $markers
##    Length     Class      Mode 
##         5 character character 
## 
## $tree
##                           Length Class         Mode     
## predictions               80     factor        numeric  
## num.trees                  1     -none-        numeric  
## num.independent.variables  1     -none-        numeric  
## mtry                       1     -none-        numeric  
## min.node.size              1     -none-        numeric  
## variable.importance       31     -none-        numeric  
## prediction.error           1     -none-        numeric  
## forest                    10     ranger.forest list     
## confusion.matrix           4     table         numeric  
## splitrule                  1     -none-        character
## treetype                   1     -none-        character
## call                       5     -none-        call     
## importance.mode            1     -none-        character
## num.samples                1     -none-        numeric  
## replace                    1     -none-        logical  
## 
## $cluster
##    Length     Class      Mode 
##         1 character character 
## 
## $df_importances
##    importance           pvalue          marker         
##  Min.   :-0.17503   Min.   :0.0000   Length:31         
##  1st Qu.: 0.02163   1st Qu.:0.0000   Class :character  
##  Median : 0.19782   Median :0.0000   Mode  :character  
##  Mean   : 0.38617   Mean   :0.2258                     
##  3rd Qu.: 0.57540   3rd Qu.:0.4000                     
##  Max.   : 1.82855   Max.   :1.0000
```

## Now applying that function to all the clusters


```r
classif_info <- lapply(unique(pbmc_df$ident), 
       function(x) find_cluster_classification_forest(
         cluster_id = x, 
         class_df = pbmc_df, 
         markers = rownames(df_importances)))
```

```
## [1] "Cluster: 0"
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```
## [1] "Cluster: 1"
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```
## [1] "Cluster: 2"
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```
## [1] "Cluster: 3"
```

```
## Warning in ranger::importance_pvalues(forest): Only few negative importance
## values found, inaccurate p-values. Consider the 'altmann' approach.
```

```r
# This name conversion step carried out by `plyr::mapvalues`
# is required because `FeaturePlot` uses the
# original seurat object which still has the old names
for (i in classif_info) {
   i$markers %>% 
    plyr::mapvalues(
      from = newnames,
      to = original_names, 
      warn_missing = FALSE) %>%
    {
       for (marker in .) {
            FeaturePlot(
                object = pbmc_small,
                marker,
                cols.use = viridis::viridis(7),
                reduction.use = "tsne")
        }
   }
}
```

![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-1.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-2.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-3.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-4.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-5.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-6.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-7.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-8.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-9.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-10.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-11.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-12.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-13.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-14.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-15.png)<!-- -->![](tree_based_variable_selection_files/figure-html/unnamed-chunk-10-16.png)<!-- -->

