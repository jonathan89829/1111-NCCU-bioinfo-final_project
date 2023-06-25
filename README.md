# Expression Divergence of Chemosensory Genes between Drosophila sechellia and Its Sibling Species and Its Implications for Host Shift

### Members  
*陳羽暉, 111753156  
*邱鉅樺, 111753166  
*山本遙人, 111753170  


### Demo 
You might provide an example commend or few commends to reproduce your analysis, i.e., the following R script

```R
Rscript final_code.r 

#sample code###
# To show the heatmap of Or genes
heatmap(heat_or_scale)
# To show the Tree of Or genes (boostrap)
plot(or_bp)
# To show the the deg table 
deg_or
```

## Folder organization and its related information
idea by Noble WS (2009) [A Quick Guide to Organizing Computational Biology Projects.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424) PLoS Comput Biol 5(7): e1000424.

### docs
* Your presentation, 1112_bioinformatics_FP_studentID.ppt/pptx/pdf (i.e.,1112_bioinformatics_FP_556688.ppt), by **01.12**
* Any related document for the project
  * software user guide (pvclust, limma)

### data
* Source
  * FPKM tables of D.sec, D.sim, and D.mel  
  downloded from Supplementary data on https://academic.oup.com/gbe/article/7/10/2843/2465881#supplementary-data  
* Format
  * .csv (csv format)
* Size
  * size of each csv file is 1~3 kb

### code
* Which packages do you use? 
  * original packages in the paper (pvclust)
  * additional packages you found (limma) 
* Analysis steps  
  * Use final_code.r to reproduce figures from the paper and analize FPKM

### results
* Which part of the paper do you reproduce?
  * Heat maps of expression profiles of Or genes, Ir genes, and Obp genes. 
  * Cluster dendrogram with P values of Or genes, Ir genes, and Obp genes.(Boostrapped Tree)
* Any improvement or change by your package?
  * using limma packeage instead of deseq2 (Change)

## References
* Packages you use  
  * pvclust  
    https://cran.r-project.org/web/packages/pvclust/index.html  
  * limma  
    https://bioconductor.org/packages/release/bioc/html/limma.html  
* Related publications
  * Shiao MS, et al. 2013. Transcriptional profiling of adult Drosophila antennae by high-throughput sequencing. Zool Stud. 52:42–51.  
    https://zoologicalstudies.springeropen.com/articles/10.1186/1810-522X-52-42    
    (This paper provides gene information and FPKM of D.mel)
