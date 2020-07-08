### 0.缘起
在[xena]([https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443](https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
)上面可以看到，TCGA和GTex、Target数据库的测序数据已经被重新计算整合到了一起，提供了各种格式的文件。
![](https://upload-images.jianshu.io/upload_images/9475888-91b6e7a9c44fe70b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
这里上游分析使用的是RSEM，而不是featurecout，导致得到的数据并不是标准的count值，是非整数。
### 1.RSEM 对应的差异分析包是EBSeq
RSEM (RNA-Seq by Expectation-Maximization)

关于它的下游分析，官网建议使用的R包是EBSeq：
[EBSeq](http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html)

但市面上公认最好的差异分析R包是DESeq2,limma,edgeR。有没有办法将RSEM的counts矩阵交给三大R包来处理呢？

### 2.使用limma和edgeR
这个问题刚好是关于TCGA的RSEM数据，
https://support.bioconductor.org/p/63981/#64004
limma的作者亲自回答了：
>The RSEM expected counts from the TCGA project will work fine with either limma-voom or edgeR. However, with such a large number of samples, limma-voom is easily the best choice from a computational point of view.

limma-voom是更好的选择。
关于expected_count和norm_count在这里也有讨论，即edgeR只能用expected，vomm理论上可以使用norm_count（只是可以不是必须）。

#### 3.使用Deseq2
https://support.bioconductor.org/p/94003/#94028

![](https://upload-images.jianshu.io/upload_images/9475888-cdca7fcb72ad98cf.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

作者说最好的办法是用tximport进行转换，其次就是四舍五入对矩阵进行取整，然后用 DESeqDataSetFromMatrix去导入即可。

### 4.tximport是Deseq2作者写的R包
[tximport](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)：将各种上游分析软件得到的数据转换给三大R包使用。
其中讲了如何将ERR格式的数据导入R，并生成矩阵。DESeq2 和edgeR除了需要count矩阵，还需要一个length矩阵，而limma则是需要经过一步 "scaledTPM" 或"lengthScaledTPM"转换，需要另外一个矩阵来做校正。在xena中的数据我们只能拿到一个count矩阵，因此这个暂时用不上，但不妨碍它是个好东西。


