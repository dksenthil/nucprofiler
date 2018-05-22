## ----eval=TRUE-----------------------------------------------------------
library(DNAshapeR)
#pred <- getShape("SingleSeqsample.fa")
pred <- getShape("GCGTA.fasta")

## ----fig.width=7, fig.height=7, fig.align='center', eval=TRUE------------
plotShape(pred$MGW)
dev.copy(jpeg,filename="plot_avg_MGW_500.jpg");
dev.off ();
