#
# R program to call differentially expressed genes in metastasis vs CTC
# Only those have the same cluster-of-origin are compared
# This program should be run under scXXX-XXX/ directory
# Note that this program takes long time to run
#
suppressPackageStartupMessages({
library(Seurat);
library(DEsingle);
library(BiocParallel);
});

load( "PAM.with.Kim.reanalysis.rds" );
sce$seurat_clusters = sce$RNA_snn_res.0.4;

allcells = colnames( sce );

## Set the parameters and register the back-end to be used
param = MulticoreParam( workers=8, progressbar=T );
register( param );

for ( cluster in 0:16 ) {
	here = subset( sce, seurat_clusters == cluster );

	PAM = subset( here, patient=="sc27" | patient=="sc29" );
	CTR = subset( here, patient!="sc27" & patient!="sc29" );

	## prepare data for DEsingle
	counts = cbind( CTR@assays[["RNA"]]@counts, PAM@assays[["RNA"]]@counts );
	group  = factor( c(rep(1,ncol(CTR)), rep(2,ncol(PAM))) );

	## Detecting the DE genes in parallelization with 32 cores,
	## this step uses lots of memory (~80g) and time
	DE.r = DEsingle( counts=counts, group=group, parallel=T, BPPARAM=param );

	DE.c = DEtype(results=DE.r, threshold=0.05);
	write.table( DE.c, file=paste0("DEsingle.c", cluster), sep="\t", quote=F );

	## calculate log-normalized fold change
	fcfix = 0.1;
	logfc = log2((DE.c$norm_total_mean_2+fcfix)/(DE.c$norm_total_mean_1+fcfix));
	DE.k  = data.frame( DE.c$norm_total_mean_1, DE.c$norm_total_mean_2,
					   logfc, DE.c$pvalue.adj.FDR, DE.c$Type );
	rownames(DE.k) = rownames( DE.c );
	colnames(DE.k) = c("Exp1", "Exp2", "Log2FC", "FDR", "Type");

	DE.s = subset( DE.k, FDR < 0.05 & (Log2FC > 1 | Log2FC < -1) );
	write.table( DE.s, file=paste0("DEsingle.c", cluster, ".sig"), sep="\t", quote=F );
}

