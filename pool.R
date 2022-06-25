#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for pool scRNA-seq data
#
suppressPackageStartupMessages({
	library(Seurat);library(plyr);library(dplyr);library(cowplot);
	library(ggplot2);library(gridExtra);library(ggthemes);library(patchwork);
	library(gridExtra);library(forcats);library(SingleR);
});
options( stringsAsFactors=F );

sid="PAM.with.Kim";
outfileName=paste0(sid, '.pool.pdf');
pdf( outfileName, width=10, height=8 );

allsids = read.table( "PAM.sample.info", head=T, row.names=1 );
samples = rownames( allsids ); 
valid   = list();
sc.list = list();
to.merge= c();

## load files and add sampleID to each barcode
for(i in 1:length(samples)) {
	if( file.exists( paste0("source/", samples[i], ".h5") ) ) {
		counts = Read10X_h5( paste0("source/", samples[i], ".h5") );
	} else {
		counts = read.csv( paste0("source/", samples[i], ".csv.gz"), row.names=1 );
	}
	colnames(counts) = paste0( samples[i], "_", colnames(counts) );
	sce = CreateSeuratObject(counts, project=samples[i]);
	SingleR = read.table( paste0("source/", samples[i], ".SingleR.result"));
	valid[[i]] = paste0(samples[i], "_", SingleR[,1]);
	sc.list[[i]] = subset(sce, cells=valid[[i]]);

	if( i != 1 ) {
		to.merge = c( to.merge, sc.list[[i]] );
	}
}

## this is Seurat's merge
sce = merge(sc.list[[1]], y=to.merge);

## add a label of patient
sce$patient = sce$orig.ident;

## regress out unwanted genes
all.gene = rownames( sce );
RP.genes = grep(rownames(sce), pattern="^RP[SL]", value=T, invert=F, ignore.case=T);
MT.genes = grep(rownames(sce), pattern="^MT-", value=T, invert=F, ignore.case=T);
HB.genes = grep(rownames(sce), pattern="^HB[AB][12]?$", value=T, invert=F, ignore.case=T);

HIST.genes = grep(rownames(sce), pattern="^HIST", value=T, invert=F, ignore.case=T);
cc.genes = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes);

keep.genes = all.gene[!(all.gene %in% c(RP.genes, MT.genes, HB.genes, HIST.genes, cc.genes))];
sce = subset(sce, features=keep.genes);

## Perform cell cycle scoring, will plot this afterwards
sce = CellCycleScoring(sce, g2m.features=cc.genes.updated.2019$g2m.genes,
				s.features=cc.genes.updated.2019$s.genes);

## use Harmony package
library( harmony );
sce = NormalizeData(sce, normalization.method='LogNormalize', scale.factor=10000, verbose=F);
sce = FindVariableFeatures(sce, selection.method="vst", verbose=F);
sce = ScaleData( sce, verbose=F );
sce = RunPCA(sce, npcs=50, verbose=F);

sce = RunHarmony(sce, "patient", plot_convergence=T);
sce = FindNeighbors(sce, reduction="harmony", dims=1:30, verbose=F);
sce = FindClusters(sce, resolution=c(0.4, 0.5), verbose=F);
#identity( sce );

sce = RunUMAP(sce, reduction="harmony", dims=1:30, verbose=F);
sce = RunTSNE(sce, reduction="harmony", dims=1:30, verbose=F);

DimPlot(sce, reduction="tsne", group.by="patient", pt.size=0.1);
DimPlot(sce, reduction="tsne", group.by="RNA_snn_res.0.4", label=T, pt.size=0.1);
DimPlot(sce, reduction="tsne", group.by="RNA_snn_res.0.5", label=T, pt.size=0.1);

for( i in 1:length(samples) ) {
	per.cell = DimPlot(sce, reduction="tsne", cells=valid[[i]])+ggtitle(samples[i]);
	grid.arrange( per.cell );
}

## plot phase for cells, which could be used to evaluate the performance of integration
DimPlot(sce, reduction="tsne", group.by="Phase");

## normal analysis
sce = SetIdent(sce,value='RNA_snn_res.0.5');
sce$seurat_clusters = sce$RNA_snn_res.0.5;

## find cluster markers
cm = FindAllMarkers(sce, logfc.threshold=0.2, min.pct=0.25, only.pos=T);
cm = cm[with(cm,order(cluster, p_val_adj, -avg_log2FC)),];
#divide cluster_markers by cluster.
#cm_list = split(cm, cm$cluster);

## plot top markers
source( "/mnt/sunkun/scRNA.library/bubblePlot.R" );
top10 = cm %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC);
ngene = length(top10$gene);
part1 = round(ngene * 0.5);
part2 = ngene - round(ngene * 0.5);
# wt could be avg_logFC or -p_val_adj
DoHeatmap(sce, features=head(top10$gene, n=part1)) + NoLegend();
qqq(sce, features=unique(top10$gene)[1 : part1],
	cols=c('#08519c','#b30000'), group.by='seurat_clusters',
	dot.scale=4, scale=F, col.min=0, col.max=3, scale.min=0, scale.max=100)+
xlab('Cluster ID')+
ggtitle('Top10 marker genes per cluster, part1')+
theme(axis.text.y=element_text(size=8));

DoHeatmap(sce, features=tail(top10$gene, n=part2)) + NoLegend();
qqq(sce, features=unique(top10$gene)[part1+1 : ngene],
	cols=c('#08519c','#b30000'),
	group.by='seurat_clusters', dot.scale=4., scale=F,
	col.min=0, col.max=3, scale.min=0,scale.max=100)+
xlab('Cluster ID')+
ggtitle('Top10 marker genes per cluster, part2')+
theme(axis.text.y=element_text(size = 8));

## plot known markers
## by default, UMAP results are shown
source("lung.markers.R");
VlnPlot(sce, features=general, pt.size=0, ncol=2)
VlnPlot(sce, features=epithelium.major, pt.size=0, ncol=2)
VlnPlot(sce, features=epithelium.minor, pt.size=0, ncol=2)
VlnPlot(sce, features=epithelium.AT, pt.size=0, ncol=2)
VlnPlot(sce, features=endothelium, pt.size=0, ncol=2)
VlnPlot(sce, features=stroma.1, pt.size=0, ncol=2)
VlnPlot(sce, features=stroma.2, pt.size=0, ncol=2)
VlnPlot(sce, features=immune.B.NK, pt.size=0, ncol=2)
VlnPlot(sce, features=immune.T, pt.size=0, ncol=2)
VlnPlot(sce, features=immune.T.subtype, pt.size=0, ncol=2)
VlnPlot(sce, features=immune.mono.macro, pt.size=0, ncol=2)
VlnPlot(sce, features=immune.dendritic, pt.size=0, ncol=2)
VlnPlot(sce, features=immune.granulocyte, pt.size=0, ncol=2)

## plot known genes
FeaturePlot(sce, reduction="tsne",
			features=c("MARCO", "ASPM", "TK1", "FAM111A"));

## automated annotation using SingleR
load("ref_Human_all.RData");
sce_for_SingleR = GetAssayData(sce, slot="data");
pred = SingleR(test=sce_for_SingleR, ref=ref_Human_all,
			   labels=ref_Human_all$label.main, method="cluster",
			   clusters=sce@meta.data$seurat_clusters,
			   assay.type.test="logcounts", assay.type.ref="logcounts");
celltype = data.frame(ClusterID=rownames(pred),
					  celltype=pred$labels, stringsAsFactors=F);
sce@meta.data$SingleR = celltype[match(sce@meta.data$seurat_clusters,
									 celltype$ClusterID), 'celltype'];
TSNEPlot(sce, group.by="SingleR", label=T, label.size=5);
write.table(sce$SingleR, file=paste0(sid, ".SingleR.result"),
			quote=F, sep='\t', row.names=T, col.names=F);

## save project
save(sce, file=paste0(sid, ".rds"));
write.table(cm, file=paste0(sid, ".cluster.markers"),
			quote=F, sep='\t', row.names=F, col.names=T);
