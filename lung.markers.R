
## For each category, keep the markers less than 12
## otherwise the figure will skew

## Travaglini et al. human lung cell atlas
general = c(
"PTPRC",									## CD45, immune cell
"EPCAM",									## epithelial
"PECAM1", "VWF",							## endothelial
"COL3A1",									## fibroblast
"MME"										## stroma
);

###########################################################
epithelium.major = c(
"FOXJ1", "TUBB1", "TP73", "CCDC78","TPPP3",	#Ciliated
"KRT5", "KRT14", "TP63", "DAPL1",			#Basal
"CYP2F1", "SCGB3A2", "CCKAR", "SCGB1A1"		#Club
);

epithelium.minor = c(
"MUC5B", "MUC5AC", "SPDEF",				#Goblet
										#MUC5B is also a marker for Mucous
"PRR4", "LPO", "LTF",					#Serous
"CFTR", "FOXI1", "ASCL3",				#Ionocyte
"CALCA", "CHGA", "ASCL1",				#Neuroendocrine
"DCLK1", "ASCL2"						#Tuft Cell
);

epithelium.AT = c(
"AGER", "PDPN", "CLIC5",					#Alveolar Type 1
"SFTPB", "SFTPD", "MUC1", "ETV5", "LAMP3"	#Alveolar Type 2, SFTPC removed
);

###########################################################
endothelium = c(
#"VWF",									##kind of general
"GJA5", "BMX",							#Artery Cell
"ACKR1",								#Vein Cell
"CA4",									#Capillary Cell
"PROX1", "PDPN"							#Lymphatic Cell
);

###########################################################
stroma.1 = c(
"CNN1", "ACTA2", "TAGLN",				#Smooth Muscle
"RGS5",									#Vascular Smooth Muscle
"LGR6",									#Airway Smooth Muscle
"COL1A1", "PDGFRA",						#Fibroblast
"ELN", "ACTA2",							#Myofibroblast
"PLIN2", "APOE"							#Lipofibroblast
);

stroma.2 = c(
"CSPG4", "TRPC6", "PDGFRB",				#Pericyte
"MSLN", "UPK3B", "WT1"					#Mesothelial Cell
);

###########################################################
immune.B.NK = c(
"CD79A", "CD24", "MS4A1", "CD19",		#B
"CD27", "SLAMF7", "IGHG4",				#Plasma Cell (with CD79A)
"GNLY", "KLRD1", "NKG7"#, "TYROBP"		#NK Cell (CD3-)
);

immune.T = c(
"CD3D", "CD3E",							#T
"CD4", "IL7R",							#CD4+ T
"CD8A",									#CD8+ T
"FCER1G", "TYROBP"						#NK T Cell
);

immune.T.subtype = c(
"COTL1", "LDHB",						#CD4+ Mem/Eff Cell, S100A4
"CCR7", "LEF1",							#CD4+ Naive T Cell
"GZMK", "DUSP2",						#CD8+ Mem/Eff T Cell
"GZMH", "GZMB"#,						#CD8+ Naive T Cell
);

immune.granulocyte = c(
"CD163", "AIF1",						## general myeloid marker
"S100A8", "S100A9", "IFITM2", "FCGR3B",	#Neutrophil	
"MS4A2", "CPA3", "TPSAB1",				#Basophil & Mast cell
"TPSB2",								#Mast Cell, very similar to Basophil
"SIGLEC8"								#Eosinophil
);

immune.dendritic = c(
"FCER1A", "CST3",						#Dendritic Cell
"LILRB4", "IRF8", "LILRA4",				#Plasmacytoid Dendritic Cell
"HLA-DRA",								#Myeloid Dendritic Cell
"CLEC9A", "LAMP3", "CLEC10A",			#Myeloid Dendritic Cell 1
"CD1C", "PLD4"							#Myeloid Dendritic Cell 2	
);

immune.mono.macro = c(
"CD14", "FCN1", "LYZ",					#classical CD14 Monocyte
"S100A8",								#intermediate Monocyte
"FCGR3A", "MS4A7",						#non-classical CD16 Monocyte
"MARCO", "MSR1", "MRC1", "CD68",		#macrophage
"NRGN", "PPBP"#, "PF4"					#megakaryocyte/Platelet, OST4
);

immune.RedBlood = c(
"HBA1", "HBA2", "HBB", "HBD", "HBE1",
"HBG1", "HBG2", "HBM", "HBQ1","HBZ"
);

