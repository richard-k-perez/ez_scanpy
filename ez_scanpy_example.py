# Simple analysis script to demonstrate proper code syntax
import ez_scanpy as ezsc


###################################################
# Remove doublets and create initial .h5ad object #
###################################################
h5csvpath = '/ye/yelabstore2/10x.lupus/eqtls/demux.v2/v2.batches.h5.csv'
no_norm = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/SLEcrossX_nonorm.h5ad'
#ezsc.remove_doublets(h5csvpath, no_norm)


#########################################
# Basic processing and remove platelets #
#########################################
processed = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/SLEcrossX_processed.h5ad'
#ezsc.basicprocessing_noplatelets(no_norm, processed)


###################################################
# Basic processing and without removing platelets #
###################################################
processed_plat = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/SLEcrossX_processed_with_platelets.h5ad'
#ezsc.basic_processing(no_norm, processed_plat)


##################
# Basic analysis #
##################
#ezsc.basic_analysis(processed)
#ezsc.basic_analysis(processed_plat)


##########################################
# Assign cell identity                   #
# Gene markers for cell type populations #
##########################################
annotation = {'CD8+ T': ['CD8A', 'CD8B', 'CD3D'], 'CD4+ T': ['CD3D', 'ANK3', 'IL7R'],
              'CD14 Mono': ['CD14', 'LYZ', 'CCL2', 'S100A9'], 'CD16 Mono': ['FCGR3A', 'MS4A7', 'VMO1'],
              'B': ['CD19', 'MS4A1', 'CD79A'], 'NK': ['KLRF1', 'GNLY', 'NKG7'],
              'cDC': ['HLA.DQA1', 'FCER1A', 'GPR183'], 'pDC': ['CLEC4C', 'TSPAN13', 'IGJ'],
              'MK': ['PPBP', 'GNG11']}
#ezsc.cell_identity(processed, annotation)  # Must run 'Basic analysis' first
#ezsc.cell_identity(processed_plat, annotation)  # Must run 'Basic analysis' first


#################################
#     Subpopulation analysis    #
# Must run basic analysis first #
#################################
Mpath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/Mono.h5ad'
Bpath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/Bcell.h5ad'
T4path = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/T4cell.h5ad'
T8path = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/T8cell.h5ad'
NKpath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/NK.h5ad'
DCpath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/DC.h5ad'
ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed, no_norm_filepath=no_norm, savepath=Mpath)
ezsc.subpopulation_analysis(cell_type_IDs=['CD4+ T'], filepath=processed, no_norm_filepath=no_norm, savepath=T4path)
ezsc.subpopulation_analysis(cell_type_IDs=['CD8+ T'], filepath=processed, no_norm_filepath=no_norm, savepath=T8path)
ezsc.subpopulation_analysis(cell_type_IDs=['NK'], filepath=processed, no_norm_filepath=no_norm, savepath=NKpath)
ezsc.subpopulation_analysis(cell_type_IDs=['B'], filepath=processed, no_norm_filepath=no_norm, savepath=Bpath)
ezsc.subpopulation_analysis(cell_type_IDs=['cDC', 'pDC'], filepath=processed, no_norm_filepath=no_norm, savepath=DCpath)


##########################################
# Subpopulation analysis regress out IFN #
##########################################
Mpath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/Mono.regressIFN.h5ad'
IFN = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/IFN.csv'
ezsc.subpopulation_analysis_regressIFN(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed, no_norm_filepath=no_norm, IFNpath=IFN, savepath=Mpath)