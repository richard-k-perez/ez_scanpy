# Simple analysis script to demonstrate proper code syntax
import ez_scanpy as ezsc


###################################################
# Remove doublets and create initial .h5ad object #
###################################################
h5csvpath = '/ye/yelabstore2/10x.lupus/eqtls/demux.v2/v2.batches.h5.csv'
no_norm = 'SLEcrossX_nonorm.h5ad'
ezsc.remove_doublets(h5csvpath, no_norm)


#########################################################
# Remove doublets and create initial .h5ad object flare #
#########################################################
no_norm_flare = 'SLE_flare_nonorm.h5ad'
flarepaths = [
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_1/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_2/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_3/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_4/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_1/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_2/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_3/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_4/outs/filtered_gene_bc_matrices_h5.h5']

flaregroups = ['3-1', '3-2', '3-3', '3-4', '4-1', '4-2', '4-3', '4-4']
demuxletpaths = [
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_1/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_2/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_3/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_4/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_1/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_2/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_3/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_4/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best']

ezsc.remove_doublets_flare(flarepaths=flarepaths, flaregps=flaregroups, bestfilepath=demuxletpaths, savepath=no_norm_flare)



#########################################
# Basic processing and remove platelets #
#########################################
processed = 'SLEcrossX_processed.h5ad'
ezsc.basicprocessing_noplatelets(no_norm, processed)


###################################################
# Basic processing and with platelets             #
###################################################
processed_plat = 'SLEcrossX_processed_with_platelets.h5ad'
processed_plat_flare = 'SLEcrossX_processed_with_platelets_flare.h5ad'
ezsc.basic_processing(no_norm, processed_plat)
ezsc.basic_processing(no_norm_flare, processed_plat_flare)

##################
# Basic analysis #
##################
ezsc.basic_analysis(processed)
ezsc.basic_analysis(processed_plat)
ezsc.basic_analysis(processed_plat_flare)


##########################################
# Assign cell identity                   #
# Gene markers for cell type populations #
##########################################
annotation = {'CD8+ T': ['CD8A', 'CD8B', 'CD3D'], 'CD4+ T': ['CD3D', 'ANK3', 'IL7R'],
              'CD14 Mono': ['CD14', 'LYZ', 'CCL2', 'S100A9'], 'CD16 Mono': ['FCGR3A', 'MS4A7', 'VMO1'],
              'B': ['CD19', 'MS4A1', 'CD79A'], 'NK': ['KLRF1', 'GNLY', 'NKG7'],
              'cDC': ['HLA.DQA1', 'FCER1A', 'GPR183'], 'pDC': ['CLEC4C', 'TSPAN13', 'IGJ'],
              'MK': ['PPBP', 'GNG11']}
ezsc.cell_identity(processed, annotation)  # Must run 'Basic analysis' first
ezsc.cell_identity(processed_plat, annotation)  # Must run 'Basic analysis' first
ezsc.cell_identity(processed_plat_flare, annotation)  # Must run 'Basic analysis' first


#################################
#     Subpopulation analysis    #
# Must run basic analysis first #
#################################
Mpath = 'Mono.h5ad'
Bpath = 'Bcell.h5ad'
T4path = 'T4cell.h5ad'
T8path = 'T8cell.h5ad'
NKpath = 'NK.h5ad'
DCpath = 'DC.h5ad'
MKpath = 'MK.h5ad'
ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=Mpath)
ezsc.subpopulation_analysis(cell_type_IDs=['CD4+ T'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=T4path)
ezsc.subpopulation_analysis(cell_type_IDs=['CD8+ T'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=T8path)
ezsc.subpopulation_analysis(cell_type_IDs=['NK'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=NKpath)
ezsc.subpopulation_analysis(cell_type_IDs=['B'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=Bpath)
ezsc.subpopulation_analysis(cell_type_IDs=['cDC', 'pDC'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=DCpath)
ezsc.subpopulation_analysis(cell_type_IDs=['MK'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=MKpath)


#################################
# Subpopulation analysis flare  #
# Must run basic analysis first #
#################################
Mpath = 'Mono_flare.h5ad'
Bpath = 'Bcell_flare.h5ad'
T4path = 'T4cell_flare.h5ad'
T8path = 'T8cell_flare.h5ad'
NKpath = 'NK_flare.h5ad'
DCpath = 'DC_flare.h5ad'
MKpath = 'MK_flare.h5ad'
ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=Mpath)
ezsc.subpopulation_analysis(cell_type_IDs=['CD4+ T'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=T4path)
ezsc.subpopulation_analysis(cell_type_IDs=['CD8+ T'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=T8path)
ezsc.subpopulation_analysis(cell_type_IDs=['NK'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=NKpath)
ezsc.subpopulation_analysis(cell_type_IDs=['B'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=Bpath)
ezsc.subpopulation_analysis(cell_type_IDs=['cDC', 'pDC'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=DCpath)
ezsc.subpopulation_analysis(cell_type_IDs=['MK'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=MKpath)


##########################################
# Subpopulation analysis regress out IFN #
##########################################
Mpath = 'Mono_regressIFN.h5ad'
ezsc.subpopulation_analysis_regressIFN(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=Mpath)


##########################################
# Subpopulation analysis remove IFN      #
##########################################
Mpath = 'Mono_remove_IFNdiffexp.h5ad'
ezsc.subpopulation_analysis_removeIFNGenes(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed_plat, no_norm_filepath=no_norm, savepath=Mpath)
Mpath = 'Mono_flare_remove_IFNdiffexp.h5ad'
ezsc.subpopulation_analysis_removeIFNGenes(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=Mpath)


##########################################
# Subpopulation analysis by dz covariate #
##########################################
Mpath_healthy = 'Mono_healthy.h5ad'
Mpath_sle = 'Mono_sle.h5ad'
ezsc.subpopulation_analysis_covariate(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], covariate='disease_cov', covariate_grp= 'healthy', filepath=processed_plat, no_norm_filepath=no_norm, savepath=Mpath_healthy)
ezsc.subpopulation_analysis_covariate(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], covariate='disease_cov', covariate_grp= 'sle', filepath=processed_plat, no_norm_filepath=no_norm, savepath=Mpath_sle)