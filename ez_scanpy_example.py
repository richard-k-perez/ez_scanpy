# Simple analysis script to demonstrate proper code syntax
import ez_scanpy as ezsc


###################################################
# Remove doublets and create initial .h5ad object #
###################################################
h5csvpath = '/ye/yelabstore2/10x.lupus/eqtls/demux.v2/v2.batches.h5.csv'
no_norm = 'SLEcrossX_nonorm.h5ad'
#ezsc.remove_doublets(h5csvpath, no_norm)


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

flaregroups = ['flare-3-1', 'flare-3-2', 'flare-3-3', 'flare-3-4', 'flare-4-1', 'flare-4-2', 'flare-4-3', 'flare-4-4']
demuxletpaths = [
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_1/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_2/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_3/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_3/AbFlare-3_4/demux_out_raw_oldvcfs/demuxlet.raw.flare3samp.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_1/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_2/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_3/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best',
        '/ye/yelabstore2/george/Sasha_Flare_Studies/10xcount/run_4/AbFlare-4_4/demux_out_raw_oldvcfs/demuxlet.raw.flare4samps.best']

#ezsc.remove_doublets_flare(flarepaths=flarepaths, flaregps=flaregroups, bestfilepath=demuxletpaths, savepath=no_norm_flare)


##########################################################
# Remove doublets and create initial .h5ad object Immvar #
##########################################################
no_norm_immvar = 'immvar_nonorm.h5ad'
immvarpaths = [
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/count_AH7TNHDMXX_YE_8-30_1/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/count_AH7TNHDMXX_YE_8-30_2/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/count_AH7TNHDMXX_YE_8-30_3/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/count_AH7TNHDMXX_YE_8-30_4/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/count_AHCM2CDMXX_YE_0831_1/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/count_AHCM2CDMXX_YE_0831_2/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/count_AHCM2CDMXX_YE_0831_3/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/count_AHCM2CDMXX_YE_0831_4/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/count_BH7YT2DMXX_YE_0907_1/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/count_BH7YT2DMXX_YE_0907_2/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/count_BH7YT2DMXX_YE_0907_3/outs/filtered_gene_bc_matrices_h5.h5',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/count_BH7YT2DMXX_YE_0907_4/outs/filtered_gene_bc_matrices_h5.h5'
]
immvargroups = [
        'YE_8-30_1','YE_8-30_2','YE_8-30_3','YE_8-30_4',
        'YE_0831_1','YE_0831_2','YE_0831_3','YE_0831_4',
        'YE_0907_1','YE_0907_2','YE_0907_3','YE_0907_4',
]
immvardemuxletpaths = [
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/cramore_AH7TNHDMXX_YE_8-30_1/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/cramore_AH7TNHDMXX_YE_8-30_2/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/cramore_AH7TNHDMXX_YE_8-30_3/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AH7TNHDMXX/cramore_AH7TNHDMXX_YE_8-30_4/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/cramore_AHCM2CDMXX_YE_0831_1/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/cramore_AHCM2CDMXX_YE_0831_2/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/cramore_AHCM2CDMXX_YE_0831_3/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/AHCM2CDMXX/cramore_AHCM2CDMXX_YE_0831_4/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/cramore_BH7YT2DMXX_YE_0907_1/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/cramore_BH7YT2DMXX_YE_0907_2/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/cramore_BH7YT2DMXX_YE_0907_3/.0.1.best',
        '/ye/yelabstore2/gracieg/scimmvar/BH7YT2DMXX/cramore_BH7YT2DMXX_YE_0907_4/.0.1.best'
        ]
#ezsc.remove_doublets_immvar(paths=immvarpaths, gps=immvargroups, bestfilepath=immvardemuxletpaths, savepath=no_norm_immvar)


#########################################################
# Remove doublets and create initial .h5ad object CML #
#########################################################
no_norm_CML = 'CML_nonorm.h5ad'
CMLpaths = [
        "/ye/yelabstore2/cml/10x_SCS_hg19/112117_KJ_1/112117_KJ_1_run3/outs/raw_gene_bc_matrices_h5.h5",
        "/ye/yelabstore2/cml/10x_SCS_hg19/112117_KJ_2/112117_KJ_2_run3/outs/raw_gene_bc_matrices_h5.h5",
        "/ye/yelabstore2/cml/10x_SCS_hg19/112917_KJ_1/112917_KJ_1_run3/outs/raw_gene_bc_matrices_h5.h5",
        "/ye/yelabstore2/cml/10x_SCS_hg19/112917_KJ_2/112917_KJ_2_run3/outs/raw_gene_bc_matrices_h5.h5"]

CMLgroups = ['CML-117KJ-1', 'CML-117KJ-2', 'CML-917KJ-1', 'CML-917KJ-2']
demuxletpaths = [
        "/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/CML/10x_112117_KJ_1_demuxlet_run3.best",
        "/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/CML/10x_112117_KJ_2_demuxlet_run3.best",
        "/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/CML/10x_112917_KJ_1_demuxlet_run3.best",
        "/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/CML/10x_112917_KJ_2_demuxlet_run3.best",]

#ezsc.remove_doublets_cml(flarepaths=CMLpaths, flaregps=CMLgroups, bestfilepath=demuxletpaths, savepath=no_norm_CML)


#########################################
#    Combine processing of cohorts      #
#########################################
savepath_CLUES_ImmVar = 'CLUESImmVar_nonorm.h5ad'
#ezsc.combine_studies_cml(crossXpath=no_norm, flarepath=no_norm_flare, cmlpath=no_norm_CML, savepath=savepath_combined_cml)
#ezsc.combine_CLUES_Immvar(crossXpath=no_norm, immvarpath=no_norm_immvar, savepath=savepath_CLUES_ImmVar)


#########################################
#    Remove anndata atrributes          #
#########################################
# In keep, insert all attributes to retain
#ezsc.clean_anndata(savepath_CLUES_ImmVar, clean=['var', 'obs'],keep={'batch_cov', 'disease_cov', 'ind_cov', 'pop_cov', 'well', 'study_cov'})


#########################################
# Basic processing and remove platelets #
#########################################
#processed = 'SLEcrossX_processed.h5ad'
processed_CLUESImmVar = 'CLUESImmVar_processed.h5ad'
#ezsc.basicprocessing_noplatelets(savepath_CLUES_ImmVar, processed_CLUESImmVar)

###################################################
# Basic processing and with platelets             #
###################################################
processed_plat_combined_wCML = 'SLE_combined_processed_with_platelets_andCML.h5ad'
#processed_plat = 'SLEcrossX_processed_with_platelets.h5ad'
#processed_plat_flare = 'SLE_flare_processed_with_platelets.h5ad'
#ezsc.basic_processing(savepath_combined_cml, processed_plat_combined_wCML)
#ezsc.basic_processing(savepath_combined, processed_plat_combined)
#ezsc.basic_processing(no_norm, processed)
#ezsc.basic_processing(no_norm, processed_plat)
#ezsc.basic_processing(no_norm_flare, processed_plat_flare)

##################
# Basic analysis #
##################
#ezsc.basic_analysis(processed_plat_combined_wCML)
#ezsc.basic_analysis(processed_plat_combined)
#ezsc.basic_analysis(processed_plat)
#ezsc.basic_analysis(processed_CLUESImmVar)
#ezsc.basic_analysis(processed_plat_flare)

##################
#  RNA Velocity  #
##################
loompath =  '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/velocyto/loom/crossx.loom'
refh5adpath =  '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/SLEcrossX_processed.h5ad'
rnavelosavepath = 'crossx_rnavelo_nonorm.h5ad'
rnaveloprocessed = 'crossx_rnavelo_processed.h5ad'
#ezsc.add_covariates_2_rnavelocity(loompath, refh5adpath, rnavelosavepath)
#ezsc.basic_processing_rnavelocity(rnavelosavepath, rnaveloprocessed)

##########################################
# Assign cell identity                   #
# Gene markers for cell type populations #
##########################################
annotation = {'CD8+ T': ['CD8A', 'CD8B', 'CD3D'], 'CD4+ T': ['CD3D', 'ANK3', 'IL7R'],
              'CD14 Mo': ['CD14', 'LYZ', 'CCL2', 'S100A9'], 'CD16 Mo': ['FCGR3A', 'MS4A7', 'VMO1'],
              'B': ['CD19', 'MS4A1', 'CD79A'], 'NK': ['KLRF1', 'GNLY', 'NKG7'],
              'cDC': ['HLA.DQA1', 'FCER1A', 'GPR183'], 'pDC': ['CLEC4C', 'TSPAN13', 'IGJ'],
              'MK': ['PPBP', 'GNG11']}
ezsc.cell_identity(processed_CLUESImmVar, annotation)  # Must run 'Basic analysis' first
#ezsc.cell_identity(processed_plat_combined, annotation)  # Must run 'Basic analysis' first
#ezsc.cell_identity(processed_plat, annotation)  # Must run 'Basic analysis' first
#ezsc.cell_identity(processed, annotation)  # Must run 'Basic analysis' first
#ezsc.cell_identity(processed_plat_flare, annotation)  # Must run 'Basic analysis' first


#################################
#     Subpopulation analysis    #
# Must run basic analysis first #
#################################
Mpath = 'CLUESImmVarMono.h5ad'
MDCpath = 'CLUESImmVarMonoDC.h5ad'
Bpath = 'CLUESImmVarBcell.h5ad'
T4path = 'CLUESImmVarT4cell.h5ad'
T8path = 'CLUESImmVarT8cell.h5ad'
Tpath = 'CLUESImmVarTcell.h5ad'
NKpath = 'CLUESImmVarNK.h5ad'
DCpath = 'CLUESImmVarDC.h5ad'
MKpath = 'CLUESImmVarMK.h5ad'

ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mo', 'CD16 Mo', 'cDC', 'pDC'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=MDCpath)
ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mo', 'CD16 Mo'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=Mpath)
ezsc.subpopulation_analysis(cell_type_IDs=['CD4+ T'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=T4path)
ezsc.subpopulation_analysis(cell_type_IDs=['CD8+ T'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=T8path)
ezsc.subpopulation_analysis(cell_type_IDs=['CD8+ T', 'CD4+ T'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=Tpath)
ezsc.subpopulation_analysis(cell_type_IDs=['NK'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=NKpath)
ezsc.subpopulation_analysis(cell_type_IDs=['B'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=Bpath)
ezsc.subpopulation_analysis(cell_type_IDs=['cDC', 'pDC'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=DCpath)
ezsc.subpopulation_analysis(cell_type_IDs=['MK'], filepath=processed_CLUESImmVar, no_norm_filepath=savepath_CLUES_ImmVar, savepath=MKpath)

'''
#################################
# Subpopulation analysis flare  #
# Must run basic analysis first #
#################################
Mpath = 'Mono_flare.h5ad'
MDCpath = 'MonoDC_flare.h5ad'
Bpath = 'Bcell_flare.h5ad'
T4path = 'T4cell_flare.h5ad'
T8path = 'T8cell_flare.h5ad'
NKpath = 'NK_flare.h5ad'
DCpath = 'DC_flare.h5ad'
MKpath = 'MK_flare.h5ad'
ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mono', 'CD16 Mono', 'cDC', 'pDC'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=MDCpath)
ezsc.subpopulation_analysis(cell_type_IDs=['CD14 Mono', 'CD16 Mono'], filepath=processed_plat_flare, no_norm_filepath=no_norm_flare, savepath=Mpath)
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
'''
