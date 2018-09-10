# This module contains a set of functions often used in the Ye lab at UCSF when processing single cell RNA sequencing data.
# Created by Richard Perez


def remove_doublets(h5csvpath, savepath):
    import numpy as np
    import pandas as pd
    import scanpy.api as sc
    import logging
    import doubletdetection

    ##################
    # Configure file #
    ##################
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    logging.basicConfig(level=logging.INFO)

    ###########################################################
    # Compile demuxlet results and remove undetected doublets #
    ###########################################################
    batches = pd.read_csv(h5csvpath)
    i = 0
    path = batches['V3'][i]
    full_batch = pd.read_csv(batches['V1'][i] + '/' + 'covariates.csv')
    adata = sc.read_10x_h5(path, 'hg19')
    adata = adata[pd.match(full_batch['barcode'], adata.obs.index), :]
    adata.obs['disease_cov'] = full_batch['disease'].tolist()
    adata.obs['ct_cov'] = full_batch['cell'].tolist()
    adata.obs['pop_cov'] = full_batch['pop'].tolist()
    adata.obs['ind_cov'] = full_batch['ind'].tolist()
    adata.obs['batch_cov'] = [batches['V2'][i]] * len(full_batch)
    adata.obs['well'] = full_batch['well'].tolist()
    adata_singlet = adata[full_batch[full_batch['multiplets'] == 'singlet'].index.tolist(), :]
    logging.info(str('Data structure details: ' + str(adata_singlet)))
    clf = doubletdetection.BoostClassifier()
    # raw_counts is a cells by genes count matrix
    labels = clf.fit(adata_singlet.X).predict()
    f2, tsne_coords, clusters = doubletdetection.plot.tsne(adata_singlet.X, labels, random_state=1, save='tsne_test.pdf', show=False)

    logging.info(str(labels))
    nan_index = np.isnan(labels)
    labels[nan_index] = 1
    logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
    labels = labels - 1
    logging.info(str(labels))
    labels = np.abs(labels) > 0
    logging.info(str(labels))

    adata_singlet = adata_singlet[labels]
    logging.info(str('Data structure details: ' + str(adata_singlet)))
    f = doubletdetection.plot.convergence(clf, save='convergence_test.pdf', show=False)

    for i in range(1, len(batches)):
        logging.info(str('Batch ' + str(i)))
        path = batches['V3'][i]
        full_batch = pd.read_csv(batches['V1'][i] + '/' + 'covariates.csv')
        bdata = sc.read_10x_h5(path, 'hg19')
        bdata = bdata[pd.match(full_batch['barcode'], bdata.obs.index), :]
        bdata.obs['disease_cov'] = full_batch['disease'].tolist()
        bdata.obs['ct_cov'] = full_batch['cell'].tolist()
        bdata.obs['pop_cov'] = full_batch['pop'].tolist()
        bdata.obs['ind_cov'] = full_batch['ind'].tolist()
        bdata.obs['well'] = full_batch['well'].tolist()
        bdata.obs['batch_cov'] = [batches['V2'][i]] * len(full_batch)
        bdata_singlet = bdata[full_batch[full_batch['multiplets'] == 'singlet'].index.tolist(), :]

        logging.info(str('Data structure details: ' + str(bdata_singlet)))
        clf = doubletdetection.BoostClassifier()
        # raw_counts is a cells by genes count matrix
        labels = clf.fit(bdata_singlet.X).predict()
        f2, tsne_coords, clusters = doubletdetection.plot.tsne(bdata_singlet.X, labels, random_state=1, save=str(str(i) + 'b_tsne_test.pdf'), show=False)

        logging.info(str(labels))
        nan_index = np.isnan(labels)
        labels[nan_index] = 1
        logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
        labels = labels - 1
        logging.info(str(labels))
        labels = np.abs(labels) > 0
        logging.info(str(labels))

        bdata_singlet = bdata_singlet[labels]
        logging.info(str('Data structure details: ' + str(bdata_singlet)))
        f = doubletdetection.plot.convergence(clf, save=str(str(i) + 'b_convergence_test.pdf'), show=False)
        adata_singlet = adata_singlet.concatenate(bdata_singlet)
    logging.info('Saving compiled demuxlet results')
    logging.info('Making hard coded change to batch 1.10 demographics')
    adata_singlet.obs['disease_cov'][adata_singlet.obs['batch_cov'] == 'lupus1.10'] = 'sle'
    # Back the file to this directory for easy memory access
    adata_singlet.filename = 'backed.h5ad'
    adata_singlet.write(savepath)


def basic_processing(filepath, savepath):
    import numpy as np
    import scanpy.api as sc
    import logging

    ##################
    # Configure file #
    ##################
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    logging.basicConfig(level=logging.INFO)

    ####################
    # Basic processing #
    ####################
    adata = sc.read(filepath)
    adata.obs['batchcore'] = adata.obs['batch_cov'].astype('category')
    adata.var_names_make_unique()
    logging.info(str('Data structure details: ' + str(adata)))
    # Extract list of genes
    genelist = adata.var_names.tolist()
    # Find mitochondrial genes
    mito_genes_names = [gn for gn in genelist if gn.startswith('MT-')]
    logging.info(str('Mito genes: ' + str(mito_genes_names)))
    # Find indices of mitochondrial genes
    mito_genes = [genelist.index(gn) for gn in mito_genes_names]
    # For each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    # Add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
    logging.info('Filtering cells')
    # Filter cells that have more than 2500 genes or more than 10% of counts coming from mitochondrial genes.
    # These are likely outliers.
    adata = adata[adata.obs['percent_mito'] < 0.10]
    logging.info(str('Data structure details: ' + str(adata)))
    # Filter cells with abnormally low gene counts, high gene counts.
    sc.pp.filter_cells(adata, min_genes=100)
    logging.info(str('Data structure details: ' + str(adata)))
    sc.pp.filter_cells(adata, max_genes=2500)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Saving raw and raw counts')
    adata.uns['raw_counts'] = adata.X
    adata.raw = sc.pp.log1p(adata, copy=True)
    logging.info('Normalizing total counts to 10,000')
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    logging.info('Filtering genes')
    filter_result = sc.pp.filter_genes_dispersion(adata.X, log=True)
    adata = adata[:, filter_result.gene_subset]
    logging.info('Log transforming data')
    sc.pp.log1p(adata)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Regressing out variance within total nUMIs and % mitochondrial UMIs')
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
    logging.info('Batch correcting by regressing out batch variance')
    sc.pp.regress_out(adata, ['batchcore'])
    logging.info('Scaling expression data')
    sc.pp.scale(adata, max_value=10)
    logging.info(str('Data structure details: ' + str(adata)))
    # Unique list of individuals
    people = np.unique(adata.obs['ind_cov'].values.tolist())
    # Allocate space for total PMBCs per individual.
    total_pbmcs = dict.fromkeys(people)
    for p in people:
        total_pbmcs[p] = len(adata.obs_names[adata.obs['ind_cov'] == p])
    adata.uns['total_pbmcs'] = total_pbmcs
    logging.info('Saving processed data')
    adata.write(savepath)


def basicprocessing_noplatelets(filepath, savepath):
    import numpy as np
    import scanpy.api as sc
    from scipy.sparse import csr_matrix
    import logging

    ##################
    # Configure file #
    ##################
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    logging.basicConfig(level=logging.INFO)

    ####################
    # Basic processing #
    ####################
    adata = sc.read(filepath)
    adata.obs['batchcore'] = adata.obs['batch_cov'].astype('category')
    adata.var_names_make_unique()
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Removing platelet contamination and Megakaryocytes.')
    mat = csr_matrix(adata.X)
    mat = mat[:, adata.var_names.isin(['PF4', 'PPBP', 'SDPR', 'GNG11'])].todense()
    mat = np.sum(mat, axis=1)
    logging.info(str(np.shape(mat)))
    # Remove platelet contaminated cells from processing.
    adata = adata[np.ravel(mat <= 0)]
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Removing Erythrocytes.')
    mat = csr_matrix(adata.X)
    mat = mat[:, adata.var_names.isin(['HBB'])].todense()
    adata = adata[np.ravel(mat <= 1)]
    logging.info(str('Data structure details: ' + str(adata)))
    # Extract list of genes
    genelist = adata.var_names.tolist()
    # Find mitochondrial genes
    mito_genes_names = [gn for gn in genelist if gn.startswith('MT-')]
    logging.info(str('Mito genes: ' + str(mito_genes_names)))
    # Find indices of mitochondrial genes
    mito_genes = [genelist.index(gn) for gn in mito_genes_names]
    # For each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    # Add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
    logging.info('Filtering cells')
    # Filter cells that have more than 2500 genes or more than 10% of counts coming from mitochondrial genes.
    # These are likely outliers.
    adata = adata[adata.obs['percent_mito'] < 0.10]
    logging.info(str('Data structure details: ' + str(adata)))
    # Filter cells with abnormally low gene counts, high gene counts.
    sc.pp.filter_cells(adata, min_genes=100)
    logging.info(str('Data structure details: ' + str(adata)))
    sc.pp.filter_cells(adata, max_genes=2500)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Saving raw and raw counts')
    adata.uns['raw_counts'] = adata.X
    adata.raw = sc.pp.log1p(adata, copy=True)
    logging.info('Normalizing total counts to 10,000')
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    logging.info('Filtering genes')
    filter_result = sc.pp.filter_genes_dispersion(adata.X, log=True)
    adata = adata[:, filter_result.gene_subset]
    logging.info('Log transforming data')
    sc.pp.log1p(adata)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Regressing out variance within total nUMIs and % mitochondrial UMIs')
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
    logging.info('Batch correcting by regressing out batch variance')
    sc.pp.regress_out(adata, ['batchcore'])
    logging.info('Scaling expression data')
    sc.pp.scale(adata, max_value=10)
    logging.info(str('Data structure details: ' + str(adata)))
    # Unique list of individuals
    people = np.unique(adata.obs['ind_cov'].values.tolist())
    # Allocate space for total PMBCs per individual.
    total_pbmcs = dict.fromkeys(people)
    for p in people:
        total_pbmcs[p] = len(adata.obs_names[adata.obs['ind_cov'] == p])
    adata.uns['total_pbmcs'] = total_pbmcs
    logging.info('Saving processed data')
    adata.write(savepath)


def cell_identity(filepath, annotation):
    import numpy as np
    import scanpy.api as sc
    import logging

    ##################
    # Configure file #
    ##################
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    logging.basicConfig(level=logging.INFO)
    adata = sc.read(filepath)

    ######################
    # Identify cell type #
    ######################
    # Differentially expressed genes
    Genes = []
    for ct in annotation.keys():
        Genes.append(annotation[ct])
    Genes = [item for sublist in Genes for item in sublist]
    # Determine cell types
    ct_cov = np.asarray(adata.obs['disease_cov'].values.tolist())
    for ii in np.unique(adata.obs['louvain'].tolist()):
        currntsum = -1
        for ct in annotation.keys():
            ctsum = []
            for Gn in Genes:
                if Gn in annotation[ct]:
                    ctsum.append(np.nanmean(adata.X[:, adata.var_names.isin([Gn])][adata.obs['louvain'] == ii]))
                else:
                    ctsum.append(
                        np.nanmean(adata.X[:, adata.var_names.isin([Gn])][adata.obs['louvain'] == ii]) * -1)
            ctsum = np.nanmean(np.asarray(ctsum))
            # Assign cell type
            if ctsum > currntsum:
                ct_cov[adata.obs['louvain'] == ii] = ct
                currntsum = ctsum
            else:
                continue
    adata.obs['ct_cov'] = ct_cov
    adata.obs['ct_cov'] = adata.obs['ct_cov'].astype('category')
    adata.write(filepath)


def basic_analysis(filepath):
    import scanpy.api as sc
    import logging

    ##################
    # Configure file #
    ##################
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    logging.basicConfig(level=logging.INFO)
    adata = sc.read(filepath)

    #######################
    # Louvain and friends #
    #######################
    # Set parameters
    intialization = 1
    n_components = 20
    resolution = 1
    # Run louvain clustering on theoretical future gene expression per cell
    logging.info('Estimating louvain cluster identities for gene expression values.')
    sc.pp.pca(adata, random_state=intialization)
    logging.info('PCA complete.')
    sc.pp.neighbors(adata)
    logging.info('KNN complete.')
    sc.tl.diffmap(adata)
    logging.info('diffmap complete.')
    sc.tl.louvain(adata, random_state=intialization, resolution=resolution)
    logging.info('Louvain complete.')
    sc.tl.paga(adata)
    sc.pl.paga(adata, random_state=intialization, show=False)
    logging.info('paga complete.')
    sc.tl.umap(adata, random_state=intialization, init_pos='paga')
    logging.info('UMAP complete.')
    # First PC for ordering of cells in the umap
    adata.obs['ordering_UMAP'] = sc.pp.pca(adata.obsm['X_umap'], n_comps=1, copy=True)
    logging.info('UMAP ordering complete.')
    sc.tl.rank_genes_groups(adata, groupby='louvain', method='wilcoxon')
    logging.info('Ranked genes complete.')
    adata.write(filepath)
    logging.info('Basic analysis complete.')


def getIFNscore(adata):
    import pandas as pd
    import scanpy.api as sc
    ifn_genes = pd.read_csv('/ye/yelabstore/10x.lupus/batch4.analysis/ifn.lupus.crow.etal.txt', header=None)
    genes = adata.var_names[adata.var_names.isin(ifn_genes[0])].tolist()
    adata_ifn = adata[:, genes]
    sc.tl.pca(adata_ifn)
    adata.obs['IFN'] = adata_ifn.obsm['X_pca'][:, 1]
    return adata

def subpopulation_analysis(cell_type_IDs, filepath, no_norm_filepath, savepath):
    import scanpy.api as sc
    import logging
    import numpy as np

    ####################
    #  Configure file  #
    ####################
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    # Configure scanpy settings
    sc.settings.verbosity = 2
    sc.settings.autoshow = False

    ####################
    # Basic processing #
    ####################
    adata = sc.read(filepath)
    logging.info(str('Processed data structure details: ' + str(adata)))
    # Extract estimated louvain groups from PBMC clustering.
    louvain = adata.obs['louvain'].tolist()
    ct_cov = adata.obs['ct_cov'].tolist()
    cellID = adata.obs_names.tolist()
    adata = sc.read(no_norm_filepath)
    logging.info(str('No normalization data structure details: ' + str(adata)))
    adata = adata[adata.obs_names.isin(cellID)]
    adata.obs['louvain'] = louvain
    adata.obs['ct_cov'] = ct_cov
    logging.info(str('New data structure details: ' + str(adata)))
    bdata = adata[adata.obs['ct_cov'].isin(cell_type_IDs)].copy(savepath)
    bdata.write(savepath)

    ######################
    # Further processing #
    ######################
    adata = sc.read(savepath)
    adata.obs['batchcore'] = adata.obs['batch_cov'].astype('category')
    adata.var_names_make_unique()
    # Extract list of genes
    genelist = adata.var_names.tolist()
    # Find mitochondrial genes
    mito_genes_names = [gn for gn in genelist if gn.startswith('MT-')]
    logging.info(str('Mito genes: ' + str(mito_genes_names)))
    # Find indices of mitochondrial genes
    mito_genes = [genelist.index(gn) for gn in mito_genes_names]
    # For each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    # Add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
    logging.info('Saving raw and raw counts')
    adata.uns['raw_counts'] = adata.X
    adata.raw = sc.pp.log1p(adata, copy=True)
    logging.info('Normalizing total counts to 10,000')
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    logging.info('Filtering genes')
    filter_result = sc.pp.filter_genes_dispersion(adata.X, max_mean=5, min_disp=0.25, min_mean=0.05, log=True)
    adata = adata[:, filter_result.gene_subset]
    logging.info('Log transforming data')
    sc.pp.log1p(adata)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Regressing out variance within total nUMIs and % mitochondrial UMIs')
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
    logging.info('Batch correcting by regressing out batch variance')
    sc.pp.regress_out(adata, ['batchcore'])
    logging.info('Scaling expression data')
    sc.pp.scale(adata, max_value=10)
    logging.info('Basic processing complete.')
    adata.write(savepath)
    # Basic analysis
    basic_analysis(savepath)


def subpopulation_analysis_regressIFN(cell_type_IDs, filepath, no_norm_filepath, IFNpath, savepath):
    import scanpy.api as sc
    import logging
    import numpy as np

    ####################
    #  Configure file  #
    ####################
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    # Configure scanpy settings
    sc.settings.verbosity = 2
    sc.settings.autoshow = False

    ####################
    # Basic processing #
    ####################
    adata = sc.read(filepath)
    logging.info(str('Processed data structure details: ' + str(adata)))
    # Extract estimated louvain groups from PBMC clustering.
    louvain = adata.obs['louvain'].tolist()
    ct_cov = adata.obs['ct_cov'].tolist()
    cellID = adata.obs_names.tolist()
    adata = sc.read(no_norm_filepath)
    logging.info(str('No normalization data structure details: ' + str(adata)))
    adata = adata[adata.obs_names.isin(cellID)]
    adata.obs['louvain'] = louvain
    adata.obs['ct_cov'] = ct_cov
    logging.info(str('New data structure details: ' + str(adata)))
    bdata = adata[adata.obs['ct_cov'].isin(cell_type_IDs)].copy(savepath)
    bdata.write(savepath)

    ######################
    # Further processing #
    ######################
    adata = sc.read(savepath)
    adata.obs['batchcore'] = adata.obs['batch_cov'].astype('category')
    adata.var_names_make_unique()
    # Extract list of genes
    genelist = adata.var_names.tolist()
    # Find mitochondrial genes
    mito_genes_names = [gn for gn in genelist if gn.startswith('MT-')]
    logging.info(str('Mito genes: ' + str(mito_genes_names)))
    # Find indices of mitochondrial genes
    mito_genes = [genelist.index(gn) for gn in mito_genes_names]
    # For each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    # Add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
    logging.info('Saving raw and raw counts')
    adata.uns['raw_counts'] = adata.X
    adata.raw = sc.pp.log1p(adata, copy=True)
    logging.info('Normalizing total counts to 10,000')
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    logging.info('Filtering genes')
    filter_result = sc.pp.filter_genes_dispersion(adata.X, max_mean=5, min_disp=0.25, min_mean=0.05, log=True)
    adata = adata[:, filter_result.gene_subset]
    logging.info('Log transforming data')
    sc.pp.log1p(adata)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Regressing out variance within total nUMIs and % mitochondrial UMIs')
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
    logging.info('Batch correcting by regressing out batch variance')
    sc.pp.regress_out(adata, ['batchcore'])
    logging.info('Regressing out IFN signature')
    getIFNscore(adata)  # Should regress out batch, n_counts and percent_mito first
    sc.pp.regress_out(adata, ['IFN'])
    logging.info('Scaling expression data')
    sc.pp.scale(adata, max_value=10)
    adata.write(savepath)
    logging.info('Basic processing complete.')
    # Basic analysis
    basic_analysis(savepath)

