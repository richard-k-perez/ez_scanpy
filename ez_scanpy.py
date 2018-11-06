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
    labels = clf.fit(adata_singlet.X).predict(voter_thresh=0.50)
    f2, tsne_coords, clusters = doubletdetection.plot.tsne(adata_singlet.X, labels, random_state=1, save='tsne_test.pdf', show=False)
    f3 = doubletdetection.plot.threshold(clf, save='threshold_test.pdf', show=False, p_step=6)
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
        labels = clf.fit(bdata_singlet.X).predict(voter_thresh=0.50)
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
    adata.write_csvs(savepath.split('.')[0])


def remove_doublets_flare(flarepaths, flaregps, bestfilepath, savepath):
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


    ############################
    # Compile demuxlet results #
    ############################
    for ii in range(len(flarepaths)):
        best = pd.read_csv(bestfilepath[ii], delimiter='\t')
        singlets_ind = []
        singlets_bc = []
        for jj in range(len(best['BEST'].values)):
            if 'SNG' in best['BEST'].values[jj]:
                singlets_bc.append(best['BARCODE'].values[jj])
                singlets_ind.append(best['BEST'].values[jj])

        if ii == 0:
            adata = sc.read_10x_h5(flarepaths[ii], 'hg19')
            adata = adata[adata.obs_names.isin(singlets_bc)]
            logging.info(str(adata))
            singlets_bc = pd.Series(singlets_bc)
            singlets_ind = pd.Series(singlets_ind)
            logging.info(str(singlets_bc.isin(adata.obs_names.tolist()).values))
            singlets_ind = singlets_ind[singlets_bc.isin(adata.obs_names.tolist()).values]
            adata.obs['ind_cov'] = singlets_ind.values
            adata.obs['well'] = np.repeat(flaregps[ii], len(adata.obs_names))
            adata.obs['batch_cov'] = np.repeat(flaregps[ii][:-2], len(adata.obs_names))
            logging.info(str('Adata.X: ' + str(np.shape(adata.X))))

            # Independent Doublet Detection
            clf = doubletdetection.BoostClassifier()
            # raw_counts is a cells by genes count matrix
            labels = clf.fit(adata.X).predict(voter_thresh=0.50)
            f2, tsne_coords, clusters = doubletdetection.plot.tsne(adata.X, labels, random_state=1,
                                                                   save='tsne_test.pdf', show=False)
            logging.info(str(labels))
            nan_index = np.isnan(labels)
            labels[nan_index] = 1
            logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
            labels = labels - 1
            logging.info(str(labels))
            labels = np.abs(labels) > 0
            logging.info(str(labels))
            adata = adata[labels]
            logging.info(str('Data structure details: ' + str(adata)))
            f = doubletdetection.plot.convergence(clf, save='convergence_test.pdf', show=False)

        else:
            bdata = sc.read_10x_h5(flarepaths[ii], 'hg19')
            bdata = bdata[bdata.obs_names.isin(singlets_bc)]
            logging.info(str(bdata))
            singlets_bc = pd.Series(singlets_bc)
            singlets_ind = pd.Series(singlets_ind)
            singlets_ind = singlets_ind[singlets_bc.isin(bdata.obs_names.tolist()).values]
            bdata.obs['ind_cov'] = singlets_ind.values
            bdata.obs['well'] = np.repeat(flaregps[ii], len(adata.obs_names))
            bdata.obs['batch_cov'] = np.repeat(flaregps[ii][:-2], len(adata.obs_names))

            # Independent Doublet Detection
            clf = doubletdetection.BoostClassifier()
            # raw_counts is a cells by genes count matrix
            labels = clf.fit(bdata.X).predict(voter_thresh=0.50)
            f2, tsne_coords, clusters = doubletdetection.plot.tsne(bdata.X, labels, random_state=1,
                                                                   save=str(str(ii) + 'b_tsne_test.pdf'), show=False)
            logging.info(str(labels))
            nan_index = np.isnan(labels)
            labels[nan_index] = 1
            logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
            labels = labels - 1
            logging.info(str(labels))
            labels = np.abs(labels) > 0
            logging.info(str(labels))
            bdata = bdata[labels]
            logging.info(str('Data structure details: ' + str(bdata)))
            f = doubletdetection.plot.convergence(clf, save=str(str(ii) + 'b_convergence_test.pdf'), show=False)

            # Concatenate files.
            adata = adata.concatenate(bdata)

    # Add disease state labels to observations
    IDpath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/Flare/all_IDs.csv'
    IDmap = pd.read_csv(IDpath)
    disease_cov = []
    ind_cov = []
    for ii in range(len(adata.obs_names)):
        tmpmap = IDmap[IDmap['Run'].isin([adata.obs['batch_cov'].values[ii][0]])]
        for jj in range(len(tmpmap['Grep for this1'].values)):
            if tmpmap['Grep for this1'].values[jj] in adata.obs['ind_cov'].values[ii]:
                disease_cov.append(tmpmap['Type'].values[jj])
                ind_cov.append(tmpmap['Grep for this1'].values[jj])
    adata.obs['disease_cov'] = disease_cov
    # Replace verbose names with succinct names
    adata.obs['ind_cov'] = ind_cov
    logging.info(str(adata))
    logging.info('Saving compiled demuxlet results')
    adata.write(savepath)
    adata.write_csvs(savepath.split('.')[0])


def remove_doublets_cml(flarepaths, flaregps, bestfilepath, savepath):
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


    ############################
    # Compile demuxlet results #
    ############################
    for ii in range(len(flarepaths)):
        best = pd.read_csv(bestfilepath[ii], delimiter='\t')
        singlets_ind = []
        singlets_bc = []
        for jj in range(len(best['BEST'].values)):
            if 'SNG' in best['BEST'].values[jj]:
                singlets_bc.append(best['BARCODE'].values[jj])
                singlets_ind.append(best['BEST'].values[jj])

        if ii == 0:
            adata = sc.read_10x_h5(flarepaths[ii], 'hg19')
            adata = adata[adata.obs_names.isin(singlets_bc)]
            logging.info(str(adata))
            singlets_bc = pd.Series(singlets_bc)
            singlets_ind = pd.Series(singlets_ind)
            logging.info(str(singlets_bc.isin(adata.obs_names.tolist()).values))
            singlets_ind = singlets_ind[singlets_bc.isin(adata.obs_names.tolist()).values]
            adata.obs['ind_cov'] = singlets_ind.values
            adata.obs['batch_cov'] = np.repeat(flaregps[ii], len(adata.obs_names))
            logging.info(str('Adata.X: ' + str(np.shape(adata.X))))
            '''
            # Independent Doublet Detection
            clf = doubletdetection.BoostClassifier()
            # raw_counts is a cells by genes count matrix
            labels = clf.fit(adata.X).predict(voter_thresh=0.50)
            f2, tsne_coords, clusters = doubletdetection.plot.tsne(adata.X, labels, random_state=1,
                                                                   save='tsne_test.pdf', show=False)
            logging.info(str(labels))
            nan_index = np.isnan(labels)
            labels[nan_index] = 1
            logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
            labels = labels - 1
            logging.info(str(labels))
            labels = np.abs(labels) > 0
            logging.info(str(labels))
            adata = adata[labels]
            logging.info(str('Data structure details: ' + str(adata)))
            f = doubletdetection.plot.convergence(clf, save='convergence_test.pdf', show=False)
            '''
        else:
            bdata = sc.read_10x_h5(flarepaths[ii], 'hg19')
            bdata = bdata[bdata.obs_names.isin(singlets_bc)]
            logging.info(str(bdata))
            singlets_bc = pd.Series(singlets_bc)
            singlets_ind = pd.Series(singlets_ind)
            singlets_ind = singlets_ind[singlets_bc.isin(bdata.obs_names.tolist()).values]
            bdata.obs['ind_cov'] = singlets_ind.values
            bdata.obs['batch_cov'] = np.repeat(flaregps[ii], len(bdata.obs_names))
            '''
            # Independent Doublet Detection
            clf = doubletdetection.BoostClassifier()
            # raw_counts is a cells by genes count matrix
            labels = clf.fit(bdata.X).predict(voter_thresh=0.50)
            f2, tsne_coords, clusters = doubletdetection.plot.tsne(bdata.X, labels, random_state=1,
                                                                   save=str(str(ii) + 'b_tsne_test.pdf'), show=False)
            logging.info(str(labels))
            nan_index = np.isnan(labels)
            labels[nan_index] = 1
            logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
            labels = labels - 1
            logging.info(str(labels))
            labels = np.abs(labels) > 0
            logging.info(str(labels))
            bdata = bdata[labels]
            logging.info(str('Data structure details: ' + str(bdata)))
            f = doubletdetection.plot.convergence(clf, save=str(str(ii) + 'b_convergence_test.pdf'), show=False)
            '''
            # Concatenate files.
            adata = adata.concatenate(bdata)

    logging.info(str(adata))
    logging.info('Saving compiled demuxlet results')
    adata.write(savepath)
    adata.write_csvs(savepath.split('.')[0])


def remove_doublets_immvar(paths, gps, bestfilepath, savepath):
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


    ############################
    # Compile demuxlet results #
    ############################
    for ii in range(len(paths)):
        best = pd.read_csv(bestfilepath[ii], delimiter='\t')
        singlets_ind = []
        singlets_bc = []
        for jj in range(len(best['BEST'].values)):
            if 'SNG' in best['BEST'].values[jj]:
                singlets_bc.append(best['BARCODE'].values[jj])
                singlets_ind.append(best['BEST'].values[jj].split('-')[1])

        if ii == 0:
            adata = sc.read_10x_h5(paths[ii], 'hg19')
            adata = adata[adata.obs_names.isin(singlets_bc)]
            logging.info(str(adata))
            singlets_bc = pd.Series(singlets_bc)
            singlets_ind = pd.Series(singlets_ind)
            logging.info(str(singlets_bc.isin(adata.obs_names.tolist()).values))
            singlets_ind = singlets_ind[singlets_bc.isin(adata.obs_names.tolist()).values]
            adata.obs['ind_cov'] = singlets_ind.values
            adata.obs['well'] = np.repeat(gps[ii], len(adata.obs_names))
            adata.obs['batch_cov'] = np.repeat(gps[ii][:-2], len(adata.obs_names))
            logging.info(str('Adata.X: ' + str(np.shape(adata.X))))

            # Independent Doublet Detection
            clf = doubletdetection.BoostClassifier()
            # raw_counts is a cells by genes count matrix
            labels = clf.fit(adata.X).predict(voter_thresh=0.50)
            f2, tsne_coords, clusters = doubletdetection.plot.tsne(adata.X, labels, random_state=1,
                                                                   save='immvar_tsne_test.pdf', show=False)
            logging.info(str(labels))
            nan_index = np.isnan(labels)
            labels[nan_index] = 1
            logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
            labels = labels - 1
            logging.info(str(labels))
            labels = np.abs(labels) > 0
            logging.info(str(labels))
            adata = adata[labels]
            logging.info(str('Data structure details: ' + str(adata)))
            f = doubletdetection.plot.convergence(clf, save='immvar_convergence_test.pdf', show=False)
        else:
            bdata = sc.read_10x_h5(paths[ii], 'hg19')
            bdata = bdata[bdata.obs_names.isin(singlets_bc)]
            logging.info(str(bdata))
            singlets_bc = pd.Series(singlets_bc)
            singlets_ind = pd.Series(singlets_ind)
            singlets_ind = singlets_ind[singlets_bc.isin(bdata.obs_names.tolist()).values]
            bdata.obs['ind_cov'] = singlets_ind.values
            bdata.obs['well'] = np.repeat(gps[ii], len(bdata.obs_names))
            bdata.obs['batch_cov'] = np.repeat(gps[ii][:-2], len(bdata.obs_names))

            # Independent Doublet Detection
            clf = doubletdetection.BoostClassifier()
            # raw_counts is a cells by genes count matrix
            labels = clf.fit(bdata.X).predict(voter_thresh=0.50)
            logging.info(str(labels))
            nan_index = np.isnan(labels)
            labels[nan_index] = 1
            logging.info(str('Number of Doublets: ' + str(np.sum(labels))))
            labels = labels - 1
            logging.info(str(labels))
            labels = np.abs(labels) > 0
            logging.info(str(labels))
            bdata = bdata[labels]
            logging.info(str('Data structure details: ' + str(bdata)))
            # Concatenate files.
            adata = adata.concatenate(bdata)

    logging.info(str(adata))
    logging.info('Saving compiled demuxlet results')
    adata.write(savepath)


def combine_studies(crossXpath, flarepath, savepath):
    import scanpy.api as sc
    import logging
    adatacrossX = sc.read(crossXpath)
    adatacrossX.obs['study_cov'] = adatacrossX.obs['batch_cov'].astype('object')
    adatacrossX.obs['study_cov'] = 'Cross-Sectional'
    logging.info(str('CrossX structure details: ' + str(adatacrossX)))
    adataflare = sc.read(flarepath)
    adataflare.obs['study_cov'] = adataflare.obs['batch_cov'].astype('object')
    adataflare.obs['study_cov'] = 'Flare'
    logging.info(str('Flare structure details: ' + str(adataflare)))
    adata = adatacrossX.concatenate(adataflare)
    adata.obs['disease_cov_simple'] = adata.obs['disease_cov'].astype('object')
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['healthy', 'ImmVar', 'UCSF Control', 'UCSF Healthy'])] = 'healthy'
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['sle', 'Treated', 'Treated Pair'])] = 'sle_treated'
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['Untreated', 'Untreated Pair'])] = 'sle_flaring'
    logging.info(str('New combined structure details: ' + str(adata)))
    adata.write(savepath)


def clean_anndata(filepath, clean, keep):
    import numpy as np
    import scanpy.api as sc
    """Deletes attributes not needed.
    adapted from scvelo.pp.cleanup
    Copyright (c) 2018, Theis Lab
    All rights reserved.
    Arguments
    ---------
    data: :class:`~anndata.AnnData`
        Annotated data matrix.
    clean: `str` or list of `str` (default: `layers`)
        Which attributes to consider for freeing memory.
    keep: `str` or list of `str` (default: `['spliced', unspliced']`)
        Which attributes to keep.
    Returns
    -------
    Returns or updates `adata` with selection of attributes kept.
    """
    adata = sc.read(filepath)
    logging.info(str('Structure details: ' + str(adata)))
    logging.info(np.unique(adata.obs['well'].tolist()))
    logging.info(np.unique(adata.obs['batch_cov'].tolist()))
    logging.info(np.unique(adata.obs['ind_cov'].tolist()))
    logging.info(np.unique(adata.obs['pop_cov'].tolist()))

    if any(['obs' in clean, 'all' in clean]):
        for key in list(adata.obs.keys()):
            if key not in keep: del adata.obs[key]

    if any(['var' in clean, 'all' in clean]):
        for key in list(adata.var.keys()):
            if key not in keep: del adata.var[key]

    if any(['uns' in clean, 'all' in clean]):
        for key in list(adata.uns.keys()):
            if key not in keep: del adata.uns[key]

    if any(['layers' in clean, 'all' in clean]):
        for key in list(adata.layers.keys()):  # remove layers that are not needed
            if key not in keep: del adata.layers[key]
    logging.info(str('Structure details: ' + str(adata)))
    adata.write(filepath)


def combine_CLUES_Immvar(crossXpath, immvarpath, savepath):
    import scanpy.api as sc
    import logging
    import numpy as np
    adatacrossX = sc.read(crossXpath)
    adatacrossX.obs['study_cov'] = adatacrossX.obs['batch_cov'].astype('object')
    adatacrossX.obs['study_cov'] = 'CLUES'
    logging.info(str('CrossX structure details: ' + str(adatacrossX)))
    adataimmvar = sc.read(immvarpath)
    adataimmvar.obs['study_cov'] = adataimmvar.obs['batch_cov'].astype('object')
    adataimmvar.obs['disease_cov'] = adataimmvar.obs['batch_cov'].astype('object')
    adataimmvar.obs['pop_cov'] = 'WHITE'
    adataimmvar.obs['study_cov'] = 'ImmVar'
    adataimmvar.obs['disease_cov'] = 'healthy'
    logging.info(str('ImmVar structure details: ' + str(adataimmvar)))
    adata = adatacrossX.concatenate(adataimmvar)
    adata.obs['disease_cov_simple'] = adata.obs['disease_cov'].astype('object')
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['healthy', 'ImmVar', 'UCSF Control', 'UCSF Healthy'])] = 'healthy'
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['sle', 'Treated', 'Treated Pair'])] = 'sle_treated'
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['Untreated', 'Untreated Pair'])] = 'sle_flaring'
    # Reduce floating point accuracy for space saving.
    adata.X = adata.X.astype(np.float16, copy=False)
    logging.info(str('New combined structure details: ' + str(adata)))
    logging.info(np.unique(adata.obs['well'].tolist()))
    logging.info(np.unique(adata.obs['batch_cov'].tolist()))
    logging.info(np.unique(adata.obs['ind_cov'].tolist()))
    logging.info(np.unique(adata.obs['pop_cov'].tolist()))
    logging.info(str('Remove variable annotations'))
    for key in list(adata.var.keys()):
        del adata.var[key]
    logging.info(str('Remove redundant batch observation'))
    del adata.obs['batch']
    logging.info(str('New combined structure details: ' + str(adata)))
    adata.write(savepath)


def combine_studies_cml(crossXpath, flarepath, cmlpath, savepath):
    import scanpy.api as sc
    import logging
    adatacrossX = sc.read(crossXpath)
    adatacrossX.obs['study_cov'] = adatacrossX.obs['batch_cov'].astype('object')
    adatacrossX.obs['study_cov'] = 'Cross-Sectional SLE'
    logging.info(str('CrossX structure details: ' + str(adatacrossX)))
    adataflare = sc.read(flarepath)
    adataflare.obs['study_cov'] = adataflare.obs['batch_cov'].astype('object')
    adataflare.obs['study_cov'] = 'Flare SLE'
    logging.info(str('Flare structure details: ' + str(adataflare)))
    adata = adatacrossX.concatenate(adataflare)
    adata.obs['disease_cov_simple'] = adata.obs['disease_cov'].astype('object')
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['healthy', 'ImmVar', 'UCSF Control', 'UCSF Healthy'])] = 'healthy'
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['sle', 'Treated', 'Treated Pair'])] = 'sle_treated'
    adata.obs['disease_cov_simple'][adata.obs['disease_cov_simple'].isin(['Untreated', 'Untreated Pair'])] = 'sle_flaring'
    logging.info(str('New combined structure details: ' + str(adata)))

    adatacml = sc.read(cmlpath)
    adatacml.obs['study_cov'] = adatacml.obs['batch_cov'].astype('object')
    adatacml.obs['study_cov'] = 'CML'
    logging.info(str('CML structure details: ' + str(adatacml)))
    adata = adata.concatenate(adatacml)
    logging.info(str('New combined structure details: ' + str(adata)))
    adata.write(savepath)

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
    adata.obs['well'] = adata.obs['well'].astype('category')
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
    logging.info('Saving raw')
    #adata.uns['raw_counts'] = adata.X
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
    sc.pp.regress_out(adata, ['well'])
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
    adata.obs['well'] = adata.obs['well'].astype('category')
    adata.var_names_make_unique()
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info('Removing platelet contamination and Megakaryocytes.')
    mat = csr_matrix(adata.X)
    mat = mat[:, adata.var_names.isin(['PF4', 'PPBP', 'SDPR', 'GNG11'])].todense()
    # Regress out platelet signature by allocating genes as observations
    #logging.info(str('PPBP: ' + str(np.asarray(mat[:, 0], dtype=int))))
    #adata.obs['PPBP'] = np.asarray(mat[:, 0], dtype=int)
    #adata.obs['PF4'] = np.asarray(mat[:, 1], dtype=int)
    #adata.obs['SDPR'] = np.asarray(mat[:, 2], dtype=int)
    #adata.obs['GNG11'] = np.asarray(mat[:, 3], dtype=int)
    # Sum across genes
    mat = np.sum(mat, axis=1)
    logging.info(str(np.shape(mat)))
    # Remove platelet contaminated cells from processing.
    adata = adata[np.ravel(mat <= 2) | np.ravel(mat >= 25)]
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
    logging.info('Saving raw')
    #adata.uns['raw_counts'] = adata.X
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
    sc.pp.regress_out(adata, ['well'])
    # Regress out platelet signature
    #logging.info('Regressing out platelet signature')
    #sc.pp.regress_out(adata, ['PF4', 'PPBP', 'SDPR', 'GNG11'])
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


def basic_processing_rnavelocity(filepath, savepath):
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
    adata.obs['batchcore'] = adata.obs['batchcore'].astype('category')
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


def add_covariates_2_rnavelocity(loompath, refh5adpath, savepath):
    import numpy as np
    import scanpy.api as sc
    import logging
    import scvelo as scv
    import pandas as pd

    ##################
    # Configure file #
    ##################
    sc.settings.verbosity = 2
    sc.settings.autoshow = False
    logging.basicConfig(level=logging.INFO)

    adata = scv.read(loompath)
    logging.info(str('loom data structure details: ' + str(adata)))
    rnavelobarcodes = adata.obs_names.tolist()
    batch = []
    BARCODES = []
    for barcode in rnavelobarcodes:
        tmp = barcode.split(':')[0]
        if tmp[:5] != 'YE110':
            tmp = tmp.replace('-', '_', 1)
        batch.append(tmp)
        newbarcode = barcode.split(':')[1]
        newbarcode = newbarcode[:-1]
        BARCODES.append(newbarcode)
    adata.obs['barcodes'] = BARCODES
    adata.obs['batchcore'] = batch
    logging.info(str('Updated data structure details: ' + str(adata)))

    refadata = sc.read(refh5adpath, cache=True)
    logging.info(str('Reference data structure details: ' + str(refadata)))

    DZ = []
    IND = []
    WELL = []
    CT = []
    LOUVAIN = []
    for ii in range(len(BARCODES)):
        #logging.info(str('ii: ' + str(ii)))
        #logging.info(str('batch: ' + str(batch[ii])))
        refbarcodes = refadata.obs['well'][refadata.obs['well']==batch[ii]].index.tolist()

        for jj in range(len(refbarcodes)):
            refbarcodes[jj] = refbarcodes[jj].split('-')[0]

        dz = refadata.obs['disease_cov'][refadata.obs['well']==batch[ii]].values.tolist()
        ind = refadata.obs['ind_cov'][refadata.obs['well']==batch[ii]].values.tolist()
        well = refadata.obs['well'][refadata.obs['well']==batch[ii]].values.tolist()
        ct = refadata.obs['ct_cov'][refadata.obs['well']==batch[ii]].values.tolist()
        louvain = refadata.obs['louvain'][refadata.obs['well']==batch[ii]].values.tolist()

        #logging.info(str('barcode: ' + str(BARCODES[ii])))
        #logging.info(str('barcode: ' + str(refbarcodes.index(BARCODES[ii]))))
        #logging.info(str('barcode: ' + str(dz[refbarcodes.index(BARCODES[ii])])))

        DZ.append(dz[refbarcodes.index(BARCODES[ii])])
        IND.append(ind[refbarcodes.index(BARCODES[ii])])
        WELL.append(well[refbarcodes.index(BARCODES[ii])])
        CT.append(ct[refbarcodes.index(BARCODES[ii])])
        LOUVAIN.append(louvain[refbarcodes.index(BARCODES[ii])])

    adata.obs['disease_cov'] = DZ
    adata.obs['ind_cov'] = IND
    adata.obs['ct_cov'] = CT
    adata.obs['well'] = WELL
    adata.obs['louvain'] = LOUVAIN
    logging.info(str('Updated data structure details: ' + str(adata)))
    logging.info(str('Saving...'))
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
    #sc.tl.diffmap(adata)
    #logging.info('diffmap complete.')
    sc.tl.louvain(adata, random_state=intialization, resolution=resolution)
    logging.info('Louvain complete.')
    sc.tl.paga(adata)
    sc.pl.paga(adata, random_state=intialization, show=False)
    logging.info('paga complete.')
    sc.tl.umap(adata, random_state=intialization, init_pos='paga')
    logging.info('UMAP complete.')
    #sc.tl.tsne(adata, n_pcs=10)
    #logging.info('TSNE complete.')
    # First PC for ordering of cells in the umap
    #adata.obs['ordering_UMAP'] = sc.pp.pca(adata.obsm['X_umap'], n_comps=1, copy=True)
    #logging.info('UMAP ordering complete.')
    #sc.tl.rank_genes_groups(adata, groupby='louvain', method='wilcoxon')
    #logging.info('Processing PHATE')
    #sc.tl.phate(adata, n_pca=10)
    #logging.info('Ranked genes complete.')
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
    logging.info(str('ct_cov: ' + str(bdata.obs['ct_cov'])))
    logging.info(str('bdata structure details: ' + str(bdata)))
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

def subpopulation_analysis_regressIFN(cell_type_IDs, filepath, no_norm_filepath, savepath):
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

def subpopulation_analysis_covariate(cell_type_IDs, covariate, covariate_grp, filepath, no_norm_filepath, savepath):
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
    bdata = bdata[bdata.obs[covariate] == covariate_grp]
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


def subpopulation_analysis_removeIFNGenes(cell_type_IDs, filepath, no_norm_filepath, savepath):
    import scanpy.api as sc
    import logging
    import numpy as np
    import pandas as pd

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
    logging.info(str('Data structure before removing genes: ' + str(adata)))
    data = pd.read_csv('/ye/yelabstore/10x.lupus/batch2.analysis/v2.diff.exp/cd14.sig.txt', sep='\t')
    data = data['featureData.symbol'].values.tolist()
    logging.info(str('Removing the following genes: ' + str(data)))
    Genes = pd.Series(adata.var_names.tolist())
    adata = adata[:, Genes.isin(data) == False]
    logging.info(str('Data structure after removing genes: ' + str(adata)))
    adata.write(savepath)
    # Basic analysis
    basic_analysis(savepath)


def gene_heatmap():
    import numpy as np
    import scanpy.api as sc
    import logging
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    ####################
    #  Configure file  #
    ####################
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    # Configure scanpy settings
    sc.settings.verbosity = 2
    sc.settings.autoshow = False

    allcellspath = '/ye/yelabstore2/10x.lupus/process.scRNAseq/Richard/PlottingScripts/ez_scanpy/Mono_flare.h5ad'

    adata = sc.read(allcellspath)
    logging.info(str('Data structure details: ' + str(adata)))
    logging.info(str('ct_cov: ' + str(adata.obs['ct_cov'])))


    # Remove doublets and platelets contaminated cells
    adata = adata[adata.obs['louvain'].isin(['0', '1', '2', '3'])]
    unique_louvain = np.unique(adata.obs['louvain'].values)

    Diff_genes = []
    n_genes = 50
    for ii in range(len(unique_louvain)):
        if ii == 0:
            Diff_genes = set(
                adata.var_names[np.flipud(np.argsort(np.mean(adata.X[adata.obs['louvain'] == str(ii)], axis=0)))][
                :n_genes])
        else:
            Diff_genes.update(set(
                adata.var_names[np.flipud(np.argsort(np.mean(adata.X[adata.obs['louvain'] == str(ii)], axis=0)))][
                :n_genes]))

    # Diff_genes = adata.var_names.tolist()
    Diff_genes = list(Diff_genes)
    logging.info(str(Diff_genes))
    # Assign each cell to a louvain group
    # palette = sns.color_palette("husl", len(np.unique(adata.obs['louvain'].tolist())))
    colors = [adata.uns['louvain_colors'][int(a)] for a in adata.obs['louvain'].tolist()]
    louvain = pd.Series(colors, index=adata.obs['ordering_UMAP'].tolist())

    # Generate Gene expression plot
    for Gn in Diff_genes:
        logging.info(str(Gn))
        if Gn == Diff_genes[0]:
            mat = np.reshape(adata.X[:, adata.var_names.isin([Gn])], (-1, 1))
        else:
            mat = np.concatenate((mat, np.reshape(adata.X[:, adata.var_names.isin([Gn])], (-1, 1))), axis=1)
    logging.info(str(np.shape(mat)))
    clustermat = pd.DataFrame(data=mat, index=adata.obs['ordering_UMAP'].tolist(), columns=Diff_genes)
    clustermat = clustermat.sort_index(axis=0)
    clustermat = clustermat.transpose()
    clustermat = clustermat.rolling(window=int(np.round(len(adata.obs_names.tolist()) * 0.05)), min_periods=1,
                                    axis=1).mean()
    min = clustermat.min(axis=1)
    max = clustermat.max(axis=1)
    clustermat = clustermat.subtract(min, axis=0)
    clustermat = clustermat.div(max - min, axis=0)
    sns.set(font_scale=0.25)
    sns.clustermap(clustermat, row_cluster=True, col_cluster=False, yticklabels=True, col_colors=louvain)
    plt.savefig('clustermap.png', dpi=1000)
