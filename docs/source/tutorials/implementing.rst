Implementing SAKURA
========================================

In this tutorial, we go through the basic steps for implementing SAKURA, using a dataset of ~5k
Peripheral Blood Mononuclear Cells (PBMC) freely available from 10XGenomics_.
The processed data can be found here_.

.. _10XGenomics: https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x

.. _here: https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x

Data preprocessing
------------------------------------------
We first include a data preprocessing pipeline using R and the Seurat package, which starts by reading in the data.
The Read10X_h5() function in Seurat reads in the HDF5 file that contains single-cell RNA sequencing (scRNA-seq) data.
We use the count matrix to initialize a Seurat object with quality control parameters:

.. code-block:: r

    library(Seurat)
    library(Matrix)
    library(tidyverse)
    library(ggplot2)

    data <- Read10X_h5(filename = "./pbmc5k/raw/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5")
    # Retains genes detected in at least 3 cells, cells with at least 200 detected genes
    seurat_object <- CreateSeuratObject(counts=data, project="pbmc5k", min.cells=3, min.features=200)

Cell Type Filtering
------------------------------------------
Next, we read in the cell type information and filter a few cell types in the PBMC_5k dataset
for simplicity, and then integrate corresponding cell type information into Seurat object metadata.

.. code-block:: r

    cell_type <- read.csv("./pbmc5k/raw/cell_types.csv",row.names = 1,header=TRUE)
    cur_cell <- rownames(seurat_object@meta.data)
    cur_celltype <- cell_type[cur_cell, ]
    seurat_object@meta.data <- cbind(seurat_object@meta.data, cur_celltype)
    cell_types_to_remove <- c("glial cell", "hematopoietic cell", "hematopoietic precursor cell", "mast cell", "stem cell")
    cells_to_keep <- !seurat_object@meta.data$coarse_cell_type %in% cell_types_to_remove
    seurat_filtered <- seurat_object[, cells_to_keep]

.. note::
    From 10x, cell type annotation (https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-cell-annotation-algorithm) is currently in beta and relies on the Chan Zuckerberg CELL by GENE reference.
    As this reference is community-driven it may not cover all tissue types, which could impact your results.

Normalization and Feature Selection
-----------------------------------
We then perform minimal pre-processing of the data, including log-transformation of the raw UMI counts using NormalizeData(),
selecting the 10,000 genes with highest variance using FindVariableFeatures(), and scaling the data based on the high variance genes using ScaleData().

.. code-block:: r

    seurat_filtered <- seurat_filtered %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 10000)
    seurat_filtered <- seurat_filtered %>% ScaleData(features = VariableFeatures(seurat_filtered))
    seurat.hv10k <- seurat_filtered[VariableFeatures(seurat_filtered),]

Data Export
-----------
After subsetting 10k highly variable genes, we export key data components as input data of SAKURA:

  - ``genenames_hv10k.csv``: List of 10k highly variable gene names
  - ``cell_names.csv``: Filtered cell barcodes
  - ``pheno_df.csv``: Cell metadata including cell type annotations
  - ``lognorm_hv10k.mtx``: Normalized expression matrix in Matrix Market format

.. code-block:: r

    write.csv(rownames(seurat.hv10k),"./tests/pbmc5k/processed/genenames_hv10k.csv")
    write.csv(colnames(seurat.hv10k),"./tests/pbmc5k/processed/cell_names.csv")
    write.csv(seurat.hv10k@meta.data,"./tests/pbmc5k/processed/pheno_df.csv")
    writeMM(seurat.hv10k@assays$RNA$data,"./tests/pbmc5k/processed/lognorm_hv10k.mtx")

Users can prepare the pre-processed example dataset as input for SAKURA or use their own datasets.
Training and testing with demo configuration on this dataset will cost less than 1 hour on the Intel Xeon CPU E5-2630 v2 @ 2.60GHz,
which can be substantially improved with GPUs.

Running SAKURA with the example dataset
-----------------------------------------
SAKURA uses a comprehensive JSON configuration file to control all aspects of the training process.
Below we break down each section of the example configuration files:

    - ``./pbmc5k/config.json``
    - ``./pbmc5k/signature_config.json``

Basic Configuration
,,,,,,,,,,,,,,,,,,,,,,

.. code-block:: json

  {
    "remarks": "",
    "log_path": "./test/pbmc5k/log/",
    "reproducible": "True",
    "rnd_seed": 3407,

**Parameters:**

    - ``remarks``: User comments or notes about this configuration
    - ``log_path``: Directory path for storing training logs and outputs
    - ``reproducible``: When set to "True", ensures reproducible results using the specified random seed
    - ``rnd_seed``: Random seed (3407) for reproducible random number generation

Dataset Configuration
,,,,,,,,,,,,,,,,,,,,,,
Related API: :class:`sakura.dataset`

.. code-block:: json

  "dataset": {
    "type": "rna_count_sparse",
    "gene_expr_MM_path": "./tests/pbmc5k/processed/lognorm_hv10k.mtx",
    "gene_name_csv_path": "./tests/pbmc5k/processed/genenames_hv10k.csv",
    "cell_name_csv_path": "./tests/pbmc5k/processed/cell_names.csv",
    "pheno_csv_path": "./tests/pbmc5k/processed/pheno_df.csv",
    "pheno_df_dtype": {
      "Batch": "string",
      "Main_cluster_name": "string"
    },
    "pheno_df_na_filter": "False",
    "expr_mat_pre_slice": "False",
    "signature_config_path": "./tests/pbmc5k/processed/signature_config.json",
    "selected_signature": ["cd8"]
  },

**Parameters:**

    - ``type``: Data format type ("rna_count_sparse" for PBMC5k sparse RNA count matrices)
    - ``gene_expr_MM_path``: Path to normalized gene expression matrix in Matrix Market format
    - ``gene_name_csv_path``: Path to CSV file containing gene names
    - ``cell_name_csv_path``: Path to CSV file containing cell barcodes/identifiers
    - ``pheno_csv_path``: Path to CSV file containing cell phenotype metadata
    - ``pheno_df_dtype``: Data types for phenotype columns (Batch and cell type annotations as strings)
    - ``pheno_df_na_filter``: Whether to filter out NA values in phenotype data
    - ``expr_mat_pre_slice``: Whether the expression matrix is pre-sliced
    - ``signature_config_path``: Path to signature configuration file defining signature learning task
    - ``selected_signature``: Name list of signatures to use (here, we use 'CD8A', 'CD8B' and name them 'cd8')

.. note::
    Similarly, users can include optional phenotype learning task configuration JSON file with
    ``pheno_meta_path`` and ``selected_pheno``. See `Signature Configuration <./pbmc5k/signature_config.json>` for more details.

Hardware and Data Splitting
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
Related API: :class:`sakura.utils.data_splitter.DataSplitter`

.. code-block:: json

  "device": "cpu",
  "overall_train_test_split": {
    "type": "auto",
    "train_dec": 5,
    "seed": 3407
  },

**Parameters:**

    - ``device``: Computation device ("cpu" for CPU-only, "cuda" for GPU acceleration)
    - ``overall_train_test_split``: Configuration for splitting data into training and testing sets
    - ``type``: "auto" for automatic splitting
    - ``train_dec``: Training set ratio denominator (5 = 80% training, 20% testing)
    - ``seed``: Random seed for reproducible data splitting

Model Architecture
,,,,,,,,,,,,,,,,,,,,,,
Related API: structure settings - :class:`sakura.models.extractor.Extractor` and loss/regularization settings - :class:`sakura.model_controllers.extractor_controller.ExtractorController`

.. code-block:: json

  "main_latent": {
    "encoder_neurons": 200,
    "decoder_neurons": 200,
    "latent_dim": 50,
    "loss": {
      "L2": {
        "type": "MSE",
        "init_weight": 1.0,
        "progressive_const": 1.0,
        "progressive_start_epoch": 100
      },
      "L1": {
        "type": "L1",
        "init_weight": 1.0,
        "progressive_const": 1.0,
        "progressive_start_epoch": 100
      }
    },
    "regularization": {
      "uniform_shape_unsupervised": {
        "type": "SW2_uniform",
        "init_weight": 0.0001,
        "progressive_const": 1.01,
        "progressive_start_epoch": 1,
        "max_weight": 1.0,
        "SW2_num_projections": 50,
        "uniform_low": -10,
        "uniform_high": 10
      }
    }
  },

**Parameters:**

    - ``encoder_neurons``: Number of neurons in encoder hidden layers (200)
    - ``decoder_neurons``: Number of neurons in decoder hidden layers (200)
    - ``latent_dim``: Dimensionality of the latent space (50)
    - ``loss``: Reconstruction loss configuration
    - ``L2``: Mean Squared Error loss with progressive weighting
    - ``L1``: L1 loss for robust reconstruction
    - ``regularization``: Regularization terms to shape the latent space
    - ``uniform_shape_unsupervised``: Sliced Wasserstein distance regularization to enforce uniform distribution


Optimizer Configuration
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
Related API: structure settings - :func:`sakura.model_controllers.extractor_controller.ExtractorController.setup_optimize()`

.. code-block:: json

  "optimizer": {
    "type": "RMSProp",
    "RMSProp_lr": 0.001,
    "RMSProp_alpha": 0.9
  },

**Parameters:**

    - ``type``: Optimization algorithm ("RMSProp")
    - ``RMSProp_lr``: Learning rate (0.001)
    - ``RMSProp_alpha``: Smoothing constant for RMSProp (0.9)

Training Pipeline (Story)
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
Related API: :func:`sakura.sakuraAE.sakuraAE.train_hybrid()`
The ``story`` section defines the complete training workflow:

.. code-block:: json

  "story": [
    {
      "action": "train_hybrid",
      "ticks": 5000,
      "hybrid_mode": "interleave",

**Training Strategy:**

    - ``action``: Training mode ("train_hybrid" for overleave reconstruction and signature training)
    - ``ticks``: Total number of training epochs (5000)
    - ``hybrid_mode``: "interleave" alternates between different training objectives

Split Configurations
~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

  "split_configs": {
    "main_lat_reconstruct": {
      "use_split": "overall_train",
      "batch_size": 100,
      "train_main_latent": "True",
      "train_pheno": "False",
      "train_signature": "False"
    },
    "cd8_focused": {
      "use_split": "overall_train",
      "batch_size": 100,
      "train_main_latent": "False",
      "train_pheno": "False",
      "train_signature": "True",
      "selected_signature": {
        "cd4cd8": {
          "loss": "*",
          "regularization": "*"
        }
      }
    }
  },

**Training Splits:**
- ``main_lat_reconstruct``: Main autoencoder reconstruction training
  - ``batch_size``: 100 cells per batch
  - Focuses on learning latent representations
- ``cd8_focused``: Signature-guided training
  - Uses CD4/CD8 signature to guide latent space organization
  - Applies signature-specific losses and regularizations

Monitoring and Checkpointing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

  "prog_loss_weight_mode": "epoch_end",
  "perform_test": "True",
  "test_segment": 500,
  "make_logs": "True",
  "log_prefix": "hybrid",
  "perform_checkpoint": "True",
  "checkpoint_segment": 500,


.. _signature_config.json:

Signature Configuration
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
testing


