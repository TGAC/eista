# nf-core/eista: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Primary analysis](#primary-analysis)
  - [Vizgen post-processing analysis](#vizgen-post-processing-analysis) - Cell segmentation and quantification
- [Secondary analysis](#secondary-analysis)
  - [QC & cell filtering](#qc--cell-filtering) - Cell filtering and QC on raw data and filtered data
  - [Clustering analysis](#clustering-analysis) - Single-cell clustering analysis
  - [Spatial statistics analysis](#spatial-statistics-analysis) - Single-cell spatial statistics analysis
- [Pipeline reporting](#pipeline-reporting)
  - [Analysis report](#analysis-report) - Single-ell Analysis Report
  - [MultiQC](#multiqc) - Aggregate report describing results and QC for tools registered in nf-core
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution


## Primary analysis

### <u>Vizgen post-processing analysis</u>

**Output directory: `results/vizgen`**
- `<sample_name>/`: Contains output files from Vizgen post-processing too VPT.
  - `counts_unfiltered/`
    - `cell_by_gene.csv`: The raw count of transcripts of each targeted gene in each cell.
    - `cellpose_micron_space.parquet`: The primary output of the segmentation. This table contains the EntityIDs for the cells and the geometries in units of microns.
    - `detected_transcripts.csv`: This file has the information about the detected genes and their locations.
    - `cell_metadata.csv`: The cell metadata file has annotation about the location, size, and shape of each cell that can be used to identify cell neighbors, sort cells into cell types, filter low quality cells, etc.
    - `sum_signals.csv`: The sum signals file has information about the brightness of each mosaic tiff image within each cell.
    - `updated_*.vzg`: An updated vzg file with new entity geometries and a new Cell-by-Entity matrix. The file is used for visualization in the MERSCOPE Desktop Vizualizer software.
- `mtx_conversions/`
  - `<sample_name>/`
    - `sample_name_st_matrix.h5ad`: AnnData object file for this sample.
  - `combined_st_matrix.h5ad`: AnnData object file for combined samples.


## Secondary analysis

### <u>QC & cell filtering</u>

**Output directory: `results/qc_cell_filter`**
- `sample_summary.csv`: overall summary of the single-cell count matrix
- `adata_filtered_normalized.h5ad`: AnnData object file after cell filtering and normalization
- `raw_counts/sample_*/`
  - `scatter_counts_genes_volume.png`: scatter plot shows the relationship between total read counts and the number of genes.
  - `violin*.png`: violin plots display the distribution of cells based on the number of genes, total counts, and cell volumes.
  - `histograms.png`: histogram plots display the distribution of cells based on the number of genes, total counts, and cell volumes. 
- `cell_filtering/`
  - `highly_variable_genes.png`: plot of mean expressions against dispersions of genes for highly variable genes.
  - `sample_summary_filtered.csv`: overall summary of the single-cell count matrix after cell filtering
  - `sample_*/`
    - `spatial_scatter_total_counts_genes.png`: spatial scatter plots for the number of genes, total counts.
    - `umap_total_counts_genes.png`: UMAP plots for the number of genes, total counts.
    - `violin*.png`: violin plots display the distribution of cells based on the number of genes, total counts, and cell volumes.
    - `histograms.png`: histogram plots display the distribution of cells based on the number of genes, total counts, and cell volumes.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.
    

### <u>Clustering analysis</u>

**Output directory: `results/clustering`**
- `adata_clustering.h5ad`: AnnData object file after clustering analysis.
- `sample_*/` or `group_*/`
  - `umap_leiden_res_*.png`: UMAP plots showing clustering results with differnt resoultuion settings.
  - `spatial_scatter_leiden_res_*.png`: spatial scatter plots for the corresponding resolutions of clustering.
- `resolution_*/`
  - `prop_leiden_res_*.png`: plot showing a stacked bar chart that presents the proportions of clusters across samples/groups.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>Spatial statistics analysis</u>

**Output directory: `results/spatialstats`**
- `sample_*/` or `group_*/`
  - `centrality_scores_leiden_res_*.png`: a plot showing 3 centrality scores include: Closeness centrality, Degree centrality, and Clustering coefficient.
  - `spatial_scatter_closeness_high_*.png`: spatial scatter plots highlight clusters/groups with high closeness centrality
  - `spatial_scatter_closeness_low_*.png`: spatial scatter plots highlight clusters/groups with low closeness centrality
  - `spatial_scatter_degree_high_*.png`: spatial scatter plots highlight clusters/groups with high degree centrality
  - `spatial_scatter_degree_low_*.png`: spatial scatter plots highlight clusters/groups with low degree centrality
  - `spatial_scatter_clustering_high_*.png`: spatial scatter plots highlight clusters/groups with high clustering coefficient
  - `spatial_scatter_clustering_low_*.png`: spatial scatter plots highlight clusters/groups with low clustering coefficient
  - `neighbors_enrichment_leiden_res_*.png`: the heatmap shows the results of the neighborhood enrichment analysis on clusters. 
  - `autocorr_moranI.csv`: a csv file shows Moran’s I statistics for autocorrelation analysis of the top 50 and bottom 50 genes.
  - `spatial_scatter_top_*.png`: spatial scatter plots showing top 6 genes among the highest Moran's I scores.
  - `spatial_scatter_bot_*.png`: spatial scatter plots showing top 6 genes among the lowest Moran's I scores.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


## Pipeline reporting

### <u>Analysis report</u>

**Output directory: `results/report`**
- `eisca_report.html`: this is HTML report file showing all major analysis results.


### <u>MultiQC</u>

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### <u>Pipeline information</u>

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
