# Theory questions

## Q1: Samples from different patients may be processed by different technicians and/or at different time points. How could this possibly affect the single cell RNA sequencing data? What type of method is supposed to correct for these potential confounding effects?

Processing samples from different patients by different technicians or at different time points can introduce **batch effects**, which are technical variations unrelated to the biological differences between samples. These variations can affect gene expression measurements, leading to inconsistent or biased data that may obscure true biological signals.

As a result the following could for example happen

> Variability in gene expression due to differences in sample preparation, reagents, or sequencing platforms.
> Introduction of systematic noise that can complicate downstream analysis and interpretation.
> Reduced ability to compare samples accurately across different conditions or patient groups.

Batch effects are typically addressed using batch correction methods such as the following:

> Linear correction approaches: e.g. ComBat which adjusts for batch effects in a linear model framework
> Integration tools: e.g. Harmony, Seurat integration methods, or Scanorama which align datasets while preserving biological variability
> Dimensionality reduction techniques: e.g., principal component analysis (PCA) or canonical correlation analysis (CCA), which can identify and correct batch-related variation.

## Q2: What are the two main categories of methods for RNA deconvolution that exist? List 2 advantages and 2 disadvantages for each of these categories.

The two main categories of RNA deconvolution methods are reference-based and reference-free methods. Reference-based methods rely on pre-existing single-cell or sorted bulk RNA-seq data as reference profiles to estimate the proportions of cell types in mixed samples. On the other hand reference-free methods do not rely on predefined reference profiles and instead use computational algorithms to infer cell type proportions directly from bulk RNA-seq data.

### Reference-based methods

**Advantages**

> High accuracy when reference profiles are representative of the cell types in the mixed sample
> Allows identification of specific cell types with predefined markers

**Disadvantages**

> Dependency on high-quality reference data, which may not be available for all cell types or conditions
> Limited applicability if the sample contains cell types not represented in the reference

### Reference-free methods

**Advantages**

> Suitable for unknown or novel cell types, where reference data is unavailable
> Can adapt to sample-specific variations in gene expression

**Disadvantages:**

- Generally less accurate compared to reference-based methods due to the lack of predefined guidance
- Higher computational complexity and often requires strong assumptions about the number of cell types or their marker genes
