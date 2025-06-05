import scanpy as sc
from sklearn.metrics import adjusted_rand_score, v_measure_score
from sklearn.model_selection import train_test_split
import scanpy.external as sce

seed = 42

TRAIN_ADATA_PATH = "./data/cluster/train/train_adata.h5ad"
TEST_ADATA_PATH = "./data/cluster/test/test_adata.h5ad"


def load_data():
    train = sc.read_h5ad(TRAIN_ADATA_PATH)
    test = sc.read_h5ad(TEST_ADATA_PATH)

    return train, test


def train_cluster(train_adata):
    # List of unique patients
    patients = train_adata.obs["Patient"].unique()

    # Split patients into training and validation sets
    train_patients, val_patients = train_test_split(
        patients, test_size=0.3, random_state=seed
    )

    # Subset the data
    train = train_adata[train_adata.obs["Patient"].isin(train_patients)].copy()
    validation = train_adata[train_adata.obs["Patient"].isin(val_patients)].copy()

    combined_train_adata = train.concatenate(
        validation,
        batch_key="validation_split",
        batch_categories=["train", "validation"],
        index_unique=None,
    )

    combined_train_adata = combined_train_adata[
        combined_train_adata.obs.pct_counts_mt < 5, :
    ]

    sc.pp.filter_cells(combined_train_adata, max_counts=25000)
    sc.pp.filter_cells(combined_train_adata, max_genes=6000)

    combined_train_adata = base_clustering(combined_train_adata)

    return combined_train_adata


def base_clustering(combined_adata):
    # Normalization and log transformation
    sc.pp.normalize_total(combined_adata, target_sum=1e4)
    sc.pp.log1p(combined_adata)

    # HVGs
    sc.pp.highly_variable_genes(
        combined_adata, flavor="seurat_v3", n_top_genes=2000, batch_key="Sample"
    )
    combined_adata = combined_adata[:, combined_adata.var.highly_variable]

    # Scaling and PCA
    sc.pp.scale(combined_adata, max_value=10)
    sc.tl.pca(combined_adata, svd_solver="arpack", n_comps=200)

    # Harmony integration
    sce.pp.harmony_integrate(combined_adata, key="Sample")

    # Clustering
    sc.pp.neighbors(combined_adata, n_neighbors=30, use_rep="X_pca_harmony")
    sc.tl.leiden(combined_adata, resolution=0.2)

    return combined_adata


def validation_score(combined_train_adata):
    """evaluate clustering performance on validation set"""

    # Split combined data back
    val_indices = combined_train_adata.obs["validation_split"] == "validation"

    validation_subset = combined_train_adata[val_indices].copy()

    # eval cluster performance

    # True labels and predicted clusters for validation data
    true_labels = validation_subset.obs["highLevelType"]
    pred_labels = validation_subset.obs["leiden"]

    # Compute ARI
    ari = adjusted_rand_score(true_labels, pred_labels)
    print(f"Validation Adjusted Rand Index (ARI): {ari:.4f}")

    # Compute V-measure
    v_measure = v_measure_score(true_labels, pred_labels)
    print(f"Validation V-measure Score: {v_measure:.4f}")

    # Combined metric
    clustering_score = 0.5 * v_measure + 0.5 * ari
    print(f"Validation Clustering Score: {clustering_score:.4f}")


def test_cluster(test_adata, train_adata):
    # Filter training data before combining
    train_adata_filtered = train_adata.copy()

    # Quality control for training data only
    train_adata_filtered = train_adata_filtered[
        train_adata_filtered.obs.pct_counts_mt < 5, :
    ]
    sc.pp.filter_cells(train_adata_filtered, max_counts=25000)
    sc.pp.filter_cells(train_adata_filtered, max_genes=6000)

    # Combine filtered training data with test data (unfiltered)
    combined_adata = train_adata_filtered.concatenate(
        test_adata,
        batch_key="dataset",
        batch_categories=["train", "test"],
        index_unique=None,
    )

    combined_adata = base_clustering(combined_adata)

    return combined_adata


def export(combined_adata, output_path="./submission/cluster_membership.csv"):
    # Extract test data from combined AnnData object
    test_data = combined_adata[combined_adata.obs["dataset"] == "test"]

    # Adjust clusters to 1-based indexing
    test_clusters_df = test_data.obs[["leiden"]].reset_index()
    test_clusters_df.columns = [
        "index",
        "cluster",
    ]  # Rename the original index column for clarity

    # Convert categorical clusters to integers, adjust, and revert to categorical
    test_clusters_df["cluster"] = test_clusters_df["cluster"].astype(int) + 1

    test_clusters_df.to_csv(output_path, index=True)
    print(f"File saved to: {output_path}")


def pipeline():
    train, test = load_data()
    train_clustered = train_cluster(train)
    validation_score(train_clustered)
    test_clustered = test_cluster(test, train)
    export(test_clustered)
